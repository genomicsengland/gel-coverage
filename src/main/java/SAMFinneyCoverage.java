import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.Deflater;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

public class SAMFinneyCoverage {
    public static int MIN_QUALITY = 30;
    private static boolean VERBOSE_MODE = false;
    private static boolean WIG_MODE = false;
    private static boolean SKIP_DUPLICATE_READS = true;
    private static final int MAX_UNSIGNED_SHORT = 0xffff;
    private boolean STDOUT_MODE = false;
    public static int CHR_I_START = 0;
    public static int CHR_I_END = 24;

    private List<Integer> chr_sizes = new ArrayList();

    private String outfile;
    PrintStream ps = System.err;
    private File bam_file = null;

    public void set_outfile(String outfile) {
        this.outfile = outfile;
    }

    public void set_bam_file(File f) {
        bam_file = f;
    }

    public void set_verbose(boolean v) {
        VERBOSE_MODE = v;
    }

    public void set_skip_duplicates(boolean v) {
        SKIP_DUPLICATE_READS = v;
    }

    public void set_wig(boolean v) {
        WIG_MODE = v;
    }

    public void set_stdout_mode(boolean v) {
        STDOUT_MODE = v;
    }

    public void find_coverage() throws IOException {
        System.err.println("skipping duplicate reads?: " + SKIP_DUPLICATE_READS);

        SAMFileReader sfr = new SAMFileReader(bam_file);

        WorkingFile wf = null;

        if (!STDOUT_MODE) {
            if (outfile == null) outfile = bam_file.getName() + (WIG_MODE ? ".wig" : ".coverage");
            System.err.println("outfile=" + outfile);  // debug
            wf = new WorkingFile(outfile);
        }

        DataOutputStream dos = null;
        PrintStream ps = null;

        boolean allow_non_hg18 = WIG_MODE;

        if (WIG_MODE) {
            //      ps = new PrintStream(new BufferedOutputStream(new FileOutputStream(wf)));
            if (STDOUT_MODE) {
                //	ps = System.out;
                // System.out is line-buffered (SLOW)

                FileOutputStream fdout = new FileOutputStream(FileDescriptor.out);
                BufferedOutputStream bos = new BufferedOutputStream(fdout, 4096);
                ps = new PrintStream(bos, false);
            } else {
                OutputStream os = new BufferedOutputStream(new FileOutputStream(wf));

                if (outfile.indexOf(".gz") == outfile.length() - 3) {
                    System.err.println("generating GZ wig file");  // debug
                    if (false) {
                        // unfinished/abandoned
                        // unsure if this output would work with wigToBigWig anyway
                        System.err.println("using zip");  // debug
                        ZipOutputStream zos = new ZipOutputStream(os);
                        zos.setMethod(ZipOutputStream.DEFLATED);
                        zos.setLevel(Deflater.BEST_SPEED);

                        ZipEntry ze = new ZipEntry("wig");
                        ze.setMethod(ZipOutputStream.DEFLATED);

                        os = zos;
                    } else if (false) {
                        System.err.println("using experimental gzip");  // debug
                        int buffer_size = 4096;
                        os = new GZIPOutputStreamEx(os, buffer_size, Deflater.BEST_SPEED);
                        // does this help speed any?
                    } else {
                        //	    System.err.println("using standard gzip");  // debug
                        os = new GZIPOutputStream(os);
                    }
                }
                ps = new PrintStream(os);
            }
        } else {
            dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(wf)));
        }

        List<String> chr_labels = new ArrayList();
        List<Chromosome> chroms = new ArrayList();

        //
        //  find .bam names for each chr and verify sizes (genome version):
        //
        SAMFileHeader h = sfr.getFileHeader();
        SAMSequenceDictionary dict = h.getSequenceDictionary();
        boolean ref_mismatch = false;
        for (SAMSequenceRecord ssr : dict.getSequences()) {
            String ref_name = ssr.getSequenceName();
            Chromosome c = Chromosome.valueOfString(ref_name);
            if (c != null) {
                //
                // .bam reference sequence name matches a known chromosome
                //
                int c_index = c.toInt() - 1;
                int ref_len = ssr.getSequenceLength();
                chr_labels.add(ref_name);
                chroms.add(c);
                chr_sizes.add(ref_len);
                /*

                int c_index = c.toInt() - 1;
                chr_labels[c_index] = ref_name;
                chroms[c_index] = c;
                int ref_len = ssr.getSequenceLength();
                int wanted_len = chr_sizes[c_index];
                if (ref_len != wanted_len) {
                    System.err.println("chr=" + c + " bam_ref_len=" + ref_len + " hg18_chr_len=" + wanted_len);
                    ref_mismatch = true;
                    if (allow_non_hg18) chr_sizes[c_index] = ref_len;
                    // rewrite to actual
                }
                */
                // reference name the .bam file uses for this chromosome
            }
        }

        if (ref_mismatch) {
            String tag = allow_non_hg18 ? "NOTE" : "ERROR";
            System.err.print(tag + ": .bam reference sequences don't appear to be hg18");
            if (!allow_non_hg18) {
                System.err.println(", not generating coverage");  // debug
                return;
            } else {
                System.err.println("");  // debug
            }
        }

        //
        //  generate coverage
        //
        byte[] read, quals;
        int read_i, ref_i, i, end;
        int null_qual = 0;
        int qual_length_problem = 0;
        int qual_bounds_problem = 0;

        for (int chr_i = CHR_I_START; chr_i <= CHR_I_END; chr_i++) {
            int coverage_len = chr_sizes.get(chr_i);
            //      System.err.println("allocating " + coverage_len);  // debug

            int[] coverage = new int[coverage_len];
            //      Arrays.fill(coverage, 0);
            // not necessary?

            //
            //  generate coverage:
            //
            String ref_name = chr_labels.get(chr_i);
            boolean has_coverage = false;
            long record_count = 0;
            if (ref_name == null) {
                System.err.println("WTF: no .bam mappings vs. chr index " + chr_i);  // debug
            } else {
                System.err.println("CHR " + ref_name);  // debug
                if (VERBOSE_MODE) System.err.println("query=" + ref_name + " size=" + coverage_len);
                CloseableIterator<SAMRecord> iterator = sfr.queryOverlapping(ref_name, 1, coverage_len);
                SAMRecord sr;
                while (iterator.hasNext()) {
                    sr = iterator.next();
                    record_count++;
                    if (VERBOSE_MODE) {
                        System.err.print(
                                sr.getReadName() +
                                        " unmapped=" + sr.getReadUnmappedFlag() +
                                        " dup=" + sr.getDuplicateReadFlag()
                        );
                        if (!sr.getReadUnmappedFlag()) {
                            System.err.print(" pos=" + sr.getAlignmentStart() +
                                    "-" + sr.getAlignmentEnd());  // debug
                        }

                        System.err.println("");  // debug
                    }

                    if (sr.getReadUnmappedFlag()) continue;
                    if (SKIP_DUPLICATE_READS && sr.getDuplicateReadFlag()) continue;

                    has_coverage = true;

                    read = sr.getReadBases();
                    quals = sr.getBaseQualities();
                    if (quals.length == 0) {
                        if (null_qual++ == 0)
                            System.err.println("ERROR: 0-length qual array for " + sr.getReadName() + " (only warning, counts at end of run)");  // debug
                        // complain only once
                        continue;
                    } else if (read.length != quals.length) {
                        if (qual_length_problem++ == 0)
                            System.err.println("ERROR: base/qual length mismatch for " + sr.getReadName() + "(" + read.length + " vs " + quals.length + "; only warning, counts at end of run)");  // debug
                        // complain only once
                        continue;
                    }

                    List<CigarElement> cigarElements = sr.getCigar().getCigarElements();
                    for (AlignmentBlock ab : sr.getAlignmentBlocks()) {
                        read_i = ab.getReadStart() - 1;
                        ref_i = ab.getReferenceStart() - 1;

                        //	    System.err.println("ref_i="+ref_i + " cov_len=" + coverage.length + " readlen=" + read.length + " qlen=" + quals.length + " read_i=" + read_i + " blen=" + ab.getLength() + " read=" + sr.getReadName() + " alignStart=" + sr.getAlignmentStart() + " cigar=" + sr.getCigar());  // debug

                        for (i = read_i, end = read_i + ab.getLength(); i < end; i++, ref_i++) {
                            if (i >= quals.length) {
                                // read index out of quality array range
                                // very odd since we already check lengths above
                                if (qual_bounds_problem++ == 0) {
                                    // complain only once
                                    System.err.println("quality index out of range for " + sr.getReadName() + " align_start=" + sr.getAlignmentStart() + " align_end=" + sr.getAlignmentEnd() + " ref_base=" + (ref_i + 1) + " block count:" + sr.getAlignmentBlocks().size() + " block_start=" + ab.getReadStart() + " block_len=" + ab.getLength() + " read_index=" + i + " read_len=" + read.length + " qual_len=" + quals.length);  // debug
                                    System.err.println("qual array:");  // debug
                                    for (int q = 0; q < quals.length; q++) {
                                        System.err.println(q + "=" + quals[q]);  // debug
                                    }

                                }
                            } else {
                                int cigarRef = 0;
                                int absPosition = i - read_i + 1;
                                int cigarLength = cigarElements.get(cigarRef).getLength();
                                if (sr.getReadNegativeStrandFlag()) { // Negative strand
                                    cigarRef = cigarElements.size() - 1;
                                    absPosition = ab.getLength() - (i - read_i);
                                    cigarLength = sr.getCigarLength() - cigarElements.get(cigarRef).getLength();
                                    while (absPosition < cigarLength) {
                                        cigarRef -= 1;
                                        cigarLength -= cigarElements.get(cigarRef).getLength();
                                    }
                                } else { // Positive strand
                                    while (absPosition > cigarLength) {
                                        cigarRef += 1;
                                        cigarLength += cigarElements.get(cigarRef).getLength();
                                    }
                                }

                                if (quals[i] >= MIN_QUALITY && ref_i >= 0 && ref_i < coverage_len
                                        && !cigarElements.get(cigarRef).getOperator().isClipping()) {
                                    coverage[ref_i]++;

                                }
                            }
                        }
                    }
                }
                iterator.close();
            }
            if (VERBOSE_MODE) System.err.println("records returned: " + record_count);

            if (false && chr_i == 21) {
                // chr22:15830145
                System.err.println("DEBUG: adding coverage spike");  // debug
                coverage[15830144] = 49152;
                // 0xc000 (also endian test, vs 0x00c0 [192])
            }

            //
            //  write results:
            //

            if (WIG_MODE) {
                // .wig format
                if (has_coverage) {
                    ps.println("fixedStep chrom=" + chroms.get(chr_i).toString() + " start=1 step=1");
                    for (i = 0; i < coverage_len; i++) {
                        ps.println(i + "\t" + Integer.toString(coverage[i]));
                    }
                }
            } else {
                int c;
                for (i = 0; i < coverage_len; i++) {
                    c = coverage[i];
                    if (c > MAX_UNSIGNED_SHORT) c = MAX_UNSIGNED_SHORT;
                    //	System.err.println("writing " + c);  // debug

                    if (false) {
                        // big endian
                        dos.write(c >> 8);
                        dos.write(c);
                    } else {
                        // little endian
                        dos.write(c);
                        dos.write(c >> 8);
                    }
                }
            }

        }

        if (null_qual > 0) {
            System.err.println("ERROR: " + null_qual + " reads w/0-length base qualities");  // debug
        }

        if (qual_length_problem > 0) {
            System.err.println("ERROR: " + qual_length_problem + " reads w/base/quality array length mismatches");  // debug
        }

        if (qual_bounds_problem > 0) {
            System.err.println("ERROR: " + qual_bounds_problem + " instances of AlignmentBlock read index > base quality array size");  // debug
        }

        if (dos != null) dos.close();
        if (ps != null) ps.close();
        if (wf != null) wf.finish();

    }

    public static void main(String[] argv) {
        SAMFileReader.setDefaultValidationStringency(SAMFileReader.getDefaultValidationStringency().SILENT);
        // STFU

        File bam_file = null;
        SAMFinneyCoverage sc = new SAMFinneyCoverage();

        for (int i = 0; i < argv.length; i++) {
            if (argv[i].equals("-bam")) {
                bam_file = new File(argv[++i]);
            } else if (argv[i].equals("-verbose")) {
                sc.set_verbose(true);
            } else if (argv[i].equals("-wig")) {
                sc.set_wig(true);
            } else if (argv[i].equals("-of")) {
                sc.set_outfile(new String(argv[++i]));
            } else if (argv[i].equals("-stdout")) {
                sc.set_stdout_mode(true);
            } else if (argv[i].equals("-count-duplicates")) {
                sc.set_skip_duplicates(false);
            } else if (argv[i].equals("-chr-i-start")) {
                SAMFinneyCoverage.CHR_I_START = Integer.parseInt(argv[++i]);
            } else if (argv[i].equals("-chr-i-end")) {
                SAMFinneyCoverage.CHR_I_END = Integer.parseInt(argv[++i]);
            } else {
                System.err.println("error: unknown switch " + argv[i]);  // debug
                System.exit(1);
            }
        }

        String error = null;
        if (bam_file == null) {
            error = "specify -bam [file]";
        }

        if (error != null) {
            System.err.println("ERROR: " + error);  // debug
        } else if (bam_file != null) {
            try {
                sc.set_bam_file(bam_file);
                sc.find_coverage();
            } catch (Exception e) {
                System.err.println("ERROR: " + e);  // debug
                e.printStackTrace();
                System.exit(1);
            }
        }

    }
}

