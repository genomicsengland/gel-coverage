import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPOutputStream;

public class SAMFinneyCoverage {
    private int BASE_MIN_QUALITY = 30;
    private int MAPPING_MIN_QUALITY = 10;
    private static boolean VERBOSE_MODE = false;
    private static boolean SKIP_DUPLICATE_READS = true;
    private static boolean STDOUT = false;
    public static int CHR_I_START = 0;
    public static int CHR_I_END = 24;

    private List<Integer> chr_sizes = new ArrayList();

    private String outfile;
    private File bam_file = null;

    public void set_output(String outfile) {
        this.outfile = outfile + ".wig";
    }

    public void set_base_quality(int quality) {
        this.BASE_MIN_QUALITY = quality;
    }

    public void set_mapping_quality(int quality) {
        this.MAPPING_MIN_QUALITY = quality;
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

    public void set_stdout(boolean v) {
        this.STDOUT = v;
    }

    public void print_thresholds() {
        System.err.println("Minimum base quality: " + BASE_MIN_QUALITY);
        System.err.println("Minimum mapping quality: " + MAPPING_MIN_QUALITY);
        if (SKIP_DUPLICATE_READS) {
            System.err.println("Skipping duplicated reads.\n");
        } else {
            System.err.println("Considering duplicated reads.\n");
        }
    }

    public void find_coverage() throws IOException {
        SAMFileReader sfr = new SAMFileReader(bam_file);

        if (outfile == null) outfile = bam_file.getName() + ".wig";
        WorkingFile wf = null;
        FileOutputStream fos;

        if (STDOUT) {
            fos = new FileOutputStream(FileDescriptor.out);
        } else {
            wf = new WorkingFile(outfile);
            fos = new FileOutputStream(wf);
        }

        OutputStream os = new BufferedOutputStream(fos);


        /*if (outfile.indexOf(".gz") == outfile.length() - 3) {
            System.err.println("Generating GZ wig file: " + outfile);
            os = new GZIPOutputStream(os);
        }*/

        PrintStream ps = new PrintStream(os);

        List<String> chr_labels = new ArrayList<String>();
        List<Chromosome> chroms = new ArrayList<Chromosome>();

        //
        //  find .bam names for each chr and verify sizes (genome version):
        //
        SAMFileHeader h = sfr.getFileHeader();
        SAMSequenceDictionary dict = h.getSequenceDictionary();

        WorkingFile chr_wf = new WorkingFile(outfile.replace(".wig", ".chr"));
        PrintStream chr_os = new PrintStream(new BufferedOutputStream(new FileOutputStream(chr_wf)));

        for (SAMSequenceRecord ssr : dict.getSequences()) {
            String ref_name = ssr.getSequenceName();
            Chromosome c = Chromosome.valueOfString(ref_name);
            if (c != null) {
                int ref_len = ssr.getSequenceLength();
                chr_labels.add(ref_name);
                chroms.add(c);
                chr_sizes.add(ref_len);
                chr_os.println(c.toString() + "\t" + ref_len);
            }
        }
        chr_os.close();
        chr_wf.finish();

        //
        //  generate coverage
        //
        byte[] read, baseQuals;
        int read_i, ref_i, i, end;
        int null_qual = 0;
        int qual_length_problem = 0;
        int qual_bounds_problem = 0;

        for (int chr_i = CHR_I_START; chr_i <= CHR_I_END; chr_i++) {
            long startTime = System.currentTimeMillis();

            int coverage_len = chr_sizes.get(chr_i);
            int[] coverage = new int[coverage_len];

            //
            //  generate coverage:
            //
            String ref_name = chr_labels.get(chr_i);
            boolean has_coverage = false;
            long record_count = 0;
            if (ref_name == null) {
                System.err.println("WTF: no .bam mappings vs. chr index " + chr_i);  // debug
            } else {
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
                    if (sr.getMappingQuality() < MAPPING_MIN_QUALITY) continue;

                    has_coverage = true;

                    read = sr.getReadBases();
                    baseQuals = sr.getBaseQualities();

                    if (baseQuals.length == 0) {
                        if (null_qual++ == 0)
                            System.err.println("ERROR: 0-length qual array for " + sr.getReadName() + " (only warning, counts at end of run)");  // debug
                        // complain only once
                        continue;
                    } else if (read.length != baseQuals.length) {
                        if (qual_length_problem++ == 0)
                            System.err.println("ERROR: base/qual length mismatch for " + sr.getReadName() + "(" + read.length + " vs " + baseQuals.length + "; only warning, counts at end of run)");  // debug
                        // complain only once
                        continue;
                    }

                    for (AlignmentBlock ab : sr.getAlignmentBlocks()) {
                        read_i = ab.getReadStart() - 1;
                        ref_i = ab.getReferenceStart() - 1;

                        for (i = read_i, end = read_i + ab.getLength(); i < end; i++, ref_i++) {
                            if (i >= baseQuals.length) {
                                // read index out of quality array range
                                // very odd since we already check lengths above
                                if (qual_bounds_problem++ == 0) {
                                    // complain only once
                                    System.err.println("quality index out of range for " + sr.getReadName() + " align_start=" + sr.getAlignmentStart() + " align_end=" + sr.getAlignmentEnd() + " ref_base=" + (ref_i + 1) + " block count:" + sr.getAlignmentBlocks().size() + " block_start=" + ab.getReadStart() + " block_len=" + ab.getLength() + " read_index=" + i + " read_len=" + read.length + " qual_len=" + baseQuals.length);  // debug
                                    System.err.println("qual array:");  // debug
                                    for (int q = 0; q < baseQuals.length; q++) {
                                        System.err.println(q + "=" + baseQuals[q]);  // debug
                                    }

                                }
                            } else {
                                if (baseQuals[i] >= BASE_MIN_QUALITY && ref_i >= 0 && ref_i < coverage_len) {
                                    coverage[ref_i]++;
                                }
                            }
                        }
                    }
                }
                iterator.close();
            }
            if (VERBOSE_MODE) System.err.println("records returned: " + record_count);

            //
            //  write results:
            //
            // .wig format
            if (has_coverage) {
                ps.println("fixedStep chrom=" + chroms.get(chr_i).toString() + " start=1 step=1");
                for (i = 0; i < coverage_len; i++) {
                    ps.println(Integer.toString(coverage[i]));
                }
            }

            // Run some code;
            long stopTime = System.currentTimeMillis();
            System.err.println("CHR " + ref_name + " in " + (stopTime - startTime)/1000 + " seconds.");
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

        ps.close();
        if (wf != null) {
            wf.finish();
        }

    }

    public static void main(String[] argv) {
        SAMFileReader.setDefaultValidationStringency(SAMFileReader.getDefaultValidationStringency().SILENT);

        File bam_file = null;
        SAMFinneyCoverage sc = new SAMFinneyCoverage();

        for (int i = 0; i < argv.length; i++) {
            if (argv[i].equals("-bam")) {
                bam_file = new File(argv[++i]);
            } else if (argv[i].equals("-verbose")) {
                sc.set_verbose(true);
            } else if (argv[i].equals("-count-duplicates")) {
                sc.set_skip_duplicates(false);
            } else if (argv[i].equals("-base-quality")) {
                sc.set_base_quality(Integer.parseInt(argv[++i]));
            } else if (argv[i].equals("-mapping-quality")) {
                sc.set_mapping_quality(Integer.parseInt(argv[++i]));
            } else if (argv[i].equals("-output")) {
                sc.set_output(argv[++i]);
            } else if (argv[i].equals("-stdout")) {
                sc.set_stdout(true);
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
                sc.print_thresholds();
                sc.find_coverage();
            } catch (Exception e) {
                System.err.println("ERROR: " + e);  // debug
                e.printStackTrace();
                System.exit(1);
            }
        }

    }
}

