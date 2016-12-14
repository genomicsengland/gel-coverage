import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;
import org.apache.commons.cli.*;
import org.ini4j.Wini;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPOutputStream;


public class SAMFinneyCoverage {
    private int BASE_MIN_QUALITY = 30;
    private int MAPPING_MIN_QUALITY = 10;
    private boolean SKIP_DUPLICATE_READS = true;
    private static boolean VERBOSE_MODE = false;
    private static boolean STDOUT = false;
    private Map<String, Integer> chr_sizes = new HashMap<String, Integer>();

    private String outfile;
    private File bam_file = null;

    public void set_output(String outfile) {
        this.outfile = outfile;
    }

    public void set_base_quality(int quality) {
        this.BASE_MIN_QUALITY = quality;
    }

    public void set_mapping_quality(int quality) { this.MAPPING_MIN_QUALITY = quality; }

    public void set_skip_duplicate_reads(boolean skip) { this.SKIP_DUPLICATE_READS = skip; }

    public void set_bam_file(File f) {
        bam_file = f;
    }

    public void set_verbose(boolean v) {
        VERBOSE_MODE = v;
    }

    public void set_stdout(boolean v) {
        SAMFinneyCoverage.STDOUT = v;
    }

    public void print_thresholds() {
        System.err.println("Minimum base quality: " + BASE_MIN_QUALITY);
        System.err.println("Minimum mapping quality: " + MAPPING_MIN_QUALITY);
    }

    public void find_coverage() throws IOException {

        // Open the input BAM/SAM file
        SAMFileReader sfr = new SAMFileReader(bam_file);
        // Open the output WIG file
        if (outfile == null) outfile = bam_file.getName() + ".wig";
        WorkingFile wf = null;
        FileOutputStream fos;
        if (SAMFinneyCoverage.STDOUT) {
            fos = new FileOutputStream(FileDescriptor.out);
        } else {
            wf = new WorkingFile(outfile);
            fos = new FileOutputStream(wf);
        }
        OutputStream os = new BufferedOutputStream(fos);
        PrintStream ps = new PrintStream(os);
        // Read chromosomes and their length from BAM header
        List<String> chr_labels = new ArrayList<String>();
        SAMFileHeader h = sfr.getFileHeader();
        SAMSequenceDictionary dict = h.getSequenceDictionary();
        WorkingFile chr_wf = new WorkingFile(outfile.replace(".wig", ".chr"));
        PrintStream chr_os = new PrintStream(new BufferedOutputStream(new FileOutputStream(chr_wf)));
        for (SAMSequenceRecord ssr : dict.getSequences()) {
            String chromosome = ssr.getSequenceName();
            chr_labels.add(chromosome);
            Integer ref_len = new Integer(ssr.getSequenceLength());
            chr_sizes.put(chromosome, ref_len);
            chr_os.println(chromosome + "\t" + ref_len);
        }
        chr_os.close();
        chr_wf.finish();
        // Create WIG coverages file
        byte[] read, baseQuals;
        int read_i, ref_i, i, end;
        int null_qual = 0;
        int qual_length_problem = 0;
        int qual_bounds_problem = 0;
        for (String chromosome: chr_labels) {
            long startTime = System.currentTimeMillis();
            int coverage_len = chr_sizes.get(chromosome).intValue();
            int[] coverage = new int[coverage_len];
            //
            //  generate coverage:
            //
            boolean has_coverage = false;
            long record_count = 0;
            if (chromosome == null) {
                System.err.println("WTF: no .bam mappings vs. chr index " + chromosome);  // debug
            } else {
                if (VERBOSE_MODE) System.err.println("query=" + chromosome + " size=" + coverage_len);
                CloseableIterator<SAMRecord> iterator = sfr.queryOverlapping(chromosome, 1, coverage_len);
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
                ps.println("fixedStep chrom=" + chromosome + " start=1 step=1");
                for (i = 0; i < coverage_len; i++) {
                    ps.println(Integer.toString(coverage[i]));
                }
            }
            long stopTime = System.currentTimeMillis();
            System.err.println("CHR " + chromosome + " in " + (stopTime - startTime)/1000 + " seconds.");
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

    public static void main(String[] argv) throws IOException {
        SAMFileReader.setDefaultValidationStringency(SAMFileReader.getDefaultValidationStringency().SILENT);
        File bam_file = null;
        SAMFinneyCoverage sc = new SAMFinneyCoverage();
        // Defines command line
        String usage = "bam2wig -bam <BAM> -wig <WIG>";
        Options options = new Options();
        options.addOption("b", "bam", true, "The input BAM file");
        options.addOption("w", "wig", true, "The output WIG file");
        options.addOption("c", "config", true, "The configuration file");
        options.addOption("h", "help", false, "This help menu");
        options.addOption("v", "verbose", false, "Outputs lengthy logs");
        HelpFormatter formatter = new HelpFormatter();
        try {
            // Parses the command line
            CommandLineParser parser = new BasicParser();
            CommandLine cmd = parser.parse(options, argv);
            // Shows help message
            if (cmd.hasOption("h")) {
                formatter.printHelp(usage, options);
                System.exit(1);
            }
            // Sets verbose mode
            if (cmd.hasOption("v")) {
                sc.set_verbose(true);
            }
            // Reads input BAM
            if (cmd.hasOption("b")) {
                String bam_input = cmd.getOptionValue("b");
                sc.set_bam_file(new File(bam_input));
            } else {
                // why is it so difficult to implement a required parameter in Java libs!?
                System.out.println("Missing argument -bam");
                formatter.printHelp(usage, options);
                System.exit(1);
            }
            // Opens output WIG
            if (cmd.hasOption("w")) {
                String wig_output = cmd.getOptionValue("w");
                if (wig_output.equals("-")) {
                    sc.set_stdout(true);
                } else {
                    sc.set_output(wig_output);
                }
            } else {
                System.out.println("Missing argument -wig");
                formatter.printHelp(usage, options);
                System.exit(1);
            }
            // Opens configuration file
            if (cmd.hasOption("c")) {
                String config_file = cmd.getOptionValue("c");
                Wini ini = new Wini(new File(config_file));
                sc.set_base_quality(ini.get("bam2wig", "base_quality", int.class));
                sc.set_mapping_quality(ini.get("bam2wig", "mapping_quality", int.class));
                sc.set_skip_duplicate_reads(ini.get("bam2wig", "filter_duplicates", boolean.class));
            } else {
                System.out.println("Missing argument -config");
                formatter.printHelp(usage, options);
                System.exit(1);
            }
        }
        catch (ParseException e) {
            System.out.println( "Incorrect call:" + e.getMessage() );
            formatter.printHelp(usage, options);
        }
        // Runs the creation of the wig file
        try {
            sc.print_thresholds();
            sc.find_coverage();
        } catch (Exception e) {
            System.err.println("ERROR: " + e);  // debug
            e.printStackTrace();
            System.exit(1);
        }
    }
}

