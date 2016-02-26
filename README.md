## Gel Coverage Pipeline

### .wig File Generation

On Skyscape:

```bash
module load java/jdk1.8.0_45
java -jar /genomes/software/apps/gel-coverage/gel-coverage-jar-with-dependencies.jar -bam <BAM_FILE>
```

Two files will be produced <BAM_FILE>.wig <BAM_FILE>.chr (this is the chromosome size file required for bigwig generation)

You can optionally specify an output, the path must exist:

```bash
java -jar /genomes/software/apps/gel-coverage/gel-coverage-jar-with-dependencies.jar -bam <BAM_FILE> -output <OUTPUT_PREFIX>
```

and you will get <OUTPUT_PREFIX>.wig <OUTPUT_PREFIX>.chr

You can also write to stdout for piping purposes.

```bash
java -jar /genomes/software/apps/gel-coverage/gel-coverage-jar-with-dependencies.jar -bam <BAM_FILE> -stdout
```

#### Run Time

It will take roughly 1hr for a 30x rare disease/cancer germline bam, for a 75x cancer bam it'll be about 2hrs.

### Coverage Summaries

coverage_summary.py prepares two coverage summary files (histograms), one for the whole genome and one for all exons.

takes about 30 minutes

1 processor 2.5gb ram required

```bash
python coverage_summary.py --bw <BW_FILE> --genome_n <Non-N REGIONS> --output <OUTPUT_FILE>
```

```bash
coverage_summary.R
```

```bash
gc_exon_boxplots.R
```