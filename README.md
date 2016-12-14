## Gel Coverage Pipeline

### .wig File Generation

Build and compile:

```
mvn compile
mvn package
```

On Skyscape:

```bash
module load java/jdk1.8.0_45
java -jar /genomes/software/apps/gel-coverage/target/bam2wig-jar-with-dependencies.jar --bam <BAM_FILE> --wig <WIG_FILE> --config <CONFIG_FILE>
```

Two files will be produced <WIG_FILE> with extension .wig and a chromosomes file replacing the extension by .chr (this is the chromosome size file required for bigwig generation)

You can also write to stdout for piping purposes.

```bash
java -jar /genomes/software/apps/gel-coverage/bam2wig-jar-with-dependencies.jar --bam <BAM_FILE> --wig - --config <CONFIG_FILE>
```

#### Run Time

It will take roughly 1hr for a 30x rare disease/cancer germline bam, for a 75x cancer bam it'll be about 2hrs. (To be validated...)

#### Piping to generate bigWig

You can pipe the wig to generate a bigWig without having to generate the wig 1st

```bash
module load java/jdk1.8.0_45
java -jar /genomes/software/apps/gel-coverage/target/bam2wig-jar-with-dependencies.jar --bam /genomes/by_date/2014-11-17/RAREP00885/LP2000275-DNA_D09/Assembly/LP2000275-DNA_D09.bam --wig - --config <CONFIG_FILE> | /accelrys/apps/gel/toolkit/bin/linux64/ucsc/wigToBigWig stdin /genomes/analysis/rare_disease/coverage/LP2000275-DNA_D09.chr /genomes/analysis/rare_disease/coverage/LP2000275-DNA_D09.bw
```

### BigWig coverage analyser

This script calculates coverage statistics, using a bigwig as input. It has different execution modes.
   * `--panel`: This mode will calculate the coverage metrics for one panel.
   * `--gene-list`: This mode will calculate the coverage metrics for a list of genes.
   * `none of above`: This version will calculate the coverage metrics for all genes.

It will output statistics at exon, transcript, gene (by creating a union transcript), chromosome, analysis coding region
(this is panel, gene list or whole coding region) and whole genome. The output format is JSON.

#### How to use it from commandline

This script is executed in the following way:

```
bigwig_analyser --bw <bigwig.bw> --output <output.json> --config <configuration.config> --wg-regions <non_n_region.bed> --disable-exon-stats
```

**NOTE:** The configuration file contain the basic configuration to annotate, get panels information, filter transcripts and calculate stats.
This file will be found in `/genomes/resources/genomeref/...`.


#### How to use it from python


```
#Create a dictionary with the configuration
config = {
    "bw" : '/path/to/bigwig.bw',
    "panel" : None,
    "panel_version": None,
    "gene_list": None,
    "coverage_threshold": 15,
    'configuration_file': '-',
    "cellbase_species": 'hsapiens',
    "cellbase_version": 'latest',
    "cellbase_assembly": 'GRCh37/GRCh38',
    "cellbase_host": '10.5.8.201:8080/cellbase-4.5.0-rc',
    "panelapp_host": 'bioinfo.extge.co.uk/crowdsourcing/WebServices',
    "panelapp_gene_confidence": 'HighEvidence',
    "transcript_filtering_flags": 'basic',
    "transcript_filtering_biotypes": 'IG_C_gene,IG_D_gene,IG_J_gene,IG_V_gene,IG_V_gene,protein_coding,nonsense_mediated_decay,non_stop_decay,TR_C_gene,TR_D_gene,TR_J_gene,TR_V_gene',
    "exon_padding": 15,
    "wg_stats_enabled": True,
    "wg_regions": '/path/to/non_n_regions.bed',
    "exon_stats_enabled": False
}

gel_coverage_engine = GelCoverageRunner(config)
(results, bed) = gel_coverage_engine.run()
# Prints output to stdout
with codecs.open(args.output, 'w', 'utf8') as output_file:
    output_file.write(
        ujson.dumps(
            results,
            ensure_ascii=False
        )
    )
# Saves the analysed region as a BED file
bed.saveas(args.output + ".bed")
```



**NOTE:** Please note that this process is highly dependent on the reference genome, use a different assembly o version assembly will produce wrong results.

**NOTE:** When running an analysis over all genes the resulting JSON will be around 4GB, unless you add the flag --disable-exon-stats,
but in this case you will be missing the exon level statistics and the coverage gaps.

**NOTE:** When running an analysis in panel or gene list mode it might be useful to disable the whole genome statistics to improve performance, by using the flag --disable-wg-stats.

**NOTE:** Beware that the reference genome and chromosome notation (i.e.: chr prefix or not) should be the same in the input bigwig file and the bed file in wg-regions.


