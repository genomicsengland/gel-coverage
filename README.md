# GEL Coverage Pipeline

## BigWig coverage analyser

This script calculates coverage statistics, using a bigwig as input. It has different execution modes.
   * `--panel`: This mode will calculate the coverage metrics for one panel.
   * `--gene-list`: This mode will calculate the coverage metrics for a list of genes.
   * `--coding-regions`: This mode will calculate the coverage metrics for a set of genes defined in the BED file provided, while avoiding any connection to CellBase.
   * `none of above`: This version will calculate the coverage metrics for all genes.

It will output statistics at exon, transcript, gene (by creating a union transcript), chromosome, analysis coding region
(this is panel, gene list or whole coding region) and whole genome. The output format is JSON.

### How to use it from commandline

This script is executed in the following way:

```
bigwig_analyser --bw <bigwig.bw> --output <output.json> --config <configuration.config> --wg-regions <non_n_region.bed> --disable-exon-stats
```

**NOTE:** The configuration file contain the basic configuration to annotate, get panels information, filter transcripts and calculate stats.
This file will be found in `/genomes/resources/genomeref/...`.


### How to use it from python


```
#Create a dictionary with the configuration
config = {
    "bw" : '/path/to/bigwig.bw',
    "panel" : None,
    "panel_version": None,
    "gene_list": None,
    "coding_regions": None,
    "coverage_threshold": 15,
    'configuration_file': '-',
    "cellbase_species": 'hsapiens',
    "cellbase_version": 'latest',
    "cellbase_assembly": 'GRCh37/GRCh38',
    "cellbase_host": '10.5.8.201:8080/cellbase-4.5.0-rc',
    "cellbase_retries": -1,
    "panelapp_host": 'bioinfo.extge.co.uk/crowdsourcing/WebServices',
    "panelapp_gene_confidence": 'HighEvidence',
    "panelapp_retries": -1,
    "panelapp_assembly": "GRCh37",
    "transcript_filtering_flags": 'basic',
    "transcript_filtering_biotypes": 'IG_C_gene,IG_D_gene,IG_J_gene,IG_V_gene,IG_V_gene,protein_coding,nonsense_mediated_decay,non_stop_decay,TR_C_gene,TR_D_gene,TR_J_gene,TR_V_gene',
    "exon_padding": 15,
    "wg_stats_enabled": True,
    "wg_regions": '/path/to/non_n_regions.bed',
    "exon_stats_enabled": False,
    "coding_region_stats_enabled": True
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

**NOTE:** When running an analysis over all genes the resulting JSON will be around 1.5GB, unless you add the flag --disable-exon-stats,
but in this case you will be missing the exon level statistics and the coverage gaps.

**NOTE:** When running an analysis in panel or gene list mode it might be useful to disable the whole genome statistics to improve performance, by using the flag --disable-wg-stats.

**NOTE:** Beware that the reference genome and chromosome notation (i.e.: chr prefix or not) should be the same in the input bigwig file and the bed file in wg-regions.


### Different configurations

The program iterates through the bigwig file twice: the first for the analysis of the coding region (panel, gene list or
full) and the second for the analysis of the whole genome.

* To **run statistics only for a panel** from exon level up to panel level, provide a panel (`panel`) and panel
version (`panel_version`) and disable the whole genome statistics (`"wg_stats_enabled": False`), while making sure that
the coding region and the exon level statistics are enabled (`"coding_region_stats_enabled": True` and `"exon_stats_enabled": True`).
    * Execution time is of some seconds or under a second for small panels. The panel of intellectual disability v1.23 having 1232 genes took 46s.
* To **run statistics only for a gene list** from exon level up to gene list level, provide a gene list (`gene_list`) instead
of panel and panel version and use the same configuration as above.
    * Execution time depends on the size of the list
* To **run statistics only for all genes in the coding region** do not provide panel (`panel`) or gene list (`gene_list`),
disable the whole genome statistics (`"wg_stats_enabled": False`) and the exon level statistics (`"exon_stats_enabled": False`)
(the output JSON will be around 4 GB if exon stats are enabled for all genes, otherwise it is around 250 MB),
while making sure that the coding region is enabled (`"coding_region_stats_enabled": True`).
    * Execution time is over 3 hours
* To **run only whole genome statistics** enable `"wg_stats_enabled": True` and disable the coding region statistics
(`"coding_region_stats_enabled": False`). The whole genome analysis might be used in combination with a bed file defining
the region to analyse (e.g.: non N regions) that is to be passed in parameter `"wg_regions": '/path/to/non_n_regions.bed'`.
This `wg_regions` can be used to calculate coverage over very specific regions, for instance Cosmic variants if they are set in
a BED file.
    * Execution time is around 1 hour

Any combination, of the previous should generate a single JSON with all the information.


### Dependencies and connection retries

The coverage module depends on two external systems: CellBase and PanelApp. CellBase is used to retrieve the genes to be analysed and the precise coordinates of each genomic region. PanelApp is used only in the panel mode to retrieve those genes belonging to a given panel. If any of these systems is down the analysis cannot run.
An exponential backoff policy will work whenever the connection to any of these two systems fails. The maximum number of retries can be configured by using the parameters `cellbase_retries` and `panelapp_retries`. If the value is -1 infinite retries will apply. These parameters are not available from the command line, but from the configuration file.

### Avoiding connection to CellBase

The connection to CellBase can be avoided by using the parameter `--coding-regions`. The BED file provided must follow the format:
```
17	67323193	67323242	ABCA5|ENST00000392676|exon1	0.68	-
17	67310455	67310571	ABCA5|ENST00000392676|exon2	0.37607	-
17	67309233	67309437	ABCA5|ENST00000392676|exon3	0.30732	-
17	67305403	67305564	ABCA5|ENST00000392676|exon4	0.33333	-
17	67304421	67304509	ABCA5|ENST00000392676|exon5	0.44944	-
```

**IMPORTANT**: the BED file must be sorted by transcript (not by gene!), otherwise coverage statistics will be erroneous!

### Exploring the results

The results of the coverage analysis are in JSON format which is intended to be machine readable, but not human readable. To explore the results the `jq` tool may ne useful (https://stedolan.github.io/jq/).

* To just print the whole JSON in a nice tabulated format use: `jq '' your.json`
* To retrieve a specific field from the json query for it like this: `jq '.results.whole_genome.stats.uneveness' your.json`

Much more complex queries can be done. See https://robots.thoughtbot.com/jq-is-sed-for-json for further insight.

**NOTE 1**: Beware that the JSONs for the coding region with coverage detail at exon level might be quite big and jq will be slow as it loads all data in memory.

**NOTE 2**: jq is installed in `bio-pp9-01` at `/genomes/software/apps/jq-1.5`


## Bed Maker

The Bed Maker builds from a list of genes the BED file defining the coding regions that the Bigwig Analyser expects.

### How to use it from commandline

```
bed_maker --config resources/bed_maker.config --output resources/test/deleteme.bed --gene-list SCN2A,SPTAN1,PLCB1,SLC25A22,SCN8A,STXBP1,PNKP
```

### How to use it from Python

```
#Create a dictionary with the configuration
config = {
    # Sets parameters from CLI
    "gene_list": None,
    "chr_prefix": False,
    "log_level": 10,
    "transcript_filtering_flags": "basic",
    "transcript_filtering_biotypes": "IG_C_gene,IG_D_gene,IG_J_gene,IG_V_gene,IG_V_gene,protein_coding,"
                                     "nonsense_mediated_decay,non_stop_decay,TR_C_gene,TR_D_gene,TR_J_gene,"
                                     "TR_V_gene",
    "cellbase_species": "hsapiens",
    "cellbase_version": "latest",
    "cellbase_assembly": "grch37",
    "cellbase_host": "bio-test-cellbase-haproxy-01.gel.zone/cellbase",
    "cellbase_retries": -1,
}

bed_maker = BedMaker(config)
bed = bed_maker.run()
# Saves the analysed region as a BED file
bed.saveas("my.bed")
```

## Mocked data

Generate mocked results for the coding region analysis using this:
```
mocked_bigwig_analyser --output-folder mocked --config resources/bigwig_analyser.GRCh38.config --coding-regions resources/coding_regions.GRCh38.bed
```

This script can be used also with panels using `--panel` and `--panel-version` or with gene lists using `--gene-list`.
Generate multiple results in the same call using `--number-results`.