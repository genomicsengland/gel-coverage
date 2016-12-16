# GEL Coverage Pipeline

## BigWig coverage analyser

This script calculates coverage statistics, using a bigwig as input. It has different execution modes.
   * `--panel`: This mode will calculate the coverage metrics for one panel.
   * `--gene-list`: This mode will calculate the coverage metrics for a list of genes.
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
    * Execution time similar as the previous
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