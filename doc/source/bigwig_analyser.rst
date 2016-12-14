bigwig_analyser
===============

This is script calculate coverage for coding regions, using a bigwig as input. It has different execution modes.
   * `--panel`: This mode will calculate the coverage metrics for one panel.
   * `--gene-list`: This mode will calculate the coverage metrics for a list of genes.
   * `none of above`: This version will calculate the coverage metrics for all genes.

How to use it from commandline
------------------------------

This script is executed in the following way:

.. code-block:: bash

    bigwig_analyser --bw <bigwig.bw> --output <otuput.json> --config <configuration.config> --wg-regions <non_n_region.bed> --disable-exon-stats

.. note::

    The configuration file contain the basic configuration to annotate, get panels information, filter transcripts and calculate stats.
    This file will be found in `/genomes/resources/genomeref/...`.


How to use it from python
-------------------------

.. code-block:: python

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



.. note::

    Please note that this process is highly dependent on the reference genome, use a different assembly o version assembly
    will produce wrong results


