import ConfigParser
import argparse
import json
import os
from collections import defaultdict

import gelCoverage.stats.coverage_stats as coverage_stats
from gelCoverage.tools.cellbase_helper import CellbaseHelper
from gelCoverage.tools.panelapp_helper import PanelappHelper


def main():

    config_file = '../resources/exon_coverage_summary.config'

    parser = argparse.ArgumentParser(description = 'Coverage summary for exons')
    parser.add_argument('--bw', metavar='bw', help = 'This is the bigwig file')
    parser.add_argument('--panel', metavar = 'panel',
                        help = 'The panel name or identifier in PanelApp /'
                             '(see https://bioinfo.extge.co.uk/crowdsourcing/WebServices/list_panels)',
                        default = None)
    parser.add_argument('--panel-version', metavar = 'panel_version',
                        help='The panel version',
                        default = None)
    parser.add_argument('--genes', metavar = 'genes',
                        default = None,
                        help = 'Comma separated list of genes (HGNC gene symbols) to analyse. Will be masked by a panel')
    #parser.add_argument('--transcripts', metavar='transcripts',
    #                    help='Comma separated list of transcripts to analyse. Will be masked by a panel or a list of genes')
    #parser.add_argument('--canonical_transcript', metavar='canonical_transcript',
    #                    help='collapse transcripts to one cannonical transcript - only used when using panel or genes as input',
    #                    default=0)
    #parser.add_argument('--output', metavar='output',
    #                    help='Output format',
    #                    choices = ["json", "tabular"],
    #                    default = "json")
    parser.add_argument('--coverage-threshold', metavar='coverage_threshold',
                        help='The coverage threshold used to compute continuous gaps with low coverage (0 = disabled)',
                        default = 15)
    #TODO: add parameter to add flanking regions to genes
    #TODO: add parameter to return only information about canonical transcripts
    #parser.add_argument('--cnv', metavar='cnv', help='cnv vcf - so that losses can be indicated', default=0)

    args = parser.parse_args()

    # Reads configuration file
    config = ConfigParser.ConfigParser()
    config_filepath = os.path.join(os.path.dirname(os.path.realpath(__file__)), config_file)
    config.readfp(open(config_filepath))

    # Initialize CellBase helper
    cellbase_helper = CellbaseHelper(species = config.get('cellbase', 'species'),
                                     version = config.get('cellbase', 'version'),
                                     assembly = config.get('cellbase', 'assembly'),
                                     host = config.get('cellbase', 'host'),
                                     filter_flags = config.get('transcript_filtering', 'flags').split(","),
                                     filter_biotypes = config.get('transcript_filtering', 'biotypes').split(",")
                                     )

    ## Gets list of genes to analyse
    if args.panel is not None and args.panel_version is not None:
        # Initialize PanelApp helper
        panelapp_helper = PanelappHelper(host=config.get('panelapp', 'host'))
        # Get list of genes from PanelApp
        gene_list = panelapp_helper.get_gene_list(args.panel, args.panel_version,
                                                  config.get("panelapp", "gene_confidence"))
    elif args.genes is not None:
        gene_list = args.genes.split(",")
    #elif args.transcripts is not None:
    #    transcripts_list = args.transcripts.split(",")
        # Get list of genes from CellBase client
    #    gene_list = cellbase_helper.get_gene_list_from_transcripts(transcripts_list)
    else:
        # Warn the user as this will be time consuming
        print "WARNING: you are going to run a whole exome coverage analysis!"
        # Retrieve the list of all genes
        gene_list = cellbase_helper.get_all_genes()

    # Get genes annotations in BED format
    bed = cellbase_helper.make_exons_bed(gene_list)

    # Initialize results data structure
    output = defaultdict()
    output["parameters"] = defaultdict()
    parameters = output["parameters"]
    parameters["gap_coverage_threshold"] = args.coverage_threshold
    parameters["input_file"] = args.bw
    parameters["species"] = config.get('cellbase', 'species')
    parameters["assembly"] = config.get('cellbase', 'assembly')
    if args.panel is not None and args.panel_version is not None:
        parameters["panel"] = args.panel
        parameters["panel_version"] = args.panel_version
        parameters["gene_list"] = gene_list
    elif args.genes is not None:
        parameters["gene_list"] = gene_list

    output["results"] = defaultdict()
    results = output["results"]

    for interval in bed:
        # Reads data from BED entry
        chrom = interval.chrom
        start = int(interval.start)
        end = int(interval.end)
        gene, txid, exon_idx = interval.name.split("|")
        strand = interval.strand
        # TODO: truncate this to two decimal positions
        gc_content = interval.score

        # Queries the bigwig for a specific interval
        if start == end:
            end += 1
        coverages = args.bw.values(chrom, start, end)

        # Computes statistics at exon level
        exon_statistics = coverage_stats.compute_exon_level_statistics(coverages, gc_content)

        # Store results in data structure
        if txid not in output:
            output[txid] = defaultdict()
            output[txid]["chrom"] = chrom
            output[txid]["gene"] = gene
            output[txid]["strand"] = strand
            output[txid]["exons"] = defaultdict()
        output[gene][txid]["exons"][exon_idx] = {
            "statistics" : exon_statistics
        }

        # compute gaps
        if args.coverage_threshold > 0:
            gaps = coverage_stats.find_gaps(coverages, start, args.coverage_threshold)
            output[txid]["exons"][exon_idx]["gaps"] = gaps

    # add an aggregation of statistics at transcript level
    for transcript, content in output.iteritems():
        output[transcript]["statistics"] = coverage_stats.compute_transcript_level_statistics(content["exons"])

    # TODO: output results in different formats
    if args.output == "json":
        print json.dumps(output, indent = 4)
    elif args.output == "tabular":
        # TODO: flatten the dict into a table
        pass


if __name__ == '__main__':
    main()