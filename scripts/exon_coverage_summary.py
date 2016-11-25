import os
import argparse
import ConfigParser
import json
from collections import defaultdict
import numpy as np
from gelCoverage.tools.cellbase_helper import CellbaseHelper
from gelCoverage.tools.panelapp_helper import PanelappHelper


def find_gaps(coverages, start_position, coverage_threshold):
    end = start_position + len(coverages)
    open_gap = False
    current_gap = defaultdict()
    gaps = []

    for idx, value in enumerate(coverages):
        if value < coverage_threshold and not open_gap:
            open_gap = True
            current_gap["start"] = start_position + idx
        elif value >= coverage_threshold and open_gap:
            open_gap = False
            current_gap["end"] = start_position + idx - 1
            gaps.append(current_gap)
            current_gap = defaultdict()
    if open_gap:
        current_gap["end"] = end
        gaps.append(current_gap)

    return gaps


def compute_exon_level_statistics(coverages, gc_content):
    stats = defaultdict()
    stats["total_bases"] = len(coverages)
    stats["mean"] = np.mean(coverages)
    stats["median"] = np.median(coverages)
    stats["pct75"] = np.percentile(coverages, 75)
    stats["pct25"] = np.percentile(coverages, 25)
    stats["bases_lt_3x"] = np.sum(1 for x in coverages if x < 3)
    stats["bases_lt_15x"] = np.sum(1 for x in coverages if x < 15)
    stats["bases_gte_15x"] = np.sum(1 for x in coverages if x >= 15)
    stats["bases_gte_30x"] = np.sum(1 for x in coverages if x >= 30)
    stats["bases_gte_50x"] = np.sum(1 for x in coverages if x >= 50)
    stats["percent_lt_15x"] = stats["bases_lt_15x"] / stats["total_bases"]
    stats["percent_gte_15x"] = stats["bases_gte_15x"] / stats["total_bases"]
    stats["percent_gte_30x"] = stats["bases_gte_30x"] / stats["total_bases"]
    stats["percent_gte_50x"] = stats["bases_gte_50x"] / stats["total_bases"]
    stats["gc_content"] = gc_content

    return stats


def compute_transcript_level_statistics(exons):
    stats = defaultdict()
    stats["total_bases"] = np.sum([x["total_bases"] for _, x in exons.iteritems()])
    stats["mean"] = np.mean([x["mean"] for _,x in exons.iteritems()])
    stats["weighted_median"] = np.sum(
        [x["median"] * x["total_bases"] for _,x in exons.iteritems()]) / stats["total_bases"]
    stats["weighted_pct75"] = np.sum(
        [x["pct75"] * x["total_bases"] for _,x in exons.iteritems()]) / stats["total_bases"]
    stats["weighted_pct25"] = np.sum(
        [x["pct25"] * x["total_bases"] for _,x in exons.iteritems()]) / stats["total_bases"]
    stats["bases_lt_3x"] = np.sum([x["bases_lt_3x"] for _,x in exons.iteritems()])
    stats["bases_lt_15x"] = np.sum([x["bases_lt_15x"] for _,x in exons.iteritems()])
    stats["bases_gte_15x"] = np.sum([x["bases_gte_15x"] for _,x in exons.iteritems()])
    stats["bases_gte_30x"] = np.sum([x["bases_gte_30x"] for _,x in exons.iteritems()])
    stats["bases_gte_50x"] = np.sum([x["bases_gte_50x"] for _,x in exons.iteritems()])
    stats["percent_lt_15x"] = stats["bases_lt_15x"] / stats["total_bases"]
    stats["percent_gte_15x"] = stats["bases_gte_15x"] / stats["total_bases"]
    stats["percent_gte_30x"] = stats["bases_gte_30x"] / stats["total_bases"]
    stats["percent_gte_50x"] = stats["bases_gte_50x"] / stats["total_bases"]
    stats["gc_content"] = np.sum(
        [x["gc_content"] * x["total_bases"] for _, x in exons.iteritems()]) / stats["total_bases"]

    return stats


def main():

    parser = argparse.ArgumentParser(description = 'Coverage summary for exons')
    parser.add_argument('--bw', metavar='bw', help = 'This is the bigwig file')
    parser.add_argument('--panel', metavar = 'panel',
                        help = 'The panel name or identifier in PanelApp /'
                             '(see https://bioinfo.extge.co.uk/crowdsourcing/WebServices/list_panels)',
                        default = None)
    parser.add_argument('--panel-version', metavar = 'panel_version',
                        help='The panel version',
                        default = None)
    # TODO: verify if wew want this parameter in the config file
    parser.add_argument('--panel-gene-confidence-thr', metavar = 'gene_confidence_threshold',
                        help = 'The minimum gene association confidence',
                        choices = ["HighEvidence", "MediumEvidence", "LowEvidence"],
                        default = "HighEvidence")
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
    config_filepath = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../resources/exon_coverage_summary.config')
    config.readfp(open(config_filepath))

    # Initialize CellBase helper
    cellbase_helper = CellbaseHelper(species = config.get('cellbase', 'species'),
                                     version = config.get('cellbase', 'version'),
                                     assembly = config.get('cellbase', 'assembly'),
                                     host = config.get('cellbase', 'host'),
                                     filter_basic_flag = bool(config.get('gene_filtering', 'filter_basic_flag')),
                                     filter_biotypes = config.get('gene_filtering', 'filter_biotypes').split(",")
                                     )

    ## Gets list of genes to analyse
    if args.panel is not None and args.panel_version is not None:
        # Initialize PanelApp helper
        panelapp_helper = PanelappHelper(host=config.get('panelapp', 'host'))
        # Get list of genes from PanelApp
        gene_list = panelapp_helper.get_gene_list(args.panel, args.panel_version, args.gene_confidence_threshold)
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
        exon_statistics = compute_exon_level_statistics(coverages,  gc_content)

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
            gaps = find_gaps(coverages, start, args.coverage_threshold)
            output[txid]["exons"][exon_idx]["gaps"] = gaps

    # add an aggregation of statistics at transcript level
    for transcript, content in output.iteritems():
        output[transcript]["statistics"] = compute_transcript_level_statistics(content["exons"])

    # TODO: output results in different formats
    if args.output == "json":
        print json.dumps(output, indent = 4)
    elif args.output == "tabular":
        # TODO: flatten the dict into a table
        pass


if __name__ == '__main__':
    main()