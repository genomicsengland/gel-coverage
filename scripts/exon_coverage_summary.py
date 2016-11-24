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

def initialize_statistics():
    stats = defaultdict()
    stats["total_bases"] = 0
    stats["mean"] = 0
    stats["median"] = 0
    stats["pct75"] = 0
    stats["pct25"] = 0
    stats["bases_lt_3x"] = 0
    stats["bases_lt_15x"] = 0
    stats["bases_gte_15x"] = 0
    stats["bases_gte_30x"] = 0
    stats["bases_gte_50x"] = 0
    stats["percent_lt_15x"] = 0
    stats["percent_gte_15x"] = 0
    stats["percent_gte_30x"] = 0
    stats["percent_gte_50x"] = 0

    return stats

def compute_exon_level_statistics(coverages):
    stats = initialize_statistics()
    stats["total_bases"] = len(coverages)
    stats["mean"] = np.mean(coverages)
    stats["median"] = np.median(coverages)
    stats["pct75"] = np.percentile(coverages, 75)
    stats["pct25"] = np.percentile(coverages, 25)
    stats["bases_lt_3x"] = sum(1 for x in coverages if x < 3)
    stats["bases_lt_15x"] = sum(1 for x in coverages if x < 15)
    stats["bases_gte_15x"] = sum(1 for x in coverages if x >= 15)
    stats["bases_gte_30x"] = sum(1 for x in coverages if x >= 30)
    stats["bases_gte_50x"] = sum(1 for x in coverages if x >= 50)
    stats["percent_lt_15x"] = stats["bases_lt_15x"] / stats["total_bases"]
    stats["percent_gte_15x"] = stats["bases_gte_15x"] / stats["total_bases"]
    stats["percent_gte_30x"] = stats["bases_gte_30x"] / stats["total_bases"]
    stats["percent_gte_50x"] = stats["bases_gte_50x"] / stats["total_bases"]

    return stats

def compute_transcript_level_statistics(exons):
    transcript_stats = initialize_statistics()
    for _, exon in exons.iteritems():
        exon_stats = exon["statistics"]
        if __name__ == '__main__':
            if transcript_stats["total_bases"] == 0:
                transcript_stats["total_bases"] = exon_stats["total_bases"]
            else:
                transcript_stats["total_bases"] = transcript_stats["total_bases"] + exon_stats["total_bases"]
            # TODO: ...

    return transcript_stats


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
    parser.add_argument('--panel-gene-confidence-thr', metavar = 'gene_confidence_threshold',
                        help = 'The minimum gene association confidence',
                        choices = ["HighEvidence", "MediumEvidence", "LowEvidence"],
                        default = "HighEvidence")
    parser.add_argument('--genes', metavar = 'genes',
                        default = None,
                        help = 'Comma separated list of genes to analyse. Will be masked by a panel')
    parser.add_argument('--transcripts', metavar='transcripts',
                        help='Comma separated list of transcripts to analyse. Will be masked by a panel or a list of genes')
    #parser.add_argument('--canonical_transcript', metavar='canonical_transcript',
    #                    help='collapse transcripts to one cannonical transcript - only used when using panel or genes as input',
    #                    default=0)
    parser.add_argument('--output', metavar='output',
                        help='Output format',
                        choices = ["json", "tabular"],
                        default = "json")
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
                                                     filter_basic_flag = bool(config.get('transcripts', 'filter_basic_flag')),
                                                     filter_biotypes = config.get('transcripts', 'filter_biotypes').split(",")
                                                     )

    ## Gets list of genes to analyse
    if args.panel is not None and args.panel_version is not None:
        # Initialize PanelApp helper
        panelapp_helper = PanelappHelper(host=config.get('panelapp', 'host'))
        # get list of genes from PanelApp
        gene_list = panelapp_helper.get_gene_list(args.panel, args.panel_version, args.gene_confidence_threshold)
    elif args.genes is not None:
        gene_list = args.genes.split(",")
    elif args.transcripts is not None:
        transcripts_list = args.transcripts.split(",")
        # get list of genes from CellBase client
        gene_list = cellbase_helper.get_gene_list_from_transcripts(transcripts_list)
    else:
        # TODO: retrieve the list of all genes
        gene_list = None

    # Get genes annotations in BED format
    bed = cellbase_helper.make_exons_bed(gene_list)

    # Initialize results data structure
    results = defaultdict()

    for interval in bed:
        # Reads data from BED entry
        chrom = interval.chrom
        start = int(interval.start)
        end = int(interval.end)
        gene, txid, exon_idx = interval.name.split("|")
        strand = interval.strand

        # Queries the bigwig for a specific interval
        if start == end:
            end += 1
        coverages = args.bw.values(chrom, start, end)

        # Computes statistics at exon level
        exon_statistics = compute_exon_level_statistics(coverages)

        # Store results in data structure
        if txid not in results:
            results[txid] = defaultdict()
            results[txid]["chrom"] = chrom
            results[txid]["gene"] = gene
            results[txid]["strand"] = strand
            results[txid]["exons"] = defaultdict()
        results[gene][txid]["exons"][exon_idx] = {
            "statistics" : exon_statistics
        }

        # compute gaps
        if args.coverage_threshold > 0:
            gaps = find_gaps(coverages, start, args.coverage_threshold)
            results[txid]["exons"][exon_idx]["gaps"] = gaps

    # add an aggregation of statistics at transcript level
    for transcript, content in results.iteritems():
        results[transcript]["statistics"] = compute_transcript_level_statistics(content["exons"])

    # TODO: output results in different formats
    if args.output == "json":
        print json.dumps(results, indent = 4)
    elif args.output == "tabular":
        # TODO: flatten the dict into a table
        pass


if __name__ == '__main__':
    main()