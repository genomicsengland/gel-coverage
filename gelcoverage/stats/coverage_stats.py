import numpy as np
import logging
import gelcoverage.tools.bed_parser as bed_parser


def find_gaps(coverages, start_position, coverage_threshold):
    """
    Find continuous genomic positions under a given threshold coverage_threshold.
    :param coverages: list of depth of coverage values
    :param start_position: starting position of the coverages sequence
    :param coverage_threshold: the coverage threshold to determine gaps
    :return: the gaps start and end genomic coordinates in JSON-friendly format.
    Chromosome is not set as this information
    will be embedded within an exon-transcript-gene where the chromosome is available.
    """
    end = start_position + len(coverages)
    open_gap = False
    current_gap = {}
    gaps = []

    # Iterates through every coverage position
    for idx, value in enumerate(coverages):
        if value < coverage_threshold and not open_gap:
            open_gap = True
            current_gap["start"] = start_position + idx
        elif value >= coverage_threshold and open_gap:
            open_gap = False
            current_gap["end"] = start_position + idx - 1
            current_gap["length"] = current_gap["end"] - current_gap["start"] + 1
            gaps.append(current_gap)
            current_gap = {}
    # Closes the last gap when it extends until the last position
    if open_gap:
        current_gap["end"] = end
        current_gap["length"] = current_gap["end"] - current_gap["start"] + 1
        gaps.append(current_gap)

    return gaps


def compute_exon_level_statistics(coverages, gc_content):
    """
    Computes coverage and GC content statistics
    :param coverages: list of depth of coverage values
    :param gc_content: the GC content for this sequence precomputed
    :return: the coverage and GC content exon statistics in JSON-friendly format
    """
    stats = {}
    stats["total_bases"] = len(coverages)
    stats["mean"] = round(float(np.mean(coverages)), 3)
    stats["median"] = round(float(np.median(coverages)), 3)
    stats["pct75"] = round(float(np.percentile(coverages, 75)), 3)
    stats["pct25"] = round(float(np.percentile(coverages, 25)), 3)
    stats["bases_lt_3x"] = int(np.sum(1 for x in coverages if x < 3))
    stats["bases_lt_15x"] = int(np.sum(1 for x in coverages if x < 15))
    stats["bases_gte_15x"] = int(np.sum(1 for x in coverages if x >= 15))
    stats["bases_gte_30x"] = int(np.sum(1 for x in coverages if x >= 30))
    stats["bases_gte_50x"] = int(np.sum(1 for x in coverages if x >= 50))
    stats["percent_lt_15x"] = round(float(stats["bases_lt_15x"]) / stats["total_bases"], 5)
    stats["percent_gte_15x"] = round(float(stats["bases_gte_15x"]) / stats["total_bases"], 5)
    stats["percent_gte_30x"] = round(float(stats["bases_gte_30x"]) / stats["total_bases"], 5)
    stats["percent_gte_50x"] = round(float(stats["bases_gte_50x"]) / stats["total_bases"], 5)
    if gc_content is not None:  # GC content is not provided for padded exons
        stats["gc_content"] = gc_content

    return stats

def compute_transcript_level_statistics(exons):
    """
    Computes coverage and GC content statistics at gene level by aggregating the statistics at exon level.
    Median and percentiles are estimated by weighting the per-exon metric by the number of bases.
    :param exons: list of exon coverage and GC content statistics
    :return: the coverage and GC content gene statistics in JSON-friendly format
    """
    stats = {}
    exons_stats = [x["statistics"] for x in exons]
    stats["total_bases"] = int(np.sum([x["total_bases"] for x in exons_stats]))
    stats["mean"] = round(float(np.mean([x["mean"] for x in exons_stats])), 3)
    stats["weighted_median"] = round(float(np.sum(
        [x["median"] * x["total_bases"] for x in exons_stats]) / stats["total_bases"]), 3)
    stats["weighted_pct75"] = round(float(np.sum(
        [x["pct75"] * x["total_bases"] for x in exons_stats]) / stats["total_bases"]), 3)
    stats["weighted_pct25"] = round(float(np.sum(
        [x["pct25"] * x["total_bases"] for x in exons_stats]) / stats["total_bases"]), 3)
    stats["bases_lt_3x"] = int(np.sum([x["bases_lt_3x"] for x in exons_stats]))
    stats["bases_lt_15x"] = int(np.sum([x["bases_lt_15x"] for x in exons_stats]))
    stats["bases_gte_15x"] = int(np.sum([x["bases_gte_15x"] for x in exons_stats]))
    stats["bases_gte_30x"] = int(np.sum([x["bases_gte_30x"] for x in exons_stats]))
    stats["bases_gte_50x"] = int(np.sum([x["bases_gte_50x"] for x in exons_stats]))
    stats["percent_lt_15x"] = round(float(stats["bases_lt_15x"]) / stats["total_bases"], 5)
    stats["percent_gte_15x"] = round(float(stats["bases_gte_15x"]) / stats["total_bases"], 5)
    stats["percent_gte_30x"] = round(float(stats["bases_gte_30x"]) / stats["total_bases"], 5)
    stats["percent_gte_50x"] = round(float(stats["bases_gte_50x"]) / stats["total_bases"], 5)
    try:
        stats["gc_content"] = round(float(np.sum(
            [x["gc_content"] * x["total_bases"] for x in exons_stats]) / stats["total_bases"]), 5)
    except KeyError:
        # There is no GC content data to show (e.g.: the union transcript)
        pass
    return stats

def compute_panel_level_statistics(genes):
    """

    :param genes:
    :return:
    """
    stats = {}
    # TODO: compute the stats aggregated for union transcript
    genes_stats = [x["union_transcript"]["statistics"] for x in genes]
    stats["total_bases"] = int(np.sum([x["total_bases"] for x in genes_stats]))
    stats["mean"] = round(float(np.mean([x["mean"] for x in genes_stats])), 3)
    stats["weighted_median"] = round(float(np.sum(
        [x["weighted_median"] * x["total_bases"] for x in genes_stats]) / stats["total_bases"]), 3)
    stats["weighted_pct75"] = round(float(np.sum(
        [x["weighted_pct75"] * x["total_bases"] for x in genes_stats]) / stats["total_bases"]), 3)
    stats["weighted_pct25"] = round(float(np.sum(
        [x["weighted_pct25"] * x["total_bases"] for x in genes_stats]) / stats["total_bases"]), 3)
    stats["bases_lt_3x"] = int(np.sum([x["bases_lt_3x"] for x in genes_stats]))
    stats["bases_lt_15x"] = int(np.sum([x["bases_lt_15x"] for x in genes_stats]))
    stats["bases_gte_15x"] = int(np.sum([x["bases_gte_15x"] for x in genes_stats]))
    stats["bases_gte_30x"] = int(np.sum([x["bases_gte_30x"] for x in genes_stats]))
    stats["bases_gte_50x"] = int(np.sum([x["bases_gte_50x"] for x in genes_stats]))
    stats["percent_lt_15x"] = round(float(stats["bases_lt_15x"]) / stats["total_bases"], 5)
    stats["percent_gte_15x"] = round(float(stats["bases_gte_15x"]) / stats["total_bases"], 5)
    stats["percent_gte_30x"] = round(float(stats["bases_gte_30x"]) / stats["total_bases"], 5)
    stats["percent_gte_50x"] = round(float(stats["bases_gte_50x"]) / stats["total_bases"], 5)
    return stats

def compute_whole_genome_statistics(bigwig_reader, bed = None):
    """
    Iterates through the whole genome in a sliding window to obtain some metrics
    :param bigwig_reader:
    :return:
    """
    logging.info("Computing whole genome statistics...")
    stats = {}
    chunk_size = 100000
    if bed is None:
        analysis_regions = bigwig_reader.get_chromosome_lengths()
    else:
        analysis_regions = bed_parser.get_regions_dictionary(bed)
    rmsds = []
    means = []
    total_bases = []
    bases_lt_3x = []
    bases_lt_15x = []
    bases_gte_15x = []
    bases_gte_30x = []
    bases_gte_50x = []
    medians = []
    pct25s = []
    pct75s = []
    lengths = []
    for chromosome, regions in analysis_regions.iteritems():
        for (start, end) in regions:
            current_start = start
            current_end = min(current_start + chunk_size - 1, end)
            while current_end < end:
                coverages = bigwig_reader.read_bigwig_coverages(chromosome, current_start, current_end, strict=False)
                length = len(coverages)
                chunk_mean = np.mean(coverages)
                if chunk_mean == 0:
                    # As coverage values are positive values we infer that all values are zero
                    # this may speed up things for missing long regions in the bigwig file, if any
                    chunk_rmsd = 0
                else:
                    # Gets the squared root sum of squares of the deviation from the mean
                    chunk_rmsd = np.sqrt(np.sum([(x - chunk_mean) ** 2 for x in coverages]) / length)
                total_bases.append(np.sum([1 for coverage in coverages]))
                bases_lt_3x.append(np.sum([1 for coverage in coverages if coverage < 3]))
                bases_lt_15x.append(np.sum([1 for coverage in coverages if coverage < 15]))
                bases_gte_15x.append(np.sum([1 for coverage in coverages if coverage >= 15]))
                bases_gte_30x.append(np.sum([1 for coverage in coverages if coverage >= 30]))
                bases_gte_50x.append(np.sum([1 for coverage in coverages if coverage >= 50]))
                means.append(chunk_mean)
                rmsds.append(chunk_rmsd)
                medians.append(np.median(coverages))
                pct25s.append(np.percentile(coverages, 25))
                pct75s.append(np.percentile(coverages, 75))
                lengths.append(length)
                current_start = current_end + 1
                current_end = min(current_start + chunk_size - 1, end)
        stats["uneveness"] = round(float(np.median(rmsds)), 3)
        stats["mean"] = round(float(np.mean(means)), 3)
        stats["total_bases"] = int(np.sum(total_bases))
        stats["bases_lt_3x"] = int(np.sum(bases_lt_3x))
        stats["bases_lt_15x"] = int(np.sum(bases_lt_15x))
        stats["bases_gte_15x"] = int(np.sum(bases_gte_15x))
        stats["bases_gte_30x"] = int(np.sum(bases_gte_30x))
        stats["bases_gte_50x"] = int(np.sum(bases_gte_50x))
        stats["percent_lt_15x"] = round(float(stats["bases_lt_15x"]) / stats["total_bases"], 5)
        stats["percent_gte_15x"] = round(float(stats["bases_gte_15x"]) / stats["total_bases"], 5)
        stats["percent_gte_30x"] = round(float(stats["bases_gte_30x"]) / stats["total_bases"], 5)
        stats["percent_gte_50x"] = round(float(stats["bases_gte_50x"]) / stats["total_bases"], 5)
        stats["weighted_median"] = round(float(np.sum([median * length for median in medians for length in lengths]) /
                                               np.sum(lengths)), 3)
        stats["weighted_pct75"] = round(float(np.sum([pct75 * length for pct75 in pct75s for length in lengths]) /
                                               np.sum(lengths)), 3)
        stats["weighted_pct25"] = round(float(np.sum([pct25 * length for pct25 in pct25s for length in lengths]) /
                                               np.sum(lengths)), 3)
        logging.info("Whole genome statistics computed!")
        return stats


