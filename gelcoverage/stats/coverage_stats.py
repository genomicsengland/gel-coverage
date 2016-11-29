import numpy as np


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
    stats["gc_content"] = round(float(np.sum(
        [x["gc_content"] * x["total_bases"] for x in exons_stats]) / stats["total_bases"]), 5)
    return stats

def compute_panel_level_statistics(genes):
    """

    :param genes:
    :return:
    """
    stats = {}
    # TODO: compute the stats aggregated for union transcript
    genes_stats = [y["statistics"] for x in genes for y in x["transcripts"]]
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
    stats["gc_content"] = round(float(np.sum(
        [x["gc_content"] * x["total_bases"] for x in genes_stats]) / stats["total_bases"]), 5)
    return stats
