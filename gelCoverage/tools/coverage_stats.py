import numpy as np
from collections import defaultdict


def find_gaps(coverages, start_position, coverage_threshold):
    """
    Find continuous genomic positions under a given threshold coverage_threshold.
    :param coverages: list of depth of coverage values
    :param start_position: starting position of the coverages sequence
    :param coverage_threshold: the coverage threshold to determine gaps
    :return: the gaps start and end genomic coordinates in JSON-friendly format. Chromosome is not set as this information
    will be embedded within an exon-transcript-gene where the chromosome is available.
    """
    end = start_position + len(coverages)
    open_gap = False
    current_gap = defaultdict()
    gaps = []

    # Iterates through every coverage position
    for idx, value in enumerate(coverages):
        if value < coverage_threshold and not open_gap:
            open_gap = True
            current_gap["start"] = start_position + idx
        elif value >= coverage_threshold and open_gap:
            open_gap = False
            current_gap["end"] = start_position + idx - 1
            gaps.append(current_gap)
            current_gap = defaultdict()
    # Closes the last gap when it extends until the last position
    if open_gap:
        current_gap["end"] = end
        gaps.append(current_gap)

    return gaps


def compute_exon_level_statistics(coverages, gc_content):
    """
    Computes coverage and GC content statistics
    :param coverages: list of depth of coverage values
    :param gc_content: the GC content for this sequence precomputed
    :return: the coverage and GC content exon statistics in JSON-friendly format
    """
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
    stats["percent_lt_15x"] = float(stats["bases_lt_15x"] / stats["total_bases"])
    stats["percent_gte_15x"] = float(stats["bases_gte_15x"] / stats["total_bases"])
    stats["percent_gte_30x"] = float(stats["bases_gte_30x"] / stats["total_bases"])
    stats["percent_gte_50x"] = float(stats["bases_gte_50x"] / stats["total_bases"])
    stats["gc_content"] = gc_content

    return stats


def compute_transcript_level_statistics(exons):
    """
    Computes coverage and GC content statistics at gene level by aggregating the statistics at exon level.
    Median and percentiles are estimated by weighting the per-exon metric by the number of bases.
    :param exons: list of exon coverage and GC content statistics
    :return: the coverage and GC content gene statistics in JSON-friendly format
    """
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
    stats["percent_lt_15x"] = float(stats["bases_lt_15x"] / stats["total_bases"])
    stats["percent_gte_15x"] = float(stats["bases_gte_15x"] / stats["total_bases"])
    stats["percent_gte_30x"] = float(stats["bases_gte_30x"] / stats["total_bases"])
    stats["percent_gte_50x"] = float(stats["bases_gte_50x"] / stats["total_bases"])
    stats["gc_content"] = np.sum(
        [x["gc_content"] * x["total_bases"] for _, x in exons.iteritems()]) / stats["total_bases"]

    return stats