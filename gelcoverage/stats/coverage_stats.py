import numpy as np
import logging
import itertools
import operator
from gelcoverage.tools.bed_reader import BedReader
from gelcoverage.tools.bigwig_reader import BigWigReader
import gelcoverage.constants as constants


def find_gaps(coverages, start_position, coverage_threshold, gap_length_threshold):
    """
    Find continuous genomic positions under a given threshold coverage_threshold.
    :param coverages: list of depth of coverage values
    :param start_position: starting position of the coverages sequence (0-based)
    :param coverage_threshold: the coverage threshold to determine gaps
    :param gap_length_threshold: the gap length threshold, every gap under this length will be discarded
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
            current_gap[constants.GAP_START] = start_position + idx
        elif value >= coverage_threshold and open_gap:
            open_gap = False
            current_gap[constants.GAP_END] = start_position + idx
            current_gap[constants.GAP_LENGTH] = current_gap[constants.GAP_END] - current_gap[constants.GAP_START]
            if current_gap[constants.GAP_LENGTH] >= gap_length_threshold:
                gaps.append(current_gap)
            current_gap = {}
    # Closes the last gap when it extends until the last position
    if open_gap:
        current_gap[constants.GAP_END] = end
        current_gap[constants.GAP_LENGTH] = current_gap[constants.GAP_END] - current_gap[constants.GAP_START] + 1
        if current_gap[constants.GAP_LENGTH] >= gap_length_threshold:
            gaps.append(current_gap)

    return gaps


def compute_exon_level_statistics(coverages, gc_content):
    """
    Computes coverage and GC content statistics
    :param coverages: list of depth of coverage values
    :param gc_content: the GC content for this sequence precomputed
    :return: the coverage and GC content exon statistics in JSON-friendly format
    """
    stats = {
        constants.BASES: len(coverages) if coverages else 0,
        constants.AVERAGE: round(float(np.mean(coverages)), 3) if coverages else 0.0,
        constants.MEDIAN: round(float(np.median(coverages)), 3) if coverages else 0.0,
        constants.PERCENTILE75: round(float(np.percentile(coverages, 75)), 3) if coverages else 0.0,
        constants.PERCENTILE25: round(float(np.percentile(coverages, 25)), 3) if coverages else 0.0,
        constants.SD: round(float(np.std(coverages)), 3) if coverages else 0.0,
        constants.BASES_LT15X: int(np.sum(1 for x in coverages if x < 15)) if coverages else 0,
        constants.BASES_GTE15X: int(np.sum(1 for x in coverages if x >= 15)) if coverages else 0,
        constants.BASES_GTE30X: int(np.sum(1 for x in coverages if x >= 30)) if coverages else 0,
        constants.BASES_GTE50X: int(np.sum(1 for x in coverages if x >= 50) if coverages else 0)
    }
    stats[constants.LT15X] = round(float(stats[constants.BASES_LT15X]) / stats[constants.BASES], 5) \
        if stats[constants.BASES] > 0 else 0.0
    stats[constants.GTE15X] = round(float(stats[constants.BASES_GTE15X]) / stats[constants.BASES], 5) \
        if stats[constants.BASES] > 0 else 0.0
    stats[constants.GTE30X] = round(float(stats[constants.BASES_GTE30X]) / stats[constants.BASES], 5) \
        if stats[constants.BASES] > 0 else 0.0
    stats[constants.GTE50X] = round(float(stats[constants.BASES_GTE50X]) / stats[constants.BASES], 5) \
        if stats[constants.BASES] > 0 else 0.0
    if gc_content is not None:  # GC content is not provided for padded exons
        stats[constants.GC_CONTENT] = gc_content

    return stats


def compute_transcript_level_statistics(exons):
    """
    Computes coverage and GC content statistics at gene level by aggregating the statistics at exon level.
    Median and percentiles are estimated by weighting the per-exon metric by the number of bases.
    :param exons: list of exon coverage and GC content statistics
    :return: the coverage and GC content gene statistics in JSON-friendly format
    """
    exons_stats = [x[constants.STATISTICS] for x in exons]
    total_bases = int(np.sum([x[constants.BASES] for x in exons_stats])) if exons_stats else 0
    bases_lt_15x = int(np.sum([x[constants.BASES_LT15X] for x in exons_stats])) if exons_stats else 0
    bases_gte_15x = int(np.sum([x[constants.BASES_GTE15X] for x in exons_stats])) if exons_stats else 0
    bases_gte_30x = int(np.sum([x[constants.BASES_GTE30X] for x in exons_stats])) if exons_stats else 0
    bases_gte_50x = int(np.sum([x[constants.BASES_GTE50X] for x in exons_stats])) if exons_stats else 0
    stats = {
        constants.BASES: total_bases,
        constants.AVERAGE: round(float(np.mean([x[constants.AVERAGE] for x in exons_stats])), 3)
        if exons_stats else 0.0,
        constants.MEDIAN: round(float(np.sum(
            [x[constants.MEDIAN] * x[constants.BASES] for x in exons_stats])) / total_bases, 3) if exons_stats else float(0.0),
        constants.PERCENTILE25: round(float(np.sum(
            [x[constants.PERCENTILE25] * x[constants.BASES] for x in exons_stats])) / total_bases, 3)
        if exons_stats else 0.0,
        constants.PERCENTILE75: round(float(np.sum(
            [x[constants.PERCENTILE75] * x[constants.BASES] for x in exons_stats])) / total_bases, 3)
        if exons_stats else 0.0,
        constants.SD: round(float(np.sum(
            [x[constants.SD] * x[constants.BASES] for x in exons_stats])) / total_bases, 3) if exons_stats else 0.0,
        constants.LT15X: round(float(bases_lt_15x) / total_bases, 5) if total_bases > 0 else 0.0,
        constants.GTE15X: round(float(bases_gte_15x) / total_bases, 5) if total_bases > 0 else 0.0,
        constants.GTE30X: round(float(bases_gte_30x) / total_bases, 5) if total_bases > 0 else 0.0,
        constants.GTE50X: round(float(bases_gte_50x) / total_bases, 5) if total_bases > 0 else 0.0,
        constants.BASES_LT15X: bases_lt_15x,
        constants.BASES_GTE15X: bases_gte_15x,
        constants.BASES_GTE30X: bases_gte_30x,
        constants.BASES_GTE50X: bases_gte_50x
    }
    try:
        stats[constants.GC_CONTENT] = round(float(np.sum(
            [x[constants.GC_CONTENT] * x[constants.BASES] for x in exons_stats]) / total_bases), 5) \
            if exons_stats and total_bases > 0 else 0.0
    except KeyError:
        # There is no GC content data to show (e.g.: the union transcript)
        pass
    return stats


def compute_coding_region_statistics(genes):
    """

    :param genes:
    :return:
    """
    logging.info("Computing coding region statistics...")
    results = {
        constants.STATISTICS: None,
        constants.CHROMOSOMES: []
    }
    # Avoids failing when no genes have been reported (might be related with wrong BAM and/or gene list)
    if len(genes) == 0:
        return results
    # Compute the stats aggregated for union transcript
    genes_stats = [x[constants.UNION_TRANSCRIPT][constants.STATISTICS] for x in genes]
    total_bases = int(np.sum([x[constants.BASES] for x in genes_stats])) if genes_stats else 0
    bases_lt_15x = int(np.sum([x[constants.BASES_LT15X] for x in genes_stats])) if genes_stats else 0
    bases_gte_15x = int(np.sum([x[constants.BASES_GTE15X] for x in genes_stats])) if genes_stats else 0
    bases_gte_30x = int(np.sum([x[constants.BASES_GTE30X] for x in genes_stats])) if genes_stats else 0
    bases_gte_50x = int(np.sum([x[constants.BASES_GTE50X] for x in genes_stats])) if genes_stats else 0
    results[constants.STATISTICS] = {
        constants.BASES: total_bases,
        constants.AVERAGE: round(float(np.mean([x[constants.AVERAGE] for x in genes_stats])), 3) if genes_stats else 0.0,
        constants.MEDIAN: round(float(np.sum(
            [x[constants.MEDIAN] * x[constants.BASES] for x in genes_stats]) / total_bases), 3)
        if genes_stats and total_bases > 0 else 0.0,
        constants.PERCENTILE75: round(float(np.sum(
            [x[constants.PERCENTILE75] * x[constants.BASES] for x in genes_stats]) / total_bases), 3)
        if genes_stats and total_bases > 0 else 0.0,
        constants.PERCENTILE25: round(float(np.sum(
            [x[constants.PERCENTILE25] * x[constants.BASES] for x in genes_stats]) / total_bases), 3)
        if genes_stats and total_bases > 0 else 0.0,
        constants.SD: round(float(np.sum(
            [x[constants.SD] * x[constants.BASES] for x in genes_stats]) / total_bases), 3)
        if genes_stats and total_bases > 0 else 0.0,
        constants.LT15X: round(float(bases_lt_15x) / total_bases, 5) if total_bases > 0 else 0.0,
        constants.GTE15X: round(float(bases_gte_15x) / total_bases, 5) if total_bases > 0 else 0.0,
        constants.GTE30X: round(float(bases_gte_30x) / total_bases, 5) if total_bases > 0 else 0.0,
        constants.GTE50X: round(float(bases_gte_50x) / total_bases, 5) if total_bases > 0 else 0.0
    }
    # Compute the stats disaggregated by chromosome
    chr2stats = [(x[constants.CHROMOSOME], x[constants.UNION_TRANSCRIPT][constants.STATISTICS]) for x in genes]

    def groupby_chromosome(list_of_tuples):
        it = itertools.groupby(list_of_tuples, operator.itemgetter(0))
        for _chromosome, subiter in it:
            yield _chromosome, [item[1] for item in subiter]

    # Aggregates stats for all chromosomes
    chromosome_stats = dict(groupby_chromosome(chr2stats))
    autosomes_stats = []
    for chromosome, chr_stats in chromosome_stats.iteritems():
        chr_total_bases = int(np.sum([x[constants.BASES] for x in chr_stats])) if chr_stats else 0
        chr_bases_lt_15x = int(np.sum([x[constants.BASES_LT15X] for x in chr_stats])) if chr_stats else 0
        chr_bases_gte_15x = int(np.sum([x[constants.BASES_GTE15X] for x in chr_stats])) if chr_stats else 0
        chr_bases_gte_30x = int(np.sum([x[constants.BASES_GTE30X] for x in chr_stats])) if chr_stats else 0
        chr_bases_gte_50x = int(np.sum([x[constants.BASES_GTE50X] for x in chr_stats])) if chr_stats else 0
        formatted_chr_stats = {
            constants.CHROMOSOME: chromosome,
            constants.BASES: chr_total_bases,
            constants.AVERAGE: round(float(np.mean([x[constants.AVERAGE] for x in chr_stats])), 3)
            if chr_stats else 0.0,
            constants.MEDIAN: round(float(np.sum(
                [x[constants.MEDIAN] * x[constants.BASES] for x in chr_stats]) / chr_total_bases), 3)
            if chr_stats and chr_total_bases > 0 else 0.0,
            constants.PERCENTILE75: round(float(np.sum(
                [x[constants.PERCENTILE75] * x[constants.BASES] for x in chr_stats]) / chr_total_bases), 3)
            if chr_stats and chr_total_bases > 0 else 0.0,
            constants.PERCENTILE25: round(float(np.sum(
                [x[constants.PERCENTILE25] * x[constants.BASES] for x in chr_stats]) / chr_total_bases), 3)
            if chr_stats and chr_total_bases > 0 else 0.0,
            constants.SD: round(float(np.sum(
                [x[constants.SD] * x[constants.BASES] for x in chr_stats]) / chr_total_bases), 3)
            if chr_stats and chr_total_bases > 0 else 0.0,
            constants.LT15X: round(float(chr_bases_lt_15x) / chr_total_bases, 5) if chr_total_bases > 0 else 0.0,
            constants.GTE15X: round(float(chr_bases_gte_15x) / chr_total_bases, 5) if chr_total_bases > 0 else 0.0,
            constants.GTE30X: round(float(chr_bases_gte_30x) / chr_total_bases, 5) if chr_total_bases > 0 else 0.0,
            constants.GTE50X: round(float(chr_bases_gte_50x) / chr_total_bases, 5) if chr_total_bases > 0 else 0.0
        }
        results[constants.CHROMOSOMES].append(formatted_chr_stats)
        logging.info("Coding region statistics for chromosome %s computed!" % chromosome)
        # Records stats for autosome
        if chromosome in constants.AUTOSOME_IDS:
            autosomes_stats.append(formatted_chr_stats)
    # Aggregates stats for autosomes
    autosomes_total_bases = int(np.sum([x[constants.BASES] for x in autosomes_stats])) if autosomes_stats else 0
    autosomes_chr_stats = {
        constants.CHROMOSOME: constants.AUTOSOMES,
        constants.BASES: autosomes_total_bases,
        constants.AVERAGE: round(float(np.mean([x[constants.AVERAGE] for x in autosomes_stats])), 3)
        if autosomes_stats else 0.0,
        constants.MEDIAN: round(float(np.sum(
            [x[constants.MEDIAN] * x[constants.BASES] for x in autosomes_stats]) / autosomes_total_bases), 3)
        if autosomes_stats and autosomes_total_bases > 0 else 0.0,
        constants.PERCENTILE75: round(float(np.sum(
            [x[constants.PERCENTILE75] * x[constants.BASES] for x in autosomes_stats]) / autosomes_total_bases), 3)
        if autosomes_stats and autosomes_total_bases > 0 else 0.0,
        constants.PERCENTILE25: round(float(np.sum(
            [x[constants.PERCENTILE25] * x[constants.BASES] for x in autosomes_stats]) / autosomes_total_bases), 3)
        if autosomes_stats and autosomes_total_bases > 0 else 0.0,
        constants.SD: round(float(np.sum(
            [x[constants.SD] * x[constants.BASES] for x in autosomes_stats]) / autosomes_total_bases), 3)
        if autosomes_stats and autosomes_total_bases > 0 else 0.0,
        constants.LT15X: round(float(np.sum(
            [x[constants.LT15X] * x[constants.BASES] for x in autosomes_stats]) / autosomes_total_bases), 5)
        if autosomes_stats and autosomes_total_bases > 0 else 0.0,
        constants.GTE15X: round(float(np.sum(
            [x[constants.GTE15X] * x[constants.BASES] for x in autosomes_stats]) / autosomes_total_bases), 5)
        if autosomes_stats and autosomes_total_bases > 0 else 0.0,
        constants.GTE30X: round(float(np.sum(
            [x[constants.GTE30X] * x[constants.BASES] for x in autosomes_stats]) / autosomes_total_bases), 5)
        if autosomes_stats and autosomes_total_bases > 0 else 0.0,
        constants.GTE50X: round(float(np.sum(
            [x[constants.GTE50X] * x[constants.BASES] for x in autosomes_stats]) / autosomes_total_bases), 5)
        if autosomes_stats and autosomes_total_bases > 0 else 0.0
    }
    results[constants.CHROMOSOMES].append(autosomes_chr_stats)

    logging.info("Coding region statistics computed!")
    return results


def compute_whole_genome_statistics(bigwig_reader, bed_reader=None, chunk_size=100000):
    """
    Iterates through the whole genome in a sliding window to obtain some metrics
    :type chunk_size: int
    :type bigwig_reader: BigWigReader
    :type bed_reader: BedReader
    :param bigwig_reader:
    :return:
    """
    logging.info("Computing whole genome statistics...")
    results = {
        constants.STATISTICS: None,
        constants.CHROMOSOMES: []
    }
    if bed_reader.is_null_bed:
        logging.info("Running on all chromosomes defined in the bigwig.")
        analysis_regions = bigwig_reader.get_chromosome_lengths()
    else:
        logging.info("Running on the regions provided in a bed file in --wg-region.")
        analysis_regions = bed_reader.get_regions_dictionary()
    # Iterates each chromosome
    chr_stats = {}
    autosomes_stats = []
    for chromosome, regions in analysis_regions.iteritems():
        chr_stats[chromosome] = {
            constants.RMSD: [],
            constants.AVERAGE: [],
            constants.BASES: [],
            constants.BASES_LT15X: [],
            constants.BASES_GTE15X: [],
            constants.BASES_GTE30X: [],
            constants.BASES_GTE50X: [],
            constants.MEDIAN: [],
            constants.PERCENTILE25: [],
            constants.PERCENTILE75: [],
            constants.SD: []
        }
        # Iterates intervals in chunks of fixed size and stores the stats for each chunk
        for (start, end) in regions:
            logging.debug("Analysing region %s:%s-%s" % (chromosome, start, end))
            current_start = start
            current_end = min(current_start + chunk_size, end)
            while current_start < current_end:
                logging.debug("Analysing chunk %s:%s-%s" % (chromosome,
                                                            current_start,
                                                            current_end))
                coverages = bigwig_reader.read_bigwig_coverages(chromosome, current_start, current_end, strict=False)
                length = current_end - current_start
                chunk_mean = np.mean(coverages)
                if chunk_mean == 0:
                    # As coverage values are positive values we infer that all values are zero
                    # this may speed up things for missing long regions in the bigwig file, if any
                    chunk_rmsd = 0
                else:
                    # Gets the squared root sum of squares of the deviation from the mean
                    chunk_rmsd = np.sqrt((np.sum([(x - chunk_mean) ** 2 for x in coverages]) if coverages else 0)
                                         / length)
                chr_stats[chromosome][constants.BASES].append(length)
                chr_stats[chromosome][constants.BASES_LT15X].append(
                    np.sum(1 for coverage in coverages if coverage < 15) if coverages else 0
                )
                chr_stats[chromosome][constants.BASES_GTE15X].append(
                    np.sum(1 for coverage in coverages if coverage >= 15) if coverages else 0
                )
                chr_stats[chromosome][constants.BASES_GTE30X].append(
                    np.sum(1 for coverage in coverages if coverage >= 30) if coverages else 0
                )
                chr_stats[chromosome][constants.BASES_GTE50X].append(
                    np.sum(1 for coverage in coverages if coverage >= 50) if coverages else 0
                )
                chr_stats[chromosome][constants.AVERAGE].append(chunk_mean)
                chr_stats[chromosome][constants.RMSD].append(chunk_rmsd)
                chr_stats[chromosome][constants.MEDIAN].append(np.median(coverages))
                chr_stats[chromosome][constants.PERCENTILE25].append(np.percentile(coverages, 25) if coverages else 0.0)
                chr_stats[chromosome][constants.PERCENTILE75].append(np.percentile(coverages, 75) if coverages else 0.0)
                chr_stats[chromosome][constants.SD].append(np.std(coverages) if coverages else 0)
                current_start = current_end
                current_end = min(current_start + chunk_size, end)
        # Set the statistics per chromosome
        chr_total_bases = np.sum(chr_stats[chromosome][constants.BASES])
        chr_bases_lt_15x = np.sum(chr_stats[chromosome][constants.BASES_LT15X])
        chr_bases_gte_15x = np.sum(chr_stats[chromosome][constants.BASES_GTE15X])
        chr_bases_gte_30x = np.sum(chr_stats[chromosome][constants.BASES_GTE30X])
        chr_bases_gte_50x = np.sum(chr_stats[chromosome][constants.BASES_GTE50X])
        formatted_chr_stats = {
            constants.CHROMOSOME: chromosome,
            constants.RMSD: round(float(np.median(chr_stats[chromosome][constants.RMSD])), 3)
            if chr_stats[chromosome][constants.RMSD] else 0.0,
            constants.AVERAGE: round(float(np.mean(chr_stats[chromosome][constants.AVERAGE])), 3)
            if chr_stats[chromosome][constants.AVERAGE] else 0.0,
            constants.BASES: int(chr_total_bases),
            constants.LT15X: round(float(chr_bases_lt_15x) / chr_total_bases, 5) if chr_total_bases > 0 else 0.0,
            constants.GTE15X: round(float(chr_bases_gte_15x) / chr_total_bases, 5) if chr_total_bases > 0 else 0.0,
            constants.GTE30X: round(float(chr_bases_gte_30x) / chr_total_bases, 5) if chr_total_bases > 0 else 0.0,
            constants.GTE50X: round(float(chr_bases_gte_50x) / chr_total_bases, 5) if chr_total_bases > 0 else 0.0,
            constants.MEDIAN: round(float(np.sum([median * weight for median, weight in
                                                  zip(chr_stats[chromosome][constants.MEDIAN],
                                                      chr_stats[chromosome][constants.BASES])])) / chr_total_bases, 3)
            if chr_total_bases > 0 and chr_stats[chromosome][constants.MEDIAN] and
               chr_stats[chromosome][constants.BASES] else 0.0,
            constants.PERCENTILE75: round(float(np.sum([pct75 * weight for pct75, weight in
                                                        zip(chr_stats[chromosome][constants.PERCENTILE75],
                                                            chr_stats[chromosome][constants.BASES])])) /
                                          chr_total_bases, 3)
            if chr_total_bases > 0 and chr_stats[chromosome][constants.PERCENTILE75] and
               chr_stats[chromosome][constants.BASES] else 0.0,
            constants.PERCENTILE25: round(float(np.sum([pct25 * weight
                                         for pct25, weight in zip(chr_stats[chromosome][constants.PERCENTILE25],
                                             chr_stats[chromosome][constants.BASES])])) /
                           chr_total_bases, 3)
            if chr_total_bases > 0 and chr_stats[chromosome][constants.PERCENTILE25] and
               chr_stats[chromosome][constants.BASES] else 0.0,
            constants.SD: round(float(np.sum([sd * weight for sd, weight in zip(chr_stats[chromosome][constants.SD],
                                                                            chr_stats[chromosome][constants.BASES])])) /
                                          chr_total_bases, 3)
            if chr_total_bases > 0 and chr_stats[chromosome][constants.SD] and
               chr_stats[chromosome][constants.BASES] else 0.0,
        }
        results[constants.CHROMOSOMES].append(formatted_chr_stats)
        if chromosome in constants.AUTOSOME_IDS:
            autosomes_stats.append(formatted_chr_stats)
        logging.info("Whole genome statistics for chromosome %s computed!" % chromosome)

    # Aggregates statistics for the whole genome (important to do before addings autosomes!)
    formatted_chr_stats = aggregate_chromosomes(results[constants.CHROMOSOMES])
    results[constants.STATISTICS] = formatted_chr_stats
    logging.info("Aggregated whole genome statistics!")

    # Aggregates statistics for the autosomes
    formatted_autosomes_stats = aggregate_chromosomes(autosomes_stats)
    formatted_autosomes_stats[constants.CHROMOSOME] = constants.AUTOSOMES
    results[constants.CHROMOSOMES].append(formatted_autosomes_stats)
    logging.info("Aggregated whole genome statistics for autosomes!")

    logging.info("Whole genome statistics computed!")
    return results


def aggregate_chromosomes(stats):
    """
    Aggregates information from several chromosomes
    :param stats: the list of chromosome stats
    :type stats: list
    :return: the aggregated stats
    """
    total_bases = int(np.sum([x[constants.BASES] for x in stats])) if stats else 0
    return {
        constants.BASES: total_bases,
        constants.AVERAGE: round(float(np.mean([x[constants.AVERAGE] for x in stats])), 3) if stats else 0.0,
        constants.MEDIAN: round(float(np.sum(
            [x[constants.MEDIAN] * x[constants.BASES] for x in stats]) / total_bases), 3)
        if total_bases > 0  and stats else 0.0,
        constants.PERCENTILE75: round(float(np.sum(
            [x[constants.PERCENTILE75] * x[constants.BASES] for x in stats]) / total_bases), 3)
        if total_bases > 0 and stats else 0.0,
        constants.PERCENTILE25: round(float(np.sum(
            [x[constants.PERCENTILE25] * x[constants.BASES] for x in stats]) / total_bases), 3)
        if total_bases > 0 and stats else 0.0,
        constants.SD: round(float(np.sum(
            [x[constants.SD] * x[constants.BASES] for x in stats]) / total_bases), 3)
        if total_bases > 0 and stats else 0.0,
        constants.LT15X: round(float(np.sum(
            [x[constants.LT15X] * x[constants.BASES] for x in stats]) / total_bases), 5)
        if total_bases > 0 and stats else 0.0,
        constants.GTE15X: round(float(np.sum(
            [x[constants.GTE15X] * x[constants.BASES] for x in stats]) / total_bases), 5)
        if total_bases > 0 and stats else 0.0,
        constants.GTE30X: round(float(np.sum(
            [x[constants.GTE30X] * x[constants.BASES] for x in stats]) / total_bases), 5)
        if total_bases > 0 and stats else 0.0,
        constants.GTE50X: round(float(np.sum(
            [x[constants.GTE50X] * x[constants.BASES] for x in stats]) / total_bases), 5)
        if total_bases > 0 and stats else 0.0,
        constants.RMSD: round(float(np.sum(
            [x[constants.RMSD] * x[constants.BASES] for x in stats]) / total_bases), 3)
        if total_bases > 0 and stats else 0.0
    }
