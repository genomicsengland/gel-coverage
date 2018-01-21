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
            length = current_gap[constants.GAP_END] - current_gap[constants.GAP_START]
            if length >= gap_length_threshold:
                gaps.append(current_gap)
            current_gap = {}
    # Closes the last gap when it extends until the last position
    if open_gap:
        current_gap[constants.GAP_END] = end
        length = current_gap[constants.GAP_END] - current_gap[constants.GAP_START] + 1
        if length >= gap_length_threshold:
            gaps.append(current_gap)

    return gaps


def compute_exon_level_statistics(coverages, gc_content):
    """
    Computes coverage and GC content statistics
    :param coverages: list of depth of coverage values
    :param gc_content: the GC content for this sequence precomputed
    :return: the coverage and GC content exon statistics in JSON-friendly format
    """
    bases = 0
    bases_lt15 = 0
    bases_gte15_lt30 = 0
    bases_gte30_lt50 = 0
    bases_gte50 = 0
    accum_coverage = 0
    for coverage in coverages:
        bases += 1
        accum_coverage += coverage
        # NOTE: this is optmised for performance, not for code simplicity or beauty
        if coverage >= 50:
            bases_gte50 += 1
        elif coverage >= 30:
            bases_gte30_lt50 += 1
        elif coverage >= 15:
            bases_gte15_lt30 += 1
        else:
            bases_lt15 += 1

    bases_gte15 = bases_gte15_lt30 + bases_gte30_lt50 + bases_gte50
    bases_gte30 = bases_gte30_lt50 + bases_gte50
    if bases > 0:
        mean_coverage = round(float(accum_coverage) / bases, 3)
        lt_15x = round(float(bases_lt15) / bases, 5)
        gte_15x = round(float(bases_gte15) / bases, 5)
        gte_30x = round(float(bases_gte30) / bases, 5)
        gte_50x = round(float(bases_gte50) / bases, 5)
    else:
        mean_coverage = 0.0
        lt_15x = 0.0
        gte_15x = 0.0
        gte_30x = 0.0
        gte_50x = 0.0
    if getattr(coverages, "size", 0) > 0:   # gets the size of the numpy array, if None returns 0
        median_coverage = round(float(np.median(coverages)), 3)
        pct75 = round(float(np.percentile(coverages, 75)), 3)
        pct25 = round(float(np.percentile(coverages, 25)), 3)
        sd = round(float(np.std(coverages)), 3)
    else:
        median_coverage = 0.0
        pct75 = 0.0
        pct25 = 0.0
        sd = 0.0

    stats = {
        constants.BASES: bases,
        constants.AVERAGE: mean_coverage,
        constants.MEDIAN: median_coverage,
        constants.PERCENTILE75: pct75,
        constants.PERCENTILE25: pct25,
        constants.SD: sd,
        constants.BASES_LT15X: bases_lt15,
        constants.BASES_GTE15X: bases_gte15,
        constants.BASES_GTE30X: bases_gte30,
        constants.BASES_GTE50X: bases_gte50,
        constants.LT15X: lt_15x,
        constants.GTE15X: gte_15x,
        constants.GTE30X: gte_30x,
        constants.GTE50X: gte_50x
    }

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
    exon_count = 0
    bases = 0
    bases_lt15 = 0
    bases_gte15 = 0
    bases_gte30 = 0
    bases_gte50 = 0
    accum_mean = 0.0
    accum_weighted_median = 0.0
    accum_weighted_pct75 = 0.0
    accum_weighted_pct25 = 0.0
    accum_weighted_sd = 0.0
    accum_weighted_gc_content = 0.0

    # NOTE: these variables are created to avoid a global lookup inside the loop
    constant_statistics = constants.STATISTICS
    constant_bases = constants.BASES
    constant_bases_lt15x = constants.BASES_LT15X
    constant_bases_gte15x = constants.BASES_GTE15X
    constant_bases_gte30x = constants.BASES_GTE30X
    constant_bases_gte50x = constants.BASES_GTE50X
    constant_average = constants.AVERAGE
    constant_median = constants.MEDIAN
    constant_pct75 = constants.PERCENTILE75
    constant_pct25 = constants.PERCENTILE25
    constant_sd = constants.SD
    constant_gc_content = constants.GC_CONTENT

    for exon in exons:
        exon_stats = exon[constant_statistics]
        if exon_stats:
            exon_count += 1
            bases += exon_stats[constant_bases]
            bases_lt15 += exon_stats[constant_bases_lt15x]
            bases_gte15 += exon_stats[constant_bases_gte15x]
            bases_gte30 += exon_stats[constant_bases_gte30x]
            bases_gte50 += exon_stats[constant_bases_gte50x]
            accum_mean += exon_stats[constant_average]
            accum_weighted_median += exon_stats[constant_median] * exon_stats[constant_bases]
            accum_weighted_pct75 += exon_stats[constant_pct75] * exon_stats[constant_bases]
            accum_weighted_pct25 += exon_stats[constant_pct25] * exon_stats[constant_bases]
            accum_weighted_sd += exon_stats[constant_sd] * exon_stats[constant_bases]
            if constants.GC_CONTENT in exon_stats:
                accum_weighted_gc_content += exon_stats[constant_gc_content] * exon_stats[constant_bases]

    mean = round(float(accum_mean) / exon_count, 3) if exon_count > 0 else 0.0  # the mean of the means
    if bases > 0:
        lt_15x = round(float(bases_lt15) / bases, 5)
        gte_15x = round(float(bases_gte15) / bases, 5)
        gte_30x = round(float(bases_gte30) / bases, 5)
        gte_50x = round(float(bases_gte50) / bases, 5)
        median = round(float(accum_weighted_median) / bases, 3)
        pct75 = round(float(accum_weighted_pct75) / bases, 3)
        pct25 = round(float(accum_weighted_pct25) / bases, 3)
        sd = round(float(accum_weighted_sd) / bases, 3)
        gc_content = round(float(accum_weighted_gc_content) / bases, 5)
    else:
        lt_15x = 0.0
        gte_15x = 0.0
        gte_30x = 0.0
        gte_50x = 0.0
        median = 0.0
        pct75 = 0.0
        pct25 = 0.0
        sd = 0.0
        gc_content = 0.0

    stats = {
        constant_bases: bases,
        constant_average: mean,
        constant_median: median,
        constant_pct25: pct25,
        constant_pct75: pct75,
        constant_sd: sd,
        constants.LT15X: lt_15x,
        constants.GTE15X: gte_15x,
        constants.GTE30X: gte_30x,
        constants.GTE50X: gte_50x,
        constant_bases_lt15x: bases_lt15,
        constant_bases_gte15x: bases_gte15,
        constant_bases_gte30x: bases_gte30,
        constant_bases_gte50x: bases_gte50
    }
    if gc_content > 0.0:
        stats[constants.GC_CONTENT] = gc_content

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

    gene_count = 0
    bases = 0
    bases_lt15 = 0
    bases_gte15 = 0
    bases_gte30 = 0
    bases_gte50 = 0
    accum_mean = 0.0
    accum_weighted_median = 0.0
    accum_weighted_pct75 = 0.0
    accum_weighted_pct25 = 0.0
    accum_weighted_sd = 0.0
    accum_weighted_gc_content = 0.0
    gene_count_by_chromosome = {}
    bases_by_chromosome = {}
    bases_lt15_by_chromosome = {}
    bases_gte15_by_chromosome = {}
    bases_gte30_by_chromosome = {}
    bases_gte50_by_chromosome = {}
    accum_mean_by_chromosome = {}
    accum_weighted_median_by_chromosome = {}
    accum_weighted_pct75_by_chromosome = {}
    accum_weighted_pct25_by_chromosome = {}
    accum_weighted_sd_by_chromosome = {}
    accum_weighted_gc_content_by_chromosome = {}
    observed_chromosomes = []
    gene_count_autosomes = 0
    bases_autosomes = 0
    bases_lt15_autosomes = 0
    bases_gte15_autosomes = 0
    bases_gte30_autosomes = 0
    bases_gte50_autosomes = 0
    accum_mean_autosomes = 0.0
    accum_weighted_median_autosomes = 0.0
    accum_weighted_pct75_autosomes = 0.0
    accum_weighted_pct25_autosomes = 0.0
    accum_weighted_sd_autosomes = 0.0
    accum_weighted_gc_content_autosomes = 0.0

    # NOTE: these variables are created to avoid a global lookup inside the loop
    constant_union_transcript = constants.UNION_TRANSCRIPT
    constant_statistics = constants.STATISTICS
    constant_bases = constants.BASES
    constant_bases_lt15x = constants.BASES_LT15X
    constant_bases_gte15x = constants.BASES_GTE15X
    constant_bases_gte30x = constants.BASES_GTE30X
    constant_bases_gte50x = constants.BASES_GTE50X
    constant_average = constants.AVERAGE
    constant_median = constants.MEDIAN
    constant_pct75 = constants.PERCENTILE75
    constant_pct25 = constants.PERCENTILE25
    constant_sd = constants.SD
    constant_gc_content = constants.GC_CONTENT
    constant_chromosome = constants.CHROMOSOME

    for gene in genes:
        gene_stats = gene[constant_union_transcript][constant_statistics]
        if gene_stats:

            # aggregates data for all genes
            gene_count += 1
            bases += gene_stats[constant_bases]
            bases_lt15 += gene_stats[constant_bases_lt15x]
            bases_gte15 += gene_stats[constant_bases_gte15x]
            bases_gte30 += gene_stats[constant_bases_gte30x]
            bases_gte50 += gene_stats[constant_bases_gte50x]
            accum_mean += gene_stats[constant_average]
            accum_weighted_median += gene_stats[constant_median] * gene_stats[constant_bases]
            accum_weighted_pct75 += gene_stats[constant_pct75] * gene_stats[constant_bases]
            accum_weighted_pct25 += gene_stats[constant_pct25] * gene_stats[constant_bases]
            accum_weighted_sd += gene_stats[constant_sd] * gene_stats[constant_bases]
            if constants.GC_CONTENT in gene_stats:
                accum_weighted_gc_content += gene_stats[constant_gc_content] * gene_stats[constant_bases]

            # aggregates data by chromosome
            chromosome = gene[constant_chromosome]
            if chromosome in observed_chromosomes:
                gene_count_by_chromosome[chromosome] += 1
                bases_by_chromosome[chromosome] += gene_stats[constant_bases]
                bases_lt15_by_chromosome[chromosome] += gene_stats[constant_bases_lt15x]
                bases_gte15_by_chromosome[chromosome] += gene_stats[constant_bases_gte15x]
                bases_gte30_by_chromosome[chromosome] += gene_stats[constant_bases_gte30x]
                bases_gte50_by_chromosome[chromosome] += gene_stats[constant_bases_gte50x]
                accum_mean_by_chromosome[chromosome] += gene_stats[constant_average]
                accum_weighted_median_by_chromosome[chromosome] += gene_stats[constant_median] * gene_stats[constant_bases]
                accum_weighted_pct75_by_chromosome[chromosome] += gene_stats[constant_pct75] * gene_stats[constant_bases]
                accum_weighted_pct25_by_chromosome[chromosome] += gene_stats[constant_pct25] * gene_stats[constant_bases]
                accum_weighted_sd_by_chromosome[chromosome] += gene_stats[constant_sd] * gene_stats[constant_bases]
                if constant_gc_content in gene_stats:
                    accum_weighted_gc_content_by_chromosome[chromosome] += gene_stats[constant_gc_content] * gene_stats[constant_bases]
            else:
                observed_chromosomes.append(chromosome)
                gene_count_by_chromosome[chromosome] = 1
                bases_by_chromosome[chromosome] = gene_stats[constant_bases]
                bases_lt15_by_chromosome[chromosome] = gene_stats[constant_bases_lt15x]
                bases_gte15_by_chromosome[chromosome] = gene_stats[constant_bases_gte15x]
                bases_gte30_by_chromosome[chromosome] = gene_stats[constant_bases_gte30x]
                bases_gte50_by_chromosome[chromosome] = gene_stats[constant_bases_gte50x]
                accum_mean_by_chromosome[chromosome] = gene_stats[constant_average]
                accum_weighted_median_by_chromosome[chromosome] = gene_stats[constant_median] * gene_stats[constant_bases]
                accum_weighted_pct75_by_chromosome[chromosome] = gene_stats[constant_pct75] * gene_stats[constant_bases]
                accum_weighted_pct25_by_chromosome[chromosome] = gene_stats[constant_pct25] * gene_stats[constant_bases]
                accum_weighted_sd_by_chromosome[chromosome] = gene_stats[constant_sd] * gene_stats[constant_bases]
                if constant_gc_content in gene_stats:
                    accum_weighted_gc_content_by_chromosome[chromosome] = gene_stats[constant_gc_content] * gene_stats[constant_bases]

            # aggregates data for autosomes
            if chromosome in constants.AUTOSOME_IDS:
                gene_count_autosomes += 1
                bases_autosomes += gene_stats[constant_bases]
                bases_lt15_autosomes += gene_stats[constant_bases_lt15x]
                bases_gte15_autosomes += gene_stats[constant_bases_gte15x]
                bases_gte30_autosomes += gene_stats[constant_bases_gte30x]
                bases_gte50_autosomes += gene_stats[constant_bases_gte50x]
                accum_mean_autosomes += gene_stats[constant_average]
                accum_weighted_median_autosomes += gene_stats[constant_median] * gene_stats[constant_bases]
                accum_weighted_pct75_autosomes += gene_stats[constant_pct75] * gene_stats[constant_bases]
                accum_weighted_pct25_autosomes += gene_stats[constant_pct25] * gene_stats[constant_bases]
                accum_weighted_sd_autosomes += gene_stats[constant_sd] * gene_stats[constant_bases]
                if constants.GC_CONTENT in gene_stats:
                    accum_weighted_gc_content_autosomes += gene_stats[constant_gc_content] * gene_stats[constant_bases]

    # computes statistics aggregated across the whole coding region
    mean = round(float(accum_mean) / gene_count, 3) if gene_count > 0 else 0.0  # the mean of the means
    lt_15x = round(float(bases_lt15) / bases, 5) if bases > 0 else 0.0
    gte_15x = round(float(bases_gte15) / bases, 5) if bases > 0 else 0.0
    gte_30x = round(float(bases_gte30) / bases, 5) if bases > 0 else 0.0
    gte_50x = round(float(bases_gte50) / bases, 5) if bases > 0 else 0.0
    median = round(float(accum_weighted_median) / bases, 3) if bases > 0 else float(0.0)
    pct75 = round(float(accum_weighted_pct75) / bases, 3) if bases > 0 else float(0.0)
    pct25 = round(float(accum_weighted_pct25) / bases, 3) if bases > 0 else float(0.0)
    sd = round(float(accum_weighted_sd) / bases, 3) if bases > 0 else float(0.0)
    gc_content = round(float(accum_weighted_gc_content) / bases, 5) if bases > 0 else float(0.0)
    stats = {
        constant_bases: bases,
        constant_average: mean,
        constant_median: median,
        constant_pct25: pct25,
        constant_pct75: pct75,
        constant_sd: sd,
        constants.LT15X: lt_15x,
        constants.GTE15X: gte_15x,
        constants.GTE30X: gte_30x,
        constants.GTE50X: gte_50x
    }
    if gc_content > 0.0:
        stats[constants.GC_CONTENT] = gc_content
    results[constants.STATISTICS] = stats

    # computes statistics aggregated by chromosome
    for chromosome in observed_chromosomes:
        mean_by_chromosome = round(float(accum_mean_by_chromosome[chromosome]) / gene_count_by_chromosome[chromosome], 3) \
            if gene_count_by_chromosome[chromosome] > 0 else 0.0  # the mean of the means
        lt_15x_by_chromosome = round(float(bases_lt15_by_chromosome[chromosome]) / bases_by_chromosome[chromosome], 5) \
            if bases_by_chromosome[chromosome] > 0 else 0.0
        gte_15x_by_chromosome = round(float(bases_gte15_by_chromosome[chromosome]) / bases_by_chromosome[chromosome], 5) \
            if bases_by_chromosome[chromosome] > 0 else 0.0
        gte_30x_by_chromosome = round(float(bases_gte30_by_chromosome[chromosome]) / bases_by_chromosome[chromosome], 5) \
            if bases_by_chromosome[chromosome] > 0 else 0.0
        gte_50x_by_chromosome = round(float(bases_gte50_by_chromosome[chromosome]) / bases_by_chromosome[chromosome], 5) \
            if bases_by_chromosome[chromosome] > 0 else 0.0
        median_by_chromosome = round(float(accum_weighted_median_by_chromosome[chromosome]) / bases_by_chromosome[chromosome], 3) \
            if bases_by_chromosome[chromosome] > 0 else float(0.0)
        pct75_by_chromosome = round(float(accum_weighted_pct75_by_chromosome[chromosome]) / bases_by_chromosome[chromosome], 3) \
            if bases_by_chromosome[chromosome] > 0 else float(0.0)
        pct25_by_chromosome = round(float(accum_weighted_pct25_by_chromosome[chromosome]) / bases_by_chromosome[chromosome], 3) \
            if bases_by_chromosome[chromosome] > 0 else float(0.0)
        sd_by_chromosome = round(float(accum_weighted_sd_by_chromosome[chromosome]) / bases_by_chromosome[chromosome], 3) \
            if bases_by_chromosome[chromosome] > 0 else float(0.0)
        gc_content_by_chromosome = round(
            float(accum_weighted_gc_content_by_chromosome.get(chromosome, 0.0)) / bases_by_chromosome[chromosome], 5) \
            if bases_by_chromosome[chromosome] > 0 else float(0.0)
        stats = {
            constant_chromosome: chromosome,
            constant_bases: bases_by_chromosome[chromosome],
            constant_average: mean_by_chromosome,
            constant_median: median_by_chromosome,
            constant_pct25: pct25_by_chromosome,
            constant_pct75: pct75_by_chromosome,
            constant_sd: sd_by_chromosome,
            constants.LT15X: lt_15x_by_chromosome,
            constants.GTE15X: gte_15x_by_chromosome,
            constants.GTE30X: gte_30x_by_chromosome,
            constants.GTE50X: gte_50x_by_chromosome
        }
        if gc_content_by_chromosome > 0.0:
            stats[constants.GC_CONTENT] = gc_content_by_chromosome
        results[constants.CHROMOSOMES].append(stats)

    # computes statistics for the autosomes
    mean_autosomes = round(float(accum_mean_autosomes) / gene_count_autosomes, 3) if gene_count_autosomes > 0 else 0.0
    lt_15x_autosomes = round(float(bases_lt15_autosomes) / bases_autosomes, 5) if bases_autosomes > 0 else 0.0
    gte_15x_autosomes = round(float(bases_gte15_autosomes) / bases_autosomes, 5) if bases_autosomes > 0 else 0.0
    gte_30x_autosomes = round(float(bases_gte30_autosomes) / bases_autosomes, 5) if bases_autosomes > 0 else 0.0
    gte_50x_autosomes = round(float(bases_gte50_autosomes) / bases_autosomes, 5) if bases_autosomes > 0 else 0.0
    median_autosomes = round(float(accum_weighted_median_autosomes) / bases_autosomes, 3) \
        if bases_autosomes > 0 else float(0.0)
    pct75_autosomes = round(float(accum_weighted_pct75_autosomes) / bases_autosomes, 3) \
        if bases_autosomes > 0 else float(0.0)
    pct25_autosomes = round(float(accum_weighted_pct25_autosomes) / bases_autosomes, 3) \
        if bases_autosomes > 0 else float(0.0)
    sd_autosomes = round(float(accum_weighted_sd_autosomes) / bases_autosomes, 3) if bases_autosomes > 0 else float(0.0)
    gc_content_autosomes = round(float(accum_weighted_gc_content_autosomes) / bases_autosomes, 5) \
        if bases_autosomes > 0 else float(0.0)
    stats = {
        constants.CHROMOSOME: constants.AUTOSOMES,
        constant_bases: bases_autosomes,
        constant_average: mean_autosomes,
        constant_median: median_autosomes,
        constant_pct25: pct25_autosomes,
        constant_pct75: pct75_autosomes,
        constant_sd: sd_autosomes,
        constants.LT15X: lt_15x_autosomes,
        constants.GTE15X: gte_15x_autosomes,
        constants.GTE30X: gte_30x_autosomes,
        constants.GTE50X: gte_50x_autosomes
    }
    if gc_content_autosomes > 0.0:
        stats[constants.GC_CONTENT] = gc_content_autosomes
    results[constants.CHROMOSOMES].append(stats)

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

    # NOTE: these variables are created to avoid a global lookup inside the loop
    constant_statistics = constants.STATISTICS
    constant_bases = constants.BASES
    constant_lt15x = constants.LT15X
    constant_gte15x = constants.GTE15X
    constant_gte30x = constants.GTE30X
    constant_gte50x = constants.GTE50X
    constant_average = constants.AVERAGE
    constant_median = constants.MEDIAN
    constant_pct75 = constants.PERCENTILE75
    constant_pct25 = constants.PERCENTILE25
    constant_sd = constants.SD
    constant_chromosome = constants.CHROMOSOME
    constant_chromosomes = constants.CHROMOSOMES
    constant_rmsd = constants.RMSD

    # Iterates each chromosome
    autosomes_stats = []
    for chromosome, regions in analysis_regions.iteritems():
        logging.info("Analysing chromosome %s" % chromosome)
        count_chunks = 0
        chr_list_rmsd = []
        chr_bases = 0
        chr_bases_lt15x = 0
        chr_bases_gte15x = 0
        chr_bases_gte30x = 0
        chr_bases_gte50x = 0
        chr_accum_mean = 0.0
        chr_accum_median = 0.0
        chr_accum_pct75 = 0.0
        chr_accum_pct25 = 0.0
        chr_accum_sd = 0.0

        # Iterates intervals in chunks of fixed size and stores the stats for each chunk
        for (start, end) in regions:
            logging.debug("Analysing region %s:%s-%s" % (chromosome, start, end))
            current_start = start
            current_end = min(current_start + chunk_size, end)
            # iterates over chunks
            while current_start < current_end:
                logging.debug("Analysing chunk %s:%s-%s" % (chromosome, current_start, current_end))
                # read coverage data
                coverages = bigwig_reader.read_bigwig_coverages(chromosome, current_start, current_end, strict=False)
                length = current_end - current_start
                if coverages is None:
                    chunk_bases = length
                    chunk_bases_lt15 = length
                    chunk_bases_gte15_lt30 = 0
                    chunk_bases_gte30_lt50 = 0
                    chunk_bases_gte50 = 0
                    chunk_mean = 0.0
                    chunk_weighted_median = 0.0
                    chunk_weighted_pct75 = 0.0
                    chunk_weighted_pct25 = 0.0
                    chunk_weighted_sd = 0.0
                    chunk_rmsd = 0
                else:
                    chunk_bases = 0
                    chunk_bases_lt15 = 0
                    chunk_bases_gte15_lt30 = 0
                    chunk_bases_gte30_lt50 = 0
                    chunk_bases_gte50 = 0
                    chunk_accum_coverage = 0

                    # iterate chunk
                    for coverage in coverages:
                        chunk_bases += 1
                        chunk_accum_coverage += coverage
                        # NOTE: this is optmised for performance, not for code simplicity or beauty
                        if coverage >= 50:
                            chunk_bases_gte50 += 1
                        elif coverage >= 30:
                            chunk_bases_gte30_lt50 += 1
                        elif coverage >= 15:
                            chunk_bases_gte15_lt30 += 1
                        else:
                            chunk_bases_lt15 += 1

                    # compute statistics for the chunk
                    chunk_mean = round(float(chunk_accum_coverage) / chunk_bases, 3) if chunk_bases > 0 else 0.0
                    if getattr(coverages, "size", 0) > 0:  # gets the size of the numpy array, if None returns 0
                        chunk_weighted_median = round((float(np.median(coverages))) * chunk_bases, 3)
                        chunk_weighted_pct75 = round((float(np.percentile(coverages, 75))) * chunk_bases, 3)
                        chunk_weighted_pct25 = round((float(np.percentile(coverages, 25))) * chunk_bases, 3)
                        chunk_weighted_sd = round(float(np.std(coverages)) * chunk_bases, 3)
                        if chunk_mean == 0:
                            # As coverage values are positive values we infer that all values are zero
                            # this may speed up things for missing long regions in the bigwig file, if any
                            chunk_rmsd = 0
                        else:
                            # Gets the squared root sum of squares of the deviation from the mean
                            chunk_rmsd = np.sqrt((np.sum([(x - chunk_mean) ** 2 for x in coverages])) / length)
                    else:
                        chunk_weighted_median = 0.0
                        chunk_weighted_pct75 = 0.0
                        chunk_weighted_pct25 = 0.0
                        chunk_weighted_sd = 0.0
                        chunk_rmsd = 0

                count_chunks += 1
                chr_bases += chunk_bases
                chr_bases_lt15x += chunk_bases_lt15
                chr_bases_gte15x += chunk_bases_gte15_lt30
                chr_bases_gte15x += chunk_bases_gte30_lt50
                chr_bases_gte15x += chunk_bases_gte50
                chr_bases_gte30x += chunk_bases_gte30_lt50
                chr_bases_gte30x += chunk_bases_gte50
                chr_bases_gte50x += chunk_bases_gte50
                chr_accum_mean += chunk_mean
                chr_accum_median += chunk_weighted_median
                chr_accum_pct25 += chunk_weighted_pct25
                chr_accum_pct75 += chunk_weighted_pct75
                chr_accum_sd += chunk_weighted_sd
                chr_list_rmsd.append(chunk_rmsd)

                # prepares coordinates for next chunk
                current_start = current_end
                current_end = min(current_start + chunk_size, end)

        # Computes the statistics per chromosome
        if chr_bases > 0:
            chr_lt15x = round(float(chr_bases_lt15x) / chr_bases, 5)
            chr_gte15x = round(float(chr_bases_gte15x) / chr_bases, 5)
            chr_gte30x = round(float(chr_bases_gte30x) / chr_bases, 5)
            chr_gte50x = round(float(chr_bases_gte50x) / chr_bases, 5)
            chr_median = round(float(chr_accum_median) / chr_bases, 3)
            chr_pct25 = round(float(chr_accum_pct25) / chr_bases, 3)
            chr_pct75 = round(float(chr_accum_pct75) / chr_bases, 3)
            chr_sd = round(float(chr_accum_sd) / chr_bases, 3)
        else:
            chr_lt15x = 0.0
            chr_gte15x = 0.0
            chr_gte30x = 0.0
            chr_gte50x = 0.0
            chr_median = 0.0
            chr_pct25 = 0.0
            chr_pct75 = 0.0
            chr_sd = 0.0
        chromosome_stats = {
            constant_chromosome: chromosome,
            constant_rmsd: round(float(np.median(chr_list_rmsd)), 3) if chr_list_rmsd else 0.0,
            constant_average: round(float(chr_accum_mean) / count_chunks, 3) if count_chunks > 0 else 0.0,
            constant_bases: chr_bases,
            constant_lt15x: chr_lt15x,
            constant_gte15x: chr_gte15x,
            constant_gte30x: chr_gte30x,
            constant_gte50x: chr_gte50x,
            constant_median: chr_median,
            constant_pct25: chr_pct25,
            constant_pct75: chr_pct75,
            constant_sd: chr_sd,
        }
        results[constant_chromosomes].append(chromosome_stats)
        if chromosome in constants.AUTOSOME_IDS:
            autosomes_stats.append(chromosome_stats)
        logging.info("Whole genome statistics for chromosome %s computed!" % chromosome)

    # Aggregates statistics for the whole genome (important to do before addings autosomes!)
    chromosome_stats = aggregate_chromosomes(results[constant_chromosomes])
    results[constant_statistics] = chromosome_stats
    logging.info("Aggregated whole genome statistics!")

    # Aggregates statistics for the autosomes
    aggregated_autosomes_stats = aggregate_chromosomes(autosomes_stats)
    aggregated_autosomes_stats[constants.CHROMOSOME] = constants.AUTOSOMES
    results[constant_chromosomes].append(aggregated_autosomes_stats)
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
