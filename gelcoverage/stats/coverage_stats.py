import numpy as np
import logging
import itertools
import operator
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
            current_gap["s"] = start_position + idx
        elif value >= coverage_threshold and open_gap:
            open_gap = False
            current_gap["e"] = start_position + idx - 1
            current_gap["l"] = current_gap["e"] - current_gap["s"] + 1
            gaps.append(current_gap)
            current_gap = {}
    # Closes the last gap when it extends until the last position
    if open_gap:
        current_gap["e"] = end
        current_gap["l"] = current_gap["e"] - current_gap["s"] + 1
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
    stats["bases"] = len(coverages)
    stats["avg"] = round(float(np.mean(coverages)), 3)
    stats["med"] = round(float(np.median(coverages)), 3)
    stats["pct75"] = round(float(np.percentile(coverages, 75)), 3)
    stats["pct25"] = round(float(np.percentile(coverages, 25)), 3)
    stats["bases_lt_15x"] = int(np.sum(1 for x in coverages if x < 15))
    stats["bases_gte_15x"] = int(np.sum(1 for x in coverages if x >= 15))
    stats["bases_gte_30x"] = int(np.sum(1 for x in coverages if x >= 30))
    stats["bases_gte_50x"] = int(np.sum(1 for x in coverages if x >= 50))
    stats["%<15x"] = round(float(stats["bases_lt_15x"]) / stats["bases"], 5)
    stats["%>=15x"] = round(float(stats["bases_gte_15x"]) / stats["bases"], 5)
    stats["%>=30x"] = round(float(stats["bases_gte_30x"]) / stats["bases"], 5)
    stats["%>=50x"] = round(float(stats["bases_gte_50x"]) / stats["bases"], 5)
    if gc_content is not None:  # GC content is not provided for padded exons
        stats["gc"] = gc_content

    return stats

def compute_transcript_level_statistics(exons):
    """
    Computes coverage and GC content statistics at gene level by aggregating the statistics at exon level.
    Median and percentiles are estimated by weighting the per-exon metric by the number of bases.
    :param exons: list of exon coverage and GC content statistics
    :return: the coverage and GC content gene statistics in JSON-friendly format
    """
    exons_stats = [x["stats"] for x in exons]
    total_bases = int(np.sum([x["bases"] for x in exons_stats]))
    bases_lt_15x = int(np.sum([x["bases_lt_15x"] for x in exons_stats]))
    bases_gte_15x = int(np.sum([x["bases_gte_15x"] for x in exons_stats]))
    bases_gte_30x = int(np.sum([x["bases_gte_30x"] for x in exons_stats]))
    bases_gte_50x = int(np.sum([x["bases_gte_50x"] for x in exons_stats]))
    stats = {
        "bases": total_bases,
        "avg": round(float(np.mean([x["avg"] for x in exons_stats])), 3),
        "med": round(float(np.sum(
            [x["med"] * x["bases"] for x in exons_stats])) / total_bases, 3),
        "pct25": round(float(np.sum(
            [x["pct25"] * x["bases"] for x in exons_stats])) / total_bases, 3),
        "pct75": round(float(np.sum(
            [x["pct75"] * x["bases"] for x in exons_stats])) / total_bases, 3),
        "%<15x" : round(float(bases_lt_15x) / total_bases, 5),
        "%>=15x": round(float(bases_gte_15x) / total_bases, 5),
        "%>=30x": round(float(bases_gte_30x) / total_bases, 5),
        "%>=50x": round(float(bases_gte_50x) / total_bases, 5),
        "bases_lt_15x": bases_lt_15x,
        "bases_gte_15x": bases_gte_15x,
        "bases_gte_30x": bases_gte_30x,
        "bases_gte_50x": bases_gte_50x
    }
    try:
        stats["gc"] = round(float(np.sum(
            [x["gc"] * x["bases"] for x in exons_stats]) / total_bases), 5)
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
        "stats": None,
        "chrs": []
    }
    # Compute the stats aggregated for union transcript
    genes_stats = [x["union_tr"]["stats"] for x in genes]
    total_bases = int(np.sum([x["bases"] for x in genes_stats]))
    bases_lt_15x = int(np.sum([x["bases_lt_15x"] for x in genes_stats]))
    bases_gte_15x = int(np.sum([x["bases_gte_15x"] for x in genes_stats]))
    bases_gte_30x = int(np.sum([x["bases_gte_30x"] for x in genes_stats]))
    bases_gte_50x  = int(np.sum([x["bases_gte_50x"] for x in genes_stats]))
    results["stats"] = {
        "bases" : total_bases,
        "avg" : round(float(np.mean([x["avg"] for x in genes_stats])), 3),
        "med" : round(float(np.sum(
            [x["med"] * x["bases"] for x in genes_stats]) / total_bases), 3),
        "pct75" : round(float(np.sum(
            [x["pct75"] * x["bases"] for x in genes_stats]) / total_bases), 3),
        "pct25" : round(float(np.sum(
            [x["pct25"] * x["bases"] for x in genes_stats]) / total_bases), 3),
        "%<15x" : round(float(bases_lt_15x) / total_bases, 5),
        "%>=15x" : round(float(bases_gte_15x) / total_bases, 5),
        "%>=30x" : round(float(bases_gte_30x) / total_bases, 5),
        "%>=50x" : round(float(bases_gte_50x) / total_bases, 5)
    }
    # Compute the stats disaggregated by chromosome
    chr2stats = [(x["chr"], x["union_tr"]["stats"]) for x in genes]
    def groupby_chromosome(list_of_tuples):
        it = itertools.groupby(list_of_tuples, operator.itemgetter(0))
        for chromosome, subiter in it:
            yield chromosome, [item[1] for item in subiter]
    chromosome_stats = dict(groupby_chromosome(chr2stats))
    for chromosome, chr_stats in chromosome_stats.iteritems():
        chr_total_bases = int(np.sum([x["bases"] for x in chr_stats]))
        chr_bases_lt_15x = int(np.sum([x["bases_lt_15x"] for x in chr_stats]))
        chr_bases_gte_15x = int(np.sum([x["bases_gte_15x"] for x in chr_stats]))
        chr_bases_gte_30x = int(np.sum([x["bases_gte_30x"] for x in chr_stats]))
        chr_bases_gte_50x = int(np.sum([x["bases_gte_50x"] for x in chr_stats]))
        results["chrs"].append(
            {
                "chr" : chromosome,
                "bases": chr_total_bases,
                "avg": round(float(np.mean([x["avg"] for x in chr_stats])), 3),
                "med": round(float(np.sum(
                    [x["med"] * x["bases"] for x in chr_stats]) / chr_total_bases), 3),
                "pct75": round(float(np.sum(
                    [x["pct75"] * x["bases"] for x in chr_stats]) / chr_total_bases), 3),
                "pct25": round(float(np.sum(
                    [x["pct25"] * x["bases"] for x in chr_stats]) / chr_total_bases), 3),
                "%<15x": round(float(chr_bases_lt_15x) / chr_total_bases, 5),
                "%>=15x": round(float(chr_bases_gte_15x) / chr_total_bases, 5),
                "%>=30x": round(float(chr_bases_gte_30x) / chr_total_bases, 5),
                "%>=50x": round(float(chr_bases_gte_50x) / chr_total_bases, 5)
            }
        )
        logging.info("Coding region statistics for chromosome %s computed!" % chromosome)
    logging.info("Coding region statistics computed!")
    return results

def compute_whole_genome_statistics(bigwig_reader, bed = None, chunk_size = 100000):
    """
    Iterates through the whole genome in a sliding window to obtain some metrics
    :param bigwig_reader:
    :return:
    """
    logging.info("Computing whole genome statistics...")
    results = {
        "stats" : None,
        "chrs" : []
    }
    if bed is None:
        logging.info("Running on all chromosomes defined in the bigwig.")
        analysis_regions = bigwig_reader.get_chromosome_lengths()
    else:
        logging.info("Running on the regions provided in a bed file in --wg-region.")
        analysis_regions = bed_parser.get_regions_dictionary(bed)
    # Iterates each chromosome
    chr_stats = {}
    for chromosome, regions in analysis_regions.iteritems():
        chr_stats[chromosome] = {
            "rmsd" : [],
            "avg" : [],
            "bases" : [],
            "bases_lt_15x" : [],
            "bases_gte_15x": [],
            "bases_gte_30x": [],
            "bases_gte_50x": [],
            "med" : [],
            "pct25" : [],
            "pct75" : []
        }
        # Iterates intervals in chunks of fixed size and stores the stats for each chunk
        for (start, end) in regions:
            logging.debug("Analysing region %s:%s-%s" % (chromosome, start, end))
            current_start = start
            current_end = min(current_start + chunk_size - 1, end)
            while current_end < end:
                logging.debug("Analysing chunk %s:%s-%s" % (chromosome,
                                                            current_start,
                                                            current_end))
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
                chr_stats[chromosome]["bases"].append(length)
                chr_stats[chromosome]["bases_lt_15x"].append(
                    np.sum(1 for coverage in coverages if coverage < 15)
                )
                chr_stats[chromosome]["bases_gte_15x"].append(
                    np.sum(1 for coverage in coverages if coverage >= 15)
                )
                chr_stats[chromosome]["bases_gte_30x"].append(
                    np.sum(1 for coverage in coverages if coverage >= 30)
                )
                chr_stats[chromosome]["bases_gte_50x"].append(
                    np.sum(1 for coverage in coverages if coverage >= 50)
                )
                chr_stats[chromosome]["avg"].append(chunk_mean)
                chr_stats[chromosome]["rmsd"].append(chunk_rmsd)
                chr_stats[chromosome]["med"].append(np.median(coverages))
                chr_stats[chromosome]["pct25"].append(np.percentile(coverages, 25))
                chr_stats[chromosome]["pct75"].append(np.percentile(coverages, 75))
                current_start = current_end + 1
                current_end = min(current_start + chunk_size - 1, end)
        # Set the statistics per chromosome
        chr_total_bases = np.sum(chr_stats[chromosome]["bases"])
        chr_bases_lt_15x = np.sum(chr_stats[chromosome]["bases_lt_15x"])
        chr_bases_gte_15x = np.sum(chr_stats[chromosome]["bases_gte_15x"])
        chr_bases_gte_30x = np.sum(chr_stats[chromosome]["bases_gte_30x"])
        chr_bases_gte_50x = np.sum(chr_stats[chromosome]["bases_gte_50x"])
        results["chrs"].append({
            "chr" : chromosome,
            "uneveness" : round(float(np.median(chr_stats[chromosome]["rmsd"])), 3),
            "avg" : round(float(np.mean(chr_stats[chromosome]["avg"])), 3),
            "bases" : int(chr_total_bases),
            "%<15x" : round(float(chr_bases_lt_15x) / chr_total_bases, 5),
            "%>=15x" : round(float(chr_bases_gte_15x) / chr_total_bases, 5),
            "%>=30x" : round(float(chr_bases_gte_30x) / chr_total_bases, 5),
            "%>=50x" : round(float(chr_bases_gte_50x) / chr_total_bases, 5),
            "med" : round(float(np.sum([median * weight
                                            for median, weight in
                                            zip(chr_stats[chromosome]["med"],
                                                chr_stats[chromosome]["bases"])])) /
                             chr_total_bases, 3),
            "pct75" : round(float(np.sum([pct75 * weight
                                           for pct75, weight in
                                           zip(chr_stats[chromosome]["pct75"],
                                               chr_stats[chromosome]["bases"])])) /
                            chr_total_bases, 3),
            "pct25" : round(float(np.sum([pct25 * weight
                                          for pct25, weight in
                                          zip(chr_stats[chromosome]["pct25"],
                                              chr_stats[chromosome]["bases"])])) /
                            chr_total_bases, 3)
        })
        logging.info("Whole genome statistics for chromosome %s computed!" % chromosome)
    # Set the statistics for the whole genome
    def aggregate_chromosomes(dictionary, field):
        return sum([dictionary[chromosome][field] for chromosome in dictionary.keys()], [])
    total_bases = int(np.sum(aggregate_chromosomes(chr_stats, "bases")))
    bases_lt_15x = int(np.sum(aggregate_chromosomes(chr_stats, "bases_lt_15x")))
    bases_gte_15x = int(np.sum(aggregate_chromosomes(chr_stats, "bases_gte_15x")))
    bases_gte_30x = int(np.sum(aggregate_chromosomes(chr_stats, "bases_gte_30x")))
    bases_gte_50x = int(np.sum(aggregate_chromosomes(chr_stats, "bases_gte_50x")))
    results["stats"] = {
        "uneveness" : round(float(np.median(aggregate_chromosomes(chr_stats, "rmsd"))), 3),
        "avg" : round(float(np.mean(aggregate_chromosomes(chr_stats, "avg"))), 3),
        "bases" : total_bases,
        "%<15x" : round(float(bases_lt_15x) / total_bases, 5),
        "%>=15x": round(float(bases_gte_15x) / total_bases, 5),
        "%>=30x": round(float(bases_gte_30x) / total_bases, 5),
        "%>=50x": round(float(bases_gte_50x) / total_bases, 5),
        "med": round(float(np.sum(
            [median * length for median, length in zip(aggregate_chromosomes(chr_stats, "med"),
                                                       aggregate_chromosomes(chr_stats, "bases"))])) /
                        np.sum(aggregate_chromosomes(chr_stats, "bases")), 3),
        "pct75": round(float(np.sum(
            [pct75 * length for pct75, length in zip(aggregate_chromosomes(chr_stats, "pct75"),
                                                     aggregate_chromosomes(chr_stats, "bases"))])) /
                       np.sum(aggregate_chromosomes(chr_stats, "bases")), 3),
        "pct25": round(float(np.sum(
            [pct25 * length for pct25, length in zip(aggregate_chromosomes(chr_stats, "pct25"),
                                                     aggregate_chromosomes(chr_stats, "bases"))])) /
                       np.sum(aggregate_chromosomes(chr_stats, "bases")), 3)
    }
    logging.info("Whole genome statistics computed!")
    return results


