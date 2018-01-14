import pyBigWig
import logging


class UncoveredIntervalException(Exception):

    pass


class BigWigReader:

    def __init__(self, bigwig):
        # Opens the bigwig file for reading
        logging.debug("Creating the BigWig reader...")
        self.bigwig = bigwig
        self.reader = pyBigWig.open(self.bigwig)
        self.chromosomes = self.reader.chroms().keys()
        self.has_chr_prefix = any([chromosome.startswith("chr") for chromosome in self.chromosomes])
        self.reported_unexisting_chr = []
        logging.debug("BigWig reader created!")

    def read_bigwig_coverages(self, chromosome, start, end, strict=True):
        """
        Reads the coverage values in a region from the bigwig file
        :param chromosome: the region chromosome
        :param start: the start position
        :param end: the end position
        :param strict: when true if a position not present in the bigwig is queried it will raise an error
        :return: the sequence of coverages (one integer per position)
        """
        logging.debug("Queries bigwig for %s:%s-%s in %s mode" %
                      (chromosome, str(start), str(end), "strict" if strict else "non strict"))
        if start == end:  # do we really need this?
            end += 1
        # Read from the bigwig file
        try:
            coverages = self.reader.values(chromosome, start, end, numpy=True)
        except RuntimeError:
            if strict:
                # When the queried interval is not present in the bigwig file it returns all 0s coverage and logs it
                # logging.debug("Missing interval in bigwig %s:%s-%s" % (chromosome, start, end))
                raise UncoveredIntervalException()
            else:
                # if chromosome not in self.reported_unexisting_chr:
                #     logging.debug("Missing chromosome in bigwig %s" % chromosome)
                #     self.reported_unexisting_chr.append(chromosome)
                coverages = None
        return coverages

    def get_chromosome_lengths(self):
        """
        get chromosome lengths from header of bigWig file
        :return: list of chromosomes and start and end positions
        """
        logging.debug("Getting chromosomes and their lenght from BigWig header")
        return {chromosome: [(0, length)] for (chromosome, length) in self.reader.chroms().iteritems()}
