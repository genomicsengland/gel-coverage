import pyBigWig
import logging


class UncoveredIntervalException(Exception):
    pass

class BigWigReader:

    def __init__(self, bigwig):
        # Opens the bigwig file for reading
        self.bigwig = bigwig
        self.reader = pyBigWig.open(self.bigwig)
        self.reported_unexisting_chr = []

    def read_bigwig_coverages(self, chromosome, start, end, strict = True):
        """
        Reads the coverage values in a region from the bigwig file
        :param chromosome: the region chromosome
        :param start: the start position
        :param end: the end position
        :strict: when true if a position not present in the bigwig is queried it will raise an error
        :return: the sequence of coverages (one integer per position)
        """
        logging.debug("Queries bigwig for %s:%s-%s in %s mode" %
                      (chromosome, str(start), str(end), "strict" if strict else "non strict"))
        # Queries the bigwig for a specific interval
        if start == end:  # do we really need this?
            end += 1
        # Converts chromosome identifiers to the expected format with "chr" prefix
        if type(chromosome) == int or not chromosome.startswith("chr"):
            chromosome = "chr" + str(chromosome)
        # Read from the bigwig file
        # TODO: why our bigwig has "chr" prefix? BAMs don't.
        # ANSWER: this is the bigwig generation pipeline that is adding the chr prefix
        if strict:
            try:
                coverages = self.reader.values(chromosome, start, end)
            except RuntimeError, e:
                # When the queried interval is not present in the bigwig file it returns all 0s coverage and logs it
                logging.warn("Missing interval in bigwig %s:%s-%s" % (chromosome, start, end))
                raise UncoveredIntervalException()
        else:
            # By using the function intervals we make sure that we are not querying coordinates not present in the
            # bigwig file
            try:
                intervals = self.reader.intervals(chromosome, start, end)
                coverages = [coverage for _, _, coverage in intervals]
            except RuntimeError, e:
                if chromosome not in self.reported_unexisting_chr:
                    logging.warn("Missing chromosome in bigwig %s" % (chromosome))
                    self.reported_unexisting_chr.append(chromosome)
                coverages = [0] * (end- start + 1)
                #raise UncoveredIntervalException()  # keeps going
        return coverages

    def get_chromosome_lengths(self):
        """
        get chromosome lengths from header of bigWig file

        :param bw: pyBigWig file object
        :param format: specify "dict" to return a dict instead of a bed recognisable format
        :return: list of chromosomes and start and end positions
        """
        return {chromosome : [(0, length)] for (chromosome, length) in self.reader.chroms().iteritems()}
