import pyBigWig
import logging


class BigWigReader:

    def __init__(self, bigwig):
        # Opens the bigwig file for reading
        self.bigwig = bigwig
        self.reader = pyBigWig.open(self.bigwig)

    def read_bigwig_coverages(self, chromosome, start, end, strict = True):
        """
        Reads the coverage values in a region from the bigwig file
        :param chromosome: the region chromosome
        :param start: the start position
        :param end: the end position
        :strict: when true if a position not present in the bigwig is queried it will raise an error
        :return: the sequence of coverages (one integer per position)
        """
        # Queries the bigwig for a specific interval
        if start == end:  # do we really need this?
            end += 1
        # Converts chromosome identifiers to the expected format with "chr" prefix
        if type(chromosome) == int or not chromosome.startswith("chr"):
            chromosome = "chr" + str(chromosome)
        # Read from the bigwig file
        try:
            # TODO: why our bigwig has "chr" prefix? BAMs don't
            if strict:
                coverages = self.reader.values(chromosome, start, end)
            else:
                # By using the function intervals we make sure that we are not querying coordinates not present in the
                # bigwig file
                intervals = self.reader.intervals(chromosome, start, end)
                coverages = [coverage for _, _, coverage in intervals]
        except RuntimeError, e:
            # TODO: deal with this errors
            logging.error("Querying for unexisting interval %s:%s-%s" % (chromosome, start, end))
            raise e
        return coverages

    def get_chromosome_lengths(self, format):
        """
        get chromosome lengths from header of bigWig file

        :param bw: pyBigWig file object
        :param format: specify "dict" to return a dict instead of a bed recognisable format
        :return: list of chromosomes and length in a bed recognisable format
        """
        return self.reader.chroms()
