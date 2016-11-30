import pyBigWig
import logging


class BigWigReader:

    def __init__(self, bigwig):
        # Opens the bigwig file for reading
        self.bigwig = bigwig
        self.reader = pyBigWig.open(self.bigwig)

    def read_bigwig_coverages(self, chromosome, start, end):
        """
        Reads the coverage values in a region from the bigwig file
        :param chromosome: the region chromosome
        :param start: the start position
        :param end: the end position
        :return: the sequence of coverages (one integer per position)
        """
        # Queries the bigwig for a specific interval
        if start == end:  # do we really need this?
            end += 1
        # Read from the bigwig file
        try:
            # TODO: why our bigwig has "chr" prefix? BAMs don't
            coverages = self.reader.values("chr" + str(chromosome), start, end)
        except RuntimeError, e:
            # TODO: deal with this errors
            logging.error("Querying for unexisting interval %s:%s-%s" % (chromosome, start, end))
            raise e
        return coverages