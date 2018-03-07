import pybedtools
import logging


class BedInterval:

    def __init__(self, bed_row):
        self.bedline = bed_row.rstrip('\n').split('\t')
        try:
            self.chrom = self.bedline[0]
            self.start = self.bedline[1]
            self.end = self.bedline[2]
            self.name = self.bedline[3]
            self.score = self.bedline[4]
            self.strand = self.bedline[5]
        except IndexError, e:
            logging.info("Expected something like "
                         "'chr17   69327052        69327101        ABCA5|ENST00000392676|exon1     0.68    -'")
            logging.error("Bad BED row '{}'".format(bed_row))
            raise e


class BedReader:

    def __init__(self, bed):
        # Opens the bigwig file for reading
        logging.debug("Creating the Bed reader...")
        self.bed = bed
        self.is_null_bed = self.bed is None or self.bed == ""
        if not self.is_null_bed:
            self.chromosomes = self.__get_chromosomes()
            self.has_chr_prefix = any([chromosome.startswith("chr") for chromosome in self.chromosomes])
        logging.debug("Bed reader created!")

    def __get_chromosomes(self):
        """
        Gets the list of chromosomes present in the BED file
        :return:
        """
        logging.debug("Reading chromosomes from BED...")
        reader = pybedtools.BedTool(self.bed)
        chromosomes = []
        for interval in reader:
            chromosome = interval.chrom
            if chromosome not in chromosomes:
                chromosomes.append(chromosome)
        return chromosomes

    def get_regions_dictionary(self):
        """
        Builds the regions in a bed file into a python dictionary structure
        :return: the dict structure
        """
        logging.debug("Loading Bed regions into a dict...")
        reader = pybedtools.BedTool(self.bed)
        regions = {}
        for interval in reader:
            chromosome = interval.chrom
            start = int(interval.start)
            end = int(interval.end)
            if chromosome not in regions:
                regions[chromosome] = []
            regions[chromosome].append((start, end))
        return regions
