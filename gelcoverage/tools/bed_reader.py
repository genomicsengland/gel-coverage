import pybedtools


class BedReader:

    def __init__(self, bed):
        # Opens the bigwig file for reading
        self.bed = bed
        self.is_null_bed = self.bed is None or self.bed == ""
        if not self.is_null_bed:
            self.chromosomes = self.__get_chromosomes()
            self.has_chr_prefix = any([chromosome.startswith("chr") for chromosome in self.chromosomes])

    def __get_chromosomes(self):
        """
        Gets the list of chromosomes present in the BED file
        :return:
        """
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
        :param bed: the bed file
        :return: the dict structure
        """
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