import pybedtools


def get_regions_dictionary(bed):
    """
    Builds the regions in a bed file into a python dictionary structure
    :param bed: the bed file
    :return: the dict structure
    """
    intervals = pybedtools.BedTool(bed)
    regions = {}
    for interval in intervals:
        chromosome = interval.chrom
        start = int(interval.start)
        end = int(interval.end)
        if chromosome not in regions:
            regions[chromosome] = []
        regions[chromosome].append((start, end))
    return regions