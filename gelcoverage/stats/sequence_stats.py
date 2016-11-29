import re


def compute_gc_content(sequence):
    """
    Calculates gc content of a DNA sequence
    :param sequence: a string of ATGC's
    :return: the gc fraction of the region
    """
    gc_count = 0
    total_base_count = 0
    gc_count += len(re.findall("[GC]", sequence))
    total_base_count += len(re.findall("[GCTA]", sequence))
    gc_fraction = round(float(gc_count) / total_base_count, 5)
    return gc_fraction
