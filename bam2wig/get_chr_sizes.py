import pysam
from collections import defaultdict
import argparse

__author__ = 'mparker'


def lt_helper(a):
    """
    helper to sort chromosomes properly

    :param a: sort object
    :return:
    """
    try:
        return int(a)
    except ValueError:
        if a == "X":
            return 24
        elif a == "Y":
            return 25
        elif a == "MT" or a == "M":
            return 26
        else:
            return 27


def __lt__(a, b):
    """
    run the chromosome sort helper

    :param a:
    :param b:
    :return: proper sorted chromosomes
    """
    return cmp(lt_helper(a), lt_helper(b))


def get_chr_sizes(path_to_bam):
    """
    read bam header and produce dict of chromosome lengths

    :param path_to_bam: path to the bam filr
    :return: dict of chromosomes and their reference length
    """
    samfile = pysam.AlignmentFile(path_to_bam, "rb")
    chromosomes = defaultdict()
    for chr_details in samfile.header["SQ"]:
        reference_name = chr_details["SN"]
        reference_length = chr_details["LN"]
        chromosomes[reference_name] = {}
        chromosomes[reference_name]["reference_length"] = reference_length
    return chromosomes


def main():
    parser = argparse.ArgumentParser(description='Script to get chromosome sizes from bam header for tools that need a \
                                                 chr size file')
    parser.add_argument('--bam', metavar='bam', default=False,
                        help='bam file')
    parser.add_argument('--output', metavar='output', default=False,
                        help='output file')
    args = parser.parse_args()
    sizes = get_chr_sizes(args.bam)
    output = open(args.output,"w")
    for chr in sorted(sizes,cmp=__lt__):
        length = sizes[chr]["reference_length"]
        line =  "%s\t%s\n" % (str(chr), str(length))
        output.write(line)
    output.close()

if __name__ == '__main__':
    main()
