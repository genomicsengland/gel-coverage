import pysam
from collections import defaultdict
import argparse
from coverage_summary import __lt__



def get_chr_sizes(path_to_bam):

    samfile = pysam.AlignmentFile(path_to_bam, "rb")

    chromosomes = defaultdict()
    for chr_details in samfile.header["SQ"]:
        reference_name = chr_details["SN"]
        reference_length = chr_details["LN"]
        chromosomes[reference_name] = {}
        chromosomes[reference_name]["reference_length"] = reference_length

    return chromosomes


def main():
    parser = argparse.ArgumentParser(description='Sam/Bam file filter')
    parser.add_argument('--bam', metavar='bam', default=False,
                        help='bam file')
    parser.add_argument('--output', metavar='output', default=False,
                        help='output file')

    args = parser.parse_args()

    sizes = get_chr_sizes(args.bam)
    output = open(args.output,"w")

    for chr in sorted(sizes,cmp=__lt__):
        length = sizes[chr]["reference_length"]
        line = "chr" + str(chr) + "\t" + str(length) + "\n"
        line = line.replace("chrMT","chrM")
        output.write(line)

    output.close()


if __name__ == '__main__':
    main()
