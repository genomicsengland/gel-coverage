import pybedtools
import logging
import subprocess
import argparse
from collections import defaultdict
import tqdm
import pyBigWig
import pandas
import json
import urllib2



def intervals(n_file,other_chr):
    ns = pybedtools.BedTool(n_file)
    ns_chrs=list()
    for feature in ns:
        ns_chrs.append(feature.chrom)
    all = pybedtools.BedTool(other_chr)

    missing=list()
    for feature in all:
        all_chr=feature.chrom
        if all_chr.replace("chr","") not in ns_chrs:
            missing.append((feature.chrom.replace("chr",""),feature.start,feature.end))

    to_cat=pybedtools.BedTool(missing)
    all = ns.cat(to_cat,postmerge=False)

    return all

def lt_helper(a):
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
    return cmp(lt_helper(a), lt_helper(b))

def chromosome_sizes(chr_sizes):
    chomosomes = dict()
    with open(chr_sizes) as fp:
        for line in fp:
            chrom,length=line.split("\t")
            chomosomes[str(chrom.replace("chr",""))] = int(length)
    return chomosomes


def other_chrs(bw):
    chromosomes=list()
    for i in range(1,24):
        chromosome = "chr"+str(i)
        if i == 23:
            chromosome = "chrX"
        elif i == 24:
            chromosome = "chrY"
        length=int(bw.chroms(chromosome))
        chromosomes.append((chromosome,1,100))

    return chromosomes

def coverage_counts_wg(n_file,bw):
    bw = pyBigWig.open( bw )

    other = other_chrs(bw)
    ns = intervals(n_file,other)
    return generic_coverage(bw,ns)


def coverage_counts_exon(bw):
    bw = pyBigWig.open( bw )
    exons=make_exons_bed()
    return generic_coverage(bw,exons)

def generic_coverage(bw,bed):
    result = defaultdict()

    for interval in bed:
        print interval
        chrom = interval.chrom
        start = interval.start
        end = interval.end

        cov = bw.values("chr"+str(chrom), start, end)
        result[chrom]=defaultdict(lambda: 0)
        for coverage in cov:
            result[chrom][int(coverage)]+=1

    order = sorted(result,cmp=__lt__)
    data = pandas.DataFrame.from_dict(result,orient="columns").fillna(0)
    sorted_data = data[order]

    sorted_data["Total"] = sorted_data.sum(axis=1)
    sorted_data.index.name = 'Coverage'
    return sorted_data


def make_exons_bed():
    url="http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/latest/hsapiens/feature/gene/all?include=chromosome,transcripts.exons.start,transcripts.exons.end"
    exons=json.load(urllib2.urlopen(url))["response"][0]["result"]
    all_exons=list()
    for result in exons:
        for transcript in result["transcripts"]:
            for exons in transcript["exons"]:
                all_exons.append((result["chromosome"],exons["start"],exons["end"]))
    bed = pybedtools.BedTool(all_exons)
    bed = bed.sort()
    bed = bed.merge()
    return bed


def main():
    parser = argparse.ArgumentParser(description='Sam/Bam file filter')
    parser.add_argument('--genome_n', metavar='Ns_in_genome', type=str, default="/accelrys/apps/gel/genomeref/data/human/human_g1k_v37_NonN_Regions.CHR.bed",
                        help='a bed file of the non-N regions in genome')
    # parser.add_argument('--chr_sizes', metavar='chr_sizes', default=False,
    #                     help='chromosome size file')
    parser.add_argument('--bw', metavar='bw', default=False,
                        help='bw file')
    parser.add_argument('--output', metavar='output', default=False,
                        help='output file')
    parser.add_argument('--debug', action='store_true', default=False,
                        help='If it is present then debugging is turned on')


    args = parser.parse_args()

    logger = logging.getLogger(__name__)

    level = logging.INFO
    if args.debug:
        level = logging.DEBUG


   # coverage_counts_exon(exons,args.bw)


    logging.basicConfig(level=level, format='%(asctime)s %(levelname)s %(name)s %(message)s')
    output = open(args.output, 'w')
    result = coverage_counts_wg(args.genome_n,args.bw)
    result.to_csv(output, sep='\t')
    output.close()

    print coverage_counts_exon(args.bw)

    command = "Rscript"
    args = [args.output]
    Rcommand = [command, "coverage_summary.R"] + args
    print Rcommand
    x = subprocess.check_output(Rcommand, universal_newlines=True)



if __name__ == '__main__':
    main()
