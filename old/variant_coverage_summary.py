import vcf
import argparse
import json
import pyBigWig
from collections import defaultdict
import numpy as np


def process_coverage(bw_file, vcf_file):
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    bw = pyBigWig.open(bw_file)
    count = 0
    coverage = defaultdict()
    gene = None
    for record in vcf_reader:
        chr = str(record.CHROM)
        chr = chr.replace("MT", "M")
        start = record.POS
        alt = str(record.ALT[0])
        ref = str(record.REF[0])
        if "GENE" in record.INFO:
            gene = record.INFO["GENE"]
        else:
            # could get annotation here with cellbase
            pass
        if chr not in coverage:
            coverage[chr] = defaultdict()
        if start not in coverage[chr]:
            coverage[chr][start] = defaultdict()
        if ref not in coverage[chr][start]:
            coverage[chr][start][ref] = defaultdict()
        if alt not in coverage[chr][start][ref]:
            coverage[chr][start][ref][alt] = defaultdict()

        # vcf is 1 based - -1 from the position when getting coverage
        cov = bw.values("chr" + str(chr), start - 1, start)
        coverage[chr][start][ref][alt]["gene"] = gene
        coverage[chr][start][ref][alt]["coverage"] = cov[0]

    covs = list()
    result = defaultdict()
    result["summary"] = defaultdict()
    result["details"] = defaultdict()
    result["details"]["lt15x"] = defaultdict()
    result["details"]["lt30x"] = defaultdict()
    result["details"]["lt40x"] = defaultdict()
    for chr in coverage:
        for start in coverage[chr]:
            for ref in coverage[chr][start]:
                for alt in coverage[chr][start][ref]:
                    cov = coverage[chr][start][ref][alt]["coverage"]
                    gene = coverage[chr][start][ref][alt]["gene"]
                    id = str(chr) + "." + str(start) + "." + ref + ">" + alt + "." + str(gene)
                    if cov < 15:
                        result["details"]["lt15x"][id] = cov
                    if cov < 30:
                        result["details"]["lt30x"][id] = cov
                    if cov < 40:
                        result["details"]["lt40x"][id] = cov

                    covs.append(cov)

    result["summary"]["std"] = np.std(covs)
    result["summary"]["mean"] = np.mean(covs)
    result["summary"]["max"] = np.max(covs)
    result["summary"]["min"] = np.min(covs)
    result["summary"]["median"] = np.median(covs)
    result["summary"]["pct25"] = np.percentile(covs, 25)
    result["summary"]["pct75"] = np.percentile(covs, 75)
    result["summary"]["n"] = len(covs)
    result["summary"]["lt15x"] = len(filter(lambda x: x < 15, covs))
    result["summary"]["lt30x"] = len(filter(lambda x: x < 30, covs))
    result["summary"]["lt50x"] = len(filter(lambda x: x < 50, covs))

    print(json.dumps(result, indent=4))


def main():
    parser = argparse.ArgumentParser(
        description='Gets coverage stats for sites in a vcf file using a bigWig, output is a summary and a list of sites that don\'t meet 15,30 or 50x')
    parser.add_argument('--vcf', metavar='vcf', type=str,
                        help='vcf file defining sites you would like to know their coverage')
    parser.add_argument('--bw', metavar='bw', default=False,
                        help='bw file')

    args = parser.parse_args()

    process_coverage(args.bw, args.vcf)


if __name__ == '__main__':
    main()
