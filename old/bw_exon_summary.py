import argparse
import json
import pyBigWig
from collections import defaultdict
from itertools import chain


def calculate_percent(count, total):
    total = int(total)
    if total == 0:
        perc = "ERROR TOTAL=" + str(total)
    else:
        perc = (count / float(total)) * 100

    if perc > 100:
        return "ERROR_OVER_100"
    else:
        return perc


def collapse_transcripts(collapse, exons, gene, chr):
    if gene not in collapse:
        collapse[gene] = {}
    if "all_exons" not in collapse[gene]:
        collapse[gene]["all_exons"] = list()
    if "result" not in collapse[gene]:
        collapse[gene]["result"] = list()
    if "chr" not in collapse[gene]:
        collapse[gene]["chr"] = chr

    collapse[gene]["all_exons"].append(exons)
    for tx_exons in collapse[gene]["all_exons"]:
        for exons in tx_exons:
            for exon_rank in exons:
                start_stop = exons[str(exon_rank)]
                start, stop = start_stop.split("-")
                exon_range = [int(start), int(stop)]
                collapse[gene]["result"].append(exon_range)
        collapse[gene]["all_exons"].remove(tx_exons)

    return collapse


flatten = chain.from_iterable

LEFT, RIGHT = 1, -1


def join_ranges(data, offset=0):
    data = sorted(flatten(((start, LEFT), (stop + offset, RIGHT))
                          for start, stop in data))
    c = 0
    for value, label in data:
        if c == 0:
            x = value
        c += label
        if c == 0:
            yield x, value - offset


parser = argparse.ArgumentParser(description='Coverage summary for exons')
parser.add_argument('--bw', metavar='bw', help='This is the bw file')
parser.add_argument('--chr2gene', metavar='chr2gene', help='This is the chr2gene file')
parser.add_argument('--genes', metavar='genes', help='This is a list of genes or you can have a list of transcripts')
parser.add_argument('--transcripts', metavar='transcripts',
                    help='This is a list of transcripts or you can have a list of genes')
parser.add_argument('--genes2tx', metavar='genes2tx', help='This is a file listing genes to transcripts')
parser.add_argument('--collapse', metavar='flatten',
                    help='collapse transcripts to one cannonical transcript - only used when using genes as input',
                    default=0)
parser.add_argument('--cnv', metavar='cnv', help='cnv vcf - so that losses can be indicated', default=0)

args = parser.parse_args()

with open(args.chr2gene) as data_file:
    chr2gene = json.load(data_file)

bw = pyBigWig.open(args.bw)

# get transcripts for genes
transcripts = list()
genes = list()
if args.genes is not None:
    transcripts = list()
    genesx = defaultdict(list)
    with open(args.genes2tx, 'r') as handle:
        lines = handle.readlines()
        for line in lines:
            eGeneId, eTxID, gene = line.split("\t")
            genesx[gene.rstrip()].append(eTxID)
    with open(args.genes, 'r') as handle:
        lines = handle.readlines()
        for line in lines:
            for transcript in genesx[line.rstrip()]:
                    transcripts.append(transcript)
            if args.collapse != 0:
                genes.append(line.rstrip())

if args.transcripts is not None:
    with open(args.transcripts) as f:
        transcripts = f.read().splitlines()


if args.collapse != 0:
    queries = genes
    collapse = defaultdict()
    genes_collapse = list()
    for transcript in transcripts:
        exons = chr2gene[transcript]["exons"]["coding_exon_regions"]
        eGeneID = chr2gene[transcript]["eGeneID"]
        gene = chr2gene[transcript]["gene"]
        chr = chr2gene[transcript]["chr"]
        genes_collapse.append(gene)
        collapse = collapse_transcripts(collapse, exons, gene, chr)

    chr2gene=defaultdict()
    for gene in genes_collapse:
        final = list(join_ranges(collapse[gene]["result"]))
        chr2gene[gene] = defaultdict()[gene] = defaultdict()
        chr2gene[gene]["exons"] = defaultdict()
        chr2gene[gene]["gene"] = gene
        chr2gene[gene]["eGeneID"] = gene
        chr2gene[gene]["strand"] = "NA"
        chr2gene[gene]["chr"] = collapse[gene]["chr"]
        chr2gene[gene]["exons"]["flattened"] = list()
        for count, exon_range in enumerate(final):
            exon = "-".join(str(x) for x in exon_range)
            exon_dict = defaultdict()
            exon_dict[count] = exon
            chr2gene[gene]["exons"]["flattened"].append(exon_dict)

    scope = "flattened"
else:
    queries = transcripts
    scope = "coding_exon_regions"

exon_counts = list()
for query in queries:
    if query in chr2gene:
        exons = chr2gene[query]["exons"][scope]
        exon_counts.append(len(exons))

max = max(exon_counts)



header = ["gene", "txID", "chr", "strand", "totalBases", "bases_lt_3x", "bases_lt_15x",
          "bases_gte_15x", "bases_gte_30x", "bases_gte_50x","percent_lt_15x","percent_gte_15x",
          "percent_gte_30x","percent_gte_50x"]

for head in range(1, max + 1):
    header.append("Exon" + str(head))

print("\t".join(header))


if args.collapse !=0:
        max-=1

for query in queries:
    exon_output = defaultdict()
    for head in range(1, max + 1):
        exon_output[int(head)] = "-"
    if query in chr2gene:
        chr = chr2gene[query]["chr"]
        strand = chr2gene[query]["strand"]
        gene = chr2gene[query]["gene"]
        eGeneID = chr2gene[query]["eGeneID"]
        chr = "chr" + str(chr)
        exons = chr2gene[query]["exons"][scope]
        exon_count = len(exons)
        total_bases_lt_3x = list()
        total_bases_lt_15x = list()
        total_bases_gte_15x = list()
        total_bases_gte_30x = list()
        total_bases_gte_50x = list()
        total_bases = 0
        exon_coverage = list()
        if int(exon_count) > 0:
            for exon_dict in exons:
                for exon_number in exon_dict:
                    start, end = exon_dict[exon_number].split("-")
                    exon_total_bases = int(end) - int(start)
                    total_bases += exon_total_bases
                    try:
                        result = bw.values(str(chr), int(start), int(end))
                    except:
                        pass

                    exon_bases_lt_3x = list()
                    exon_bases_lt_15x = list()
                    exon_bases_gte_15x = list()
                    exon_bases_gte_30x = list()
                    exon_bases_gte_50x = list()

                    for coverage in result:
                        coverage = int(coverage)
                        if coverage < 3:
                            exon_bases_lt_3x.append(coverage)
                            total_bases_lt_3x.append(coverage)
                        if coverage < 15:
                            exon_bases_lt_15x.append(coverage)
                            total_bases_lt_15x.append(coverage)
                        if coverage >= 15:
                            exon_bases_gte_15x.append(coverage)
                            total_bases_gte_15x.append(coverage)
                        if coverage >= 30:
                            exon_bases_gte_30x.append(coverage)
                            total_bases_gte_30x.append(coverage)
                        if coverage >= 50:
                            exon_bases_gte_50x.append(coverage)
                            total_bases_gte_50x.append(coverage)

                exon_perc_lt_15x = calculate_percent(len(exon_bases_gte_15x), exon_total_bases)
                exon_perc_gte_15x = calculate_percent(len(exon_bases_gte_15x), exon_total_bases)
                exon_perc_gte_30x = calculate_percent(len(exon_bases_gte_30x), exon_total_bases)
                exon_perc_gte_50x = calculate_percent(len(exon_bases_gte_50x), exon_total_bases)
                exon_output[int(exon_number)] = str(len(exon_bases_lt_3x)) + "/" + str(exon_total_bases)
                # exon_coverage.append(str(len(exon_bases_lt_3x))+"/"+str(exon_total_bases))

            total_perc_lt_15x = calculate_percent(len(total_bases_lt_15x), total_bases)
            total_perc_gte_15x = calculate_percent(len(total_bases_gte_15x), total_bases)
            total_perc_gte_30x = calculate_percent(len(total_bases_gte_30x), total_bases)
            total_perc_gte_50x = calculate_percent(len(total_bases_gte_50x), total_bases)
        else:
            pass


        out = [gene, query, chr, strand, total_bases, len(total_bases_lt_3x), len(total_bases_lt_15x),
               len(total_bases_gte_15x), len(total_bases_gte_30x), len(total_bases_gte_50x),
               total_perc_lt_15x,total_perc_gte_15x,total_perc_gte_30x,total_perc_gte_50x,
               "\t".join(str(value) for exon, value in exon_output.items())]
        if total_bases > 0:
            print("\t".join(str(x) for x in out))
