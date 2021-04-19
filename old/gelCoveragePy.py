import pybedtools
from collections import defaultdict
import tqdm
import numpy as np
import pyBigWig
import pandas
import json
import re
import urllib2


class Coverage:
    def __init__(self):
        pass

    def intervals(self, n_file, other_chr):
        """
        this is a method to make up for the lack of non-autosomal chromosomes in the nonN's bed file, fill in chrX
        and chrY in the bed file for calculating wg coverage

        :param n_file: nonN regions of the genome
        :param other_chr: all chromosome in bw file
        :return: bed file with all regions
        """
        ns = pybedtools.BedTool(n_file)
        ns_chrs = list()
        for feature in ns:
            ns_chrs.append(feature.chrom)
        all = pybedtools.BedTool(other_chr)

        missing = list()
        for feature in all:
            all_chr = feature.chrom
            if all_chr.replace("chr", "") not in ns_chrs:
                missing.append((feature.chrom.replace("chr", ""), feature.start, feature.end))

        to_cat = pybedtools.BedTool(missing)
        all = ns.cat(to_cat, postmerge=False)

        return all

    @staticmethod
    def lt_helper(a):
        """
        helper to sort chromosmes properly

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

    @staticmethod
    def __lt__(a, b):
        """
        run the chromosome sort helper

        :param a:
        :param b:
        :return: proper sorted chromosomes
        """
        return cmp(Coverage.lt_helper(a), Coverage.lt_helper(b))

    def get_chr_lengths(self, bw, format):
        """
        get chromosome lengths from header of bigWig file

        :param bw: pyBigWig file object
        :param format: specify "dict" to return a dict instead of a bed recognisable format
        :return: list of chromosomes and length in a bed recognisable format
        """
        chromosomes = list()
        chromosomes_dict = defaultdict()
        for i in range(1, 25):
            chromosome = "chr" + str(i)
            if i == 23:
                chromosome = "chrX"
            elif i == 24:
                chromosome = "chrY"
            length = int(bw.chroms(chromosome))
            chromosomes.append((chromosome, 1, length))
            chromosomes_dict[chromosome] = length
        if format == "dict":
            return chromosomes_dict
        else:
            return chromosomes

    def coverage_counts_wg(self, n_file, bw):
        """
        get coverage summary counts in non-N regions of the genome

        :param n_file: non-N regions of the genome bed file
        :param bw: path to bigWig file
        :return: dataframe of exon coverage summary
        """
        bw = pyBigWig.open(bw)
        other = self.get_chr_lengths(bw, "bed")
        ns = self.intervals(n_file, other)
        return self.generic_coverage(bw, ns)

    def coverage_counts_exon(self, bw, bed):
        """
        get coverage summary counts in exonic regions

        :param bw: path to a bigWig file
        :param bed: bed file object of regions
        :return: dataframe of exon coverage summary
        """
        bw = pyBigWig.open(bw)
        exons = self.flatten_bed(bed)
        return self.generic_coverage(bw, exons)

    def exon_gc_coverage(self, bw, bed, prefix='chr'):
        """
        this is a specific method to calculate coverage in all exonic regions with gc content

        :param prefix: this prefix will be appended to the
        :param bw: path to bigWig file
        :param bed: bed file object of regions
        :return: list of results (NOT a bed file compatible object)
        """
        bw = pyBigWig.open(bw)

        total_exons = bed.count()
        exon_count = 0
        temp_store = list()
        headers = ["chrm", "start", "end", "id", "exon", "gc", "strand", "cov"]
        for interval in bed:
            chrom = interval.chrom
            start = interval.start
            end = interval.end
            if start == end:
                end += 1

            if prefix + str(chrom) in self.get_chr_lengths(bw, "dict"):
                cov = bw.values(prefix + str(chrom), start, end)
                cov_mean = np.mean(cov)
                id = interval.name
                gene, txid, exon_raw = id.split("|")
                exon = exon_raw.replace("exon", "")
                temp_store.append([chrom, start, end, interval.name, exon, interval.score, interval.strand, cov_mean])

            exon_count += 1
            if exon_count % 100000 == 0:
                print("Done " + str(exon_count) + " exons out of " + str(total_exons))
        return pandas.DataFrame(temp_store, columns=headers)

    def generic_coverage(self, bw, bed):
        """
        using pyBigWig calculates number of bases at each coverage level in regions of a bed file

        :param bw: pyBigWig bw file object
        :param bed: bed file of regions to calculate the summary over
        :return: sorted_data: pandas dataframe of coverage summary - rows coverage, columns chromosomes
        """
        result = defaultdict()

        for i in self.get_chr_lengths(bw, "dict"):
            result[i.replace("chr", "")] = defaultdict(lambda: 0)

        for interval in bed:
            chrom = interval.chrom
            start = interval.start
            end = interval.end
            if start == end:
                end += 1
            if "chr" + str(chrom) in self.get_chr_lengths(bw, "dict"):
                cov = bw.values("chr" + str(chrom), start, end)
                for coverage in cov:
                    result[chrom][int(coverage)] += 1

        order = sorted(result, cmp=self.__lt__)
        data = pandas.DataFrame.from_dict(result, orient="columns").fillna(0)
        sorted_data = data[order]

        sorted_data["allchrs"] = sorted_data.sum(axis=1)
        col_list = list(sorted_data)
        col_list.remove("X")
        col_list.remove("Y")
        col_list.remove("allchrs")
        sorted_data["autosomes"] = sorted_data[col_list].sum(axis=1)
        sorted_data.index.name = 'Coverage'
        return sorted_data

    def gc_content(self, sequence):
        """
        Calculates gc content of a DNA sequence

        :param sequence: a string of ATGC's
        :return: the gc fraction of the region
        """
        gcCount = 0
        totalBaseCount = 0
        gcCount += len(re.findall("[GC]", sequence))
        totalBaseCount += len(re.findall("[GCTA]", sequence))
        gcFraction = float(gcCount) / totalBaseCount
        return (gcFraction)

    def make_exons_bed(self, assembly='GRCh37'):
        """
        Gets all exons from cellbase and makes a bed - also calculates gc content, returns a valid bed with gc in the
        score column

        :return: pybedtools object
        """
        url = 'http://10.5.8.201:8080/cellbase-4.5.0-rc/webservices/rest/latest/hsapiens/feature/gene/search?limit=1' \
              '&include=name&assembly={assembly}'.format(assembly=assembly)

        numTotalResults = json.load(urllib2.urlopen(url))["response"][0]["numTotalResults"]
        # numTotalResults=100
        all_exons = list()

        # TODO while doing this calculate coding exons too

        print("making exon bed file... <3min")
        for i in tqdm.tqdm(xrange(0, numTotalResults, 5000)):
            url = "http://10.5.8.201:8080/cellbase-4.5.0-rc/webservices/rest/latest/hsapiens/feature/gene/search" \
                  "?include=name,chromosome,transcripts.exons.start," \
                  "transcripts.exons.exonNumber," \
                  "transcripts.id,transcripts.strand," \
                  "transcripts.exons.end,transcripts.exons.sequence," \
                  "exonNumber&limit=5000&skip={i}&assembly={assembly}".format(i=str(i), assembly=assembly)

            exons = json.load(urllib2.urlopen(url))["response"][0]["result"]
            for result in exons:
                gene = result["name"]
                for transcript in result["transcripts"]:
                    txid = transcript["id"]
                    strand = transcript["strand"]
                    for exons in transcript["exons"]:
                        if "sequence" in exons:
                            gc = self.gc_content(exons["sequence"])
                            exonNumber = exons["exonNumber"]
                            row_id = gene + "|" + txid + "|exon" + str(exonNumber)
                            all_exons.append(
                                (result["chromosome"], exons["start"], exons["end"], row_id, str(gc), strand))
        bed = pybedtools.BedTool(all_exons)
        return bed

    def flatten_bed(self, bed):
        """
        flatten a bed file, i.e. collapse overlapping regions

        :param bed: bed file object
        :return: bed: sorted bed
        """
        bed = bed.sort()
        bed = bed.merge()
        return bed
