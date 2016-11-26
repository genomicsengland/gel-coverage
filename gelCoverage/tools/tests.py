import unittest
from gelCoverage.tools.cellbase_helper import CellbaseHelper
from gelCoverage.tools.panelapp_helper import PanelappHelper
import gelCoverage.tools.coverage_stats as coverage_stats
import gelCoverage.tools.sequence_stats as sequence_stats
import random
import collections
import numpy


class CellbaseHelperTests(unittest.TestCase):

    def setUp(self):
        self.species = "hsapiens"
        self.version = "latest"
        self.assembly = "GRCh37"
        self.host = "10.5.8.201:8080/cellbase-4.5.0-rc"
        self.filter_basic_flag = ["basic"]
        self.filter_biotypes = ["IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "IG_V_gene", "protein_coding",
                           "nonsense_mediated_decay", "non_stop_decay", "TR_C_gene",
                           "TR_D_gene", "TR_J_gene", "TR_V_gene"]
        self.cellbase_helper = CellbaseHelper(
            self.species, self.version, self.assembly, self.host,
            self.filter_basic_flag, self.filter_biotypes)

    def test1(self):
        """
        Tests get_all_genes()
        :return:
        """
        genes = self.cellbase_helper.get_all_genes()
        self.assertEqual(type(genes), list)
        self.assertEqual(len(genes), 20760)
        print "%s genes were returned" % str(len(genes))
        print "10 first results: %s..." % ",".join(genes[1:10])

    def test1_1(self):
        """
        Tests get_all_genes() without filtering enabled
        :return:
        """
        genes = self.cellbase_helper.get_all_genes(filter = False)
        self.assertEqual(type(genes), list)
        self.assertEqual(len(genes), 57905)
        print "%s genes were returned" % str(len(genes))
        print "10 first results: %s..." % ",".join(genes[1:10])

    def test2(self):
        """
        Tests make_exons_bed() for a list of 3 genes
        :return:
        """
        gene_list = ["CFTR", "BRCA1", "BRCA2",
                     "IGHE" # gene with 3 transcripts flagged as basic and 1 not flagged as basic
                     ]
        bed = self.cellbase_helper.make_exons_bed(gene_list)
        self.assertIsNotNone(bed)
        ighe_transcripts = []
        for interval in bed:
            chr = interval.chrom
            self.assertEqual(type(chr), unicode)
            gene, txid, exon_idx = interval.name.split("|")
            self.assertEqual(type(gene), unicode)
            self.assertEqual(type(txid), unicode)
            self.assertEqual(type(exon_idx), unicode)
            strand = interval.strand
            self.assertEqual(type(strand), unicode)
            start = int(interval.start)
            self.assertEqual(type(start), int)
            end = int(interval.end)
            self.assertEqual(type(end), int)
            if gene == "BRCA1":
                self.assertEqual(chr, "17")
                self.assertEqual(strand, "-")
                gene_start = 41196312
                gene_end = 41322262
                self.assertGreaterEqual(start, gene_start, "Transcript %s exon %s out of bounds: %s < %s" %
                                   (txid, exon_idx, str(start), str(gene_start)))
                self.assertLessEqual(end, gene_end, "Transcript %s exon %s out of bounds: %s > %s" %
                                 (txid, exon_idx, str(end), str(gene_end)))
            elif gene == "BRCA2":
                self.assertEqual(chr, "13")
                self.assertEqual(strand, "+")
                gene_start = 32889611
                gene_end = 32974403
                self.assertGreaterEqual(start, gene_start, "Transcript %s exon %s out of bounds: %s < %s" %
                                   (txid, exon_idx, str(start), str(gene_start)))
                self.assertLessEqual(end, gene_end, "Transcript %s exon %s out of bounds: %s > %s" %
                                 (txid, exon_idx, str(end), str(gene_end)))
            elif gene == "CFTR":
                self.assertEqual(chr, "7")
                self.assertEqual(strand, "+")
                gene_start = 117105838
                gene_end = 117356025
                self.assertGreaterEqual(start, gene_start, "Transcript %s exon %s out of bounds: %s < %s" %
                                   (txid, exon_idx, str(start), str(gene_start)))
                self.assertLessEqual(end, gene_end, "Transcript %s exon %s out of bounds: %s > %s" %
                                 (txid, exon_idx, str(end), str(gene_end)))
            elif gene == "IGHE":
                ighe_transcripts.append(txid)
            else:
                self.assertTrue(False, "Unexpected information in the BED file for gene %s" % gene)
        self.assertEqual(len(set(ighe_transcripts)), 2) # checks that non basic flagged transcript has been filtered out

    def test2_1(self):
        """
        Tests make_exons_bed() for a list of 3 genes not filtering
        :return:
        """
        gene_list = ["CFTR", "BRCA1", "BRCA2",
                     "IGHE"  # gene with 3 transcripts flagged as basic and 1 not flagged as basic
                     ]
        bed = self.cellbase_helper.make_exons_bed(gene_list, filter = False)
        self.assertIsNotNone(bed)
        ighe_transcripts = []
        for interval in bed:
            chr = interval.chrom
            self.assertEqual(type(chr), unicode)
            gene, txid, exon_idx = interval.name.split("|")
            self.assertEqual(type(gene), unicode)
            self.assertEqual(type(txid), unicode)
            self.assertEqual(type(exon_idx), unicode)
            strand = interval.strand
            self.assertEqual(type(strand), unicode)
            start = int(interval.start)
            self.assertEqual(type(start), int)
            end = int(interval.end)
            self.assertEqual(type(end), int)
            if gene == "BRCA1":
                self.assertEqual(chr, "17")
                self.assertEqual(strand, "-")
                gene_start = 41196312
                gene_end = 41322262
                self.assertGreaterEqual(start, gene_start, "Transcript %s exon %s out of bounds: %s < %s" %
                                   (txid, exon_idx, str(start), str(gene_start)))
                self.assertLessEqual(end, gene_end, "Transcript %s exon %s out of bounds: %s > %s" %
                                 (txid, exon_idx, str(end), str(gene_end)))
            elif gene == "BRCA2":
                self.assertEqual(chr, "13")
                self.assertEqual(strand, "+")
                gene_start = 32889611
                gene_end = 32974403
                self.assertGreaterEqual(start, gene_start, "Transcript %s exon %s out of bounds: %s < %s" %
                                   (txid, exon_idx, str(start), str(gene_start)))
                self.assertLessEqual(end, gene_end, "Transcript %s exon %s out of bounds: %s > %s" %
                                 (txid, exon_idx, str(end), str(gene_end)))
            elif gene == "CFTR":
                self.assertEqual(chr, "7")
                self.assertEqual(strand, "+")
                gene_start = 117105838
                gene_end = 117356025
                self.assertGreaterEqual(start, gene_start, "Transcript %s exon %s out of bounds: %s < %s" %
                                   (txid, exon_idx, str(start), str(gene_start)))
                self.assertLessEqual(end, gene_end, "Transcript %s exon %s out of bounds: %s > %s" %
                                 (txid, exon_idx, str(end), str(gene_end)))
            elif gene == "IGHE":
                ighe_transcripts.append(txid)
            else:
                self.assertTrue(False, "Unexpected information in the BED file for gene %s" % gene)
        self.assertEqual(len(set(ighe_transcripts)), 3)  # checks that non basic flagged transcript has not been filtered out

    def test3(self):
        """
        Tests make_exons_bed() with an empty gene list
        :return:
        """
        gene_list = []
        try:
            self.cellbase_helper.make_exons_bed(gene_list)
            self.assertTrue(False, "Function did not fail with an empty list")
        except SystemError:
            pass

    def test4(self):
        """
        Tests make_exons_bed() with a None gene list
        :return:
        """
        gene_list = None
        try:
            self.cellbase_helper.make_exons_bed(gene_list)
            self.assertTrue(False, "Function did not fail with a None list")
        except SystemError:
            pass

    def test5(self):
        """
        Tests make_exons_bed() with a non existing gene in the list
        :return:
        """
        # TODO: this test might need fixing if we check for the existence of genes
        gene_list = ["BRCA1", "BRCA2", "non_existing_gene"]
        bed = self.cellbase_helper.make_exons_bed(gene_list)
        self.assertIsNotNone(bed)
        for interval in bed:
            chr = interval.chrom
            self.assertEqual(type(chr), unicode)
            gene, txid, exon_idx = interval.name.split("|")
            if gene not in ["BRCA1", "BRCA2"]:
                self.assertTrue(False, "BED contains an unexpected gene %s" % gene)

    def test6(self):
        """
        Tests make_exons_bed() with a list that should return no results
        :return:
        """
        # TODO: this test might need fixing if we check for the existence of genes
        gene_list = ["non_existing_gene"]
        bed = self.cellbase_helper.make_exons_bed(gene_list)
        self.assertIsNotNone(bed)


class PanelappHelperTests(unittest.TestCase):

    def setUp(self):
        host = "bioinfo.extge.co.uk/crowdsourcing/WebServices"
        self.panelapp_helper = PanelappHelper(host)

    def test1(self):
        """
        Tests querying of an existing panel filtered by HighEvidence
        :return:
        """
        panel_name = "Epileptic encephalopathy"
        panel_version = "0.2"
        gene_confidence_threshold = "HighEvidence"
        gene_list = self.panelapp_helper.get_gene_list(panel = panel_name,
                                           panel_version = panel_version,
                                           gene_confidence_threshold = gene_confidence_threshold)
        self.assertEqual(len(gene_list), 7)
        print "%s genes were returned" % str(len(gene_list))
        print "10 first results: %s..." % ",".join(gene_list[1:min(len(gene_list), 10)])

    def test2(self):
        """
        Tests querying of an existing panel filtered by LowEvidence
        :return:
        """
        panel_name = "Epileptic encephalopathy"
        panel_version = "0.2"
        gene_confidence_threshold = "LowEvidence"
        gene_list = self.panelapp_helper.get_gene_list(panel = panel_name,
                                           panel_version = panel_version,
                                           gene_confidence_threshold = gene_confidence_threshold)
        self.assertEqual(len(gene_list), 117)
        print "%s genes were returned" % str(len(gene_list))
        print "10 first results: %s..." % ",".join(gene_list[1:min(len(gene_list), 10)])

    def test3(self):
        """
        Tests querying of an existing panel filtered by LowEvidence
        :return:
        """
        panel_name = "Epileptic encephalopathy"
        panel_version = "0.2"
        gene_confidence_threshold = ["LowEvidence", "HighEvidence"]
        gene_list = self.panelapp_helper.get_gene_list(panel = panel_name,
                                           panel_version = panel_version,
                                           gene_confidence_threshold = gene_confidence_threshold)
        self.assertEqual(len(gene_list), 124)
        print "%s genes were returned" % str(len(gene_list))
        print "10 first results: %s..." % ",".join(gene_list[1:min(len(gene_list), 10)])

    def test3(self):
        """
        Tests querying of an unexisting panel
        :return:
        """
        panel_name = "unexisting_disease"
        panel_version = "0.2"
        gene_confidence_threshold = ["LowEvidence", "HighEvidence"]
        try:
            gene_list = self.panelapp_helper.get_gene_list(panel = panel_name,
                                               panel_version = panel_version,
                                               gene_confidence_threshold = gene_confidence_threshold)
            self.assertTrue(False, "Function did not fail with an unexisting panel")
        except SystemError:
            pass


class CoverageStatsTests(unittest.TestCase):

    def setUp(self):
        # Creates 1000 random coverage values between 0 and 400
        self.coverages = [random.randint(a=30, b=400) for p in range(0,1000)]
        self.coverages[100:150] = [random.randint(a=0, b=29) for p in range(0,50)]
        self.coverages[300:310] = [random.randint(a=0, b=29) for p in range(0, 10)]
        self.coverages[340:345] = [random.randint(a=0, b=29) for p in range(0, 5)]

        self.coverages2 = [random.randint(a=0, b=400) for p in range(0, 2000)]
        self.coverages3 = [random.randint(a=0, b=400) for p in range(0, 3000)]

        self.start_position = 1234
        self.coverage_threshold = 30
        self.gc_content = 0.59

    def test1(self):
        """
        Tests find_gaps(coverages, start_position, coverage_threshold)
        :return:
        """
        gaps = coverage_stats.find_gaps(self.coverages, self.start_position, self.coverage_threshold)
        self.assertEqual(type(gaps), list)
        self.assertEqual(len(gaps), 3)
        for gap in gaps:
            self.assertEqual(type(gap), collections.defaultdict)
            self.assertTrue("start" in gap)
            self.assertTrue("end" in gap)
            self.assertEqual(type(gap["start"]), int)
            self.assertEqual(type(gap["end"]), int)
            print "Found a gap at %s-%s" % (str(gap["start"]), str(gap["end"]))

    def test2(self):
        """
        compute_exon_level_statistics(coverages, gc_content)
        :return:
        """
        stats = coverage_stats.compute_exon_level_statistics(self.coverages, self.gc_content)
        self.assertEqual(type(stats), collections.defaultdict)
        self.assertTrue("total_bases" in stats)
        self.assertTrue("mean" in stats)
        self.assertTrue("median" in stats)
        self.assertTrue("pct75" in stats)
        self.assertTrue("pct25" in stats)
        self.assertTrue("bases_lt_3x" in stats)
        self.assertTrue("bases_lt_15x" in stats)
        self.assertTrue("bases_gte_15x" in stats)
        self.assertTrue("bases_gte_30x" in stats)
        self.assertTrue("bases_gte_50x" in stats)
        self.assertTrue("percent_lt_15x" in stats)
        self.assertTrue("percent_gte_15x" in stats)
        self.assertTrue("percent_gte_30x" in stats)
        self.assertTrue("percent_gte_50x" in stats)
        self.assertTrue("gc_content" in stats)
        self.assertEqual(type(stats["mean"]), numpy.float64)
        self.assertEqual(type(stats["median"]), numpy.float64)
        self.assertEqual(type(stats["pct75"]), numpy.float64)
        self.assertEqual(type(stats["total_bases"]), int)
        self.assertEqual(type(stats["pct25"]), numpy.float64)
        self.assertEqual(type(stats["bases_lt_3x"]), int)
        self.assertEqual(type(stats["bases_lt_15x"]), int)
        self.assertEqual(type(stats["bases_gte_15x"]), int)
        self.assertEqual(type(stats["bases_gte_30x"]), int)
        self.assertEqual(type(stats["bases_gte_50x"]), int)
        self.assertEqual(type(stats["percent_lt_15x"]), float)
        self.assertEqual(type(stats["percent_gte_15x"]), float)
        self.assertEqual(type(stats["percent_gte_30x"]), float)
        self.assertEqual(type(stats["percent_gte_50x"]), float)
        self.assertEqual(type(stats["gc_content"]), float)

    def test3(self):
        """
        compute_transcript_level_statistics(exons)
        :return:
        """
        exon1 = coverage_stats.compute_exon_level_statistics(self.coverages, self.gc_content)
        exon2 = coverage_stats.compute_exon_level_statistics(self.coverages2, self.gc_content)
        exon3 = coverage_stats.compute_exon_level_statistics(self.coverages3, self.gc_content)
        exons = {1:exon1, 2:exon2, 3:exon3}
        stats = coverage_stats.compute_transcript_level_statistics(exons)
        self.assertEqual(type(stats), collections.defaultdict)
        self.assertTrue("total_bases" in stats)
        self.assertTrue("mean" in stats)
        self.assertTrue("weighted_median" in stats)
        self.assertTrue("weighted_pct75" in stats)
        self.assertTrue("weighted_pct25" in stats)
        self.assertTrue("bases_lt_3x" in stats)
        self.assertTrue("bases_lt_15x" in stats)
        self.assertTrue("bases_gte_15x" in stats)
        self.assertTrue("bases_gte_30x" in stats)
        self.assertTrue("bases_gte_50x" in stats)
        self.assertTrue("percent_lt_15x" in stats)
        self.assertTrue("percent_gte_15x" in stats)
        self.assertTrue("percent_gte_30x" in stats)
        self.assertTrue("percent_gte_50x" in stats)
        self.assertTrue("gc_content" in stats)
        self.assertEqual(type(stats["mean"]), numpy.float64)
        self.assertEqual(type(stats["weighted_median"]), numpy.float64)
        self.assertEqual(type(stats["weighted_pct75"]), numpy.float64)
        self.assertEqual(type(stats["total_bases"]), numpy.int64)
        self.assertEqual(type(stats["weighted_pct25"]), numpy.float64)
        self.assertEqual(type(stats["bases_lt_3x"]), numpy.int64)
        self.assertEqual(type(stats["bases_lt_15x"]), numpy.int64)
        self.assertEqual(type(stats["bases_gte_15x"]), numpy.int64)
        self.assertEqual(type(stats["bases_gte_30x"]), numpy.int64)
        self.assertEqual(type(stats["bases_gte_50x"]), numpy.int64)
        self.assertEqual(type(stats["percent_lt_15x"]), float)
        self.assertEqual(type(stats["percent_gte_15x"]), float)
        self.assertEqual(type(stats["percent_gte_30x"]), float)
        self.assertEqual(type(stats["percent_gte_50x"]), float)
        self.assertEqual(type(stats["gc_content"]), numpy.float64)


class SequenceStatsTests(unittest.TestCase):

    def setUp(self):
        # Creates random sequence with a GC content close 0.6
        self.sequence = numpy.random.choice(["G", "C", "A", "T"], 10000, p=[0.3, 0.3, 0.2, 0.2])

    def test1(self):
        """
        Tests find_gaps(coverages, start_position, coverage_threshold)
        :return:
        """
        gc_content = sequence_stats.compute_gc_content(self.sequence)
        print "Found a GC content of %s" % gc_content
        self.assertEqual(type(gc_content), float)
        self.assertTrue(gc_content <= 0.65 and gc_content >= 0.55)

if __name__ == '__main__':
    unittest.main()
