import unittest
from cellbase_helper import CellbaseHelper
from panelapp_helper import PanelappHelper


class CellbaseHelperTests(unittest.TestCase):

    def setUp(self):
        self.species = "hsapiens"
        self.version = "latest"
        self.assembly = "GRCh37"
        self.host = "10.5.8.201:8080/cellbase-4.5.0-rc"
        self.filter_basic_flag = 1
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
        gene_list = ["CFTR", "BRCA1", "BRCA2"]
        bed = self.cellbase_helper.make_exons_bed(gene_list)
        self.assertIsNotNone(bed)
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
            else:
                self.assertTrue(False, "Unexpected information in the BED file for gene %s" % gene)

    def test2_1(self):
        """
        Tests make_exons_bed() for a list of 3 genes not filtering
        :return:
        """
        gene_list = ["CFTR", "BRCA1", "BRCA2"]
        bed = self.cellbase_helper.make_exons_bed(gene_list, filter = False)
        self.assertIsNotNone(bed)
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
            else:
                self.assertTrue(False, "Unexpected information in the BED file for gene %s" % gene)

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


if __name__ == '__main__':
    unittest.main()
