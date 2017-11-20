import unittest
import logging
import requests
import urllib2
from gelcoverage.tools.cellbase_helper import CellbaseHelper
from gelcoverage.tools.panelapp_helper import PanelappHelper
import gelcoverage.tools.backoff_retrier as backoff_retrier


PANELAPP_HOST = "bio-test-panelapp.gel.zone/WebServices"
ASSEMBLY = "GRCh37"
SPECIES = "hsapiens"
CELLBASE_VERSION = "latest"
CELLBASE_HOST = "bio-test-cellbase-tomcat-01.gel.zone:8080/cellbase"
FILTER_BASIC_FLAG = ["basic"]
FILTER_BIOTYPES = ["IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "IG_V_gene", "protein_coding",
                   "nonsense_mediated_decay", "non_stop_decay", "TR_C_gene",
                   "TR_D_gene", "TR_J_gene", "TR_V_gene"]


class CellbaseHelperTests(unittest.TestCase):

    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)
        self.retries = -1
        self.cellbase_helper = CellbaseHelper(
            SPECIES, CELLBASE_VERSION, ASSEMBLY, CELLBASE_HOST, self.retries,
            FILTER_BASIC_FLAG, FILTER_BIOTYPES
        )

    def test1(self):
        """
        Tests get_all_genes()
        :return:
        """
        genes = self.cellbase_helper.get_all_gene_names()
        self.assertEqual(type(genes), list)
        self.assertEqual(len(genes), 20760)
        print "%s genes were returned" % str(len(genes))
        print "10 first results: %s..." % ",".join(genes[1:10])

    def test1_1(self):
        """
        Tests get_all_genes() without filtering enabled
        :return:
        """
        genes = self.cellbase_helper.get_all_gene_names(_filter=False)
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
        #bed.saveas("/home/priesgo/test.bed")

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
        bed = self.cellbase_helper.make_exons_bed(gene_list, _filter=False)
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

    def test7(self):
        """
        Test panel with no transcript passing filters.
        :return:
        """
        "5763f2ea8f620350a1996048", "1.0", "HighEvidence"


class CellbaseHelperRetriesTests(unittest.TestCase):

    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)
        self.retries = -1
        self.cellbase_helper = CellbaseHelper(
            SPECIES, CELLBASE_VERSION, ASSEMBLY, CELLBASE_HOST, self.retries,
            FILTER_BASIC_FLAG, FILTER_BIOTYPES
        )
        # overwrites wrapped CB search
        self.count_failures = 0
        def simulate_connection_failures(*args, **kwargs):
            if self.count_failures < 3:
                self.count_failures += 1
                raise requests.exceptions.RequestException("Faked exception")
            return self.cellbase_helper.cellbase_gene_client.search(*args, **kwargs)

        self.cellbase_helper.cellbase_search = backoff_retrier.wrapper(simulate_connection_failures, self.retries)

    def test1(self):
        """
        Tests get_all_genes()
        :return:
        """
        genes = self.cellbase_helper.get_all_gene_names()
        self.assertEqual(type(genes), list)
        self.assertEqual(len(genes), 20760)
        self.assertEqual(self.count_failures, 3)
        print "%s genes were returned" % str(len(genes))
        print "10 first results: %s..." % ",".join(genes[1:10])

    def test1_1(self):
        """
        Tests get_all_genes() without filtering enabled
        :return:
        """
        genes = self.cellbase_helper.get_all_gene_names(_filter=False)
        self.assertEqual(type(genes), list)
        self.assertEqual(len(genes), 57905)
        self.assertEqual(self.count_failures, 3)
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
        #bed.saveas("/home/priesgo/test.bed")
        self.assertEqual(self.count_failures, 3)

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
        bed = self.cellbase_helper.make_exons_bed(gene_list, _filter=False)
        self.assertIsNotNone(bed)
        self.assertEqual(self.count_failures, 3)
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
            self.assertEqual(self.count_failures, 3)
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
            self.assertEqual(self.count_failures, 3)
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
        self.assertEqual(self.count_failures, 3)
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
        self.assertEqual(self.count_failures, 3)

    def test7(self):
        """
        Test panel with no transcript passing filters.
        :return:
        """
        "5763f2ea8f620350a1996048", "1.0", "HighEvidence"


class CellbaseHelperRetriesTests2(unittest.TestCase):

    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)
        self.retries = -1
        self.cellbase_helper = CellbaseHelper(
            SPECIES, CELLBASE_VERSION, ASSEMBLY, CELLBASE_HOST, self.retries,
            FILTER_BASIC_FLAG, FILTER_BIOTYPES
        )
        # overwrites wrapped CB search
        self.count_failures = 0
        def simulate_connection_failures(*args, **kwargs):
            if self.count_failures < 3:
                self.count_failures += 1
                raise requests.exceptions.ConnectionError("Faked exception")
            return self.cellbase_helper.cellbase_gene_client.search(*args, **kwargs)

        self.cellbase_helper.cellbase_search = backoff_retrier.wrapper(simulate_connection_failures, self.retries)

    def test1(self):
        """
        Tests get_all_genes()
        :return:
        """
        genes = self.cellbase_helper.get_all_gene_names()
        self.assertEqual(type(genes), list)
        self.assertEqual(len(genes), 20760)
        self.assertEqual(self.count_failures, 3)
        print "%s genes were returned" % str(len(genes))
        print "10 first results: %s..." % ",".join(genes[1:10])


class CellbaseHelperRetriesTests3(unittest.TestCase):

    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)
        self.retries = -1
        self.cellbase_helper = CellbaseHelper(
            SPECIES, CELLBASE_VERSION, ASSEMBLY, CELLBASE_HOST, self.retries,
            FILTER_BASIC_FLAG, FILTER_BIOTYPES
        )
        # overwrites wrapped CB search
        self.count_failures = 0
        def simulate_connection_failures(*args, **kwargs):
            if self.count_failures < 3:
                self.count_failures += 1
                raise ValueError("Faked exception")
            return self.cellbase_helper.cellbase_gene_client.search(*args, **kwargs)

        self.cellbase_helper.cellbase_search = backoff_retrier.wrapper(simulate_connection_failures, self.retries)

    def test1(self):
        """
        Tests get_all_genes()
        :return:
        """
        caught = False
        try:
            genes = self.cellbase_helper.get_all_gene_names()
            self.assertTrue(False)
        except ValueError:
            caught = True
            self.assertEqual(self.count_failures, 1)
        self.assertTrue(caught)

class CellbaseHelperRetriesTests4(unittest.TestCase):

    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)
        self.retries = 3
        self.cellbase_helper = CellbaseHelper(
            SPECIES, CELLBASE_VERSION, ASSEMBLY, CELLBASE_HOST, self.retries,
            FILTER_BASIC_FLAG, FILTER_BIOTYPES
        )
        # overwrites wrapped CB search
        self.count_failures = 0
        def simulate_connection_failures(*args, **kwargs):
            self.count_failures += 1
            raise requests.exceptions.ConnectionError("Faked exception")

        self.cellbase_helper.cellbase_search = backoff_retrier.wrapper(simulate_connection_failures, self.retries)

    def test1(self):
        """
        Tests get_all_genes()
        :return:
        """
        try:
            genes = self.cellbase_helper.get_all_gene_names()
            self.assertTrue(False)
        except requests.exceptions.ConnectionError:
            caught = True
            self.assertEqual(self.count_failures, 4)
        self.assertTrue(caught)


class PanelappHelperTests(unittest.TestCase):

    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)
        retries = 3
        assembly = "GRCh37"
        self.panelapp_helper = PanelappHelper(PANELAPP_HOST, retries, assembly)
        self.panel_name = "Adult solid tumours"
        self.panel_version = "0.19"

    def test1(self):
        """
        Tests querying of an existing panel filtered by HighEvidence
        :return:
        """
        gene_confidence_threshold = "HighEvidence"
        gene_list = self.panelapp_helper.get_gene_list(panel = self.panel_name,
                                           panel_version = self.panel_version,
                                           gene_confidence_threshold = gene_confidence_threshold)
        self.assertEqual(len(gene_list), 54)
        print "%s genes were returned" % str(len(gene_list))
        print "10 first results: %s..." % ",".join(gene_list[1:min(len(gene_list), 10)])

    def test2(self):
        """
        Tests querying of an existing panel filtered by LowEvidence
        :return:
        """
        gene_confidence_threshold = "LowEvidence"
        gene_list = self.panelapp_helper.get_gene_list(panel = self.panel_name,
                                           panel_version = self.panel_version,
                                           gene_confidence_threshold = gene_confidence_threshold)
        self.assertEqual(len(gene_list), 0)
        print "%s genes were returned" % str(len(gene_list))
        print "10 first results: %s..." % ",".join(gene_list[1:min(len(gene_list), 10)])

    def test3(self):
        """
        Tests querying of an existing panel filtered by LowEvidence
        :return:
        """
        gene_confidence_threshold = ["LowEvidence", "HighEvidence"]
        gene_list = self.panelapp_helper.get_gene_list(panel = self.panel_name,
                                           panel_version = self.panel_version,
                                           gene_confidence_threshold = gene_confidence_threshold)
        self.assertEqual(len(gene_list), 54)
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


class PanelappHelperConnectionRetriesTests(unittest.TestCase):


    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)

        retries = 3
        assembly = "GRCh37"
        self.panelapp_helper = PanelappHelper(PANELAPP_HOST, retries, assembly)
        self.count_failures = 0

        def simulate_connection_failures(*args, **kwargs):
            if self.count_failures < 3:
                self.count_failures += 1
                raise urllib2.URLError("Faked exception")
            return urllib2.urlopen(*args, **kwargs)

        self.panelapp_helper.urlopen = backoff_retrier.wrapper(simulate_connection_failures, retries)
        self.panel_name = "Adult solid tumours"
        self.panel_version = "0.19"

    def test1(self):
        """
        Tests querying of an existing panel filtered by HighEvidence
        :return:
        """
        gene_confidence_threshold = "HighEvidence"
        gene_list = self.panelapp_helper.get_gene_list(panel = self.panel_name,
                                           panel_version = self.panel_version,
                                           gene_confidence_threshold = gene_confidence_threshold)
        self.assertEqual(self.count_failures, 3)
        self.assertEqual(len(gene_list), 54)
        print "%s genes were returned" % str(len(gene_list))
        print "10 first results: %s..." % ",".join(gene_list[1:min(len(gene_list), 10)])

    def test2(self):
        """
        Tests querying of an existing panel filtered by LowEvidence
        :return:
        """
        gene_confidence_threshold = "LowEvidence"
        gene_list = self.panelapp_helper.get_gene_list(panel = self.panel_name,
                                           panel_version = self.panel_version,
                                           gene_confidence_threshold = gene_confidence_threshold)
        self.assertEqual(self.count_failures, 3)
        self.assertEqual(len(gene_list), 0)
        print "%s genes were returned" % str(len(gene_list))
        print "10 first results: %s..." % ",".join(gene_list[1:min(len(gene_list), 10)])

    def test3(self):
        """
        Tests querying of an existing panel filtered by LowEvidence
        :return:
        """
        gene_confidence_threshold = ["LowEvidence", "HighEvidence"]
        gene_list = self.panelapp_helper.get_gene_list(panel = self.panel_name,
                                           panel_version = self.panel_version,
                                           gene_confidence_threshold = gene_confidence_threshold)
        self.assertEqual(self.count_failures, 3)
        self.assertEqual(len(gene_list), 54)
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
            self.assertEqual(self.count_failures, 3)
            self.assertTrue(False, "Function did not fail with an unexisting panel")
        except SystemError:
            pass


if __name__ == '__main__':
    unittest.main()
