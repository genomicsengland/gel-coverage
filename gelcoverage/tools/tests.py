import json
import unittest
import logging
import requests
import urllib.request, urllib.error, urllib.parse
import os

from httpretty import httpretty

from gelcoverage.tools.cellbase_helper import CellbaseHelper
from gelcoverage.tools.panelapp_helper import PanelappHelper
import gelcoverage.tools.backoff_retrier as backoff_retrier

# TODO: mock all REST calls and fix tests
# panelapp.genomicsengland.co.uk/WebServices
PANELAPP_HOST = os.environ.get("PANELAPP_URL")
ASSEMBLY = "GRCh37"
SPECIES = "hsapiens"
CELLBASE_VERSION = "latest"
# https://bio-uat-cellbase.gel.zone/cellbase
CELLBASE_HOST = os.environ.get("CELLBASE_URL")
FILTER_BASIC_FLAG = ["basic"]
FILTER_BIOTYPES = [
    "IG_C_gene",
    "IG_D_gene",
    "IG_J_gene",
    "IG_V_gene",
    "IG_V_gene",
    "protein_coding",
    "nonsense_mediated_decay",
    "non_stop_decay",
    "TR_C_gene",
    "TR_D_gene",
    "TR_J_gene",
    "TR_V_gene",
]


class CellbaseHelperTests(unittest.TestCase):
    _RESOURCES_DIRECTORY = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "resources", "test"
    )

    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)
        self.retries = -1
        self.cellbase_helper = CellbaseHelper(
            SPECIES,
            CELLBASE_VERSION,
            ASSEMBLY,
            CELLBASE_HOST,
            self.retries,
            FILTER_BASIC_FLAG,
            FILTER_BIOTYPES,
        )

    def test1(self):
        """
        Tests get_all_genes()
        :return:
        """
        genes = self.cellbase_helper.get_all_gene_names()
        self.assertEqual(type(genes), list)
        self.assertEqual(len(genes), 20760)

    def test1_1(self):
        """
        Tests get_all_genes() without filtering enabled
        :return:
        """
        genes = self.cellbase_helper.get_all_gene_names(_filter=False)
        self.assertEqual(type(genes), list)
        self.assertEqual(len(genes), 57905)

    def test2(self):
        """
        Tests make_exons_bed() for a list of 3 genes
        :return:
        """
        gene_list = [
            "CFTR",
            "BRCA1",
            "BRCA2",
            "IGHE",  # gene with 3 transcripts flagged as basic and 1 not flagged as basic
        ]
        bed = self.cellbase_helper.make_exons_bed(gene_list)
        # bed.saveas("/home/priesgo/test.bed")

        self.assertIsNotNone(bed)
        ighe_transcripts = []
        for interval in bed:
            chr = interval.chrom
            self.assertEqual(type(chr), str)
            gene, txid, exon_idx = interval.name.split("|")
            self.assertEqual(type(gene), str)
            self.assertEqual(type(txid), str)
            self.assertEqual(type(exon_idx), str)
            strand = interval.strand
            self.assertEqual(type(strand), str)
            start = int(interval.start)
            self.assertEqual(type(start), int)
            end = int(interval.end)
            self.assertEqual(type(end), int)
            if gene == "BRCA1":
                self.assertEqual(chr, "17")
                self.assertEqual(strand, "-")
                gene_start = 41196312
                gene_end = 41322262
                self.assertGreaterEqual(
                    start,
                    gene_start,
                    "Transcript %s exon %s out of bounds: %s < %s"
                    % (txid, exon_idx, str(start), str(gene_start)),
                )
                self.assertLessEqual(
                    end,
                    gene_end,
                    "Transcript %s exon %s out of bounds: %s > %s"
                    % (txid, exon_idx, str(end), str(gene_end)),
                )
            elif gene == "BRCA2":
                self.assertEqual(chr, "13")
                self.assertEqual(strand, "+")
                gene_start = 32889611
                gene_end = 32974403
                self.assertGreaterEqual(
                    start,
                    gene_start,
                    "Transcript %s exon %s out of bounds: %s < %s"
                    % (txid, exon_idx, str(start), str(gene_start)),
                )
                self.assertLessEqual(
                    end,
                    gene_end,
                    "Transcript %s exon %s out of bounds: %s > %s"
                    % (txid, exon_idx, str(end), str(gene_end)),
                )
            elif gene == "CFTR":
                self.assertEqual(chr, "7")
                self.assertEqual(strand, "+")
                gene_start = 117105838
                gene_end = 117356025
                self.assertGreaterEqual(
                    start,
                    gene_start,
                    "Transcript %s exon %s out of bounds: %s < %s"
                    % (txid, exon_idx, str(start), str(gene_start)),
                )
                self.assertLessEqual(
                    end,
                    gene_end,
                    "Transcript %s exon %s out of bounds: %s > %s"
                    % (txid, exon_idx, str(end), str(gene_end)),
                )
            elif gene == "IGHE":
                ighe_transcripts.append(txid)
            else:
                self.assertTrue(
                    False, "Unexpected information in the BED file for gene %s" % gene
                )
        self.assertEqual(
            len(set(ighe_transcripts)), 2
        )  # checks that non basic flagged transcript has been filtered out

    def test2_1(self):
        """
        Tests make_exons_bed() for a list of 3 genes not filtering
        :return:
        """
        gene_list = [
            "CFTR",
            "BRCA1",
            "BRCA2",
            "IGHE",  # gene with 3 transcripts flagged as basic and 1 not flagged as basic
        ]
        bed = self.cellbase_helper.make_exons_bed(gene_list, _filter=False)
        self.assertIsNotNone(bed)
        ighe_transcripts = []
        for interval in bed:
            chr = interval.chrom
            self.assertEqual(type(chr), str)
            gene, txid, exon_idx = interval.name.split("|")
            self.assertEqual(type(gene), str)
            self.assertEqual(type(txid), str)
            self.assertEqual(type(exon_idx), str)
            strand = interval.strand
            self.assertEqual(type(strand), str)
            start = int(interval.start)
            self.assertEqual(type(start), int)
            end = int(interval.end)
            self.assertEqual(type(end), int)
            if gene == "BRCA1":
                self.assertEqual(chr, "17")
                self.assertEqual(strand, "-")
                gene_start = 41196312
                gene_end = 41322262
                self.assertGreaterEqual(
                    start,
                    gene_start,
                    "Transcript %s exon %s out of bounds: %s < %s"
                    % (txid, exon_idx, str(start), str(gene_start)),
                )
                self.assertLessEqual(
                    end,
                    gene_end,
                    "Transcript %s exon %s out of bounds: %s > %s"
                    % (txid, exon_idx, str(end), str(gene_end)),
                )
            elif gene == "BRCA2":
                self.assertEqual(chr, "13")
                self.assertEqual(strand, "+")
                gene_start = 32889611
                gene_end = 32974403
                self.assertGreaterEqual(
                    start,
                    gene_start,
                    "Transcript %s exon %s out of bounds: %s < %s"
                    % (txid, exon_idx, str(start), str(gene_start)),
                )
                self.assertLessEqual(
                    end,
                    gene_end,
                    "Transcript %s exon %s out of bounds: %s > %s"
                    % (txid, exon_idx, str(end), str(gene_end)),
                )
            elif gene == "CFTR":
                self.assertEqual(chr, "7")
                self.assertEqual(strand, "+")
                gene_start = 117105838
                gene_end = 117356025
                self.assertGreaterEqual(
                    start,
                    gene_start,
                    "Transcript %s exon %s out of bounds: %s < %s"
                    % (txid, exon_idx, str(start), str(gene_start)),
                )
                self.assertLessEqual(
                    end,
                    gene_end,
                    "Transcript %s exon %s out of bounds: %s > %s"
                    % (txid, exon_idx, str(end), str(gene_end)),
                )
            elif gene == "IGHE":
                ighe_transcripts.append(txid)
            else:
                self.assertTrue(
                    False, "Unexpected information in the BED file for gene %s" % gene
                )
        self.assertEqual(
            len(set(ighe_transcripts)), 3
        )  # checks that non basic flagged transcript has not been filtered out

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
            self.assertEqual(type(chr), str)
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

    def test8(self):
        """
        Tests __get_gene_info no CellBase errors
        :return:
        """

        httpretty.reset()
        httpretty.enable()

        # Mock CellBase alive query
        httpretty.register_uri(
            httpretty.HEAD, "http://cellbase-fake/cellbase", status=302, body=""
        )

        # Mock CellBase gene query
        httpretty.register_uri(
            httpretty.GET,
            "http://cellbase-fake/cellbase/webservices/rest/v4/hsapiens/feature/gene/search?"
            "assembly=GRCh38&name=BRCA2,ADSL,APOE&transcripts.biotype=IG_C_gene,IG_D_gene,IG_J_gene,"
            "IG_V_gene,IG_V_gene,protein_coding,nonsense_mediated_decay,non_stop_decay,TR_C_gene,"
            "TR_D_gene,TR_J_gene,TR_V_gene&skip=0&limit=1000&include=name,chromosome,"
            "transcripts.exons.start,transcripts.exons.exonNumber,transcripts.id,transcripts.strand,"
            "transcripts.exons.end,transcripts.exons.sequence,exonNumber,"
            "transcripts.annotationFlags",
            status=200,
            body=self._load_json_string(
                os.path.join(
                    self._RESOURCES_DIRECTORY, "cellbase_brca2adslapoe_response.json"
                )
            ),
            content_type="application/json",
        )

        cellbase_helper = CellbaseHelper(
            "hsapiens",
            "v4",
            "GRCh38",
            "http://cellbase-fake/cellbase",
            3,
            ["basic"],
            [
                "IG_C_gene",
                "IG_D_gene",
                "IG_J_gene",
                "IG_V_gene",
                "IG_V_gene",
                "protein_coding",
                "nonsense_mediated_decay",
                "non_stop_decay",
                "TR_C_gene",
                "TR_D_gene",
                "TR_J_gene",
                "TR_V_gene",
            ],
        )
        bed = cellbase_helper.make_exons_bed(["BRCA2", "ADSL", "APOE"], _filter=True)
        # fdw = open("/tmp/test8.tsv", "w")
        # for feature in bed.features():
        #     fdw.write("\t".join([feature.chrom, str(feature.start), str(feature.end), feature.name]) + "\n")
        # fdw.close()

        expected_tuple_set = self._load_tuple_set(
            os.path.join(self._RESOURCES_DIRECTORY, "test8.tsv")
        )
        self.assertEqual(len(expected_tuple_set), bed.count())

        # Check all expected intervals are properly retrieved
        for feature in bed.features():
            self.assertTrue(
                (feature.chrom, str(feature.start), str(feature.end), feature.name)
                in expected_tuple_set
            )

    def test9(self):
        """
        Tests __get_gene_info CellBase returns empty results: resulting bed must be empty - it is not
        make_exons_bed's responsibility to fail if this happens, as the list of genes and filters passed
        to the method can well make CellBase filters to return empty results. This was created to deal
        with the CellBase silent issue in which re-initialising connection to mongo would lead to empty
        results without any further warning or errors. It is InterpretationPipeline's responsibility to
        fail if empty coverage info is returned.
        :return:
        """

        httpretty.reset()
        httpretty.enable()

        # Mock CellBase alive query
        httpretty.register_uri(
            httpretty.HEAD, "http://cellbase-fake/cellbase", status=302, body=""
        )

        # Mock CellBase gene query
        httpretty.register_uri(
            httpretty.GET,
            "http://cellbase-fake/cellbase/webservices/rest/v4/hsapiens/feature/gene/search?"
            "assembly=GRCh38&name=BRCA2,ADSL,APOE&transcripts.biotype=IG_C_gene,IG_D_gene,IG_J_gene,"
            "IG_V_gene,IG_V_gene,protein_coding,nonsense_mediated_decay,non_stop_decay,TR_C_gene,"
            "TR_D_gene,TR_J_gene,TR_V_gene&skip=0&limit=1000&include=name,chromosome,"
            "transcripts.exons.start,transcripts.exons.exonNumber,transcripts.id,transcripts.strand,"
            "transcripts.exons.end,transcripts.exons.sequence,exonNumber,"
            "transcripts.annotationFlags",
            status=200,
            body=self._load_json_string(
                os.path.join(self._RESOURCES_DIRECTORY, "cellbase_empty_response.json")
            ),
            content_type="application/json",
        )

        cellbase_helper = CellbaseHelper(
            "hsapiens",
            "v4",
            "GRCh38",
            "http://cellbase-fake/cellbase",
            3,
            ["basic"],
            [
                "IG_C_gene",
                "IG_D_gene",
                "IG_J_gene",
                "IG_V_gene",
                "IG_V_gene",
                "protein_coding",
                "nonsense_mediated_decay",
                "non_stop_decay",
                "TR_C_gene",
                "TR_D_gene",
                "TR_J_gene",
                "TR_V_gene",
            ],
        )
        bed = cellbase_helper.make_exons_bed(["BRCA2", "ADSL", "APOE"], _filter=True)
        self.assertEqual(0, bed.count())

    @staticmethod
    def _load_json_string(filename):
        fd = open(filename)
        dictionary = json.load(fd)
        fd.close()

        return json.dumps(dictionary)

    @staticmethod
    def _load_tuple_set(filename):
        tuple_set = set()
        fd = open(filename, "r")
        for line in fd:
            tuple_set.add(tuple(line.rstrip().split("\t")))
        fd.close()

        return tuple_set


class CellbaseHelperRetriesTests(unittest.TestCase):
    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)
        self.retries = -1
        self.cellbase_helper = CellbaseHelper(
            SPECIES,
            CELLBASE_VERSION,
            ASSEMBLY,
            CELLBASE_HOST,
            self.retries,
            FILTER_BASIC_FLAG,
            FILTER_BIOTYPES,
        )
        # overwrites wrapped CB search
        self.count_failures = 0

        def simulate_connection_failures(*args, **kwargs):
            if self.count_failures < 3:
                self.count_failures += 1
                raise requests.exceptions.RequestException("Faked exception")
            return self.cellbase_helper.cellbase_gene_client.search(*args, **kwargs)

        self.cellbase_helper.cellbase_search = backoff_retrier.wrapper(
            simulate_connection_failures, self.retries
        )

    def test1(self):
        """
        Tests get_all_genes()
        :return:
        """
        genes = self.cellbase_helper.get_all_gene_names()
        self.assertEqual(type(genes), list)
        self.assertEqual(len(genes), 20760)
        self.assertEqual(self.count_failures, 3)

    def test1_1(self):
        """
        Tests get_all_genes() without filtering enabled
        :return:
        """
        genes = self.cellbase_helper.get_all_gene_names(_filter=False)
        self.assertEqual(type(genes), list)
        self.assertEqual(len(genes), 57905)
        self.assertEqual(self.count_failures, 3)

    def test2(self):
        """
        Tests make_exons_bed() for a list of 3 genes
        :return:
        """
        gene_list = [
            "CFTR",
            "BRCA1",
            "BRCA2",
            "IGHE",  # gene with 3 transcripts flagged as basic and 1 not flagged as basic
        ]
        bed = self.cellbase_helper.make_exons_bed(gene_list)
        # bed.saveas("/home/priesgo/test.bed")
        self.assertEqual(self.count_failures, 3)

        self.assertIsNotNone(bed)
        ighe_transcripts = []
        for interval in bed:
            chr = interval.chrom
            self.assertEqual(type(chr), str)
            gene, txid, exon_idx = interval.name.split("|")
            self.assertEqual(type(gene), str)
            self.assertEqual(type(txid), str)
            self.assertEqual(type(exon_idx), str)
            strand = interval.strand
            self.assertEqual(type(strand), str)
            start = int(interval.start)
            self.assertEqual(type(start), int)
            end = int(interval.end)
            self.assertEqual(type(end), int)
            if gene == "BRCA1":
                self.assertEqual(chr, "17")
                self.assertEqual(strand, "-")
                gene_start = 41196312
                gene_end = 41322262
                self.assertGreaterEqual(
                    start,
                    gene_start,
                    "Transcript %s exon %s out of bounds: %s < %s"
                    % (txid, exon_idx, str(start), str(gene_start)),
                )
                self.assertLessEqual(
                    end,
                    gene_end,
                    "Transcript %s exon %s out of bounds: %s > %s"
                    % (txid, exon_idx, str(end), str(gene_end)),
                )
            elif gene == "BRCA2":
                self.assertEqual(chr, "13")
                self.assertEqual(strand, "+")
                gene_start = 32889611
                gene_end = 32974403
                self.assertGreaterEqual(
                    start,
                    gene_start,
                    "Transcript %s exon %s out of bounds: %s < %s"
                    % (txid, exon_idx, str(start), str(gene_start)),
                )
                self.assertLessEqual(
                    end,
                    gene_end,
                    "Transcript %s exon %s out of bounds: %s > %s"
                    % (txid, exon_idx, str(end), str(gene_end)),
                )
            elif gene == "CFTR":
                self.assertEqual(chr, "7")
                self.assertEqual(strand, "+")
                gene_start = 117105838
                gene_end = 117356025
                self.assertGreaterEqual(
                    start,
                    gene_start,
                    "Transcript %s exon %s out of bounds: %s < %s"
                    % (txid, exon_idx, str(start), str(gene_start)),
                )
                self.assertLessEqual(
                    end,
                    gene_end,
                    "Transcript %s exon %s out of bounds: %s > %s"
                    % (txid, exon_idx, str(end), str(gene_end)),
                )
            elif gene == "IGHE":
                ighe_transcripts.append(txid)
            else:
                self.assertTrue(
                    False, "Unexpected information in the BED file for gene %s" % gene
                )
        self.assertEqual(
            len(set(ighe_transcripts)), 2
        )  # checks that non basic flagged transcript has been filtered out

    def test2_1(self):
        """
        Tests make_exons_bed() for a list of 3 genes not filtering
        :return:
        """
        gene_list = [
            "CFTR",
            "BRCA1",
            "BRCA2",
            "IGHE",  # gene with 3 transcripts flagged as basic and 1 not flagged as basic
        ]
        bed = self.cellbase_helper.make_exons_bed(gene_list, _filter=False)
        self.assertIsNotNone(bed)
        self.assertEqual(self.count_failures, 3)
        ighe_transcripts = []
        for interval in bed:
            chr = interval.chrom
            self.assertEqual(type(chr), str)
            gene, txid, exon_idx = interval.name.split("|")
            self.assertEqual(type(gene), str)
            self.assertEqual(type(txid), str)
            self.assertEqual(type(exon_idx), str)
            strand = interval.strand
            self.assertEqual(type(strand), str)
            start = int(interval.start)
            self.assertEqual(type(start), int)
            end = int(interval.end)
            self.assertEqual(type(end), int)
            if gene == "BRCA1":
                self.assertEqual(chr, "17")
                self.assertEqual(strand, "-")
                gene_start = 41196312
                gene_end = 41322262
                self.assertGreaterEqual(
                    start,
                    gene_start,
                    "Transcript %s exon %s out of bounds: %s < %s"
                    % (txid, exon_idx, str(start), str(gene_start)),
                )
                self.assertLessEqual(
                    end,
                    gene_end,
                    "Transcript %s exon %s out of bounds: %s > %s"
                    % (txid, exon_idx, str(end), str(gene_end)),
                )
            elif gene == "BRCA2":
                self.assertEqual(chr, "13")
                self.assertEqual(strand, "+")
                gene_start = 32889611
                gene_end = 32974403
                self.assertGreaterEqual(
                    start,
                    gene_start,
                    "Transcript %s exon %s out of bounds: %s < %s"
                    % (txid, exon_idx, str(start), str(gene_start)),
                )
                self.assertLessEqual(
                    end,
                    gene_end,
                    "Transcript %s exon %s out of bounds: %s > %s"
                    % (txid, exon_idx, str(end), str(gene_end)),
                )
            elif gene == "CFTR":
                self.assertEqual(chr, "7")
                self.assertEqual(strand, "+")
                gene_start = 117105838
                gene_end = 117356025
                self.assertGreaterEqual(
                    start,
                    gene_start,
                    "Transcript %s exon %s out of bounds: %s < %s"
                    % (txid, exon_idx, str(start), str(gene_start)),
                )
                self.assertLessEqual(
                    end,
                    gene_end,
                    "Transcript %s exon %s out of bounds: %s > %s"
                    % (txid, exon_idx, str(end), str(gene_end)),
                )
            elif gene == "IGHE":
                ighe_transcripts.append(txid)
            else:
                self.assertTrue(
                    False, "Unexpected information in the BED file for gene %s" % gene
                )
        self.assertEqual(
            len(set(ighe_transcripts)), 3
        )  # checks that non basic flagged transcript has not been filtered out

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
            self.assertEqual(type(chr), str)
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
            SPECIES,
            CELLBASE_VERSION,
            ASSEMBLY,
            CELLBASE_HOST,
            self.retries,
            FILTER_BASIC_FLAG,
            FILTER_BIOTYPES,
        )
        # overwrites wrapped CB search
        self.count_failures = 0

        def simulate_connection_failures(*args, **kwargs):
            if self.count_failures < 3:
                self.count_failures += 1
                raise requests.exceptions.ConnectionError("Faked exception")
            return self.cellbase_helper.cellbase_gene_client.search(*args, **kwargs)

        self.cellbase_helper.cellbase_search = backoff_retrier.wrapper(
            simulate_connection_failures, self.retries
        )

    def test1(self):
        """
        Tests get_all_genes()
        :return:
        """
        genes = self.cellbase_helper.get_all_gene_names()
        self.assertEqual(type(genes), list)
        self.assertEqual(len(genes), 20760)
        self.assertEqual(self.count_failures, 3)


class CellbaseHelperRetriesTests3(unittest.TestCase):
    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)
        self.retries = -1
        self.cellbase_helper = CellbaseHelper(
            SPECIES,
            CELLBASE_VERSION,
            ASSEMBLY,
            CELLBASE_HOST,
            self.retries,
            FILTER_BASIC_FLAG,
            FILTER_BIOTYPES,
        )
        # overwrites wrapped CB search
        self.count_failures = 0

        def simulate_connection_failures(*args, **kwargs):
            if self.count_failures < 3:
                self.count_failures += 1
                raise ValueError("Faked exception")
            return self.cellbase_helper.cellbase_gene_client.search(*args, **kwargs)

        self.cellbase_helper.cellbase_search = backoff_retrier.wrapper(
            simulate_connection_failures, self.retries
        )

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
            SPECIES,
            CELLBASE_VERSION,
            ASSEMBLY,
            CELLBASE_HOST,
            self.retries,
            FILTER_BASIC_FLAG,
            FILTER_BIOTYPES,
        )
        # overwrites wrapped CB search
        self.count_failures = 0

        def simulate_connection_failures(*args, **kwargs):
            self.count_failures += 1
            raise requests.exceptions.ConnectionError("Faked exception")

        self.cellbase_helper.cellbase_search = backoff_retrier.wrapper(
            simulate_connection_failures, self.retries
        )

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
        assembly = "GRCh38"
        self.panelapp_helper = PanelappHelper(PANELAPP_HOST, assembly)
        self.panel_name = "245"
        self.panel_version = "2.4"

    def test1(self):
        """
        Tests querying of an existing panel filtered by HighEvidence
        :return:
        """
        gene_confidence_threshold = 3
        gene_list = self.panelapp_helper.get_gene_list(
            panel=self.panel_name,
            panel_version=self.panel_version,
            gene_confidence_threshold=gene_confidence_threshold,
        )
        self.assertEqual(len(gene_list), 82)

    def test2(self):
        """
        Tests querying of an existing panel filtered by LowEvidence
        :return:
        """
        gene_confidence_threshold = 2
        gene_list = self.panelapp_helper.get_gene_list(
            panel=self.panel_name,
            panel_version=self.panel_version,
            gene_confidence_threshold=gene_confidence_threshold,
        )
        self.assertEqual(len(gene_list), 100)

    def test3(self):
        """
        Tests querying of an unexisting panel
        :return:
        """
        gene_confidence_threshold = 3
        try:
            gene_list = self.panelapp_helper.get_gene_list(
                panel="unexisting_disease",
                panel_version="0.99",
                gene_confidence_threshold=gene_confidence_threshold,
            )
            self.assertTrue(False, "Function did not fail with an unexisting panel")
        except SystemError:
            pass


if __name__ == "__main__":
    unittest.main()
