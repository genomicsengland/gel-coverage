import unittest
import json
from gelcoverage.runner import GelCoverageRunner


class GelCoverageRunnerTests(unittest.TestCase):

    def setUp(self):
        self.config = {
            # Sets parameters from CLI
            "bw": "../resources/test/test1.bw",  # TODO: find a testing file
            "panel": "Epileptic encephalopathy",
            "panel_version": "0.2",
            #"gene_list": args.gene_list,
            "coverage_threshold": 30,
             # Sets parameters from config file
            "cellbase_species": "hsapiens",
            "cellbase_version": "latest",
            "cellbase_assembly": "grch37",
            "cellbase_host": "10.5.8.201:8080/cellbase-4.5.0-rc",
            "panelapp_host": "bioinfo.extge.co.uk/crowdsourcing/WebServices",
            "panelapp_gene_confidence": "HighEvidence",
            "transcript_filtering_flags": "basic",
            "transcript_filtering_biotypes": "IG_C_gene,IG_D_gene,IG_J_gene,IG_V_gene,IG_V_gene,protein_coding,nonsense_mediated_decay,non_stop_decay,TR_C_gene,TR_D_gene,TR_J_gene,TR_V_gene"
        }

    def verify_parameters(self, parameters, expected_gene_list):
        self.assertEqual(type(parameters), dict)
        self.assertEqual(parameters["gap_coverage_threshold"], self.config["coverage_threshold"])
        self.assertEqual(parameters["input_file"], self.config["bw"])
        self.assertEqual(parameters["species"], self.config["cellbase_species"])
        self.assertEqual(parameters["assembly"], self.config["cellbase_assembly"])
        #self.assertEqual(parameters["panel"], self.config["panel"])
        #self.assertEqual(parameters["panel_version"], self.config["panel_version"])
        self.assertEqual(parameters["gene_list"],
                         expected_gene_list)
        self.assertEqual(parameters["transcript_filtering_flags"],
                         self.config["transcript_filtering_flags"])
        self.assertEqual(parameters["transcript_filtering_biotypes"],
                         self.config["transcript_filtering_biotypes"])

    def verify_transcript(self, transcript):
        self.assertEqual(type(transcript["id"]), unicode)
        self.assertTrue(str(transcript["id"]).startswith("ENS"))
        # TODO: test flags and biotype
        self.assertEqual(type(transcript["statistics"]), dict)
        self.assertEqual(type(transcript["statistics"]["bases_gte_15x"]), int)
        self.assertTrue(transcript["statistics"]["bases_gte_15x"] >= 0)
        self.assertEqual(type(transcript["statistics"]["bases_gte_30x"]), int)
        self.assertTrue(transcript["statistics"]["bases_gte_30x"] >= 0)
        self.assertEqual(type(transcript["statistics"]["bases_gte_50x"]), int)
        self.assertTrue(transcript["statistics"]["bases_gte_50x"] >= 0)
        self.assertEqual(type(transcript["statistics"]["bases_lt_15x"]), int)
        self.assertTrue(transcript["statistics"]["bases_lt_15x"] >= 0)
        self.assertEqual(type(transcript["statistics"]["bases_lt_3x"]), int)
        self.assertTrue(transcript["statistics"]["bases_lt_3x"] >= 0)
        self.assertEqual(type(transcript["statistics"]["gc_content"]), float)
        self.assertTrue(transcript["statistics"]["gc_content"] >= 0 and
                        transcript["statistics"]["gc_content"] <= 1)
        self.assertEqual(type(transcript["statistics"]["mean"]), float)
        self.assertTrue(transcript["statistics"]["mean"] >= 0)
        self.assertEqual(type(transcript["statistics"]["percent_gte_15x"]), float)
        self.assertTrue(transcript["statistics"]["percent_gte_15x"] >= 0 and
                        transcript["statistics"]["percent_gte_15x"] <= 1)
        self.assertTrue(transcript["statistics"]["percent_gte_30x"] >= 0 and
                        transcript["statistics"]["percent_gte_30x"] <= 1)
        self.assertTrue(transcript["statistics"]["percent_gte_50x"] >= 0 and
                        transcript["statistics"]["percent_gte_50x"] <= 1)
        self.assertTrue(transcript["statistics"]["percent_lt_15x"] >= 0 and
                        transcript["statistics"]["percent_lt_15x"] <= 1)
        self.assertEqual(type(transcript["statistics"]["total_bases"]), int)
        self.assertTrue(transcript["statistics"]["total_bases"] >= 0)
        self.assertEqual(type(transcript["statistics"]["weighted_median"]), float)
        self.assertTrue(transcript["statistics"]["weighted_median"] >= 0)
        self.assertEqual(type(transcript["statistics"]["weighted_pct75"]), float)
        self.assertTrue(transcript["statistics"]["weighted_pct75"] >= 0)
        self.assertEqual(type(transcript["statistics"]["weighted_pct25"]), float)
        self.assertTrue(transcript["statistics"]["weighted_pct25"] >= 0)

    def verify_exon(self, exon):
        self.assertEqual(type(exon), dict)
        self.assertEqual(type(exon["start"]), int)
        self.assertTrue(exon["start"] >= 0)
        self.assertEqual(type(exon["end"]), int)
        self.assertTrue(exon["end"] >= 0)
        self.assertEqual(type(exon["exon_number"]), unicode)
        self.assertTrue(str(exon["exon_number"]).startswith("exon"))
        self.assertEqual(type(exon["statistics"]), dict)
        self.assertEqual(type(exon["statistics"]["bases_gte_15x"]), int)
        self.assertTrue(exon["statistics"]["bases_gte_15x"] >= 0)
        self.assertEqual(type(exon["statistics"]["bases_gte_30x"]), int)
        self.assertTrue(exon["statistics"]["bases_gte_30x"] >= 0)
        self.assertEqual(type(exon["statistics"]["bases_gte_50x"]), int)
        self.assertTrue(exon["statistics"]["bases_gte_50x"] >= 0)
        self.assertEqual(type(exon["statistics"]["bases_lt_15x"]), int)
        self.assertTrue(exon["statistics"]["bases_lt_15x"] >= 0)
        self.assertEqual(type(exon["statistics"]["bases_lt_3x"]), int)
        self.assertTrue(exon["statistics"]["bases_lt_3x"] >= 0)
        self.assertEqual(type(exon["statistics"]["gc_content"]), float)
        self.assertTrue(exon["statistics"]["gc_content"] >= 0 and
                        exon["statistics"]["gc_content"] <= 1)
        self.assertEqual(type(exon["statistics"]["mean"]), float)
        self.assertTrue(exon["statistics"]["mean"] >= 0)
        self.assertEqual(type(exon["statistics"]["percent_gte_15x"]), float)
        self.assertTrue(exon["statistics"]["percent_gte_15x"] >= 0 and
                        exon["statistics"]["percent_gte_15x"] <= 1)
        self.assertTrue(exon["statistics"]["percent_gte_30x"] >= 0 and
                        exon["statistics"]["percent_gte_30x"] <= 1)
        self.assertTrue(exon["statistics"]["percent_gte_50x"] >= 0 and
                        exon["statistics"]["percent_gte_50x"] <= 1)
        self.assertTrue(exon["statistics"]["percent_lt_15x"] >= 0 and
                        exon["statistics"]["percent_lt_15x"] <= 1)
        self.assertEqual(type(exon["statistics"]["total_bases"]), int)
        self.assertTrue(exon["statistics"]["total_bases"] >= 0)
        self.assertEqual(type(exon["statistics"]["median"]), float)
        self.assertTrue(exon["statistics"]["median"] >= 0)
        self.assertEqual(type(exon["statistics"]["pct75"]), float)
        self.assertTrue(exon["statistics"]["pct75"] >= 0)
        self.assertEqual(type(exon["statistics"]["pct25"]), float)
        self.assertTrue(exon["statistics"]["pct25"] >= 0)

    def verify_gap(self, gap, exon):
        self.assertEqual(type(gap), dict)
        self.assertEqual(type(gap["start"]), int)
        self.assertTrue(gap["start"] >= exon["start"] and gap["start"] <= exon["end"])
        self.assertEqual(type(gap["end"]), int)
        self.assertTrue(gap["end"] >= exon["start"] and gap["start"] <= exon["end"] and
                        gap["end"] >= gap["start"])
        self.assertEqual(type(gap["length"]), int)
        self.assertTrue(gap["length"] >= 1 and gap["length"] <= gap["end"] - gap["start"] + 1)

    def test1(self):
        """
        Test 1
        :return:
        """
        expected_gene_list = [u'SCN2A', u'SPTAN1', u'PLCB1', u'SLC25A22', u'SCN8A', u'STXBP1', u'PNKP']
        self.config["bw"] = "../resources/test/test1.bw"
        self.config["panel"] = "Epileptic encephalopathy"
        self.config["panel_version"] = "0.2"
        runner = GelCoverageRunner(
            config=self.config
        )
        output = runner.run()
        self.assertEqual(type(output), dict)
        # Verify that content in parameters is correct
        self.verify_parameters(output["parameters"], expected_gene_list)
        # Verify that coverage results are correct
        self.assertEqual(type(output["results"]), list)
        # Verify every gene
        for gene in output["results"]:
            self.assertTrue(gene["name"] in expected_gene_list)
            self.assertEqual(type(gene["chromosome"]), unicode)
            # Verify every transcript
            for transcript in gene["transcripts"]:
                self.verify_transcript(transcript)
                # Verify every exon
                for exon in transcript["exons"]:
                    self.verify_exon(exon)
                    # Verify gaps
                    for gap in exon["gaps"]:
                        self.verify_gap(gap, exon)
        with open('../resources/test/output1.json', 'w') as fp:
            json.dump(output, fp)

    def test2(self):
        """
        Test 2
        :return:
        """
        self.config["panel"] = None
        self.config["panel_version"] = None
        self.config["bw"] = "../resources/test/test2.bw"
        self.config["gene_list"] = "BRCA1,BRCA2,CFTR,IGHE"
        expected_gene_list = map(
            lambda x: unicode(x),
            self.config["gene_list"].split(",")
        )
        runner = GelCoverageRunner(
            config=self.config
        )
        output = runner.run()
        self.assertEqual(type(output), dict)
        runner = GelCoverageRunner(
            config=self.config
        )
        output = runner.run()
        self.assertEqual(type(output), dict)
        # Verify that content in parameters is correct
        self.verify_parameters(output["parameters"], expected_gene_list)
        # Verify that coverage results are correct
        self.assertEqual(type(output["results"]), list)
        # Verify every gene
        for gene in output["results"]:
            self.assertTrue(gene["name"] in expected_gene_list)
            self.assertEqual(type(gene["chromosome"]), unicode)
            # Verify every transcript
            for transcript in gene["transcripts"]:
                self.verify_transcript(transcript)
                # Verify every exon
                for exon in transcript["exons"]:
                    self.verify_exon(exon)
                    # Verify gaps
                    for gap in exon["gaps"]:
                        self.verify_gap(gap, exon)
        with open('../resources/test/output2.json', 'w') as fp:
            json.dump(output, fp)