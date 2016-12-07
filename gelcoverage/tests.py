import unittest
import json
import logging
import pybedtools
from gelcoverage.runner import GelCoverageRunner
from gelcoverage.test.output_verifier import OutputVerifier


class GelCoverageRunnerTests(OutputVerifier):

    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)
        self.config = {
            # Sets parameters from CLI
            "bw": "../resources/test/test1.bw",
            "configuration_file": "../resources/exon_coverage_summary.config",
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
            "transcript_filtering_biotypes": "IG_C_gene,IG_D_gene,IG_J_gene,IG_V_gene,IG_V_gene,protein_coding,nonsense_mediated_decay,non_stop_decay,TR_C_gene,TR_D_gene,TR_J_gene,TR_V_gene",
            "wg_stats_enabled": False,
            "exon_padding": 15
        }

    def test1(self):
        """
        Test 1: panel from PanelApp
        :return:
        """
        expected_gene_list = [u'SCN2A', u'SPTAN1', u'PLCB1', u'SLC25A22', u'SCN8A', u'STXBP1', u'PNKP']
        self.config["bw"] = "../resources/test/test1.bw"
        self.config["panel"] = "Epileptic encephalopathy"
        self.config["panel_version"] = "0.2"
        self.config["exon_padding"] = 0
        runner = GelCoverageRunner(
            config=self.config
        )
        output, bed = runner.run()
        # Writes the JSON
        with open('../resources/test/sample_output_1.json', 'w') as fp:
            json.dump(output, fp)
        # Verifies the bed...
        self.assertEqual(type(bed), pybedtools.bedtool.BedTool)
        # Saves the analysed region as a BED file
        bed.saveas('../resources/test/sample_output_1.bed')
        # Runs verifications on output JSON
        self.verify_output(output, expected_gene_list)



    def test2(self):
        """
        Test 2: gene list
        :return:
        """
        self.config["panel"] = None
        self.config["panel_version"] = None
        self.config["bw"] = "../resources/test/test2.bw"
        self.config["gene_list"] = "BRCA1,BRCA2,CFTR,IGHE"
        self.config["exon_padding"] = 0
        expected_gene_list = map(
            lambda x: unicode(x),
            self.config["gene_list"].split(",")
        )
        runner = GelCoverageRunner(
            config=self.config
        )
        output, bed = runner.run()
        # Writes the JSON
        with open('../resources/test/sample_output_2.json', 'w') as fp:
            json.dump(output, fp)
        # Verifies the bed...
        self.assertEqual(type(bed), pybedtools.bedtool.BedTool)
        # Saves the analysed region as a BED file
        bed.saveas('../resources/test/sample_output_2.bed')
        # Runs verifications on output JSON
        self.verify_output(output, expected_gene_list)


    def test3(self):
        """
        Test 3: panel from PanelApp with exon padding of 15 bp
        :return:
        """
        expected_gene_list = [u'SCN2A', u'SPTAN1', u'PLCB1', u'SLC25A22', u'SCN8A', u'STXBP1', u'PNKP']
        self.config["bw"] = "../resources/test/test1.bw"
        self.config["panel"] = "Epileptic encephalopathy"
        self.config["panel_version"] = "0.2"
        self.config["exon_padding"] = 15
        runner = GelCoverageRunner(
            config=self.config
        )
        output, bed = runner.run()
        # Writes the JSON
        with open('../resources/test/sample_output_3.json', 'w') as fp:
            json.dump(output, fp)
        # Verifies the bed...
        self.assertEqual(type(bed), pybedtools.bedtool.BedTool)
        # Saves the analysed region as a BED file
        bed.saveas('../resources/test/sample_output_3.bed')
        # Runs verifications on output JSON
        self.verify_output(output, expected_gene_list)


    def test4(self):
        """
        Test 4: gene list with exon padding of 15 bp
        :return:
        """
        self.config["panel"] = None
        self.config["panel_version"] = None
        self.config["bw"] = "../resources/test/test2.bw"
        self.config["gene_list"] = "BRCA1,BRCA2,CFTR,IGHE"
        self.config["exon_padding"] = 15
        expected_gene_list = map(
            lambda x: unicode(x),
            self.config["gene_list"].split(",")
        )
        runner = GelCoverageRunner(
            config=self.config
        )
        output, bed = runner.run()
        # Saves the analysed region as a BED file
        bed.saveas('../resources/test/sample_output_4.bed')
        # Runs verifications on output JSON
        self.verify_output(output, expected_gene_list)
        # Writes the JSON
        with open('../resources/test/sample_output_4.json', 'w') as fp:
            json.dump(output, fp)
        # Verifies the bed...
        self.assertEqual(type(bed), pybedtools.bedtool.BedTool)

    def test5(self):
        """
        Test 5: tests the union transcript build up
        :return:
        """
        def create_test_exon(start, end, exon_number):
            return {
                "start": start,
                "end": end,
                "padded_start": start - self.config["exon_padding"],
                "padded_end": end + self.config["exon_padding"],
                "exon_number": exon_number
            }
        # Interval covered in the input file SCN2A: 2: 165,995,882-166,349,242
        offset = 165995882
        runner = GelCoverageRunner(
            config=self.config
        )
        # Runs union transcript with exon padding 15bp
        self.config["exon_padding"] = 15
        gene_15bp_padding = {
            "chromosome": "2",
            "name": "TEST",
            "transcripts": [
                {
                    "id": "1",
                    "exons": [
                        create_test_exon(offset + 10, offset + 15, 1),
                        create_test_exon(offset + 20, offset + 25, 2),
                        create_test_exon(offset + 30, offset + 35, 3),
                        create_test_exon(offset + 40, offset + 45, 4),
                        create_test_exon(offset + 50, offset + 55, 5),
                        create_test_exon(offset + 60, offset + 65, 6)
                    ]
                },
                {
                    "id": "2",
                    "exons": [
                        create_test_exon(offset + 10, offset + 15, 1),
                        create_test_exon(offset + 20, offset + 25, 2),
                        create_test_exon(offset + 30, offset + 35, 3),
                        create_test_exon(offset + 40, offset + 45, 4),
                        create_test_exon(offset + 50, offset + 55, 5),
                        create_test_exon(offset + 60, offset + 65, 6),
                        create_test_exon(offset + 70, offset + 75, 7)
                    ]
                },
                {
                    "id": "3",
                    "exons": [
                        create_test_exon(offset + 0, offset + 5, 1),
                        create_test_exon(offset + 10, offset + 15, 2),
                        create_test_exon(offset + 20, offset + 25, 3),
                        create_test_exon(offset + 30, offset + 35, 4),
                        create_test_exon(offset + 40, offset + 45, 5),
                        create_test_exon(offset + 50, offset + 55, 6),
                        create_test_exon(offset + 60, offset + 65, 7)
                    ]
                }
            ]
        }
        union_transcript = runner._GelCoverageRunner__create_union_transcript(gene_15bp_padding)
        self.assertEqual(len(union_transcript["exons"]), 1)
        gene_15bp_padding["union_transcript"] = union_transcript
        self.verify_union_transcript(gene_15bp_padding)
        # Runs union transcript with exon padding 0bp
        self.config["exon_padding"] = 0
        gene_0bp_padding = {
            "chromosome": "2",
            "name": "TEST",
            "transcripts": [
                {
                    "id": "1",
                    "exons": [
                        create_test_exon(offset + 10, offset + 15, 1),
                        create_test_exon(offset + 20, offset + 25, 2),
                        create_test_exon(offset + 30, offset + 35, 3),
                        create_test_exon(offset + 40, offset + 45, 4),
                        create_test_exon(offset + 50, offset + 55, 5),
                        create_test_exon(offset + 60, offset + 65, 6)
                    ]
                },
                {
                    "id": "2",
                    "exons": [
                        create_test_exon(offset + 10, offset + 15, 1),
                        create_test_exon(offset + 20, offset + 25, 2),
                        create_test_exon(offset + 30, offset + 35, 3),
                        create_test_exon(offset + 40, offset + 45, 4),
                        create_test_exon(offset + 50, offset + 55, 5),
                        create_test_exon(offset + 60, offset + 65, 6),
                        create_test_exon(offset + 70, offset + 75, 7)
                    ]
                },
                {
                    "id": "3",
                    "exons": [
                        create_test_exon(offset + 0, offset + 5, 1),
                        create_test_exon(offset + 10, offset + 15, 2),
                        create_test_exon(offset + 20, offset + 25, 3),
                        create_test_exon(offset + 30, offset + 35, 4),
                        create_test_exon(offset + 40, offset + 45, 5),
                        create_test_exon(offset + 50, offset + 55, 6),
                        create_test_exon(offset + 60, offset + 65, 7)
                    ]
                }
            ]
        }
        union_transcript = runner._GelCoverageRunner__create_union_transcript(gene_0bp_padding)
        self.assertEqual(len(union_transcript["exons"]), 8)
        gene_0bp_padding["union_transcript"] = union_transcript
        self.verify_union_transcript(gene_0bp_padding)