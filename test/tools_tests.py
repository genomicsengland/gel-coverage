import unittest
import json
import os
import logging
import pybedtools
from gelcoverage.runner import GelCoverageRunner, GelCoverageInputError, BedMaker
from .test.output_verifier import OutputVerifier
import gelcoverage.constants as constants



PANELAPP_HOST = "panelapp.genomicsengland.co.uk/WebServices"
ASSEMBLY = "GRCh37"
SPECIES = "hsapiens"
CELLBASE_VERSION = "latest"
CELLBASE_HOST = os.environ.get('CELLBASE_URL')
FILTER_BASIC_FLAG = "basic"
FILTER_BIOTYPES = "IG_C_gene,IG_D_gene,IG_J_gene,IG_V_gene,IG_V_gene,protein_coding,nonsense_mediated_decay," \
                  "non_stop_decay,TR_C_gene,TR_D_gene,TR_J_gene,TR_V_gene"


class BedMakerTests(unittest.TestCase):

    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)
        self.config = {
            # Sets parameters from CLI
            "gene_list": None,
            "chr_prefix": False,
            "log_level": 10,
            "transcript_filtering_flags": "basic",
            "transcript_filtering_biotypes": "IG_C_gene,IG_D_gene,IG_J_gene,IG_V_gene,IG_V_gene,protein_coding,"
                                             "nonsense_mediated_decay,non_stop_decay,TR_C_gene,TR_D_gene,TR_J_gene,"
                                             "TR_V_gene",
            "cellbase_species": "hsapiens",
            "cellbase_version": "latest",
            "cellbase_assembly": "grch37",
            "cellbase_host": CELLBASE_HOST,
            "cellbase_retries": -1,
        }

    def test1(self):
        """
        Generates a bed file for all genes in Cellbase
        :return:
        """
        runner = BedMaker(config=self.config)
        bed = runner.run()
        # Verifies the bed...
        self.assertEqual(type(bed), pybedtools.bedtool.BedTool)
        # Saves the analysed region as a BED file
        bed.saveas('../resources/test/sample_output_bedmaker_1.bed')
        observed_genes = self._get_genes_from_bed(bed)
        self.assertEqual(len(observed_genes), 20567)

    def test2(self):
        """
        Generates a bed file for all genes in Cellbase in assembly GRCh38
        :return:
        """
        self.config['cellbase_assembly'] = "grch38"
        runner = BedMaker(config=self.config)
        bed = runner.run()
        # Verifies the bed...
        self.assertEqual(type(bed), pybedtools.bedtool.BedTool)
        # Saves the analysed region as a BED file
        bed.saveas('../resources/test/sample_output_bedmaker_1.bed')
        observed_genes = self._get_genes_from_bed(bed)
        self.assertEqual(len(observed_genes), 20542)

    def test3(self):
        """
        Generates a bed file for some genes
        :return:
        """
        self.config['cellbase_assembly'] = "grch38"
        expected_genes = ['SCN2A', 'SPTAN1', 'PLCB1', 'SLC25A22', 'SCN8A', 'STXBP1', 'PNKP']
        self.config['gene_list'] = expected_genes
        runner = BedMaker(config=self.config)
        bed = runner.run()
        # Verifies the bed...
        self.assertEqual(type(bed), pybedtools.bedtool.BedTool)
        # Saves the analysed region as a BED file
        bed.saveas('../resources/test/sample_output_bedmaker_3.bed')
        observed_genes = self._get_genes_from_bed(bed)
        self.assertEqual(len(observed_genes), len(expected_genes))
        self.assertEqual(set(observed_genes), set(expected_genes))

    def test4(self):
        """
        Generates a bed file for some genes
        :return:
        """
        expected_genes = ['SCN2A', 'SPTAN1', 'PLCB1', 'SLC25A22', 'SCN8A', 'STXBP1', 'PNKP']
        self.config['gene_list'] = expected_genes
        runner = BedMaker(config=self.config)
        bed = runner.run()
        # Verifies the bed...
        self.assertEqual(type(bed), pybedtools.bedtool.BedTool)
        # Saves the analysed region as a BED file
        bed.saveas('../resources/test/sample_output_bedmaker_4.bed')
        observed_genes = self._get_genes_from_bed(bed)
        self.assertEqual(len(observed_genes), len(expected_genes))
        self.assertEqual(set(observed_genes), set(expected_genes))

    def test5(self):
        """
        Generates a bed file for some genes
        :return:
        """
        expected_genes = ['SCN2A', 'SPTAN1', 'PLCB1', 'SLC25A22', 'SCN8A', 'STXBP1', 'PNKP']
        self.config['gene_list'] = expected_genes
        self.config['exon_padding'] = 15
        runner = BedMaker(config=self.config)
        bed = runner.run()
        # Verifies the bed...
        self.assertEqual(type(bed), pybedtools.bedtool.BedTool)
        # Saves the analysed region as a BED file
        bed.saveas('../resources/test/sample_output_bedmaker_5.bed')
        observed_genes = self._get_genes_from_bed(bed)
        self.assertEqual(len(observed_genes), len(expected_genes))
        self.assertEqual(set(observed_genes), set(expected_genes))

    def _get_genes_from_bed(self, bed):
        genes = set()
        for interval in bed:
            # Reads data from BED entry
            chromosome, start, end, gene_name, transcript_id, \
            exon_number, strand, gc_content = GelCoverageRunner._parse_bed_interval(interval)
            genes.add(gene_name)
        return list(genes)


class GelCoverageRunnerTests(OutputVerifier):

    def setUp(self):
        logging.basicConfig(level=logging.INFO)
        self.config = {
            # Sets parameters from CLI
            "bw": "../resources/test/test1.bw",
            "configuration_file": "../resources/bigwig_analyser.config",
            "panel": "Epileptic encephalopathy",
            "panel_version": "1.2",
            "coverage_threshold": 30,
            "cellbase_species": SPECIES,
            "cellbase_version": CELLBASE_VERSION,
            "cellbase_assembly": ASSEMBLY,
            "cellbase_host": CELLBASE_HOST,
            "cellbase_retries": -1,
            "panelapp_host": PANELAPP_HOST,
            "panelapp_gene_confidence": "HighEvidence",
            "panelapp_retries": -1,
            "panelapp_assembly": ASSEMBLY,
            "transcript_filtering_flags": FILTER_BASIC_FLAG,
            "transcript_filtering_biotypes": FILTER_BIOTYPES,
            "wg_stats_enabled": False,
            "exon_padding": 15,
            "wg_regions": None,
            "exon_stats_enabled": True,
            "coding_region_stats_enabled": True
        }

    def test1(self):
        """
        Test 1: panel from PanelApp
        :return:
        """
        expected_gene_list = ['ADSL', 'ALG13', 'ARHGEF9', 'ARX', 'ATP1A3', 'ATRX', 'CDKL5', 'CHD2',
                              'CNTNAP2', 'DNM1', 'DOCK7', 'DYRK1A', 'EHMT1', 'FOXG1', 'GABRA1', 'GABRB3',
                              'GNAO1', 'GRIN1', 'GRIN2A', 'GRIN2B', 'HCN1', 'IQSEC2', 'KCNA2', 'KCNB1',
                              'KCNJ10', 'KCNQ2', 'KCNQ3', 'KCNT1', 'KIAA2022', 'KIF1BP', 'MAPK10', 'MBD5',
                              'MECP2', 'MEF2C', 'PCDH19', 'PIGA', 'PLCB1', 'PNKP', 'POLG', 'PRRT2', 'PURA',
                              'QARS', 'SCN1A', 'SCN1B', 'SCN2A', 'SCN8A', 'SETD5', 'SIK1', 'SLC12A5',
                              'SLC13A5', 'SLC16A2', 'SLC25A22', 'SLC2A1', 'SLC6A1', 'SLC9A6', 'SPTAN1',
                              'STX1B', 'STXBP1', 'SYNGAP1', 'TCF4', 'UBE2A', 'UBE3A', 'WDR45', 'WWOX', 'ZEB2']
        self.config["panel"] = "Epileptic encephalopathy"
        self.config["panel_version"] = "1.2"
        self.config["bw"] = "../resources/test/test1.bw"
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

    def test1_1(self):
        """
        Test 1: panel from PanelApp
        :return:
        """
        expected_gene_list = ['MBTPS2']
        self.config["panel"] = "568e844522c1fc1c78b67156"
        self.config["panel_version"] = "1.0"
        self.config["bw"] = "../resources/test/test1.bw"
        self.config["exon_padding"] = 0
        runner = GelCoverageRunner(
            config=self.config
        )
        output, bed = runner.run()
        # Writes the JSON
        with open('../resources/test/sample_output_1_1.json', 'w') as fp:
            json.dump(output, fp)
        # Verifies the bed...
        self.assertEqual(type(bed), pybedtools.bedtool.BedTool)
        # Saves the analysed region as a BED file
        bed.saveas('../resources/test/sample_output_1_1.bed')
        # Runs verifications on output JSON
        #self.verify_output(output, expected_gene_list)

    def test1_2(self):
        """
        Test 1: panel from PanelApp
        :return:
        """
        expected_gene_list = []
        self.config["panel"] = "5550b7bebb5a161bf644a3bc"
        self.config["panel_version"] = "1.2"
        self.config["bw"] = "../resources/test/test1.bw"
        self.config["exon_padding"] = 0
        runner = GelCoverageRunner(
            config=self.config
        )
        output, bed = runner.run()
        # Writes the JSON
        with open('../resources/test/sample_output_1_2.json', 'w') as fp:
            json.dump(output, fp)
        # Verifies the bed...
        self.assertEqual(type(bed), pybedtools.bedtool.BedTool)
        # Saves the analysed region as a BED file
        bed.saveas('../resources/test/sample_output_1_2.bed')
        # Runs verifications on output JSON
        #self.verify_output(output, expected_gene_list)

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
        expected_gene_list = [str(x) for x in self.config["gene_list"].split(",")]
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
        expected_gene_list = ['ADSL', 'ALG13', 'ARHGEF9', 'ARX', 'ATP1A3', 'ATRX', 'CDKL5', 'CHD2',
                              'CNTNAP2', 'DNM1', 'DOCK7', 'DYRK1A', 'EHMT1', 'FOXG1', 'GABRA1', 'GABRB3',
                              'GNAO1', 'GRIN1', 'GRIN2A', 'GRIN2B', 'HCN1', 'IQSEC2', 'KCNA2', 'KCNB1',
                              'KCNJ10', 'KCNQ2', 'KCNQ3', 'KCNT1', 'KIAA2022', 'KIF1BP', 'MAPK10', 'MBD5',
                              'MECP2', 'MEF2C', 'PCDH19', 'PIGA', 'PLCB1', 'PNKP', 'POLG', 'PRRT2', 'PURA',
                              'QARS', 'SCN1A', 'SCN1B', 'SCN2A', 'SCN8A', 'SETD5', 'SIK1', 'SLC12A5',
                              'SLC13A5', 'SLC16A2', 'SLC25A22', 'SLC2A1', 'SLC6A1', 'SLC9A6', 'SPTAN1',
                              'STX1B', 'STXBP1', 'SYNGAP1', 'TCF4', 'UBE2A', 'UBE3A', 'WDR45', 'WWOX', 'ZEB2']
        self.config["panel"] = "Epileptic encephalopathy"
        self.config["panel_version"] = "1.2"
        self.config["bw"] = "../resources/test/test1.bw"
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
        expected_gene_list = [str(x) for x in self.config["gene_list"].split(",")]
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
                constants.EXON_START: start,
                constants.EXON_END: end,
                constants.EXON_PADDED_START: start - self.config["exon_padding"],
                constants.EXON_PADDED_END: end + self.config["exon_padding"],
                constants.EXON: exon_number
            }
        # Interval covered in the input file SCN2A: 2: 165,995,882-166,349,242
        offset = 165995882
        runner = GelCoverageRunner(
            config=self.config
        )
        # Runs union transcript with exon padding 15bp
        self.config["exon_padding"] = 15
        gene_15bp_padding = {
            constants.CHROMOSOME: "chr2",
            constants.GENE_NAME: "TEST",
            constants.TRANSCRIPTS: [
                {
                    constants.TRANSCRIPT_ID: "1",
                    constants.EXONS: [
                        create_test_exon(offset + 10, offset + 15, 1),
                        create_test_exon(offset + 20, offset + 25, 2),
                        create_test_exon(offset + 30, offset + 35, 3),
                        create_test_exon(offset + 40, offset + 45, 4),
                        create_test_exon(offset + 50, offset + 55, 5),
                        create_test_exon(offset + 60, offset + 65, 6)
                    ]
                },
                {
                    constants.TRANSCRIPT_ID: "2",
                    constants.EXONS: [
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
                    constants.TRANSCRIPT_ID: "3",
                    constants.EXONS: [
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
        self.assertEqual(len(union_transcript[constants.EXONS]), 1)
        gene_15bp_padding[constants.UNION_TRANSCRIPT] = union_transcript
        self.verify_union_transcript(gene_15bp_padding, True)
        # Runs union transcript with exon padding 0bp
        self.config["exon_padding"] = 0
        gene_0bp_padding = {
            constants.CHROMOSOME: "chr2",
            constants.GENE_NAME: "TEST",
            constants.TRANSCRIPTS: [
                {
                    constants.TRANSCRIPT_ID: "1",
                    constants.EXONS: [
                        create_test_exon(offset + 10, offset + 15, 1),
                        create_test_exon(offset + 20, offset + 25, 2),
                        create_test_exon(offset + 30, offset + 35, 3),
                        create_test_exon(offset + 40, offset + 45, 4),
                        create_test_exon(offset + 50, offset + 55, 5),
                        create_test_exon(offset + 60, offset + 65, 6)
                    ]
                },
                {
                    constants.TRANSCRIPT_ID: "2",
                    constants.EXONS: [
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
                    constants.TRANSCRIPT_ID: "3",
                    constants.EXONS: [
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
        self.assertEqual(len(union_transcript[constants.EXONS]), 8)
        gene_0bp_padding[constants.UNION_TRANSCRIPT] = union_transcript
        self.verify_union_transcript(gene_0bp_padding, True)

    @unittest.skip("long running test")
    def test6(self):
        """
        Test 6: provided bed of nonN regions and whole genome metrics enabled
        :return:
        """
        expected_gene_list = ['ADSL', 'ALG13', 'ARHGEF9', 'ARX', 'ATP1A3', 'ATRX', 'CDKL5', 'CHD2',
                              'CNTNAP2', 'DNM1', 'DOCK7', 'DYRK1A', 'EHMT1', 'FOXG1', 'GABRA1', 'GABRB3',
                              'GNAO1', 'GRIN1', 'GRIN2A', 'GRIN2B', 'HCN1', 'IQSEC2', 'KCNA2', 'KCNB1',
                              'KCNJ10', 'KCNQ2', 'KCNQ3', 'KCNT1', 'KIAA2022', 'KIF1BP', 'MAPK10', 'MBD5',
                              'MECP2', 'MEF2C', 'PCDH19', 'PIGA', 'PLCB1', 'PNKP', 'POLG', 'PRRT2', 'PURA',
                              'QARS', 'SCN1A', 'SCN1B', 'SCN2A', 'SCN8A', 'SETD5', 'SIK1', 'SLC12A5',
                              'SLC13A5', 'SLC16A2', 'SLC25A22', 'SLC2A1', 'SLC6A1', 'SLC9A6', 'SPTAN1',
                              'STX1B', 'STXBP1', 'SYNGAP1', 'TCF4', 'UBE2A', 'UBE3A', 'WDR45', 'WWOX', 'ZEB2']
        self.config["panel"] = "Epileptic encephalopathy"
        self.config["panel_version"] = "1.2"
        self.config["bw"] = "../resources/test/test1.bw"
        self.config["exon_padding"] = 0
        self.config["wg_stats_enabled"] = True
        self.config["wg_regions"] = \
            "../resources/Homo_sapiens.GRCh37.75.dna.primary_assembly.NonN_Regions.CHR.bed"
        runner = GelCoverageRunner(
            config=self.config
        )
        output, bed = runner.run()
        # Writes the JSON
        with open('../resources/test/sample_output_6.json', 'w') as fp:
            json.dump(output, fp)
        # Verifies the bed...
        self.assertEqual(type(bed), pybedtools.bedtool.BedTool)
        # Saves the analysed region as a BED file
        bed.saveas('../resources/test/sample_output_6.bed')
        # Runs verifications on output JSON
        self.verify_output(output, expected_gene_list)

    #@unittest.skip("long running test")
    def test6_1(self):
        """
        Test 6.1: provided bed of nonN regions and whole genome metrics enabled.
        Different chromosome notations, should fail.
        :return:
        """
        expected_gene_list = ['ADSL', 'ALG13', 'ARHGEF9', 'ARX', 'ATP1A3', 'ATRX', 'CDKL5', 'CHD2',
                              'CNTNAP2', 'DNM1', 'DOCK7', 'DYRK1A', 'EHMT1', 'FOXG1', 'GABRA1', 'GABRB3',
                              'GNAO1', 'GRIN1', 'GRIN2A', 'GRIN2B', 'HCN1', 'IQSEC2', 'KCNA2', 'KCNB1',
                              'KCNJ10', 'KCNQ2', 'KCNQ3', 'KCNT1', 'KIAA2022', 'KIF1BP', 'MAPK10', 'MBD5',
                              'MECP2', 'MEF2C', 'PCDH19', 'PIGA', 'PLCB1', 'PNKP', 'POLG', 'PRRT2', 'PURA',
                              'QARS', 'SCN1A', 'SCN1B', 'SCN2A', 'SCN8A', 'SETD5', 'SIK1', 'SLC12A5',
                              'SLC13A5', 'SLC16A2', 'SLC25A22', 'SLC2A1', 'SLC6A1', 'SLC9A6', 'SPTAN1',
                              'STX1B', 'STXBP1', 'SYNGAP1', 'TCF4', 'UBE2A', 'UBE3A', 'WDR45', 'WWOX', 'ZEB2']
        self.config["panel"] = "Epileptic encephalopathy"
        self.config["panel_version"] = "1.2"
        self.config["bw"] = "../resources/test/test1.bw"
        self.config["exon_padding"] = 0
        self.config["wg_stats_enabled"] = True
        self.config["wg_regions"] = \
            "../resources/Homo_sapiens.GRCh37.75.dna.primary_assembly.NonN_Regions.CHR.prefix.bed"
        try:
            runner = GelCoverageRunner(
                config=self.config
            )
            caught_exception = False
        except GelCoverageInputError:
            caught_exception = True
        self.assertTrue(caught_exception, msg="It should have raised an exception as bed and bigwig are using"
                                              "different chromosome notations")

    def test6_2(self):
        """
        Test 6.2: provided a small bed for the whole genome analysis
        :return:
        """
        self.config["panel"] = None
        self.config["panel_version"] = None
        self.config["bw"] = "../resources/test/test1.bw"
        self.config["exon_padding"] = 0
        self.config["coding_region_stats_enabled"] = False
        self.config["wg_stats_enabled"] = True
        self.config["wg_regions"] = "../resources/test/test1.bed"
        runner = GelCoverageRunner(
            config=self.config
        )
        output, _ = runner.run()
        # Writes the JSON
        with open('../resources/test/sample_output_6_2.json', 'w') as fp:
            json.dump(output, fp)
        # Runs verifications on output JSON
        self.verify_output(output, None)

    @unittest.skip("long running test")
    def test7(self):
        """
        Test 7: whole genome metrics enabled
        :return:
        """
        expected_gene_list = ['ADSL', 'ALG13', 'ARHGEF9', 'ARX', 'ATP1A3', 'ATRX', 'CDKL5', 'CHD2',
                              'CNTNAP2', 'DNM1', 'DOCK7', 'DYRK1A', 'EHMT1', 'FOXG1', 'GABRA1', 'GABRB3',
                              'GNAO1', 'GRIN1', 'GRIN2A', 'GRIN2B', 'HCN1', 'IQSEC2', 'KCNA2', 'KCNB1',
                              'KCNJ10', 'KCNQ2', 'KCNQ3', 'KCNT1', 'KIAA2022', 'KIF1BP', 'MAPK10', 'MBD5',
                              'MECP2', 'MEF2C', 'PCDH19', 'PIGA', 'PLCB1', 'PNKP', 'POLG', 'PRRT2', 'PURA',
                              'QARS', 'SCN1A', 'SCN1B', 'SCN2A', 'SCN8A', 'SETD5', 'SIK1', 'SLC12A5',
                              'SLC13A5', 'SLC16A2', 'SLC25A22', 'SLC2A1', 'SLC6A1', 'SLC9A6', 'SPTAN1',
                              'STX1B', 'STXBP1', 'SYNGAP1', 'TCF4', 'UBE2A', 'UBE3A', 'WDR45', 'WWOX', 'ZEB2']
        self.config["panel"] = "Epileptic encephalopathy"
        self.config["panel_version"] = "1.2"
        self.config["bw"] = "../resources/test/test1.bw"
        self.config["exon_padding"] = 0
        self.config["wg_stats_enabled"] = True
        runner = GelCoverageRunner(
            config=self.config
        )
        output, bed = runner.run()
        # Writes the JSON
        with open('../resources/test/sample_output_7.json', 'w') as fp:
            json.dump(output, fp)
        # Verifies the bed...
        self.assertEqual(type(bed), pybedtools.bedtool.BedTool)
        # Saves the analysed region as a BED file
        bed.saveas('../resources/test/sample_output_7.bed')
        # Runs verifications on output JSON
        self.verify_output(output, expected_gene_list)

    def test8(self):
        """
        Test 8: exon stats disabled
        :return:
        """
        expected_gene_list = ['ADSL', 'ALG13', 'ARHGEF9', 'ARX', 'ATP1A3', 'ATRX', 'CDKL5', 'CHD2',
                              'CNTNAP2', 'DNM1', 'DOCK7', 'DYRK1A', 'EHMT1', 'FOXG1', 'GABRA1', 'GABRB3',
                              'GNAO1', 'GRIN1', 'GRIN2A', 'GRIN2B', 'HCN1', 'IQSEC2', 'KCNA2', 'KCNB1',
                              'KCNJ10', 'KCNQ2', 'KCNQ3', 'KCNT1', 'KIAA2022', 'KIF1BP', 'MAPK10', 'MBD5',
                              'MECP2', 'MEF2C', 'PCDH19', 'PIGA', 'PLCB1', 'PNKP', 'POLG', 'PRRT2', 'PURA',
                              'QARS', 'SCN1A', 'SCN1B', 'SCN2A', 'SCN8A', 'SETD5', 'SIK1', 'SLC12A5',
                              'SLC13A5', 'SLC16A2', 'SLC25A22', 'SLC2A1', 'SLC6A1', 'SLC9A6', 'SPTAN1',
                              'STX1B', 'STXBP1', 'SYNGAP1', 'TCF4', 'UBE2A', 'UBE3A', 'WDR45', 'WWOX', 'ZEB2']
        self.config["panel"] = "Epileptic encephalopathy"
        self.config["panel_version"] = "1.2"
        self.config["bw"] = "../resources/test/test1.bw"
        self.config["exon_padding"] = 0
        self.config["exon_stats_enabled"] = False
        runner = GelCoverageRunner(
            config=self.config
        )
        output, bed = runner.run()
        # Writes the JSON
        with open('../resources/test/sample_output_8.json', 'w') as fp:
            json.dump(output, fp)
        # Verifies the bed...
        self.assertEqual(type(bed), pybedtools.bedtool.BedTool)
        # Saves the analysed region as a BED file
        bed.saveas('../resources/test/sample_output_8.bed')
        # Runs verifications on output JSON
        self.verify_output(output, expected_gene_list)

    def test9(self):
        """
        Test 9: gene list containing a gene not covered by the BAM
        :return:
        """
        self.config["panel"] = None
        self.config["panel_version"] = None
        self.config["bw"] = "../resources/test/test2.bw"
        self.config["gene_list"] = "BRCA1,BRCA2,CFTR,IGHE,PTEN"
        self.config["exon_padding"] = 0
        expected_gene_list = [str(x) for x in self.config["gene_list"].split(",")]
        runner = GelCoverageRunner(
            config=self.config
        )
        output, bed = runner.run()
        # Writes the JSON
        with open('../resources/test/sample_output_9.json', 'w') as fp:
            json.dump(output, fp)
        # Verifies the bed...
        self.assertEqual(type(bed), pybedtools.bedtool.BedTool)
        # Saves the analysed region as a BED file
        bed.saveas('../resources/test/sample_output_9.bed')
        # Runs verifications on output JSON
        self.verify_output(output, expected_gene_list)
        self.assertEqual(len(output["results"]["uncovered_genes"]), 1,
                         msg="Uncovered genes should be of length 1")
        self.assertEqual(output["results"]["uncovered_genes"][0][constants.GENE_NAME], "PTEN")

    def test10(self):
        """
        Test 10: gene list containing only a gene not covered by the BAM
        :return:
        """
        self.config["panel"] = None
        self.config["panel_version"] = None
        self.config["bw"] = "../resources/test/test2.bw"
        self.config["gene_list"] = "PTEN"
        self.config["exon_padding"] = 0
        expected_gene_list = [str(x) for x in self.config["gene_list"].split(",")]
        runner = GelCoverageRunner(
            config=self.config
        )
        output, bed = runner.run()
        # Writes the JSON
        with open('../resources/test/sample_output_10.json', 'w') as fp:
            json.dump(output, fp)
        # Verifies the bed...
        self.assertEqual(type(bed), pybedtools.bedtool.BedTool)
        # Saves the analysed region as a BED file
        bed.saveas('../resources/test/sample_output_10.bed')
        # Runs verifications on output JSON
        self.assertEqual(len(output["results"]["uncovered_genes"]), 1,
                         msg="Uncovered genes should be of length 1")
        self.assertEqual(output["results"]["uncovered_genes"][0][constants.GENE_NAME], "PTEN")

    @unittest.skip("long running test")
    def test11(self):
        """
        Test 11: coding region analysis disabled
        :return:
        """
        expected_gene_list = ['ADSL', 'ALG13', 'ARHGEF9', 'ARX', 'ATP1A3', 'ATRX', 'CDKL5', 'CHD2',
                              'CNTNAP2', 'DNM1', 'DOCK7', 'DYRK1A', 'EHMT1', 'FOXG1', 'GABRA1', 'GABRB3',
                              'GNAO1', 'GRIN1', 'GRIN2A', 'GRIN2B', 'HCN1', 'IQSEC2', 'KCNA2', 'KCNB1',
                              'KCNJ10', 'KCNQ2', 'KCNQ3', 'KCNT1', 'KIAA2022', 'KIF1BP', 'MAPK10', 'MBD5',
                              'MECP2', 'MEF2C', 'PCDH19', 'PIGA', 'PLCB1', 'PNKP', 'POLG', 'PRRT2', 'PURA',
                              'QARS', 'SCN1A', 'SCN1B', 'SCN2A', 'SCN8A', 'SETD5', 'SIK1', 'SLC12A5',
                              'SLC13A5', 'SLC16A2', 'SLC25A22', 'SLC2A1', 'SLC6A1', 'SLC9A6', 'SPTAN1',
                              'STX1B', 'STXBP1', 'SYNGAP1', 'TCF4', 'UBE2A', 'UBE3A', 'WDR45', 'WWOX', 'ZEB2']
        self.config["panel"] = "Epileptic encephalopathy"
        self.config["panel_version"] = "1.2"
        self.config["bw"] = "../resources/test/test1.bw"
        self.config["exon_padding"] = 0
        self.config["wg_stats_enabled"] = True
        self.config["wg_regions"] = \
            "../resources/Homo_sapiens.GRCh37.75.dna.primary_assembly.NonN_Regions.CHR.bed"
        self.config["coding_region_stats_enabled"] = False
        runner = GelCoverageRunner(
            config=self.config
        )
        output, bed = runner.run()
        # Writes the JSON
        with open('../resources/test/sample_output_11.json', 'w') as fp:
            json.dump(output, fp)
        # Verifies the bed...
        self.assertEqual(bed, None)
        # Runs verifications on output JSON
        self.verify_output(output, expected_gene_list)

    def test12(self):
        """
        Test 12: coding region and whole genome analysis disabled
        :return:
        """
        expected_gene_list = ['ADSL', 'ALG13', 'ARHGEF9', 'ARX', 'ATP1A3', 'ATRX', 'CDKL5', 'CHD2',
                              'CNTNAP2', 'DNM1', 'DOCK7', 'DYRK1A', 'EHMT1', 'FOXG1', 'GABRA1', 'GABRB3',
                              'GNAO1', 'GRIN1', 'GRIN2A', 'GRIN2B', 'HCN1', 'IQSEC2', 'KCNA2', 'KCNB1',
                              'KCNJ10', 'KCNQ2', 'KCNQ3', 'KCNT1', 'KIAA2022', 'KIF1BP', 'MAPK10', 'MBD5',
                              'MECP2', 'MEF2C', 'PCDH19', 'PIGA', 'PLCB1', 'PNKP', 'POLG', 'PRRT2', 'PURA',
                              'QARS', 'SCN1A', 'SCN1B', 'SCN2A', 'SCN8A', 'SETD5', 'SIK1', 'SLC12A5',
                              'SLC13A5', 'SLC16A2', 'SLC25A22', 'SLC2A1', 'SLC6A1', 'SLC9A6', 'SPTAN1',
                              'STX1B', 'STXBP1', 'SYNGAP1', 'TCF4', 'UBE2A', 'UBE3A', 'WDR45', 'WWOX', 'ZEB2']
        self.config["panel"] = "Epileptic encephalopathy"
        self.config["panel_version"] = "1.2"
        self.config["bw"] = "../resources/test/test1.bw"
        self.config["exon_padding"] = 0
        self.config["wg_stats_enabled"] = False
        self.config["wg_regions"] = \
            "../resources/Homo_sapiens.GRCh37.75.dna.primary_assembly.NonN_Regions.CHR.prefix.bed"
        self.config["coding_region_stats_enabled"] = False
        runner = GelCoverageRunner(
            config=self.config
        )
        output, bed = runner.run()
        # Writes the JSON
        with open('../resources/test/sample_output_12.json', 'w') as fp:
            json.dump(output, fp)
        # Verifies the bed...
        self.assertEqual(bed, None)
        # Runs verifications on output JSON
        self.verify_output(output, expected_gene_list)

    @unittest.skip("long running test")
    def test13(self):
        """
        Test 13: panel from PanelApp with exon padding of 15 bp, the biggest panel with 1232 genes
        :return:
        """
        expected_gene_list = None  # too big to set here
        self.config["bw"] = "../resources/test/test1.bw"
        self.config["panel"] = "Intellectual disability"
        self.config["panel_version"] = "1.23"
        self.config["exon_padding"] = 15
        runner = GelCoverageRunner(
            config=self.config
        )
        output, bed = runner.run()
        # Writes the JSON
        with open('../resources/test/sample_output_13.json', 'w') as fp:
            json.dump(output, fp)
        # Verifies the bed...
        self.assertEqual(type(bed), pybedtools.bedtool.BedTool)
        # Saves the analysed region as a BED file
        bed.saveas('../resources/test/sample_output_13.bed')
        # Runs verifications on output JSON
        self.verify_output(output, expected_gene_list)

    def test14(self):
        """
        Test panel with no transcript passing filters.
        :return:
        """
        expected_gene_list = None  # too big to set here
        self.config["bw"] = "../resources/test/test1.bw"
        self.config["panel"] = "5763f2ea8f620350a1996048"
        self.config["panel_version"] = "1.0"
        self.config["exon_padding"] = 15
        runner = GelCoverageRunner(
            config=self.config
        )
        output, bed = runner.run()
        # Writes the JSON
        with open('../resources/test/sample_output_14.json', 'w') as fp:
            json.dump(output, fp)
        # Verifies the bed...
        self.assertEqual(type(bed), pybedtools.bedtool.BedTool)
        # Saves the analysed region as a BED file
        bed.saveas('../resources/test/sample_output_14.bed')
        # Runs verifications on output JSON
        self.expected_gene_list = expected_gene_list
        self.assertEqual(type(output), dict)
        # Verify that content in parameters is correct
        self._verify_dict_field(output, "parameters", dict)
        self._verify_parameters(output["parameters"])
