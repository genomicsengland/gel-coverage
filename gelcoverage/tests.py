import unittest
from gelcoverage.runner import GelCoverageRunner


class GelCoverageRunnerTests(unittest.TestCase):

    def setUp(self):
        self.config = {
            # Sets parameters from CLI
            "bw": "/any/path/file.bw",  # TODO: find a testing file
            "panel": "Epileptic encephalopathy",
            "panel_version": "0.2",
            #"gene_list": args.gene_list,
            "coverage_threshold": 15,
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
        self.runner = GelCoverageRunner(
            config= self.config
        )

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
