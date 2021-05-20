import logging
from pypanelapp.python_panel_app_client import PanelApp


class PanelappHelper:
    def __init__(self, host, assembly):
        self.assembly = assembly
        self.host = host
        self.client = PanelApp(server=self.host)

    def get_gene_list(self, panel, panel_version, gene_confidence_threshold=3):
        """
        Gets the HGNC gene names in a given panel identified by panel name and panel version.
        Also, if provided, only gets those genes having a level of evidence listed in gene_confidence_threshold.
        :param panel: the panel name
        :param panel_version: the panel version
        :param gene_confidence_threshold: the gene's level of evidence
        :return: a list of HGNC gene names
        """
        logging.debug("Getting gene list from PanelApp...")
        try:
            panel = self.client.panel_get(panel_id=panel, version=panel_version)
        except:
            raise SystemError(
                "PanelApp returned an error for the panel {}, version {}".format(
                    panel, panel_version
                )
            )
        gene_list = [
            gene.gene_data.gene_symbol
            for gene in panel.genes
            if int(gene.confidence_level) >= int(gene_confidence_threshold)
        ]
        logging.debug("Gene list obtained from PanelApp of %s!" % str(len(gene_list)))
        return gene_list
