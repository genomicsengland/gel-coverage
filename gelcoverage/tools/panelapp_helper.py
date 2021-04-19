import logging

from requests.adapters import RetryError
from pypanelapp.python_panel_app_client import PanelApp, PanelAppAPIException


class PanelappHelper:

    def __init__(self, server, assembly):
        self.assembly = assembly
        self.panelapp_client = PanelApp(server=server)

    def get_gene_list(self, panel, panel_version, gene_confidence_threshold=None):
        """
        Gets the HGNC gene names in a given panel identified by panel name and panel version.
        Also, if provided, only gets those genes having minimum evidence gene_confidence_threshold

        :param panel: the panel name
        :param panel_version: the panel version
        :param str gene_confidence_threshold: the gene's level of evidence
        :return: a list of HGNC gene names
        """

        try:
            panel_data = self.panelapp_client.panel_get(
                assembly=self.assembly,
                panel_id=panel,
                version=panel_version
            )
        except (PanelAppAPIException, RetryError) as pae:
            logging.error(
                "Encountered an error querying the PanelApp API for {}:{}".format(panel, panel_version),
                exc_info=True
            )
            raise pae

        if gene_confidence_threshold:
            gene_list = [
                gene.entity_name
                for gene in panel_data.genes
                if gene.confidence_level >= gene_confidence_threshold
            ]
        else:
            gene_list = [gene.entity_name for gene in panel_data.genes]

        logging.debug(
            "{} Genes, Confidence {}, Panel {}, Version {}".format(
                len(gene_list),
                gene_confidence_threshold,
                panel,
                panel_version
            )
        )

        return gene_list
