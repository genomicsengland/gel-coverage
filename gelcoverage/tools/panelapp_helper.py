import ujson
import urllib2
import logging

import gelcoverage.tools.backoff_retrier as backoff_retrier


class PanelappHelper:

    def __init__(self, host, retries, assembly):
        self.assembly = assembly
        self.host = host
        self.retries = retries
        self.urlopen = backoff_retrier.wrapper(urllib2.urlopen, self.retries)

    def get_gene_list(self, panel, panel_version, gene_confidence_threshold):
        """
        Gets the HGNC gene names in a given panel identified by panel name and panel version.
        Also, if provided, only gets those genes having a level of evidence listed in gene_confidence_threshold.
        :param panel: the panel name
        :param panel_version: the panel version
        :param gene_confidence_threshold: the gene's level of evidence
        :return: a list of HGNC gene names
        """
        logging.debug("Getting gene list from PanelApp...")
        url = "{host}/get_panel/{panel}/".format(
            host=self.host,
            panel=panel)

        parameters = "?version={version}&LevelOfConfidence={confidence}&assembly={assembly}".format(
            version=panel_version,
            confidence=",".join(gene_confidence_threshold) if type(gene_confidence_threshold) == list
            else gene_confidence_threshold,
            assembly=self.assembly
        )
        url = urllib2.quote(url) + parameters  # we don't want parameters quoted
        panel = ujson.load(self.urlopen("https://" + url))
        # TODO: refine error management
        if type(panel) != dict:
            raise SystemError("PanelApp returned an error for the query %s" % url)
        gene_list = [x["GeneSymbol"] for x in panel["result"]["Genes"]]
        logging.debug("Gene list obtained from PanelApp of %s!" % str(len(gene_list)))
        return gene_list
