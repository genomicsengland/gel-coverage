import json
import urllib2


class PanelappHelper:

    def __init__(self, host):
        self.host = host

    def get_gene_list(self, panel, panel_version, gene_confidence_threshold):
        # TODO: implement a call to the PanelApp webservice
        url = "{host}/get_panel/{panel}/".format(
            host=self.host,
            panel=panel)

        parameters = "?version={version}&LevelOfConfidence={confidence}".format(
            version=panel_version,
            confidence=",".join(gene_confidence_threshold) if type(gene_confidence_threshold) == list
            else gene_confidence_threshold)
        url = urllib2.quote(url) + parameters  # we don't want parameters quoted
        panel = json.load(urllib2.urlopen("https://" + url))
        # TODO: refine error management
        if type(panel) != dict:
            raise SystemError("PanelApp returned an error for the query %s" % url)
        gene_list = [x["GeneSymbol"] for x in panel["result"]["Genes"]]

        return gene_list
