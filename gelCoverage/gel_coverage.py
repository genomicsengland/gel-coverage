import pyBigWig
from collections import defaultdict

import gelCoverage.stats.coverage_stats as coverage_stats
from gelCoverage.tools.cellbase_helper import CellbaseHelper
from gelCoverage.tools.panelapp_helper import PanelappHelper


class GelCoverageEngine:

    def __init__(self, config):
        self.config = config
        # Initialize CellBase helper
        self.cellbase_helper = CellbaseHelper(
            species=config['cellbase_species'],
            version=config['cellbase_version'],
            assembly=config['cellbase_assembly'],
            host=config['cellbase_host'],
            filter_flags=config['transcript_filtering_flags'].split(","),
            filter_biotypes=config['transcript_filtering_biotypes'].split(",")
        )
        # Initialize PanelApp helper
        self.panelapp_helper = PanelappHelper(host=config['panelapp_host'])
        # Gets the list of genes to analyse
        self.gene_list = self.get_gene_list()
        # Opens the bigwig file for reading
        self.bigwig = pyBigWig.open(self.config['bw'])

    def get_gene_list(self):
        """
        Retrieves the gene list to analyse based on the input parameters.
        This list might be retrieved from PanelApp, from the parameter genes or
        just get every gene in the human genome.
        When all parameters are provided, genes are retrieved from PanelApp and
        the parameter genes is ignored.
        :return: the gene list for coverage analysis
        """
        ## Gets list of genes to analyse
        if self.config["panel"] is not None and self.config["panel_version"] is not None:
            # Get list of genes from PanelApp
            gene_list = self.panelapp_helper.get_gene_list(
                self.config["panel"],
                self.config["panel_version"],
                self.config["panelapp_gene_confidence"]
            )
        elif self.config["gene_list"] is not None:
            # Get list of genes from parameter genes
            gene_list = self.config["gene_list"].split(",")
            # elif args.transcripts is not None:
            #    transcripts_list = args.transcripts.split(",")
            # Get list of genes from CellBase client
            #    gene_list = cellbase_helper.get_gene_list_from_transcripts(transcripts_list)
        else:
            # Warn the user as this will be time consuming
            print "WARNING: you are going to run a whole exome coverage analysis!"
            # Retrieve the list of all genes
            gene_list = self.cellbase_helper.get_all_genes()
        return gene_list

    def get_parameters_output(self):
        """

        :return:
        """
        parameters = defaultdict()
        parameters["gap_coverage_threshold"] = self.config["coverage_threshold"]
        parameters["input_file"] = self.config["bw"]
        parameters["species"] = self.config['cellbase_species']
        parameters["assembly"] = self.config['cellbase_assembly']
        if self.config['panelapp_panel'] is not None and self.config['panelapp_panel_version'] is not None:
            parameters["panel"] = self.config['panelapp_panel']
            parameters["panel_version"] = self.config['panelapp_panel_version']
            parameters["gene_list"] = self.gene_list
        elif self.config['gene_list'] is not None:
            parameters["gene_list"] = self.gene_list
        # Beware that when performing analysis on the whole exome the gene list field is
        # not set. We don't want a list of 20K genes in here. Do we?
        return parameters

    def run(self):
        """

        :return:
        """
        # Get genes annotations in BED format
        bed = self.cellbase_helper.make_exons_bed(self.gene_list)
        # Initialize and fill results data structure
        output = defaultdict()
        output["parameters"] = self.get_parameters_output()
        output["results"] = defaultdict()
        results = output["results"]
        # TODO: adapt results to the new output schema
        for interval in bed:
            # Reads data from BED entry
            chrom = interval.chrom
            start = int(interval.start)
            end = int(interval.end)
            gene, txid, exon_idx = interval.name.split("|")
            strand = interval.strand
            # TODO: truncate this to two decimal positions
            gc_content = interval.score
            # Queries the bigwig for a specific interval
            if start == end:
                end += 1
            # Read from the bigwig file
            coverages = self.bigwig.values(chrom, start, end)
            # Compute statistics at exon level
            exon_statistics = coverage_stats.compute_exon_level_statistics(coverages, gc_content)
            # Store results in data structure
            if txid not in output:
                output[txid] = defaultdict()
                output[txid]["chrom"] = chrom
                output[txid]["gene"] = gene
                output[txid]["strand"] = strand
                output[txid]["exons"] = defaultdict()
            output[gene][txid]["exons"][exon_idx] = {
                "statistics": exon_statistics
            }
            # compute gaps
            if self.config['coverage_threshold'] > 0:
                gaps = coverage_stats.find_gaps(coverages, start, self.config['coverage_threshold'])
                output[txid]["exons"][exon_idx]["gaps"] = gaps
        # add an aggregation of statistics at transcript level
        for transcript, content in output.iteritems():
            output[transcript]["statistics"] = coverage_stats.compute_transcript_level_statistics(content["exons"])
        return output