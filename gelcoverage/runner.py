import logging
import pyBigWig

import gelcoverage.stats.coverage_stats as coverage_stats
from gelcoverage.tools.cellbase_helper import CellbaseHelper
from gelcoverage.tools.panelapp_helper import PanelappHelper


class GelCoverageRunner:

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
        logging.info("Gene list to analyse: %s" % ",".join(self.gene_list))
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
        # Gets list of genes to analyse
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
            logging.warning("You are about to run a whole exome coverage analysis!")
            # Retrieve the list of all genes
            gene_list = self.cellbase_helper.get_all_gene_names()
        return gene_list

    def __get_parameters_output(self):
        """
        Returns a dictionary with the formatted configuration parameters used to run the analysis.
        :return: dict
        """
        parameters = {
            "gap_coverage_threshold": self.config["coverage_threshold"],
            "input_file": self.config["bw"],
            "species": self.config['cellbase_species'],
            "assembly": self.config['cellbase_assembly'],
            "transcript_filtering_flags": self.config['transcript_filtering_flags'],
            "transcript_filtering_biotypes": self.config['transcript_filtering_biotypes']
        }
        if 'panel' in self.config and self.config['panel'] is not None \
                and 'panel_version' in self.config and self.config['panel_version'] is not None:
            parameters["panel"] = self.config['panel']
            parameters["panel_version"] = self.config['panel_version']
            parameters["panel_gene_confidence"] = self.config['panelapp_gene_confidence']
            parameters["gene_list"] = self.gene_list
        elif 'gene_list' in self.config and self.config['gene_list'] is not None:
            parameters["gene_list"] = self.gene_list
        # Beware that when performing analysis on the whole exome the gene list field is
        # not set. We don't want a list of 20K genes in here. Do we?
        return parameters

    @staticmethod
    def __initialize_gene_dict(gene_name, chromosome):
        """
        Returns the dictionary that will store gene information
        :param gene_name: the gene name
        :param chromosome: the chromosome
        :return: the dictionary
        """
        return {
            "name": gene_name,
            "chromosome": chromosome,
            "transcripts": []
        }

    @staticmethod
    def __initialize_transcript_dict(transcript_id):
        """
        Returns the dictionary that will store transcript information
        :param transcript_id: the transcript id
        :return: the dictionary
        """
        return {
            "id": transcript_id,
            "exons": []
        }

    @staticmethod
    def __initialize_exon_dict(exon_number, start, end):
        """
        Returns the dictionary that will store exon information
        :param exon_number: the exon number
        :param start: the start position
        :param end: the end position
        :return:
        """
        return {
                "exon_number": exon_number,
                "start": start,
                "end": end,
                "length": end - start + 1
            }

    @staticmethod
    def __read_bed_interval(interval):
        """
        Extracts information from a BED interval.
        :param interval: the BED interval
        :return: the extracted information in separate variables
        """
        chromosome = interval.chrom
        start = int(interval.start)
        end = int(interval.end)
        gene_name, transcript_id, exon_number = interval.name.split("|")
        strand = interval.strand
        gc_content = float(interval.score)
        return chromosome, start, end, gene_name, transcript_id, exon_number, strand, gc_content

    def __read_bigwig_coverages(self, chromosome, start, end):
        """
        Reads the coverage values in a region from the bigwig file
        :param chromosome: the region chromosome
        :param start: the start position
        :param end: the end position
        :return: the sequence of coverages (one integer per position)
        """
        # Queries the bigwig for a specific interval
        if start == end:  # do we really need this?
            end += 1
        # Read from the bigwig file
        try:
            # TODO: why our bigwig has "chr" prefix? BAMs don't
            coverages = self.bigwig.values("chr" + str(chromosome), start, end)
        except RuntimeError, e:
            # TODO: deal with this errors
            logging.error("Querying for unexisting interval %s:%s-%s" % (chromosome, start, end))
            raise e
        return coverages

    def __process_bed_file(self, bed):
        """
        Reads every interval defined in the BED file, gets the coverage information from the BigWig data, computes
        coverage statistics and formats the information in a JSON-friendly data structure.
        :param bed: the BED file handler
        :return: the data structure containing all statistics
        """
        # Checks that the BED file contains the expected information
        if bed is None or len(bed) == 0:
            logging.error("Incorrect BED file!")
            raise RuntimeError("Incorrect BED file!")
        # Initialize results data structure
        results = []
        current_gene = {}
        current_transcript = {}
        # Process every interval in the BED file
        for interval in bed:
            # Reads data from BED entry
            chromosome, start, end, gene_name, transcript_id, \
                exon_number, strand, gc_content = GelCoverageRunner.__read_bed_interval(interval)
            # Store transcript in data structure. Needs to run before gene
            try:
                if current_transcript["id"] != transcript_id:
                    # Compute transcript level statistics by aggregating stats on every exon
                    current_transcript["statistics"] = coverage_stats.compute_transcript_level_statistics(
                        current_transcript["exons"]
                    )
                    # Save previous result
                    current_gene["transcripts"].append(current_transcript)
                    logging.info("Processed transcript %s of gene %s" % (transcript_id, gene_name))
                    # Create a new data structure for new gene
                    current_transcript = GelCoverageRunner.__initialize_transcript_dict(transcript_id)
            except KeyError:
                current_transcript = GelCoverageRunner.__initialize_transcript_dict(transcript_id)
            # Store gene in data structure
            try:
                if current_gene["name"] != gene_name:
                    # Save previous result
                    results.append(current_gene)
                    # Create a new data structure for new gene
                    current_gene = GelCoverageRunner.__initialize_gene_dict(gene_name, chromosome)
            except KeyError:
                current_gene = GelCoverageRunner.__initialize_gene_dict(gene_name, chromosome)
            # Store exon in data structure
            exon = GelCoverageRunner.__initialize_exon_dict(exon_number, start, end)
            # Read from the bigwig file
            coverages = self.__read_bigwig_coverages(chromosome, start, end)
            # Compute statistics at exon level
            exon["statistics"] = coverage_stats.compute_exon_level_statistics(coverages, gc_content)
            # compute gaps
            if self.config['coverage_threshold'] > 0:
                exon["gaps"] = coverage_stats.find_gaps(coverages, start, self.config['coverage_threshold'])
            current_transcript["exons"].append(exon)
        # Adds last ocurrences
        # Compute transcript level statistics by aggregating stats on every exon
        current_transcript["statistics"] = coverage_stats.compute_transcript_level_statistics(
            current_transcript["exons"]
        )
        # Save previous result
        current_gene["transcripts"].append(current_transcript)
        logging.info("Processed transcript %s of gene %s" % (transcript_id, gene_name))
        results.append(current_gene)
        return results

    def __output(self, results):
        """
        Builds the output data structure
        :param results:
        :return:
        """
        return {
            "parameters": self.__get_parameters_output(),
            "results": results
        }

    def run(self):
        """
        PRE: we assume that the BED generated is sorted by transcript and then by gene position
        this assumption is true for BEDs generated with make_exons_bed() but might not be
        always true.
        :return:
        """
        logging.info("Starting coverage analysis")
        # Get genes annotations in BED format
        bed = self.cellbase_helper.make_exons_bed(self.gene_list)
        # Process the intervals in the BED file
        results = self.__process_bed_file(bed)
        return self.__output(results)
