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
        # Flag indicating if exon padding is enabled
        if "exon_padding" not in self.config or type(self.config["exon_padding"]) != int or \
            self.config["exon_padding"] <= 0:
            self.config["exon_padding"] = 0
        self.is_exon_padding_enabled = self.config["exon_padding"] > 0

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
            "transcript_filtering_biotypes": self.config['transcript_filtering_biotypes'],
            "exon_padding": self.config["exon_padding"]
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
    def __initialize_exon_dict(exon_number, start, end, padded_start, padded_end):
        """
        Returns the dictionary that will store exon information
        :param exon_number: the exon number
        :param start: the start position
        :param end: the end position
        :return:
        """
        exon = {
                "exon_number": exon_number,
                "start": start,
                "end": end,
                "length": end - start + 1
            }
        if padded_start != start:
            exon["padded_start"] = padded_start
            exon["padded_end"] = padded_end
        return exon

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

    def __get_union_transcript(self, gene):
        all_exons = sum([transcript["exons"] for transcript in gene["transcripts"]], [])
        all_exons.sort(key = lambda x: x["start"])
        is_exon_padding = self.config["exon_padding"] > 0
        union_exons = []
        current_exon = GelCoverageRunner.__initialize_exon_dict(
            1,
            all_exons[0]["start"],
            all_exons[0]["end"],
            all_exons[0]["padded_start"] if is_exon_padding else None,
            all_exons[0]["padded_end"] if is_exon_padding else None
        )
        exon_number = 1
        for exon in all_exons[1:]:
            start = exon["padded_start"] if is_exon_padding else exon["start"]
            end = exon["padded_end"] if is_exon_padding else exon["end"]
            if start <= end:
                # Exons overlaps, we join them
                current_exon["padded_end"] = max(
                    current_exon["padded_end"],
                    exon["padded_end"]
                )
                current_exon["end"] = max(
                    current_exon["end"],
                    exon["end"]
                )
            else:
                # Exons do not overlap, they are different exons in the union transcript
                # Update length
                current_exon["length"] = current_exon["end"] - current_exon["start"] + 1
                # Read from the bigwig file
                coverages = self.__read_bigwig_coverages(
                    gene["chromosome"],
                    current_exon["padded_start"],
                    current_exon["padded_end"])
                # Compute statistics at exon level (no GC content information)
                current_exon["statistics"] = coverage_stats.compute_exon_level_statistics(coverages, None)
                # Compute gaps
                if self.config['coverage_threshold'] > 0:
                    current_exon["gaps"] = coverage_stats.find_gaps(
                        coverages,
                        current_exon["padded_start"],
                        self.config['coverage_threshold']
                    )
                # Saves current exon and stores the next
                union_exons.append(current_exon)
                exon_number += 1
                exon["exon_number"] = "exon%s" % str(exon_number)
                current_exon = exon
        # Stores the last exon
        # Update length
        current_exon["length"] = current_exon["end"] - current_exon["start"] + 1
        # Read from the bigwig file
        coverages = self.__read_bigwig_coverages(
            gene["chromosome"],
            current_exon["padded_start"],
            current_exon["padded_end"])
        # Compute statistics at exon level (no GC content information)
        current_exon["statistics"] = coverage_stats.compute_exon_level_statistics(coverages, None)
        # Compute gaps
        if self.config['coverage_threshold'] > 0:
            current_exon["gaps"] = coverage_stats.find_gaps(
                coverages,
                current_exon["padded_start"],
                self.config['coverage_threshold']
            )
        # Saves current exon and stores the next
        union_exons.append(current_exon)
        # Create union transcript dict
        union_transcript = {
            "exons" : union_exons,
            "statistics": coverage_stats.compute_transcript_level_statistics(union_exons)
        }

        return union_transcript

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
        results = {
            "genes": []
        }
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
                    # Obtain the union transcript and compute coverage statistics
                    current_gene["union_transcript"] = self.__get_union_transcript(current_gene)
                    # Save previous result
                    results["genes"].append(current_gene)
                    # Create a new data structure for new gene
                    current_gene = GelCoverageRunner.__initialize_gene_dict(gene_name, chromosome)
            except KeyError:
                current_gene = GelCoverageRunner.__initialize_gene_dict(gene_name, chromosome)
            # Store exon in data structure
            padded_start = max(0, start - self.config["exon_padding"])
            padded_end = end + self.config["exon_padding"]  # TODO: potential problem here overflowing chromosome
            exon = GelCoverageRunner.__initialize_exon_dict(exon_number, start, end, padded_start, padded_end)
            # Read from the bigwig file
            coverages = self.__read_bigwig_coverages(chromosome, padded_start, padded_end)
            # Compute statistics at exon level
            exon["statistics"] = coverage_stats.compute_exon_level_statistics(coverages, gc_content)
            # compute gaps
            if self.config['coverage_threshold'] > 0:
                exon["gaps"] = coverage_stats.find_gaps(coverages, padded_start, self.config['coverage_threshold'])
            current_transcript["exons"].append(exon)
        # Adds last ocurrences
        # Compute transcript level statistics by aggregating stats on every exon
        current_transcript["statistics"] = coverage_stats.compute_transcript_level_statistics(
            current_transcript["exons"]
        )
        # Save previous result
        current_gene["transcripts"].append(current_transcript)
        logging.info("Processed transcript %s of gene %s" % (transcript_id, gene_name))
        # Obtain the union transcript and compute coverage statistics
        current_gene["union_transcript"] = self.__get_union_transcript(current_gene)
        results["genes"].append(current_gene)
        results["statistics"] = coverage_stats.compute_panel_level_statistics(
            results["genes"]
        )
        return results

    def __is_exon_padding_enabled(self):
        """

        :return:
        """
        return self.is_exon_padding_enabled

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
