import logging

import gelcoverage.stats.coverage_stats as coverage_stats
from gelcoverage.tools.cellbase_helper import CellbaseHelper
from gelcoverage.tools.panelapp_helper import PanelappHelper
from gelcoverage.tools.bigwig_reader import BigWigReader


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
        # Opens the bigwig reader
        self.bigwig_reader = BigWigReader(self.config['bw'])
        # Flag indicating if exon padding is enabled
        if "exon_padding" not in self.config or type(self.config["exon_padding"]) != int or \
            self.config["exon_padding"] <= 0:
            self.config["exon_padding"] = 0
        self.is_exon_padding = self.config["exon_padding"] > 0
        self.is_find_gaps_enabled = self.config["coverage_threshold"] > 0

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
            "configuration_file": self.config["configuration_file"],
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
    def __parse_bed_interval(interval):
        """
        Extracts information from a BED interval.
        :param interval: the BED interval
        :return: the extracted information in separate variables
        """
        chromosome = interval.chrom
        start = int(interval.start)
        end = int(interval.end)
        gene_name, transcript_id, exon_number = interval.name.split("|")
        exon_number = str(exon_number)
        strand = interval.strand
        gc_content = float(interval.score)
        return chromosome, start, end, gene_name, transcript_id, exon_number, strand, gc_content

    @staticmethod
    def __initialize_gene_dict(gene_name, chromosome, transcripts):
        """
        Returns the dictionary that will store gene information
        :param gene_name: the gene name
        :param chromosome: the chromosome
        :return: the basic gene data structure
        """
        return {
            "name": gene_name,
            "chromosome": chromosome,
            "transcripts": transcripts
        }

    @staticmethod
    def __initialize_transcript_dict(transcript_id, exons):
        """
        Returns the dictionary that will store transcript information
        :param transcript_id: the transcript id
        :return: the basic transcript data structure
        """
        return {
            "id": transcript_id,
            "exons": exons
        }

    @staticmethod
    def __initialize_exon_dict(exon_number, start, end, padded_start, padded_end):
        """
        Returns the dictionary that will store exon information
        :param exon_number: the exon number
        :param start: the start position
        :param end: the end position
        :return: the basic exon data structure
        """
        exon = {
                "exon_number": exon_number,
                "start": start,
                "end": end,
                "length": end - start + 1
            }
        if padded_start is not None and padded_start != start:
            exon["padded_start"] = padded_start
            exon["padded_end"] = padded_end
        return exon

    def __create_exon(self, chromosome, start, end, exon_idx, gc_content = None):
        """
        Creates an exon data structure, computes the coverage statistics and find gaps
        :param chromosome: the chromosome
        :param start: the padded start position
        :param end: the padded end position
        :param exon_idx: the exon index (two possible formats: int or str like exonN)
        :param gc_content: the GC content for the exon
        :return: the exon data structure
        """
        exon_number = "exon%s" % exon_idx if type(exon_idx) == int else exon_idx
        exon = GelCoverageRunner.__initialize_exon_dict(
            exon_number,
            start + self.config["exon_padding"] if self.is_exon_padding else start,
            end - self.config["exon_padding"] if self.is_exon_padding else end,
            start if self.is_exon_padding else None,
            end if self.is_exon_padding else None,
        )
        # Update length
        exon["length"] = end - start + 1
        # Read from the bigwig file
        coverages = self.bigwig_reader.read_bigwig_coverages(
            chromosome,
            start,
            end)
        # Compute statistics at exon level (no GC content information)
        exon["statistics"] = coverage_stats.compute_exon_level_statistics(coverages, gc_content)
        # Compute gaps
        if self.is_find_gaps_enabled:
            exon["gaps"] = coverage_stats.find_gaps(
                coverages,
                start,
                self.config['coverage_threshold']
            )
        logging.debug("Created exon %s" % exon_number)
        return exon

    def __create_transcript(self, id, exons):
        """
        Creates a transcript data structure and computes statistics
        :param id: the transcript id
        :param exons: the exons data structure
        :return: the transcript data structure
        """
        transcript = self.__initialize_transcript_dict(id, exons)
        # Compute transcript level statistics by aggregating stats on every exon
        transcript["statistics"] = coverage_stats.compute_transcript_level_statistics(
            transcript["exons"]
        )
        logging.debug("Created transcript %s" % id)
        return transcript

    def __create_union_transcript(self, gene):
        """
        Creates union transcript and computes statistics for each of the exons and the union transcript
        :param gene: the gene data structure containing all exons
        :return: the data structure for the union transcript
        """
        logging.debug("Creating union transcript for gene %s" % gene["name"])
        all_exons = sum([transcript["exons"] for transcript in gene["transcripts"]], [])
        all_exons.sort(key = lambda x: x["start"])
        union_exons = []
        first_exon = all_exons[0]
        current_start = first_exon["padded_start"] if self.is_exon_padding else first_exon["start"]
        current_end = first_exon["padded_end"] if self.is_exon_padding else first_exon["end"]
        exon_idx = 1
        for exon in all_exons[1:]:
            start = exon["padded_start"] if self.is_exon_padding else exon["start"]
            end = exon["padded_end"] if self.is_exon_padding else exon["end"]
            if start <= current_end:
                # Exons overlaps, we join them
                current_end = max(current_end, end)
            else:
                # Exons do not overlap, they are different exons in the union transcript
                current_exon = self.__create_exon(
                    gene["chromosome"],
                    current_start,
                    current_end,
                    exon_idx
                )
                # Saves current exon and stores the next
                union_exons.append(current_exon)
                exon_idx += 1
        # Stores the last exon
        last_exon = self.__create_exon(
            gene["chromosome"],
            start,
            end,
            exon_idx
        )
        # Saves current exon and stores the next
        union_exons.append(last_exon)
        # Create union transcript dict
        union_transcript = {
            "exons" : union_exons,
            "statistics": coverage_stats.compute_transcript_level_statistics(union_exons)
        }
        logging.debug("Created union transcript for gene %s" % gene["name"])
        return union_transcript

    def __create_gene(self, gene_name, chromosome, transcripts):
        """
        Creates the gene data structure and creates the union transcript
        :param gene_name: the gene name
        :param chromosome: the chromosome
        :param transcripts: the data structure for all transcripts belonging to this gene
        :return: the gene data structure
        """
        gene = self.__initialize_gene_dict(gene_name, chromosome, transcripts)
        # Calculate union transcript and compute stats
        gene["union_transcript"] = self.__create_union_transcript(gene)
        logging.debug("Created gene %s" % gene_name)
        return gene

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

        current_exons = []
        current_transcripts = []
        current_transcript_id = None
        current_gene_name = None
        current_chromosome = None
        # Process every interval in the BED file
        for interval in bed:
            # Reads data from BED entry
            chromosome, start, end, gene_name, transcript_id, \
                exon_number, strand, gc_content = GelCoverageRunner.__parse_bed_interval(interval)

            # Store transcript in data structure. Needs to run before gene
            if current_transcript_id != transcript_id:
                if current_transcript_id is None:
                    # Initialize for the first transcript
                    current_transcript_id = transcript_id
                else:
                    transcript = self.__create_transcript(current_transcript_id, current_exons)
                    current_transcripts.append(transcript)
                    logging.info("Processed transcript %s of gene %s" % (current_transcript_id, current_gene_name))
                    # Create a new data structure for new gene
                    current_transcript_id = transcript_id
                    current_exons = []
            # Store gene in data structure
            if current_gene_name != gene_name:
                if current_gene_name is None:
                    # Initialize for the first transcript
                    current_gene_name = gene_name
                    current_chromosome = chromosome
                else:
                    gene = self.__create_gene(current_gene_name, current_chromosome, current_transcripts)
                    # Save previous result
                    results["genes"].append(gene)
                    logging.info("Processed all transcripts for gene %s" % current_gene_name)
                    # Sets data for next gene
                    current_transcripts = []
                    current_gene_name = gene_name
                    current_chromosome = chromosome

            # Store exon in data structure
            exon = self.__create_exon(
                chromosome,
                start,
                end,
                exon_number,
                gc_content
            )
            current_exons.append(exon)
        # Adds last ocurrences
        transcript = self.__create_transcript(current_transcript_id, current_exons)
        current_transcripts.append(transcript)
        gene = self.__create_gene(gene_name, chromosome, current_transcripts)
        results["genes"].append(gene)
        results["statistics"] = coverage_stats.compute_panel_level_statistics(
            results["genes"]
        )
        results["whole_genome_statistics"] = coverage_stats.compute_whole_genome_statistics(
            self.bigwig_reader
        )
        return results

    def __output(self, results):
        """
        Builds the output data structure
        :param results: the coverage analysis data structure
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
        :return: the output data structure
        """
        logging.info("Starting coverage analysis")
        # Get genes annotations in BED format
        bed = self.cellbase_helper.make_exons_bed(self.gene_list)
        # Process the intervals in the BED file
        results = self.__process_bed_file(bed)
        return self.__output(results)
