import logging

import gelcoverage.stats.coverage_stats as coverage_stats
from gelcoverage.tools.cellbase_helper import CellbaseHelper
from gelcoverage.tools.panelapp_helper import PanelappHelper
from gelcoverage.tools.bigwig_reader import BigWigReader, UncoveredIntervalException
from gelcoverage.tools.bed_reader import BedReader, BedInterval
import gelcoverage.constants as constants


class GelCoverageInputError(Exception):

    pass


class GelCoverageRunner:

    def __init__(self, config):
        self.config = config
        # Run sanity checks on the configuration
        self.__config_sanity_checks()
        # Configure logs
        logging.basicConfig(level=self.config["log_level"])

        # Helper flags for the configuration
        self.is_exon_padding = self.config["exon_padding"] > 0
        self.is_find_gaps_enabled = self.config["coverage_threshold"] > 0
        self.is_wg_stats_enabled = self.config["wg_stats_enabled"]
        self.is_coding_region_stats_enabled = self.config["coding_region_stats_enabled"]
        self.is_exon_stats_enabled = self.config["exon_stats_enabled"]
        self.is_panel_analysis = True if "panel" in self.config and self.config["panel"] and \
                                 "panel_version" in self.config and self.config["panel_version"] else False
        self.is_gene_list_analysis = True if "gene_list" in self.config and self.config["gene_list"] else False

        # Attach bigwig reader
        self.bigwig_reader = BigWigReader(self.config['bw'])

        # Initialize PanelApp helper
        if self.is_panel_analysis:
            self.panelapp_helper = PanelappHelper(
                host=self.config['panelapp_host'],
                retries=self.config['panelapp_retries']
            )

        self.coding_regions = self.config['coding_regions'] if 'coding_regions' in self.config else None
        self.requires_cellbase = self.is_panel_analysis or self.is_gene_list_analysis \
                                 or (self.is_coding_region_stats_enabled and not self.coding_regions)
        if self.requires_cellbase:
            # Initialize CellBase helper only if no coding regions have been provided
            self.cellbase_helper = CellbaseHelper(
                species=self.config['cellbase_species'],
                version=self.config['cellbase_version'],
                assembly=self.config['cellbase_assembly'],
                host=self.config['cellbase_host'],
                retries=self.config['cellbase_retries'],
                filter_flags=self.config['transcript_filtering_flags'].split(","),
                filter_biotypes=self.config['transcript_filtering_biotypes'].split(",")
            )

        # Gets the list of genes to analyse
        requires_gene_list = self.is_panel_analysis or \
                             self.is_gene_list_analysis or \
                             self.is_coding_region_stats_enabled
        self.gene_list = []
        if requires_gene_list:
            self.gene_list = self.get_gene_list()

        if self.is_panel_analysis or self.is_gene_list_analysis:
            if len(self.gene_list) < 100:
                # Only prints the gene list when under 100 genes, otherwise becomes useless
                logging.info("Gene list to analyse: %s" % ",".join(self.gene_list))
            else:
                logging.info("%s genes to analyse" % str(len(self.gene_list)))

        # Opens the bigwig reader
        if self.is_wg_stats_enabled:
            self.wg_regions = self.config["wg_regions"]
            self.bed_reader = BedReader(self.wg_regions)
        self.uncovered_genes = {}
        self.__sanity_checks()

    def __config_sanity_checks(self):
        """
        Checks the input configuration is not missing any value
        :return:
        """
        # add default values to missing parameters
        if "panelapp_retries" not in self.config:
            # setting default value, infinite retries
            self.config["panelapp_retries"] = -1
        if "exon_padding" not in self.config or type(self.config["exon_padding"]) != int or \
                        self.config["exon_padding"] <= 0:
            # default value exon padding disabled
            self.config["exon_padding"] = 0
        if "log_level" not in self.config:
            # default value log level to INFO
            self.config["log_level"] = 20
        # check required parameters
        errors = []
        if "bw" not in self.config:
            errors.append("'bw' field is mising")
            # We must discuss this, it is not clear whether we always want a cellbase_helper just in
            # case or we want to get rid of it unless it is specified
        if "coverage_threshold" not in self.config:
            errors.append("'coverage_threshold' field is mising")
        if "coding_regions" not in self.config or not self.config['coding_regions']:
            if "cellbase_retries" not in self.config:
                # setting default value, infinite retries
                self.config["cellbase_retries"] = -1
            if "cellbase_species" not in self.config:
                errors.append("'cellbase_species' field is mising")
            if "cellbase_version" not in self.config:
                errors.append("'cellbase_version' field is mising")
            if "cellbase_assembly" not in self.config:
                errors.append("'cellbase_assembly' field is mising")
            if "cellbase_host" not in self.config:
                errors.append("'cellbase_host' field is mising")
        if "panelapp_host" not in self.config:
            errors.append("'panelapp_host' field is mising")
        if "panelapp_gene_confidence" not in self.config:
            errors.append("'panelapp_gene_confidence' field is mising")
        if "transcript_filtering_flags" not in self.config:
            errors.append("'transcript_filtering_flags' field is mising")
        if "transcript_filtering_biotypes" not in self.config:
            errors.append("'transcript_filtering_biotypes' field is mising")
        if "wg_stats_enabled" not in self.config:
            errors.append("'wg_stats_enabled' field is mising")
        if "exon_stats_enabled" not in self.config:
            errors.append("'exon_stats_enabled' field is mising")
        if "coding_region_stats_enabled" not in self.config:
            errors.append("'coding_region_stats_enabled' field is mising")
        # raise exceptino if necessary
        if len(errors) > 0:
            for error in errors:
                logging.error(error)
            raise GelCoverageInputError("Error in configuration data!")
        self.has_gene_list = "gene_list" in self.config and self.config["gene_list"] is not None
        self.has_coding_regions = "coding_regions" in self.config and self.config["coding_regions"] is not None
        # raise an error if both gene list and custom coding regions bed are provided
        if self.has_gene_list and self.has_coding_regions:
            raise GelCoverageInputError('Both modes gene list and custom coding region are not allowed, please'
                                        ' provide either "gene-list" or "coding-regions"')


    def __sanity_checks(self):
        """
        Checks on the configuration data to raise errors before starting long computations.
        :return:
        """
        if self.is_wg_stats_enabled:
            if not self.bed_reader.is_null_bed and self.bed_reader.has_chr_prefix != self.bigwig_reader.has_chr_prefix:
                raise GelCoverageInputError("Bed file defining whole genome analysis region and bigwig use different "
                                            "chromosome notations (i.e.: chr prefix).")
        logging.info("Sanity checks on the configuration OK!")

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
        if self.is_panel_analysis:
            # Get list of genes from PanelApp
            gene_list = self.panelapp_helper.get_gene_list(
                self.config["panel"],
                self.config["panel_version"],
                self.config["panelapp_gene_confidence"]
            )
        elif self.is_gene_list_analysis:
            # Get list of genes from parameter genes
            gene_list = self.config["gene_list"].split(",")
        elif self.has_coding_regions:
            gene_list = self.__get_gene_list_from_bed()
        else:
            # Warn the user as this will be time consuming
            logging.warning("You are about to run a whole exome coverage analysis!")
            # Retrieve the list of all genes
            gene_list = self.cellbase_helper.get_all_gene_names()
        return gene_list

    def __get_gene_list_from_bed(self):
        """
        Gets the gene list from a BED file
        :return:
        """
        gene_list = set()
        bedfile_handler = open(self.coding_regions, 'r')
        for line in bedfile_handler:
             gene_list.add(BedInterval(line).name.split('|')[0])
        bedfile_handler.close()
        return list(gene_list)

    def __get_parameters_output(self):
        """
        Returns a dictionary with the formatted configuration parameters used to run the analysis.
        :return: dict
        """
        parameters = {
            "gap_coverage_threshold": self.config["coverage_threshold"],
            "input_file": self.config["bw"],
            "transcript_filtering_flags": self.config['transcript_filtering_flags'],
            "transcript_filtering_biotypes": self.config['transcript_filtering_biotypes'],
            "exon_padding": self.config["exon_padding"],
            "panelapp_host": self.config["panelapp_host"],
            "panelapp_gene_confidence": self.config["panelapp_gene_confidence"],
            "wg_stats_enabled": self.config["wg_stats_enabled"],
            "wg_regions": self.config["wg_regions"],
            "exon_stats_enabled": self.config["exon_stats_enabled"],
            "coding_region_stats_enabled": self.config["coding_region_stats_enabled"],
        }
        if not self.has_coding_regions:
            parameters["cellbase_host"] = self.config["cellbase_host"],
            parameters["cellbase_version"] = self.config["cellbase_version"],
            parameters["species"] = self.config['cellbase_species']
            parameters["assembly"] = self.config['cellbase_assembly']
            # we don't want to add the gene list read from the bed as it might be too big
        if self.is_panel_analysis:
            parameters["panel"] = self.config['panel']
            parameters["panel_version"] = self.config['panel_version']
            parameters["panel_gene_confidence"] = self.config['panelapp_gene_confidence']
            parameters["gene_list"] = self.gene_list
        elif self.is_gene_list_analysis:
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
            constants.GENE_NAME: gene_name,
            constants.CHROMOSOME: chromosome,
            constants.TRANSCRIPTS: transcripts
        }

    @staticmethod
    def __initialize_transcript_dict(transcript_id, exons):
        """
        Returns the dictionary that will store transcript information
        :param transcript_id: the transcript id
        :return: the basic transcript data structure
        """
        return {
            constants.TRANSCRIPT_ID: transcript_id,
            constants.EXONS: exons
        }

    @staticmethod
    def __initialize_exon_dict(exon_number, start, end, padded_start, padded_end):
        """
        Returns the dictionary that will store exon information
        :param exon_number: the exon number
        :param start: the start position
        :param end: the end position
        :param padded_start: the padded start position
        :param padded_end: the padded end position
        :return: the basic exon data structure
        """
        exon = {
            constants.EXON: exon_number,
            constants.EXON_START: start,
            constants.EXON_END: end,
            constants.EXON_LENGTH: end - start + 1
        }
        if padded_start is not None and padded_start != start:
            exon[constants.EXON_PADDED_START] = padded_start
            exon[constants.EXON_PADDED_END] = padded_end
        return exon

    def __create_exon(self, chromosome, start, end, exon_idx, gc_content=None):
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
            start=start,
            end=end,
            padded_start=start - self.config["exon_padding"] if self.is_exon_padding else None,
            padded_end=end + self.config["exon_padding"] if self.is_exon_padding else None
        )
        # Update length
        exon[constants.EXON_LENGTH] = end - start + 1
        # Read from the bigwig file
        coverages = self.bigwig_reader.read_bigwig_coverages(
            chromosome,
            start,
            end)
        # Compute statistics at exon level (no GC content information)
        exon[constants.STATISTICS] = coverage_stats.compute_exon_level_statistics(coverages, gc_content)
        # Compute gaps
        if self.is_find_gaps_enabled:
            exon[constants.GAPS] = coverage_stats.find_gaps(
                coverages,
                start,
                self.config['coverage_threshold']
            )
        # logging.debug("Created exon %s" % exon_number)
        return exon

    def __create_transcript(self, _id, exons):
        """
        Creates a transcript data structure and computes statistics
        :param _id: the transcript id
        :param exons: the exons data structure
        :return: the transcript data structure
        """
        transcript = self.__initialize_transcript_dict(_id, exons)
        # Compute transcript level statistics by aggregating stats on every exon
        transcript[constants.STATISTICS] = coverage_stats.compute_transcript_level_statistics(
            transcript[constants.EXONS]
        )
        return transcript

    def __create_union_transcript(self, gene):
        """
        Creates union transcript and computes statistics for each of the exons and the union transcript
        :param gene: the gene data structure containing all exons
        :return: the data structure for the union transcript
        """
        logging.debug("Creating union transcript for gene %s" % gene[constants.GENE_NAME])
        all_exons = sum([transcript[constants.EXONS] for transcript in gene[constants.TRANSCRIPTS]], [])
        all_exons.sort(key=lambda x: x[constants.EXON_START])
        union_exons = []
        first_exon = all_exons[0]
        current_padded_end = first_exon[constants.EXON_PADDED_END] if self.is_exon_padding \
            else first_exon[constants.EXON_END]
        current_start = first_exon[constants.EXON_START]
        current_end = first_exon[constants.EXON_END]
        # TODO: consider strand when assigning exon indices
        exon_idx = 1
        for exon in all_exons[1:]:
            padded_start = exon[constants.EXON_PADDED_START] if self.is_exon_padding else exon[constants.EXON_START]
            padded_end = exon[constants.EXON_PADDED_END] if self.is_exon_padding else exon[constants.EXON_END]
            start = exon[constants.EXON_START]
            end = exon[constants.EXON_END]
            if padded_start <= current_padded_end:
                # Exons overlaps, we join them
                current_padded_end = max(current_padded_end, padded_end)
                current_end = max(current_end, end)
            else:
                # Exons do not overlap, they are different exons in the union transcript
                current_exon = self.__create_exon(
                    gene[constants.CHROMOSOME],
                    current_start,
                    current_end,
                    exon_idx
                )
                # Saves current exon and stores the next
                union_exons.append(current_exon)
                current_padded_end = padded_end
                current_start = start
                current_end = end
                exon_idx += 1
        # Stores the last exon
        last_exon = self.__create_exon(
            gene[constants.CHROMOSOME],
            current_start,
            current_end,
            exon_idx
        )
        # Saves current exon and stores the next
        union_exons.append(last_exon)
        # Create union transcript dict
        union_transcript = {
            constants.EXONS: union_exons,
            constants.STATISTICS: coverage_stats.compute_transcript_level_statistics(union_exons)
        }
        logging.debug("Built union transcript for gene %s" % gene[constants.GENE_NAME])
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
        gene[constants.UNION_TRANSCRIPT] = self.__create_union_transcript(gene)
        logging.info("Processed gene %s" % gene_name)
        return gene

    def __process_coding_region(self, bed):
        """
        Reads every interval defined in the BED file, gets the coverage information from the BigWig data, computes
        coverage statistics and formats the information in a JSON-friendly data structure.
        :param bed: the BED file handler
        :return: the data structure containing all statistics
        """
        # Checks that the BED file contains the expected information
        logging.info("Processing the information in the coverage bigwig on the bed intervals...")
        if bed is None:
            logging.error("Incorrect BED file!")
            raise RuntimeError("Incorrect BED file!")
        # Initialize results data structure
        genes = []
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
                    if len(current_exons) > 0:
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
                    if len(current_transcripts) > 0:
                        gene = self.__create_gene(current_gene_name, current_chromosome, current_transcripts)
                        # Save previous result
                        genes.append(gene)
                    # Sets data for next gene
                    current_transcripts = []
                    current_gene_name = gene_name
                    current_chromosome = chromosome
            # Store exon in data structure
            try:
                exon = self.__create_exon(
                    chromosome,
                    start,
                    end,
                    exon_number,
                    gc_content
                )
                current_exons.append(exon)
            except UncoveredIntervalException:
                self.uncovered_genes[gene_name] = chromosome
        # Adds last ocurrences
        if len(current_exons) > 0:
            transcript = self.__create_transcript(current_transcript_id, current_exons)
            current_transcripts.append(transcript)
        if len(current_transcripts) > 0:
            gene = self.__create_gene(gene_name, chromosome, current_transcripts)
            genes.append(gene)
        logging.info("Coverage bigwig coding region processed!")
        return genes

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
        results = {}
        bed = None
        if self.is_coding_region_stats_enabled:
            # Get genes annotations in BED format
            bed = self.__get_bed_for_exons()
            # Process the intervals for the coding region in the BED file
            results[constants.GENES] = self.__process_coding_region(bed)
            # Aggregate coding region statistics
            results["coding_region"] = coverage_stats.compute_coding_region_statistics(
                results[constants.GENES]
            )
            # Add uncovered genes
            results["uncovered_genes"] = [
                {constants.GENE_NAME: k, constants.CHROMOSOME: v} for k, v in self.uncovered_genes.iteritems()
            ]
            results = self.__delete_unnecessary_info_from_genes(results, constants)
        else:
            logging.info("Coding region analysis disabled")
        # Compute the whole genome statistics if enabled (this is time consuming)
        if self.is_wg_stats_enabled:
            results["whole_genome"] = coverage_stats.compute_whole_genome_statistics(
                self.bigwig_reader,
                self.bed_reader
            )
        else:
            logging.info("Whole genome analysis disabled")

        return self.__output(results), bed


    def __get_bed_for_exons(self):
        """
        Loads a bed file or creates it using cellbase connector
        :return: An iterator on a bedfile
        """
        if self.coding_regions:
            bed = self.load_bed_from_file()
        else:
            bed = self.cellbase_helper.make_exons_bed(self.gene_list,
                                                   has_chr_prefix=self.bigwig_reader.has_chr_prefix)
        return bed

    def load_bed_from_file(self):
        """
        Loads a bedfile
        :return: An iterator of BedInterval
        """
        bedfile_handler = open(self.coding_regions, 'r')
        for line in bedfile_handler:
            yield BedInterval(line)

    def __delete_unnecessary_info_from_genes(self, results, constants):
        """
        Removes unnecessary statistics of count bases at given coverage thresholds
        NOTE: the clean way would be not to store them and infer them dynamically from
        the percentage value...
        :param results:
        :param constants:
        :return:
        """
        for gene in results[constants.GENES]:
            del gene[constants.UNION_TRANSCRIPT][constants.STATISTICS][constants.BASES_LT15X]
            del gene[constants.UNION_TRANSCRIPT][constants.STATISTICS][constants.BASES_GTE15X]
            del gene[constants.UNION_TRANSCRIPT][constants.STATISTICS][constants.BASES_GTE30X]
            del gene[constants.UNION_TRANSCRIPT][constants.STATISTICS][constants.BASES_GTE50X]
            for exon in gene[constants.UNION_TRANSCRIPT][constants.EXONS]:
                del exon[constants.STATISTICS][constants.BASES_LT15X]
                del exon[constants.STATISTICS][constants.BASES_GTE15X]
                del exon[constants.STATISTICS][constants.BASES_GTE30X]
                del exon[constants.STATISTICS][constants.BASES_GTE50X]
            for transcript in gene[constants.TRANSCRIPTS]:
                del transcript[constants.STATISTICS][constants.BASES_LT15X]
                del transcript[constants.STATISTICS][constants.BASES_GTE15X]
                del transcript[constants.STATISTICS][constants.BASES_GTE30X]
                del transcript[constants.STATISTICS][constants.BASES_GTE50X]
                for exon in transcript[constants.EXONS]:
                    del exon[constants.STATISTICS][constants.BASES_LT15X]
                    del exon[constants.STATISTICS][constants.BASES_GTE15X]
                    del exon[constants.STATISTICS][constants.BASES_GTE30X]
                    del exon[constants.STATISTICS][constants.BASES_GTE50X]
        # Remove the exon statistics to save space if enabled
        if not self.is_exon_stats_enabled:
            for gene in results[constants.GENES]:
                for transcript in gene[constants.TRANSCRIPTS]:
                    del transcript[constants.EXONS]
                del gene[constants.UNION_TRANSCRIPT][constants.EXONS]
        return results
