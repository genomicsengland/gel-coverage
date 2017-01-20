import json
import unittest
import logging
import json
import gelcoverage.constants as constants


class OutputVerifier(unittest.TestCase):

    def verify_output(self, json, expected_gene_list=None):
        logging.info("Verifying JSON output...")
        self.expected_gene_list = expected_gene_list
        self.assertEqual(type(json), dict)
        # Verify that content in parameters is correct
        self.__verify_dict_field(json, "parameters", dict)
        self.__verify_parameters(json["parameters"])
        # Verify that coverage results are correct
        self.__verify_dict_field(json, "results", dict)
        if json["parameters"]["wg_stats_enabled"]:
            self.__verify_dict_field(json["results"], "whole_genome", dict)
            self.__verify_dict_field(json["results"]["whole_genome"], constants.STATISTICS, dict)
            self.__verify_dict_field(json["results"]["whole_genome"], constants.CHROMOSOMES, list)
            self.__verify_wg_stats(json["results"]["whole_genome"][constants.STATISTICS])
            for chr_stats in json["results"]["whole_genome"][constants.CHROMOSOMES]:
                self.__verify_wg_stats(chr_stats)
        if json["parameters"]["coding_region_stats_enabled"]:
            self.__verify_dict_field(json["results"], "coding_region", dict)
            self.__verify_dict_field(json["results"]["coding_region"], constants.STATISTICS, dict)
            self.__verify_dict_field(json["results"]["coding_region"], constants.CHROMOSOMES, list)
            self.__verify_panel_stats(json["results"]["coding_region"][constants.STATISTICS])
            for chr_stats in json["results"]["coding_region"][constants.CHROMOSOMES]:
                self.__verify_panel_stats(chr_stats)
            self.__verify_uncovered_genes(json["results"])
            self.__verify_genes(json["results"], json["parameters"]["exon_stats_enabled"])
        logging.info("JSON verified!")

    def __verify_dict_field(self, _dict, name, types):
        """
        Generic verification  of a field in a dictionary
        :param _dict: the dictionary
        :param name: the name of the field
        :param _type: the expected type
        :return:
        """
        if type(types) != list:
            types = [types]
        if str in types and unicode not in types:
            types.append(unicode)
        if unicode in types and str not in types:
            types.append(str)
        self.assertTrue(name in _dict, msg="Missing field '%s'" % name)
        self.assertTrue(type(_dict[name]) in types,
                        msg="Erroneous type of the field '%s': "
                            "found %s, expected any of %s" % (
                            name, str(type(_dict[name])), ",".join([str(x) for x in types])))

    def __verify_parameters(self, parameters):
        try:
            self.__verify_dict_field(parameters, "gap_coverage_threshold", int)
            self.__verify_dict_field(parameters, "input_file", str)
            self.__verify_dict_field(parameters, "configuration_file", str)
            self.__verify_dict_field(parameters, "species", str)
            self.__verify_dict_field(parameters, "assembly", str)
            self.__verify_dict_field(parameters, "transcript_filtering_flags", str)
            self.__verify_dict_field(parameters, "transcript_filtering_biotypes", str)
            self.__verify_dict_field(parameters, "cellbase_host", str)
            self.__verify_dict_field(parameters, "cellbase_version", str)
            self.__verify_dict_field(parameters, "panelapp_host", str)
            self.__verify_dict_field(parameters, "panelapp_gene_confidence", str)
            self.__verify_dict_field(parameters, "wg_stats_enabled", bool)
            self.__verify_dict_field(parameters, "wg_regions", [str, type(None)])
            self.__verify_dict_field(parameters, "exon_stats_enabled", bool)
            if self.expected_gene_list is not None:
                self.__verify_dict_field(parameters, "gene_list", list)
                self.assertEqual(parameters["gene_list"],
                                 self.expected_gene_list,
                                 msg="Gene list not matching the expected list: "
                                     "found '%s', expected '%s'" % (
                                     parameters["gene_list"],
                                     self.expected_gene_list
                                 ))
            if "panel" in parameters:
                self.__verify_dict_field(parameters, "panel", str)
                self.__verify_dict_field(parameters, "panel_version", str)
        except AssertionError, e:
            logging.error("Error verifying configuration parameters")
            logging.error(json.dumps(parameters, indent=4))
            raise e

    def __verify_uncovered_genes(self, results):

        self.__verify_dict_field(results, "uncovered_genes", list)
        observed_genes = []
        for uncovered_gene in results["uncovered_genes"]:
            self.__verify_dict_field(uncovered_gene, constants.CHROMOSOME, str)
            self.__verify_dict_field(uncovered_gene, constants.GENE_NAME, str)
            self.assertTrue(uncovered_gene[constants.GENE_NAME] not in observed_genes,
                                msg="Duplicated gene '%s'" % uncovered_gene[constants.GENE_NAME])
            observed_genes.append(uncovered_gene[constants.GENE_NAME])

    def __verify_genes(self, results, exon_stats_enabled):
        # Verify every gene
        self.__verify_dict_field(results, constants.GENES, list)
        observed_genes = []
        for gene in results[constants.GENES]:
            if self.expected_gene_list is not None:
                self.assertTrue(gene[constants.GENE_NAME] in self.expected_gene_list,
                                msg="Unexpected gene found in results '%s'" % gene[constants.GENE_NAME])
            self.__verify_dict_field(gene, constants.CHROMOSOME, str)
            # Checks that genes are not repeated
            self.assertTrue(gene[constants.GENE_NAME] not in observed_genes,
                            msg="Duplicated gene '%s'" % gene[constants.GENE_NAME])
            observed_genes.append(gene[constants.GENE_NAME])
            # Verify every transcript
            observed_transcripts = []
            for transcript in gene[constants.TRANSCRIPTS]:
                self.__verify_transcript(transcript)
                self.assertTrue(transcript[constants.TRANSCRIPT_ID] not in observed_transcripts,
                                msg="Duplicated transcript '%s'" % transcript[constants.TRANSCRIPT_ID])
                observed_transcripts.append(transcript[constants.TRANSCRIPT_ID])
                # Verify every exon
                if exon_stats_enabled:
                    observed_exons = []
                    for exon in transcript[constants.EXONS]:
                        self.__verify_exon(exon, gene[constants.GENE_NAME], transcript[constants.TRANSCRIPT_ID])
                        self.assertTrue(exon[constants.EXON] not in observed_exons,
                                        msg="Duplicated exon '%s'" % exon[constants.EXON])
                        observed_exons.append(exon[constants.EXON])
                        # Verify gaps
                        observed_gaps = []
                        for gap in exon[constants.GAPS]:
                            self.__verify_gap(gap, exon)
                            self.assertTrue((gap[constants.GAP_START], gap[constants.GAP_END]) not in observed_gaps,
                                            msg="Repeated gap '%s-%s'" % (
                                                gap[constants.GAP_START], gap[constants.GAP_END]
                                            ))
                            observed_gaps.append((gap[constants.GAP_START], gap[constants.GAP_END]))
                else:
                    self.assertTrue(constants.EXONS not in transcript)
            union_transcript = gene[constants.UNION_TRANSCRIPT]
            self.__verify_transcript_stats(union_transcript[constants.STATISTICS], has_gc=False)
            self.verify_union_transcript(gene, exon_stats_enabled)

    def __verify_transcript(self, transcript):
        self.__verify_dict_field(transcript, constants.TRANSCRIPT_ID, str)
        self.assertTrue(str(transcript[constants.TRANSCRIPT_ID]).startswith("ENS"),
                        msg="Wrong transcript id '%s'" % transcript[constants.TRANSCRIPT_ID])
        self.__verify_transcript_stats(transcript[constants.STATISTICS])

    def verify_union_transcript(self, gene, exon_stats_enabled):
        # Basic checks
        self.__verify_dict_field(gene, constants.UNION_TRANSCRIPT, dict)
        union_transcript = gene[constants.UNION_TRANSCRIPT]
        self.__verify_transcript_stats(gene[constants.UNION_TRANSCRIPT][constants.STATISTICS], has_gc=False)
        # Verifies exons
        if exon_stats_enabled:
            for ut_exon in union_transcript[constants.EXONS]:
                self.__verify_exon(ut_exon, gene[constants.GENE_NAME], constants.UNION_TRANSCRIPT, has_gc=False)
                # Verify gaps
                for gap in ut_exon[constants.GAPS]:
                    self.__verify_gap(gap, ut_exon)
            # Verifies union transcript build up
            all_exons = sum([transcript[constants.EXONS] for transcript in gene[constants.TRANSCRIPTS]], [])
            ut_exon_coordinates = [self.__get_exon_coordinates(x) for x in all_exons]
            ut_start = min([start for (start, _) in ut_exon_coordinates])
            ut_end = max([end for (_, end) in ut_exon_coordinates])
            #ut_positions = [True if (
            #    self.__is_position_overlapped(x, ut_exon)
            #    for ut_exon in union_transcript["exons"]
            #) else False for x in range(ut_start, ut_end)]
            ut_positions = []
            for x in range(ut_start, ut_end):
                found = False
                for ut_exon in union_transcript[constants.EXONS]:
                    if self.__is_position_overlapped(x, ut_exon):
                        found = True
                        break
                ut_positions.append(found)
            # Checks if exons in the union transcript are coherent with other exons

            for idx, ut_position in enumerate(ut_positions):
                real_position = ut_start + idx
                found = False
                for exon in all_exons:
                    start, end = self.__get_exon_coordinates(exon)
                    if real_position >= start and real_position <= end:
                        found = True
                        break
                if ut_position:
                    # Verify position included in the union transcript
                    self.assertTrue(found,
                                    msg="Position '%s' belonging to union transcript "
                                    "is not supported by any exon at gene %s" %
                                    (real_position, gene[constants.GENE_NAME])
                                    )
                else:
                    # Verify position not included in the union transcript
                    self.assertFalse(found,
                                     "Position '%s' not belonging to union transcript "
                                     "is present in at least one exon at gene %s" %
                                     (real_position, gene[constants.GENE_NAME])
                                     )
        else:
            self.assertTrue(constants.EXONS not in union_transcript)

    def __is_position_overlapped(self, position, exon):
        """
        Checks if the position overlaps the exon
        :param position:
        :param exon:
        :return:
        """
        start, end = self.__get_exon_coordinates(exon)
        return position >= start and position <= end

    def __get_exon_coordinates(self, exon):
        """
        Get coordinates dealing with exon padding parameter
        :param exon:
        :return: start and end positions
        """
        start = None
        end = None
        if self.__is_padding_enabled():
            start = exon[constants.EXON_PADDED_START]
            end = exon[constants.EXON_PADDED_END]
        else:
            start = exon[constants.EXON_START]
            end = exon[constants.EXON_END]
        return (start, end)

    def __is_padding_enabled(self):
        return self.config["exon_padding"] > 0

    def __verify_transcript_stats(self, stats, has_gc = True):
        try:
            self.assertEqual(type(stats), dict)
            if has_gc:
                self.__verify_dict_field(stats, constants.GC_CONTENT, float)
                self.assertTrue(stats[constants.GC_CONTENT] >= 0 and
                                stats[constants.GC_CONTENT] <= 1)
            self.__verify_dict_field(stats, constants.AVERAGE, float)
            self.assertTrue(stats[constants.AVERAGE] >= 0)
            self.__verify_dict_field(stats, constants.GTE15X, float)
            self.assertTrue(stats[constants.GTE15X] >= 0 and
                            stats[constants.GTE15X] <= 1)
            self.__verify_dict_field(stats, constants.GTE30X, float)
            self.assertTrue(stats[constants.GTE30X] >= 0 and
                            stats[constants.GTE30X] <= 1)
            self.__verify_dict_field(stats, constants.GTE50X, float)
            self.assertTrue(stats[constants.GTE50X] >= 0 and
                            stats[constants.GTE50X] <= 1)
            self.__verify_dict_field(stats, constants.LT15X, float)
            self.assertTrue(stats[constants.LT15X] >= 0 and
                            stats[constants.LT15X] <= 1)
            self.__verify_dict_field(stats, constants.BASES, [int, long])
            self.assertTrue(stats[constants.BASES] >= 0)
            self.__verify_dict_field(stats, constants.MEDIAN, float)
            self.assertTrue(stats[constants.MEDIAN] >= 0)
            self.__verify_dict_field(stats, constants.PERCENTILE75, float)
            self.assertTrue(stats[constants.PERCENTILE75] >= 0)
            self.__verify_dict_field(stats, constants.PERCENTILE25, float)
            self.assertTrue(stats[constants.PERCENTILE25] >= 0)
        except AssertionError, e:
            logging.error("Error verifying transcript statistics")
            logging.error(json.dumps(stats, indent=4))
            raise e

    def __verify_exon(self, exon, gene_name, transcript_id, has_gc = True):
        try:
            self.assertEqual(type(exon), dict)
            self.__verify_dict_field(exon, constants.EXON_START, int)
            self.assertTrue(exon[constants.EXON_START] >= 0)
            self.__verify_dict_field(exon, constants.EXON_END, int)
            self.assertTrue(exon[constants.EXON_END] >= 0)
            self.assertTrue(exon[constants.EXON_END] >= exon[constants.EXON_START], msg="End < start")
            if self.__is_padding_enabled():
                self.__verify_dict_field(exon, constants.EXON_PADDED_START, int)
                self.assertTrue(exon[constants.EXON_PADDED_START] >= 0)
                self.__verify_dict_field(exon, constants.EXON_PADDED_END, int)
                self.assertTrue(exon[constants.EXON_PADDED_END] >= 0)
                self.assertTrue(exon[constants.EXON_PADDED_END] > exon[constants.EXON_PADDED_START],
                                msg="Padded end <= padded start")
                self.assertTrue(exon[constants.EXON_START] - exon[constants.EXON_PADDED_START] == self.config["exon_padding"],
                                msg="Incorrect start coordinate padding")
                self.assertTrue(exon[constants.EXON_PADDED_END] - exon[constants.EXON_END] == self.config["exon_padding"],
                                msg="Incorrect end coordinate padding")
            self.__verify_dict_field(exon, constants.EXON, str)
            self.assertTrue(str(exon[constants.EXON]).startswith(constants.EXON),
                            msg="Exon number is not well formed")
            self.__verify_exon_statistics(exon, has_gc)
        except AssertionError, e:
            logging.error("Error verifying exon at %s:%s" % (gene_name, transcript_id))
            logging.error(json.dumps(exon, indent=4))
            raise e

    def __verify_exon_statistics(self, exon, has_gc=True):
        self.__verify_dict_field(exon, constants.STATISTICS, dict)
        statistics = exon[constants.STATISTICS]
        try:
            if has_gc:
                self.__verify_dict_field(statistics, constants.GC_CONTENT, float)
                self.assertTrue(statistics[constants.GC_CONTENT] >= 0 and
                                statistics[constants.GC_CONTENT] <= 1)
            self.__verify_dict_field(statistics, constants.AVERAGE, float)
            self.assertTrue(statistics[constants.AVERAGE] >= 0)
            self.__verify_dict_field(statistics, constants.GTE15X, float)
            self.assertTrue(statistics[constants.GTE15X] >= 0 and
                            statistics[constants.GTE15X] <= 1)
            self.assertTrue(statistics[constants.GTE30X] >= 0 and
                            statistics[constants.GTE30X] <= 1)
            self.assertTrue(statistics[constants.GTE50X] >= 0 and
                            statistics[constants.GTE50X] <= 1)
            self.assertTrue(statistics[constants.LT15X] >= 0 and
                            statistics[constants.LT15X] <= 1)
            self.__verify_dict_field(statistics, constants.BASES, [int, long])
            self.assertTrue(statistics[constants.BASES] >= 0)
            self.__verify_dict_field(statistics, constants.MEDIAN, float)
            self.assertTrue(statistics[constants.MEDIAN] >= 0)
            self.__verify_dict_field(statistics, constants.PERCENTILE75, float)
            self.assertTrue(statistics[constants.PERCENTILE75] >= 0)
            self.__verify_dict_field(statistics, constants.PERCENTILE25, float)
            self.assertTrue(statistics[constants.PERCENTILE25] >= 0)
        except AssertionError, e:
            logging.error("Error verifying exon statistics")
            logging.error(json.dumps(statistics, indent=4))
            raise e

    def __verify_gap(self, gap, exon):
        try:
            self.assertEqual(type(gap), dict)
            self.__verify_dict_field(gap, constants.GAP_START, int)
            (start, end) = self.__get_exon_coordinates(exon)
            self.assertTrue(gap[constants.GAP_START] >= start and gap[constants.GAP_START] <= end)
            self.__verify_dict_field(gap, constants.GAP_END, int)
            self.assertTrue(gap[constants.GAP_END] >= start and gap[constants.GAP_START] <= end and
                            gap[constants.GAP_END] >= gap[constants.GAP_START])
            self.__verify_dict_field(gap, constants.GAP_LENGTH, int)
            self.assertTrue(gap[constants.GAP_LENGTH] >= 1 and gap[constants.GAP_LENGTH] <=
                            gap[constants.GAP_END] - gap[constants.GAP_START] + 1)
        except AssertionError, e:
            logging.error("Error verifying gap")
            logging.error(json.dumps(gap, indent=4))
            raise e

    def __verify_panel_stats(self, statistics):
        try:
            self.__verify_dict_field(statistics, constants.AVERAGE, float)
            self.assertTrue(statistics[constants.AVERAGE] >= 0)
            self.__verify_dict_field(statistics, constants.GTE15X, float)
            self.__verify_dict_field(statistics, constants.GTE30X, float)
            self.__verify_dict_field(statistics, constants.GTE50X, float)
            self.__verify_dict_field(statistics, constants.LT15X, float)
            self.assertTrue(statistics[constants.GTE15X] >= 0 and
                            statistics[constants.GTE15X] <= 1)
            self.assertTrue(statistics[constants.GTE30X] >= 0 and
                            statistics[constants.GTE30X] <= 1)
            self.assertTrue(statistics[constants.GTE50X] >= 0 and
                            statistics[constants.GTE50X] <= 1)
            self.assertTrue(statistics[constants.LT15X] >= 0 and
                            statistics[constants.LT15X] <= 1)
            self.__verify_dict_field(statistics, constants.BASES, [int, long])
            self.assertTrue(statistics[constants.BASES] >= 0)
            self.__verify_dict_field(statistics, constants.MEDIAN, float)
            self.assertTrue(statistics[constants.MEDIAN] >= 0)
            self.__verify_dict_field(statistics, constants.PERCENTILE75, float)
            self.assertTrue(statistics[constants.PERCENTILE75] >= 0)
            self.__verify_dict_field(statistics, constants.PERCENTILE25, float)
            self.assertTrue(statistics[constants.PERCENTILE25] >= 0)
        except AssertionError, e:
            logging.error("Error panel statistics")
            logging.error(json.dumps(statistics, indent=4))
            raise e

    def __verify_wg_stats(self, statistics):
        try:
            self.__verify_dict_field(statistics, constants.AVERAGE, float)
            self.assertTrue(statistics[constants.AVERAGE] >= 0)
            self.__verify_dict_field(statistics, constants.GTE15X, float)
            self.__verify_dict_field(statistics, constants.GTE30X, float)
            self.__verify_dict_field(statistics, constants.GTE50X, float)
            self.__verify_dict_field(statistics, constants.LT15X, float)
            self.assertTrue(statistics[constants.GTE15X] >= 0 and
                            statistics[constants.GTE15X] <= 1)
            self.assertTrue(statistics[constants.GTE30X] >= 0 and
                            statistics[constants.GTE30X] <= 1)
            self.assertTrue(statistics[constants.GTE50X] >= 0 and
                            statistics[constants.GTE50X] <= 1)
            self.assertTrue(statistics[constants.LT15X] >= 0 and
                            statistics[constants.LT15X] <= 1)
            self.__verify_dict_field(statistics, constants.BASES, [int, long])
            self.assertTrue(statistics[constants.BASES] >= 0)
            self.__verify_dict_field(statistics, constants.MEDIAN, float)
            self.assertTrue(statistics[constants.MEDIAN] >= 0)
            self.__verify_dict_field(statistics, constants.PERCENTILE75, float)
            self.assertTrue(statistics[constants.PERCENTILE75] >= 0)
            self.__verify_dict_field(statistics, constants.PERCENTILE25, float)
            self.assertTrue(statistics[constants.PERCENTILE25] >= 0)
            self.__verify_dict_field(statistics, constants.RMSD, float)
            self.assertTrue(statistics[constants.RMSD] >= 0)
        except AssertionError, e:
            logging.error("Error panel statistics")
            logging.error(json.dumps(statistics, indent=4))
            raise e