import json
import unittest
import logging
import json


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
        self.__verify_dict_field(json["results"], "statistics", dict)
        self.__verify_panel_stats(json["results"]["statistics"])
        self.__verify_uncovered_genes(json["results"])
        self.__verify_genes(json["results"])
        logging.info("JSON verified!")

    def __verify_dict_field(self, _dict, name, _type):
        """
        Generic verification  of a field in a dictionary
        :param _dict: the dictionary
        :param name: the name of the field
        :param _type: the expected type
        :return:
        """
        types = [_type]
        if _type == str:
            types.append(unicode)
        if _type == unicode:
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
            self.__verify_dict_field(uncovered_gene, "chromosome", str)
            self.__verify_dict_field(uncovered_gene, "gene_name", str)
            self.assertTrue(uncovered_gene["gene_name"] not in observed_genes,
                                msg="Duplicated gene '%s'" % uncovered_gene["gene_name"])
            observed_genes.append(uncovered_gene["gene_name"])

    def __verify_genes(self, results):
        # Verify every gene
        self.__verify_dict_field(results, "genes", list)
        observed_genes = []
        for gene in results["genes"]:
            if self.expected_gene_list is not None:
                self.assertTrue(gene["name"] in self.expected_gene_list,
                                msg="Unexpected gene found in results '%s'" % gene["name"])
            self.__verify_dict_field(gene, "chromosome", str)
            # Checks that genes are not repeated
            self.assertTrue(gene["name"] not in observed_genes,
                            msg="Duplicated gene '%s'" % gene['name'])
            observed_genes.append(gene["name"])
            # Verify every transcript
            observed_transcripts = []
            for transcript in gene["transcripts"]:
                self.__verify_transcript(transcript)
                self.assertTrue(transcript["id"] not in observed_transcripts,
                                msg="Duplicated transcript '%s'" % transcript['id'])
                observed_transcripts.append(transcript["id"])
                # Verify every exon
                observed_exons = []
                for exon in transcript["exons"]:
                    self.__verify_exon(exon, gene["name"], transcript["id"])
                    self.assertTrue(exon["exon_number"] not in observed_exons,
                                    msg="Duplicated exon '%s'" % exon['exon_number'])
                    observed_exons.append(exon["exon_number"])
                    # Verify gaps
                    observed_gaps = []
                    for gap in exon["gaps"]:
                        self.__verify_gap(gap, exon)
                        self.assertTrue((gap["start"], gap["end"]) not in observed_gaps,
                                        msg="Repeated gap '%s-%s'" % (
                                            gap["start"], gap["end"]
                                        ))
                        observed_gaps.append((gap["start"], gap["end"]))
            union_transcript = gene["union_transcript"]
            self.__verify_transcript_stats(union_transcript["statistics"], has_gc=False)
            self.verify_union_transcript(gene)

    def __verify_transcript(self, transcript):
        self.__verify_dict_field(transcript, "id", str)
        self.assertTrue(str(transcript["id"]).startswith("ENS"),
                        msg="Wrong transcript id '%s'" % transcript["id"])
        self.__verify_transcript_stats(transcript["statistics"])

    def verify_union_transcript(self, gene):
        # Basic checks
        self.__verify_dict_field(gene, "union_transcript", dict)
        union_transcript = gene["union_transcript"]
        # Verifies exons
        for ut_exon in union_transcript["exons"]:
            self.__verify_exon(ut_exon, gene["name"], "union_transcript", has_gc=False)
            # Verify gaps
            for gap in ut_exon["gaps"]:
                self.__verify_gap(gap, ut_exon)
        # Verifies union transcript build up
        all_exons = sum([transcript["exons"] for transcript in gene["transcripts"]], [])
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
            for ut_exon in union_transcript["exons"]:
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
                                (real_position, gene["name"])
                                )
            else:
                # Verify position not included in the union transcript
                self.assertFalse(found,
                                 "Position '%s' not belonging to union transcript "
                                 "is present in at least one exon at gene %s" %
                                 (real_position, gene["name"])
                                 )

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
            start = exon["padded_start"]
            end = exon["padded_end"]
        else:
            start = exon["start"]
            end = exon["end"]
        return (start, end)

    def __is_padding_enabled(self):
        return self.config["exon_padding"] > 0

    def __verify_transcript_stats(self, stats, has_gc = True):
        try:
            self.assertEqual(type(stats), dict)
            self.__verify_dict_field(stats, "bases_gte_15x", int)
            self.assertTrue(stats["bases_gte_15x"] >= 0)
            self.__verify_dict_field(stats, "bases_gte_30x", int)
            self.assertTrue(stats["bases_gte_30x"] >= 0)
            self.__verify_dict_field(stats, "bases_gte_50x", int)
            self.assertTrue(stats["bases_gte_50x"] >= 0)
            self.__verify_dict_field(stats, "bases_lt_15x", int)
            self.assertTrue(stats["bases_lt_15x"] >= 0)
            self.__verify_dict_field(stats, "bases_lt_3x", int)
            self.assertTrue(stats["bases_lt_3x"] >= 0)
            if has_gc:
                self.__verify_dict_field(stats, "gc_content", float)
                self.assertTrue(stats["gc_content"] >= 0 and
                                stats["gc_content"] <= 1)
            self.__verify_dict_field(stats, "mean", float)
            self.assertTrue(stats["mean"] >= 0)
            self.__verify_dict_field(stats, "percent_gte_15x", float)
            self.assertTrue(stats["percent_gte_15x"] >= 0 and
                            stats["percent_gte_15x"] <= 1)
            self.__verify_dict_field(stats, "percent_gte_30x", float)
            self.assertTrue(stats["percent_gte_30x"] >= 0 and
                            stats["percent_gte_30x"] <= 1)
            self.__verify_dict_field(stats, "percent_gte_50x", float)
            self.assertTrue(stats["percent_gte_50x"] >= 0 and
                            stats["percent_gte_50x"] <= 1)
            self.__verify_dict_field(stats, "percent_lt_15x", float)
            self.assertTrue(stats["percent_lt_15x"] >= 0 and
                            stats["percent_lt_15x"] <= 1)
            self.__verify_dict_field(stats, "total_bases", int)
            self.assertTrue(stats["total_bases"] >= 0)
            self.__verify_dict_field(stats, "weighted_median", float)
            self.assertTrue(stats["weighted_median"] >= 0)
            self.__verify_dict_field(stats, "weighted_pct75", float)
            self.assertTrue(stats["weighted_pct75"] >= 0)
            self.__verify_dict_field(stats, "weighted_pct25", float)
            self.assertTrue(stats["weighted_pct25"] >= 0)
        except AssertionError, e:
            logging.error("Error verifying transcript statistics")
            logging.error(json.dumps(stats, indent=4))
            raise e

    def __verify_exon(self, exon, gene_name, transcript_id, has_gc = True):
        try:
            self.assertEqual(type(exon), dict)
            self.__verify_dict_field(exon, "start", int)
            self.assertTrue(exon["start"] >= 0)
            self.__verify_dict_field(exon, "end", int)
            self.assertTrue(exon["end"] >= 0)
            self.assertTrue(exon["end"] > exon["start"], msg="End <= start")
            if self.__is_padding_enabled():
                self.__verify_dict_field(exon, "padded_start", int)
                self.assertTrue(exon["padded_start"] >= 0)
                self.__verify_dict_field(exon, "padded_end", int)
                self.assertTrue(exon["padded_end"] >= 0)
                self.assertTrue(exon["padded_end"] > exon["padded_start"],
                                msg="Padded end <= padded start")
                self.assertTrue(exon["start"] - exon["padded_start"] == self.config["exon_padding"],
                                msg="Incorrect start coordinate padding")
                self.assertTrue(exon["padded_end"] - exon["end"] == self.config["exon_padding"],
                                msg="Incorrect end coordinate padding")
            self.__verify_dict_field(exon, "exon_number", str)
            self.assertTrue(str(exon["exon_number"]).startswith("exon"),
                            msg="Exon number is not well formed")
            self.__verify_exon_statistics(exon, has_gc)
        except AssertionError, e:
            logging.error("Error verifying exon at %s:%s" % (gene_name, transcript_id))
            logging.error(json.dumps(exon, indent=4))
            raise e

    def __verify_exon_statistics(self, exon, has_gc=True):
        self.__verify_dict_field(exon, "statistics", dict)
        statistics = exon["statistics"]
        try:
            self.__verify_dict_field(statistics, "bases_gte_15x", int)
            self.assertTrue(statistics["bases_gte_15x"] >= 0)
            self.__verify_dict_field(statistics, "bases_gte_30x", int)
            self.assertTrue(statistics["bases_gte_30x"] >= 0)
            self.__verify_dict_field(statistics, "bases_gte_50x", int)
            self.assertTrue(statistics["bases_gte_50x"] >= 0)
            self.__verify_dict_field(statistics, "bases_lt_15x", int)
            self.assertTrue(statistics["bases_lt_15x"] >= 0)
            self.__verify_dict_field(statistics, "bases_lt_3x", int)
            self.assertTrue(statistics["bases_lt_3x"] >= 0)
            if has_gc:
                self.__verify_dict_field(statistics, "gc_content", float)
                self.assertTrue(statistics["gc_content"] >= 0 and
                                statistics["gc_content"] <= 1)
            self.__verify_dict_field(statistics, "mean", float)
            self.assertTrue(statistics["mean"] >= 0)
            self.__verify_dict_field(statistics, "percent_gte_15x", float)
            self.assertTrue(statistics["percent_gte_15x"] >= 0 and
                            statistics["percent_gte_15x"] <= 1)
            self.assertTrue(statistics["percent_gte_30x"] >= 0 and
                            statistics["percent_gte_30x"] <= 1)
            self.assertTrue(statistics["percent_gte_50x"] >= 0 and
                            statistics["percent_gte_50x"] <= 1)
            self.assertTrue(statistics["percent_lt_15x"] >= 0 and
                            statistics["percent_lt_15x"] <= 1)
            self.__verify_dict_field(statistics, "total_bases", int)
            self.assertTrue(statistics["total_bases"] >= 0)
            self.__verify_dict_field(statistics, "median", float)
            self.assertTrue(statistics["median"] >= 0)
            self.__verify_dict_field(statistics, "pct75", float)
            self.assertTrue(statistics["pct75"] >= 0)
            self.__verify_dict_field(statistics, "pct25", float)
            self.assertTrue(statistics["pct25"] >= 0)
        except AssertionError, e:
            logging.error("Error verifying exon statistics")
            logging.error(json.dumps(statistics, indent=4))
            raise e

    def __verify_gap(self, gap, exon):
        try:
            self.assertEqual(type(gap), dict)
            self.__verify_dict_field(gap, "start", int)
            (start, end) = self.__get_exon_coordinates(exon)
            self.assertTrue(gap["start"] >= start and gap["start"] <= end)
            self.__verify_dict_field(gap, "end", int)
            self.assertTrue(gap["end"] >= start and gap["start"] <= end and
                            gap["end"] >= gap["start"])
            self.__verify_dict_field(gap, "length", int)
            self.assertTrue(gap["length"] >= 1 and gap["length"] <= gap["end"] - gap["start"] + 1)
        except AssertionError, e:
            logging.error("Error verifying gap")
            logging.error(json.dumps(gap, indent=4))
            raise e

    def __verify_panel_stats(self, statistics):
        try:
            self.__verify_dict_field(statistics, "bases_gte_15x", int)
            self.assertTrue(statistics["bases_gte_15x"] >= 0)
            self.__verify_dict_field(statistics, "bases_gte_30x", int)
            self.assertTrue(statistics["bases_gte_30x"] >= 0)
            self.__verify_dict_field(statistics, "bases_gte_50x", int)
            self.assertTrue(statistics["bases_gte_50x"] >= 0)
            self.__verify_dict_field(statistics, "bases_lt_15x", int)
            self.assertTrue(statistics["bases_lt_15x"] >= 0)
            self.__verify_dict_field(statistics, "bases_lt_3x", int)
            self.assertTrue(statistics["bases_lt_3x"] >= 0)
            self.__verify_dict_field(statistics, "mean", float)
            self.assertTrue(statistics["mean"] >= 0)
            self.__verify_dict_field(statistics, "percent_gte_15x", float)
            self.__verify_dict_field(statistics, "percent_gte_30x", float)
            self.__verify_dict_field(statistics, "percent_gte_50x", float)
            self.__verify_dict_field(statistics, "percent_lt_15x", float)
            self.assertTrue(statistics["percent_gte_15x"] >= 0 and
                            statistics["percent_gte_15x"] <= 1)
            self.assertTrue(statistics["percent_gte_30x"] >= 0 and
                            statistics["percent_gte_30x"] <= 1)
            self.assertTrue(statistics["percent_gte_50x"] >= 0 and
                            statistics["percent_gte_50x"] <= 1)
            self.assertTrue(statistics["percent_lt_15x"] >= 0 and
                            statistics["percent_lt_15x"] <= 1)
            self.__verify_dict_field(statistics, "total_bases", int)
            self.assertTrue(statistics["total_bases"] >= 0)
            self.__verify_dict_field(statistics, "weighted_median", float)
            self.assertTrue(statistics["weighted_median"] >= 0)
            self.__verify_dict_field(statistics, "weighted_pct75", float)
            self.assertTrue(statistics["weighted_pct75"] >= 0)
            self.__verify_dict_field(statistics, "weighted_pct25", float)
            self.assertTrue(statistics["weighted_pct25"] >= 0)
        except AssertionError, e:
            logging.error("Error panel statistics")
            logging.error(json.dumps(statistics, indent=4))
            raise e