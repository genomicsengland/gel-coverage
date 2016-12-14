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
        if json["parameters"]["wg_stats_enabled"]:
            self.__verify_dict_field(json["results"], "whole_genome", dict)
            self.__verify_dict_field(json["results"]["whole_genome"], "stats", dict)
            self.__verify_dict_field(json["results"]["whole_genome"], "chrs", list)
            self.__verify_wg_stats(json["results"]["whole_genome"]["stats"])
            for chr_stats in json["results"]["whole_genome"]["chrs"]:
                self.__verify_wg_stats(chr_stats)
        self.__verify_dict_field(json["results"], "coding_region", dict)
        self.__verify_dict_field(json["results"]["coding_region"], "stats", dict)
        self.__verify_dict_field(json["results"]["coding_region"], "chrs", list)
        self.__verify_panel_stats(json["results"]["coding_region"]["stats"])
        for chr_stats in json["results"]["coding_region"]["chrs"]:
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
            self.__verify_dict_field(uncovered_gene, "chr", str)
            self.__verify_dict_field(uncovered_gene, "name", str)
            self.assertTrue(uncovered_gene["name"] not in observed_genes,
                                msg="Duplicated gene '%s'" % uncovered_gene["name"])
            observed_genes.append(uncovered_gene["name"])

    def __verify_genes(self, results, exon_stats_enabled):
        # Verify every gene
        self.__verify_dict_field(results, "genes", list)
        observed_genes = []
        for gene in results["genes"]:
            if self.expected_gene_list is not None:
                self.assertTrue(gene["name"] in self.expected_gene_list,
                                msg="Unexpected gene found in results '%s'" % gene["name"])
            self.__verify_dict_field(gene, "chr", str)
            # Checks that genes are not repeated
            self.assertTrue(gene["name"] not in observed_genes,
                            msg="Duplicated gene '%s'" % gene['name'])
            observed_genes.append(gene["name"])
            # Verify every transcript
            observed_transcripts = []
            for transcript in gene["trs"]:
                self.__verify_transcript(transcript)
                self.assertTrue(transcript["id"] not in observed_transcripts,
                                msg="Duplicated transcript '%s'" % transcript['id'])
                observed_transcripts.append(transcript["id"])
                # Verify every exon
                if exon_stats_enabled:
                    observed_exons = []
                    for exon in transcript["exons"]:
                        self.__verify_exon(exon, gene["name"], transcript["id"])
                        self.assertTrue(exon["exon"] not in observed_exons,
                                        msg="Duplicated exon '%s'" % exon['exon'])
                        observed_exons.append(exon["exon"])
                        # Verify gaps
                        observed_gaps = []
                        for gap in exon["gaps"]:
                            self.__verify_gap(gap, exon)
                            self.assertTrue((gap["s"], gap["e"]) not in observed_gaps,
                                            msg="Repeated gap '%s-%s'" % (
                                                gap["s"], gap["e"]
                                            ))
                            observed_gaps.append((gap["s"], gap["e"]))
                else:
                    self.assertTrue("exons" not in transcript)
            union_transcript = gene["union_tr"]
            self.__verify_transcript_stats(union_transcript["stats"], has_gc=False)
            self.verify_union_transcript(gene, exon_stats_enabled)

    def __verify_transcript(self, transcript):
        self.__verify_dict_field(transcript, "id", str)
        self.assertTrue(str(transcript["id"]).startswith("ENS"),
                        msg="Wrong transcript id '%s'" % transcript["id"])
        self.__verify_transcript_stats(transcript["stats"])

    def verify_union_transcript(self, gene, exon_stats_enabled):
        # Basic checks
        self.__verify_dict_field(gene, "union_tr", dict)
        union_transcript = gene["union_tr"]
        self.__verify_transcript_stats(gene["union_tr"]["stats"], has_gc=False)
        # Verifies exons
        if exon_stats_enabled:
            for ut_exon in union_transcript["exons"]:
                self.__verify_exon(ut_exon, gene["name"], "union_tr", has_gc=False)
                # Verify gaps
                for gap in ut_exon["gaps"]:
                    self.__verify_gap(gap, ut_exon)
            # Verifies union transcript build up
            all_exons = sum([transcript["exons"] for transcript in gene["trs"]], [])
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
        else:
            self.assertTrue("exons" not in union_transcript)

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
            start = exon["padded_s"]
            end = exon["padded_e"]
        else:
            start = exon["s"]
            end = exon["e"]
        return (start, end)

    def __is_padding_enabled(self):
        return self.config["exon_padding"] > 0

    def __verify_transcript_stats(self, stats, has_gc = True):
        try:
            self.assertEqual(type(stats), dict)
            if has_gc:
                self.__verify_dict_field(stats, "gc", float)
                self.assertTrue(stats["gc"] >= 0 and
                                stats["gc"] <= 1)
            self.__verify_dict_field(stats, "avg", float)
            self.assertTrue(stats["avg"] >= 0)
            self.__verify_dict_field(stats, "%>=15x", float)
            self.assertTrue(stats["%>=15x"] >= 0 and
                            stats["%>=15x"] <= 1)
            self.__verify_dict_field(stats, "%>=30x", float)
            self.assertTrue(stats["%>=30x"] >= 0 and
                            stats["%>=30x"] <= 1)
            self.__verify_dict_field(stats, "%>=50x", float)
            self.assertTrue(stats["%>=50x"] >= 0 and
                            stats["%>=50x"] <= 1)
            self.__verify_dict_field(stats, "%<15x", float)
            self.assertTrue(stats["%<15x"] >= 0 and
                            stats["%<15x"] <= 1)
            self.__verify_dict_field(stats, "bases", [int, long])
            self.assertTrue(stats["bases"] >= 0)
            self.__verify_dict_field(stats, "med", float)
            self.assertTrue(stats["med"] >= 0)
            self.__verify_dict_field(stats, "pct75", float)
            self.assertTrue(stats["pct75"] >= 0)
            self.__verify_dict_field(stats, "pct25", float)
            self.assertTrue(stats["pct25"] >= 0)
        except AssertionError, e:
            logging.error("Error verifying transcript statistics")
            logging.error(json.dumps(stats, indent=4))
            raise e

    def __verify_exon(self, exon, gene_name, transcript_id, has_gc = True):
        try:
            self.assertEqual(type(exon), dict)
            self.__verify_dict_field(exon, "s", int)
            self.assertTrue(exon["s"] >= 0)
            self.__verify_dict_field(exon, "e", int)
            self.assertTrue(exon["e"] >= 0)
            self.assertTrue(exon["e"] >= exon["s"], msg="End < start")
            if self.__is_padding_enabled():
                self.__verify_dict_field(exon, "padded_s", int)
                self.assertTrue(exon["padded_s"] >= 0)
                self.__verify_dict_field(exon, "padded_e", int)
                self.assertTrue(exon["padded_e"] >= 0)
                self.assertTrue(exon["padded_e"] > exon["padded_s"],
                                msg="Padded end <= padded start")
                self.assertTrue(exon["s"] - exon["padded_s"] == self.config["exon_padding"],
                                msg="Incorrect start coordinate padding")
                self.assertTrue(exon["padded_e"] - exon["e"] == self.config["exon_padding"],
                                msg="Incorrect end coordinate padding")
            self.__verify_dict_field(exon, "exon", str)
            self.assertTrue(str(exon["exon"]).startswith("exon"),
                            msg="Exon number is not well formed")
            self.__verify_exon_statistics(exon, has_gc)
        except AssertionError, e:
            logging.error("Error verifying exon at %s:%s" % (gene_name, transcript_id))
            logging.error(json.dumps(exon, indent=4))
            raise e

    def __verify_exon_statistics(self, exon, has_gc=True):
        self.__verify_dict_field(exon, "stats", dict)
        statistics = exon["stats"]
        try:
            if has_gc:
                self.__verify_dict_field(statistics, "gc", float)
                self.assertTrue(statistics["gc"] >= 0 and
                                statistics["gc"] <= 1)
            self.__verify_dict_field(statistics, "avg", float)
            self.assertTrue(statistics["avg"] >= 0)
            self.__verify_dict_field(statistics, "%>=15x", float)
            self.assertTrue(statistics["%>=15x"] >= 0 and
                            statistics["%>=15x"] <= 1)
            self.assertTrue(statistics["%>=30x"] >= 0 and
                            statistics["%>=30x"] <= 1)
            self.assertTrue(statistics["%>=50x"] >= 0 and
                            statistics["%>=50x"] <= 1)
            self.assertTrue(statistics["%<15x"] >= 0 and
                            statistics["%<15x"] <= 1)
            self.__verify_dict_field(statistics, "bases", [int, long])
            self.assertTrue(statistics["bases"] >= 0)
            self.__verify_dict_field(statistics, "med", float)
            self.assertTrue(statistics["med"] >= 0)
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
            self.__verify_dict_field(gap, "s", int)
            (start, end) = self.__get_exon_coordinates(exon)
            self.assertTrue(gap["s"] >= start and gap["s"] <= end)
            self.__verify_dict_field(gap, "e", int)
            self.assertTrue(gap["e"] >= start and gap["s"] <= end and
                            gap["e"] >= gap["s"])
            self.__verify_dict_field(gap, "l", int)
            self.assertTrue(gap["l"] >= 1 and gap["l"] <= gap["e"] - gap["s"] + 1)
        except AssertionError, e:
            logging.error("Error verifying gap")
            logging.error(json.dumps(gap, indent=4))
            raise e

    def __verify_panel_stats(self, statistics):
        try:
            self.__verify_dict_field(statistics, "avg", float)
            self.assertTrue(statistics["avg"] >= 0)
            self.__verify_dict_field(statistics, "%>=15x", float)
            self.__verify_dict_field(statistics, "%>=30x", float)
            self.__verify_dict_field(statistics, "%>=50x", float)
            self.__verify_dict_field(statistics, "%<15x", float)
            self.assertTrue(statistics["%>=15x"] >= 0 and
                            statistics["%>=15x"] <= 1)
            self.assertTrue(statistics["%>=30x"] >= 0 and
                            statistics["%>=30x"] <= 1)
            self.assertTrue(statistics["%>=50x"] >= 0 and
                            statistics["%>=50x"] <= 1)
            self.assertTrue(statistics["%>=15x"] >= 0 and
                            statistics["%>=15x"] <= 1)
            self.__verify_dict_field(statistics, "bases", [int, long])
            self.assertTrue(statistics["bases"] >= 0)
            self.__verify_dict_field(statistics, "med", float)
            self.assertTrue(statistics["med"] >= 0)
            self.__verify_dict_field(statistics, "pct75", float)
            self.assertTrue(statistics["pct75"] >= 0)
            self.__verify_dict_field(statistics, "pct25", float)
            self.assertTrue(statistics["pct25"] >= 0)
        except AssertionError, e:
            logging.error("Error panel statistics")
            logging.error(json.dumps(statistics, indent=4))
            raise e

    def __verify_wg_stats(self, statistics):
        try:
            self.__verify_dict_field(statistics, "avg", float)
            self.assertTrue(statistics["avg"] >= 0)
            self.__verify_dict_field(statistics, "%>=15x", float)
            self.__verify_dict_field(statistics, "%>=30x", float)
            self.__verify_dict_field(statistics, "%>=50x", float)
            self.__verify_dict_field(statistics, "%<15x", float)
            self.assertTrue(statistics["%>=15x"] >= 0 and
                            statistics["%>=15x"] <= 1)
            self.assertTrue(statistics["%>=30x"] >= 0 and
                            statistics["%>=30x"] <= 1)
            self.assertTrue(statistics["%>=50x"] >= 0 and
                            statistics["%>=50x"] <= 1)
            self.assertTrue(statistics["%<15x"] >= 0 and
                            statistics["%<15x"] <= 1)
            self.__verify_dict_field(statistics, "bases", [int, long])
            self.assertTrue(statistics["bases"] >= 0)
            self.__verify_dict_field(statistics, "med", float)
            self.assertTrue(statistics["med"] >= 0)
            self.__verify_dict_field(statistics, "pct75", float)
            self.assertTrue(statistics["pct75"] >= 0)
            self.__verify_dict_field(statistics, "pct25", float)
            self.assertTrue(statistics["pct25"] >= 0)
            self.__verify_dict_field(statistics, "uneveness", float)
            self.assertTrue(statistics["uneveness"] >= 0)
        except AssertionError, e:
            logging.error("Error panel statistics")
            logging.error(json.dumps(statistics, indent=4))
            raise e