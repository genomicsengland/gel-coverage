import collections
import random
import unittest
import numpy
import logging
import time
import gelcoverage.stats.coverage_stats as coverage_stats
import gelcoverage.stats.sequence_stats as sequence_stats
import gelcoverage.constants as constants
from gelcoverage.tools.bigwig_reader import BigWigReader
from gelcoverage.tools.bed_reader import BedReader
from gelcoverage.test.output_verifier import OutputVerifier


class CoverageStatsTests(OutputVerifier):

    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)
        # Creates 1000 random coverage values between 0 and 400
        self.coverages = [random.randint(a=30, b=400) for p in range(0,1000)]
        self.coverages[100:150] = [random.randint(a=0, b=29) for p in range(0,50)]
        self.coverages[300:310] = [random.randint(a=0, b=29) for p in range(0, 10)]
        self.coverages[340:345] = [random.randint(a=0, b=29) for p in range(0, 5)]

        self.coverages2 = [random.randint(a=0, b=400) for p in range(0, 2000)]
        self.coverages3 = [random.randint(a=0, b=400) for p in range(0, 3000)]

        self.start_position = 1234
        self.coverage_threshold = 30
        self.gap_length_threshold = 1
        self.gc_content = 0.59

        self.bw = "../../resources/test/test1.bw"
        self.bigwig_reader = BigWigReader(self.bw)
        self.wg_regions = "../../resources/Homo_sapiens.GRCh37.75.dna.primary_assembly.NonN_Regions.CHR.prefix.bed"
        self.bed_reader = BedReader(self.wg_regions)

    def _check_stats(self, stats):
        self.assertEqual(type(stats), dict)
        self.assertTrue(constants.BASES in stats)
        self.assertTrue(constants.AVERAGE in stats)
        self.assertTrue(constants.MEDIAN in stats)
        self.assertTrue(constants.PERCENTILE75 in stats)
        self.assertTrue(constants.PERCENTILE25 in stats)
        self.assertTrue(constants.LT15X in stats)
        self.assertTrue(constants.GTE15X in stats)
        self.assertTrue(constants.GTE30X in stats)
        self.assertTrue(constants.GTE50X in stats)
        self.assertTrue(constants.GC_CONTENT in stats)
        self.assertEqual(type(stats[constants.AVERAGE]), float)
        self.assertEqual(type(stats[constants.MEDIAN]), float)
        self.assertEqual(type(stats[constants.PERCENTILE75]), float)
        self.assertEqual(type(stats[constants.BASES]), int)
        self.assertEqual(type(stats[constants.PERCENTILE25]), float)
        self.assertEqual(type(stats[constants.LT15X]), float)
        self.assertEqual(type(stats[constants.GTE15X]), float)
        self.assertEqual(type(stats[constants.GTE30X]), float)
        self.assertEqual(type(stats[constants.GTE50X]), float)
        self.assertEqual(type(stats[constants.GC_CONTENT]), float)

    def test1(self):
        """
        Tests find_gaps(coverages, start_position, coverage_threshold)
        :return:
        """
        gaps = coverage_stats.find_gaps(self.coverages, self.start_position, self.coverage_threshold, self.gap_length_threshold)
        self.assertEqual(type(gaps), list)
        self.assertEqual(len(gaps), 3)
        for gap in gaps:
            self.assertEqual(type(gap), dict)
            self.assertTrue(constants.GAP_START in gap)
            self.assertTrue(constants.GAP_END in gap)
            self.assertEqual(type(gap[constants.GAP_START]), int)
            self.assertEqual(type(gap[constants.GAP_END]), int)
            print "Found a gap at %s-%s" % (str(gap[constants.GAP_START]), str(gap[constants.GAP_END]))

    def test2(self):
        """
        compute_exon_level_statistics(coverages, gc_content)
        :return:
        """
        stats = coverage_stats.compute_exon_level_statistics(self.coverages, self.gc_content)
        self._check_stats(stats)

    def test3(self):
        """
        compute_transcript_level_statistics(exons)
        :return:
        """
        exon1 = coverage_stats.compute_exon_level_statistics(self.coverages, self.gc_content)
        exon2 = coverage_stats.compute_exon_level_statistics(self.coverages2, self.gc_content)
        exon3 = coverage_stats.compute_exon_level_statistics(self.coverages3, self.gc_content)
        exons = [{constants.STATISTICS:exon1}, {constants.STATISTICS:exon2}, {constants.STATISTICS:exon3}]
        stats = coverage_stats.compute_transcript_level_statistics(exons)
        self._check_stats(stats)

    def test4(self):
        """
        compute_coding_region_statistics(genes)
        :return:
        """
        coverage1_1 = [random.randint(a=0, b=400) for p in range(0, 200)]
        coverage1_2 = [random.randint(a=0, b=400) for p in range(0, 200)]
        coverage1_3 = [random.randint(a=0, b=400) for p in range(0, 200)]
        coverage2_1 = [random.randint(a=0, b=400) for p in range(0, 200)]
        coverage2_2 = [random.randint(a=0, b=400) for p in range(0, 200)]
        coverage2_3 = [random.randint(a=0, b=400) for p in range(0, 200)]
        coverage3_1 = [random.randint(a=0, b=400) for p in range(0, 200)]
        coverage3_2 = [random.randint(a=0, b=400) for p in range(0, 200)]
        coverage3_3 = [random.randint(a=0, b=400) for p in range(0, 200)]

        exon1_1 = coverage_stats.compute_exon_level_statistics(coverage1_1, 0.22)
        exon1_2 = coverage_stats.compute_exon_level_statistics(coverage1_2, 0.23)
        exon1_3 = coverage_stats.compute_exon_level_statistics(coverage1_3, 0.24)
        exon2_1 = coverage_stats.compute_exon_level_statistics(coverage2_1, 0.22)
        exon2_2 = coverage_stats.compute_exon_level_statistics(coverage2_2, 0.23)
        exon2_3 = coverage_stats.compute_exon_level_statistics(coverage2_3, 0.24)
        exon3_1 = coverage_stats.compute_exon_level_statistics(coverage3_1, 0.22)
        exon3_2 = coverage_stats.compute_exon_level_statistics(coverage3_2, 0.23)
        exon3_3 = coverage_stats.compute_exon_level_statistics(coverage3_3, 0.24)

        transcript1 = coverage_stats.compute_transcript_level_statistics([
            {constants.STATISTICS: exon1_1}, {constants.STATISTICS: exon1_2}, {constants.STATISTICS: exon1_3}
        ])
        transcript2 = coverage_stats.compute_transcript_level_statistics([
            {constants.STATISTICS: exon2_1}, {constants.STATISTICS: exon2_2}, {constants.STATISTICS: exon2_3}
        ])
        transcript3 = coverage_stats.compute_transcript_level_statistics([
            {constants.STATISTICS: exon3_1}, {constants.STATISTICS: exon3_2}, {constants.STATISTICS: exon3_3}
        ])

        gene1 = {
            constants.UNION_TRANSCRIPT: {
                constants.STATISTICS: transcript1
            },
            constants.CHROMOSOME: "chr1"
        }
        gene2 = {
            constants.UNION_TRANSCRIPT: {
                constants.STATISTICS: transcript2
            },
            constants.CHROMOSOME: "chr2"
        }
        gene3 = {
            constants.UNION_TRANSCRIPT: {
                constants.STATISTICS: transcript3
            },
            constants.CHROMOSOME: "chrX"
        }

        coding_region_stats = coverage_stats.compute_coding_region_statistics([gene1, gene2, gene3])
        stats = coding_region_stats[constants.STATISTICS]
        self._check_stats(stats)
        self.assertEqual(len(coding_region_stats[constants.CHROMOSOMES]), 4)
        for chromosome_stats in coding_region_stats[constants.CHROMOSOMES]:
            self._check_stats(chromosome_stats)

    def test5(self):
        """
        Compute whole genome statistics
        :return:
        """
        start_time_2 = time.time()
        results = coverage_stats.compute_whole_genome_statistics(self.bigwig_reader, self.bed_reader)
        execution_time_2 = time.time() - start_time_2
        print execution_time_2
        print results
        print "hey"

        # start_time_2 = time.time()
        # results = coverage_stats.compute_whole_genome_statistics(self.bigwig_reader, self.bed_reader,
        #                                                          chunk_size=1000000)
        # execution_time_2 = time.time() - start_time_2
        # print execution_time_2
        #
        # start_time_3 = time.time()
        # results = coverage_stats.compute_whole_genome_statistics(self.bigwig_reader, self.bed_reader,
        #                                                          chunk_size=10000000)
        # execution_time_3 = time.time() - start_time_3
        # print execution_time_3
        #
        # start_time_4 = time.time()
        # results = coverage_stats.compute_whole_genome_statistics(self.bigwig_reader, self.bed_reader,
        #                                                          chunk_size=10000)
        # execution_time_4 = time.time() - start_time_4
        # print execution_time_4
        #
        # start_time = time.time()
        # results = coverage_stats.compute_whole_genome_statistics(self.bigwig_reader, self.bed_reader, numpy=True)
        # execution_time = time.time() - start_time
        # print execution_time
        #
        # start_time_2 = time.time()
        # results = coverage_stats.compute_whole_genome_statistics(self.bigwig_reader, self.bed_reader,
        #                                                          chunk_size=1000000, numpy=True)
        # execution_time_2 = time.time() - start_time_2
        # print execution_time_2
        #
        # start_time_3 = time.time()
        # results = coverage_stats.compute_whole_genome_statistics(self.bigwig_reader, self.bed_reader,
        #                                                          chunk_size=10000000, numpy=True)
        # execution_time_3 = time.time() - start_time_3
        # print execution_time_3
        #
        # start_time_4 = time.time()
        # results = coverage_stats.compute_whole_genome_statistics(self.bigwig_reader, self.bed_reader,
        #                                                          chunk_size=10000, numpy=True)
        # execution_time_4 = time.time() - start_time_4
        # print execution_time_4

        self._verify_dict_field(results, constants.STATISTICS, dict)
        self._verify_dict_field(results, constants.CHROMOSOMES, list)
        self._verify_wg_stats(results[constants.STATISTICS])
        found_autosomes = False
        for chr_stats in results[constants.CHROMOSOMES]:
            self._verify_wg_stats(chr_stats)
            if chr_stats[constants.CHROMOSOME] == constants.AUTOSOMES:
                found_autosomes = True
        self.assertTrue(found_autosomes, "No aggregated stats in whole genome for autosomes")

    def test5_delete(self):
        """
        compute_exon_level_statistics(coverages, gc_content)
        :return:
        """
        import time
        start_time = time.time()
        stats = coverage_stats.compute_exon_level_statistics(self.coverages, self.gc_content)
        execution_time = time.time() - start_time
        print execution_time
        start_time = time.time()
        stats2 = coverage_stats.compute_exon_level_statistics_2(self.coverages, self.gc_content)
        execution_time = time.time() - start_time
        print execution_time
        self.assertEqual(stats, stats2)

    def test_delete_2(self):
        """
        compute_exon_level_statistics(coverages, gc_content)
        :return:
        """
        exon1 = coverage_stats.compute_exon_level_statistics(self.coverages, self.gc_content)
        exon2 = coverage_stats.compute_exon_level_statistics(self.coverages2, self.gc_content)
        exon3 = coverage_stats.compute_exon_level_statistics(self.coverages3, self.gc_content)
        exons = [{constants.STATISTICS: exon1}, {constants.STATISTICS: exon2}, {constants.STATISTICS: exon3}]
        import time
        start_time = time.time()
        stats = coverage_stats.compute_transcript_level_statistics(exons)
        execution_time = time.time() - start_time
        print execution_time
        print stats
        #start_time = time.time()
        #stats2 = coverage_stats.compute_transcript_level_statistics_2(exons)
        #execution_time_2 = time.time() - start_time
        #print execution_time
        #self.assertEqual(stats, stats2)




class SequenceStatsTests(unittest.TestCase):

    def setUp(self):
        logging.basicConfig(level=logging.DEBUG)
        # Creates random sequence with a GC content close 0.6
        self.sequence = numpy.random.choice(["G", "C", "A", "T"], 10000, p=[0.3, 0.3, 0.2, 0.2])

    def test1(self):
        """
        Tests find_gaps(coverages, start_position, coverage_threshold)
        :return:
        """
        gc_content = sequence_stats.compute_gc_content(self.sequence)
        print "Found a GC content of %s" % gc_content
        self.assertEqual(type(gc_content), float)
        self.assertTrue(gc_content <= 0.65 and gc_content >= 0.55)

if __name__ == '__main__':
    unittest.main()
