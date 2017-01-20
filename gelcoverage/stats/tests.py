import collections
import random
import unittest
import numpy
import logging
import gelcoverage.stats.coverage_stats as coverage_stats
import gelcoverage.stats.sequence_stats as sequence_stats
import gelcoverage.constants as constants


class CoverageStatsTests(unittest.TestCase):

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
        self.gc_content = 0.59

    def test1(self):
        """
        Tests find_gaps(coverages, start_position, coverage_threshold)
        :return:
        """
        gaps = coverage_stats.find_gaps(self.coverages, self.start_position, self.coverage_threshold)
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
        self.assertEqual(type(stats), dict)
        self.assertTrue(constants.BASES in stats)
        self.assertTrue(constants.AVERAGE in stats)
        self.assertTrue(constants.MEDIAN in stats)
        self.assertTrue(constants.PERCENTILE75 in stats)
        self.assertTrue(constants.PERCENTILE25 in stats)
        self.assertTrue(constants.BASES_LT15X in stats)
        self.assertTrue(constants.BASES_GTE15X in stats)
        self.assertTrue(constants.BASES_GTE30X in stats)
        self.assertTrue(constants.BASES_GTE50X in stats)
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
        self.assertEqual(type(stats[constants.BASES_LT15X]), int)
        self.assertEqual(type(stats[constants.BASES_GTE15X]), int)
        self.assertEqual(type(stats[constants.BASES_GTE30X]), int)
        self.assertEqual(type(stats[constants.BASES_GTE50X]), int)
        self.assertEqual(type(stats[constants.LT15X]), float)
        self.assertEqual(type(stats[constants.GTE15X]), float)
        self.assertEqual(type(stats[constants.GTE30X]), float)
        self.assertEqual(type(stats[constants.GTE50X]), float)
        self.assertEqual(type(stats[constants.GC_CONTENT]), float)

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
        self.assertEqual(type(stats), dict)
        self.assertTrue(constants.BASES in stats)
        self.assertTrue(constants.AVERAGE in stats)
        self.assertTrue(constants.MEDIAN in stats)
        self.assertTrue(constants.PERCENTILE75 in stats)
        self.assertTrue(constants.PERCENTILE25 in stats)
        self.assertTrue(constants.BASES_LT15X in stats)
        self.assertTrue(constants.BASES_GTE15X in stats)
        self.assertTrue(constants.BASES_GTE30X in stats)
        self.assertTrue(constants.BASES_GTE50X in stats)
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
        self.assertEqual(type(stats[constants.BASES_LT15X]), int)
        self.assertEqual(type(stats[constants.BASES_GTE15X]), int)
        self.assertEqual(type(stats[constants.BASES_GTE30X]), int)
        self.assertEqual(type(stats[constants.BASES_GTE50X]), int)
        self.assertEqual(type(stats[constants.LT15X]), float)
        self.assertEqual(type(stats[constants.GTE15X]), float)
        self.assertEqual(type(stats[constants.GTE30X]), float)
        self.assertEqual(type(stats[constants.GTE50X]), float)
        self.assertEqual(type(stats[constants.GC_CONTENT]), float)

    # TODO: test the panel level statistics

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
