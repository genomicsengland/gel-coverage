import collections
import random
import unittest
import numpy
import gelcoverage.stats.coverage_stats as coverage_stats
import gelcoverage.stats.sequence_stats as sequence_stats

class CoverageStatsTests(unittest.TestCase):

    def setUp(self):
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
            self.assertEqual(type(gap), collections.defaultdict)
            self.assertTrue("start" in gap)
            self.assertTrue("end" in gap)
            self.assertEqual(type(gap["start"]), int)
            self.assertEqual(type(gap["end"]), int)
            print "Found a gap at %s-%s" % (str(gap["start"]), str(gap["end"]))

    def test2(self):
        """
        compute_exon_level_statistics(coverages, gc_content)
        :return:
        """
        stats = coverage_stats.compute_exon_level_statistics(self.coverages, self.gc_content)
        self.assertEqual(type(stats), collections.defaultdict)
        self.assertTrue("total_bases" in stats)
        self.assertTrue("mean" in stats)
        self.assertTrue("median" in stats)
        self.assertTrue("pct75" in stats)
        self.assertTrue("pct25" in stats)
        self.assertTrue("bases_lt_3x" in stats)
        self.assertTrue("bases_lt_15x" in stats)
        self.assertTrue("bases_gte_15x" in stats)
        self.assertTrue("bases_gte_30x" in stats)
        self.assertTrue("bases_gte_50x" in stats)
        self.assertTrue("percent_lt_15x" in stats)
        self.assertTrue("percent_gte_15x" in stats)
        self.assertTrue("percent_gte_30x" in stats)
        self.assertTrue("percent_gte_50x" in stats)
        self.assertTrue("gc_content" in stats)
        self.assertEqual(type(stats["mean"]), float)
        self.assertEqual(type(stats["median"]), float)
        self.assertEqual(type(stats["pct75"]), float)
        self.assertEqual(type(stats["total_bases"]), int)
        self.assertEqual(type(stats["pct25"]), float)
        self.assertEqual(type(stats["bases_lt_3x"]), int)
        self.assertEqual(type(stats["bases_lt_15x"]), int)
        self.assertEqual(type(stats["bases_gte_15x"]), int)
        self.assertEqual(type(stats["bases_gte_30x"]), int)
        self.assertEqual(type(stats["bases_gte_50x"]), int)
        self.assertEqual(type(stats["percent_lt_15x"]), float)
        self.assertEqual(type(stats["percent_gte_15x"]), float)
        self.assertEqual(type(stats["percent_gte_30x"]), float)
        self.assertEqual(type(stats["percent_gte_50x"]), float)
        self.assertEqual(type(stats["gc_content"]), float)

    def test3(self):
        """
        compute_transcript_level_statistics(exons)
        :return:
        """
        exon1 = coverage_stats.compute_exon_level_statistics(self.coverages, self.gc_content)
        exon2 = coverage_stats.compute_exon_level_statistics(self.coverages2, self.gc_content)
        exon3 = coverage_stats.compute_exon_level_statistics(self.coverages3, self.gc_content)
        exons = {1:exon1, 2:exon2, 3:exon3}
        stats = coverage_stats.compute_transcript_level_statistics(exons)
        self.assertEqual(type(stats), collections.defaultdict)
        self.assertTrue("total_bases" in stats)
        self.assertTrue("mean" in stats)
        self.assertTrue("weighted_median" in stats)
        self.assertTrue("weighted_pct75" in stats)
        self.assertTrue("weighted_pct25" in stats)
        self.assertTrue("bases_lt_3x" in stats)
        self.assertTrue("bases_lt_15x" in stats)
        self.assertTrue("bases_gte_15x" in stats)
        self.assertTrue("bases_gte_30x" in stats)
        self.assertTrue("bases_gte_50x" in stats)
        self.assertTrue("percent_lt_15x" in stats)
        self.assertTrue("percent_gte_15x" in stats)
        self.assertTrue("percent_gte_30x" in stats)
        self.assertTrue("percent_gte_50x" in stats)
        self.assertTrue("gc_content" in stats)
        self.assertEqual(type(stats["mean"]), float)
        self.assertEqual(type(stats["weighted_median"]), float)
        self.assertEqual(type(stats["weighted_pct75"]), float)
        self.assertEqual(type(stats["total_bases"]), int)
        self.assertEqual(type(stats["weighted_pct25"]), float)
        self.assertEqual(type(stats["bases_lt_3x"]), int)
        self.assertEqual(type(stats["bases_lt_15x"]), int)
        self.assertEqual(type(stats["bases_gte_15x"]), int)
        self.assertEqual(type(stats["bases_gte_30x"]), int)
        self.assertEqual(type(stats["bases_gte_50x"]), int)
        self.assertEqual(type(stats["percent_lt_15x"]), float)
        self.assertEqual(type(stats["percent_gte_15x"]), float)
        self.assertEqual(type(stats["percent_gte_30x"]), float)
        self.assertEqual(type(stats["percent_gte_50x"]), float)
        self.assertEqual(type(stats["gc_content"]), float)


class SequenceStatsTests(unittest.TestCase):

    def setUp(self):
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
