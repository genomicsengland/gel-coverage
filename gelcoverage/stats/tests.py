import collections
import random
import unittest
import numpy
import logging
import gelcoverage.stats.coverage_stats as coverage_stats
import gelcoverage.stats.sequence_stats as sequence_stats


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
            self.assertTrue("s" in gap)
            self.assertTrue("e" in gap)
            self.assertEqual(type(gap["s"]), int)
            self.assertEqual(type(gap["e"]), int)
            print "Found a gap at %s-%s" % (str(gap["s"]), str(gap["e"]))

    def test2(self):
        """
        compute_exon_level_statistics(coverages, gc_content)
        :return:
        """
        stats = coverage_stats.compute_exon_level_statistics(self.coverages, self.gc_content)
        self.assertEqual(type(stats), dict)
        self.assertTrue("bases" in stats)
        self.assertTrue("avg" in stats)
        self.assertTrue("med" in stats)
        self.assertTrue("pct75" in stats)
        self.assertTrue("pct25" in stats)
        self.assertTrue("bases_lt_15x" in stats)
        self.assertTrue("bases_gte_15x" in stats)
        self.assertTrue("bases_gte_30x" in stats)
        self.assertTrue("bases_gte_50x" in stats)
        self.assertTrue("%<15x" in stats)
        self.assertTrue("%>=15x" in stats)
        self.assertTrue("%>=30x" in stats)
        self.assertTrue("%>=50x" in stats)
        self.assertTrue("gc" in stats)
        self.assertEqual(type(stats["avg"]), float)
        self.assertEqual(type(stats["med"]), float)
        self.assertEqual(type(stats["pct75"]), float)
        self.assertEqual(type(stats["bases"]), int)
        self.assertEqual(type(stats["pct25"]), float)
        self.assertEqual(type(stats["bases_lt_15x"]), int)
        self.assertEqual(type(stats["bases_gte_15x"]), int)
        self.assertEqual(type(stats["bases_gte_30x"]), int)
        self.assertEqual(type(stats["bases_gte_50x"]), int)
        self.assertEqual(type(stats["%<15x"]), float)
        self.assertEqual(type(stats["%>=15x"]), float)
        self.assertEqual(type(stats["%>=30x"]), float)
        self.assertEqual(type(stats["%>=50x"]), float)
        self.assertEqual(type(stats["gc"]), float)

    def test3(self):
        """
        compute_transcript_level_statistics(exons)
        :return:
        """
        exon1 = coverage_stats.compute_exon_level_statistics(self.coverages, self.gc_content)
        exon2 = coverage_stats.compute_exon_level_statistics(self.coverages2, self.gc_content)
        exon3 = coverage_stats.compute_exon_level_statistics(self.coverages3, self.gc_content)
        exons = [{"stats":exon1}, {"stats":exon2}, {"stats":exon3}]
        stats = coverage_stats.compute_transcript_level_statistics(exons)
        self.assertEqual(type(stats), dict)
        self.assertTrue("bases" in stats)
        self.assertTrue("avg" in stats)
        self.assertTrue("med" in stats)
        self.assertTrue("pct75" in stats)
        self.assertTrue("pct25" in stats)
        self.assertTrue("bases_lt_15x" in stats)
        self.assertTrue("bases_gte_15x" in stats)
        self.assertTrue("bases_gte_30x" in stats)
        self.assertTrue("bases_gte_50x" in stats)
        self.assertTrue("%<15x" in stats)
        self.assertTrue("%>=15x" in stats)
        self.assertTrue("%>=30x" in stats)
        self.assertTrue("%>=50x" in stats)
        self.assertTrue("gc" in stats)
        self.assertEqual(type(stats["avg"]), float)
        self.assertEqual(type(stats["med"]), float)
        self.assertEqual(type(stats["pct75"]), float)
        self.assertEqual(type(stats["bases"]), int)
        self.assertEqual(type(stats["pct25"]), float)
        self.assertEqual(type(stats["bases_lt_15x"]), int)
        self.assertEqual(type(stats["bases_gte_15x"]), int)
        self.assertEqual(type(stats["bases_gte_30x"]), int)
        self.assertEqual(type(stats["bases_gte_50x"]), int)
        self.assertEqual(type(stats["%<15x"]), float)
        self.assertEqual(type(stats["%>=15x"]), float)
        self.assertEqual(type(stats["%>=30x"]), float)
        self.assertEqual(type(stats["%>=50x"]), float)
        self.assertEqual(type(stats["gc"]), float)

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
