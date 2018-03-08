import factory.fuzzy
from mock import patch, Mock
from gelcoverage.runner import GelCoverageRunner
import gelcoverage.tools.bigwig_reader
import gelcoverage.constants as constants
import protocols.coverage_0_1_0
from protocols.util.factories.avro_factory import GenericFactoryAvro, FactoryAvro


class FuzzyFloatWithPrecision(factory.fuzzy.FuzzyFloat):
    """Random float within a given range."""

    def __init__(self, low, high=None, precision=3, **kwargs):

        self.precision = precision

        super(FuzzyFloatWithPrecision, self).__init__(low, high, **kwargs)

    def fuzz(self):
        return round(super(FuzzyFloatWithPrecision, self).fuzz(), self.precision)


class RegionStatisticsFactory(FactoryAvro):
    def __init__(self, *args, **kwargs):
        super(RegionStatisticsFactory, self).__init__(*args, **kwargs)

    class Meta:
        model = protocols.coverage_0_1_0.RegionStatistics

    _version = '6.0.0'

    avg = FuzzyFloatWithPrecision(0, 100)
    med = FuzzyFloatWithPrecision(0, 100)
    sd = FuzzyFloatWithPrecision(0, 100)
    gc = FuzzyFloatWithPrecision(0, 100)
    pct75 = FuzzyFloatWithPrecision(0, 100)
    pct25 = FuzzyFloatWithPrecision(0, 100)
    bases = factory.fuzzy.FuzzyInteger(0, 100)
    gte50x = FuzzyFloatWithPrecision(0, 1, precision=5)
    gte30x = FuzzyFloatWithPrecision(0, 1, precision=5)
    gte15x = FuzzyFloatWithPrecision(0, 1, precision=5)
    lt15x = FuzzyFloatWithPrecision(0, 1, precision=5)
    rmsd = FuzzyFloatWithPrecision(0, 100)
    bases_gte_50x = factory.fuzzy.FuzzyInteger(50, 200)
    bases_gte_30x = FuzzyFloatWithPrecision(50, 200)
    bases_gte_15x = FuzzyFloatWithPrecision(50, 200)
    bases_lt_15x = FuzzyFloatWithPrecision(50, 200)


class CoverageGapFactory(FactoryAvro):
    def __init__(self, *args, **kwargs):
        super(CoverageGapFactory, self).__init__(*args, **kwargs)

    class Meta:
        model = protocols.coverage_0_1_0.CoverageGap

    _version = '6.0.0'

    s = factory.fuzzy.FuzzyInteger(10000, 50000000)
    e = factory.fuzzy.FuzzyInteger(10000, 50000000)


def mocked_compute_statistics(*args):
    region_statistics_factory = GenericFactoryAvro.get_factory_avro(
        protocols.coverage_0_1_0.RegionStatistics,
        version='6.0.0'
    )
    region_statistics = region_statistics_factory.create()
    return region_statistics.toJsonDict()


def mocked_compute_coding_region_statistics(genes):
    chromosomes = {}
    for gene in genes:
        chromosome = gene[constants.CHROMOSOME]
        if chromosome not in chromosomes:
            stats = mocked_compute_statistics()
            stats[constants.CHROMOSOME] = chromosome
            chromosomes[chromosome] = stats
    chromosomes[constants.AUTOSOMES] = mocked_compute_statistics()
    results = {
        constants.STATISTICS: mocked_compute_statistics(),
        constants.CHROMOSOMES: chromosomes
    }
    return results


def mocked_find_gaps(coverages, start_position, coverage_threshold, gap_length_threshold):
    """
    Creates between 0 and 3 gaps of length between length threshold and 50bp
    :param coverages:
    :param start_position:
    :param coverage_threshold:
    :param gap_length_threshold:
    :return:
    """
    gap_dicts = []
    if coverage_threshold > 0:
        gap_factory = GenericFactoryAvro.get_factory_avro(
            protocols.coverage_0_1_0.CoverageGap,
            version='6.0.0'
        )
        gaps = gap_factory.create_batch(factory.fuzzy.FuzzyInteger(0, 3).fuzz())
        for gap in gaps:
            gap_dict = gap.toJsonDict()
            del gap_dict["l"]  # removes length from gaps
            gap_dict["e"] = gap_dict["s"] + factory.fuzzy.FuzzyInteger(
                gap_length_threshold, max(gap_length_threshold + 1, 50)).fuzz()
            gap_dicts.append(gap_dict)
    return gap_dicts


class GelCoverageMocker(GelCoverageRunner):

    def __init__(self, config):

        # mock everything that needs mocking
        GelCoverageMocker.prepare_mockers()
        super(GelCoverageMocker, self).__init__(config)

    @staticmethod
    def prepare_mockers():
        # Registers the mocker for the region statistics
        GenericFactoryAvro.register_factory(
            protocols.coverage_0_1_0.RegionStatistics, RegionStatisticsFactory, version="6.0.0"
        )

        # Registers the mocker for the region statistics
        GenericFactoryAvro.register_factory(
            protocols.coverage_0_1_0.CoverageGap, CoverageGapFactory, version="6.0.0"
        )

    @patch.object(gelcoverage.tools.bigwig_reader.BigWigReader, 'read_bigwig_coverages', lambda x, y, z, w: None)
    @patch('gelcoverage.stats.coverage_stats.compute_exon_level_statistics', mocked_compute_statistics)
    @patch('gelcoverage.stats.coverage_stats.compute_transcript_level_statistics', mocked_compute_statistics)
    @patch('gelcoverage.stats.coverage_stats.compute_coding_region_statistics', mocked_compute_coding_region_statistics)
    @patch('gelcoverage.stats.coverage_stats.find_gaps', mocked_find_gaps)
    def run(self):
        return super(GelCoverageMocker, self).run()
