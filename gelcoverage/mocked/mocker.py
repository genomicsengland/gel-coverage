import factory.fuzzy
from mock import patch, Mock
from gelcoverage.runner import GelCoverageRunner
import gelcoverage.tools.bigwig_reader
import gelcoverage.constants as constants

import protocols.coverage_0_1_0
from protocols.util.factories.avro_factory import GenericFactoryAvro, FactoryAvro


class RegionStatisticsFactory(FactoryAvro):
    def __init__(self, *args, **kwargs):
        super(RegionStatisticsFactory, self).__init__(*args, **kwargs)

    class Meta:
        model = protocols.coverage_0_1_0.RegionStatistics

    _version = '6.0.0'

    avg = factory.fuzzy.FuzzyFloat(0, 100)
    med = factory.fuzzy.FuzzyFloat(0, 100)
    sd = factory.fuzzy.FuzzyFloat(0, 100)
    gc = factory.fuzzy.FuzzyFloat(0, 100)
    pct75 = factory.fuzzy.FuzzyFloat(0, 100)
    pct25 = factory.fuzzy.FuzzyFloat(0, 100)
    bases = factory.fuzzy.FuzzyInteger(0, 100)
    gte50x = factory.fuzzy.FuzzyFloat(0, 1)
    gte30x = factory.fuzzy.FuzzyFloat(0, 1)
    gte15x = factory.fuzzy.FuzzyFloat(0, 1)
    lt15x = factory.fuzzy.FuzzyFloat(0, 1)
    rmsd = factory.fuzzy.FuzzyFloat(0, 100)
    bases_gte_50x = factory.fuzzy.FuzzyInteger(50, 200)
    bases_gte_30x = factory.fuzzy.FuzzyFloat(50, 200)
    bases_gte_15x = factory.fuzzy.FuzzyFloat(50, 200)
    bases_lt_15x = factory.fuzzy.FuzzyFloat(50, 200)


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
    gap_factory = GenericFactoryAvro.get_factory_avro(
        protocols.coverage_0_1_0.CoverageGap,
        version='6.0.0'
    )
    gaps = gap_factory.create_batch(3)
    gap_dicts = []
    for gap in gaps:
        gap_dict = gap.toJsonDict()
        del gap_dict["l"]  # removes length from gaps
        gap_dict["e"] = gap_dict["s"] + factory.fuzzy.FuzzyInteger(5, 50).fuzz()
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

    # FIXME: this mock does not work :S
    #@patch.object(gelcoverage.tools.bigwig_reader.BigWigReader, '__init_', lambda x: None)
    @patch.object(gelcoverage.tools.bigwig_reader.BigWigReader, 'read_bigwig_coverages', lambda x, y, z, w: None)
    @patch('gelcoverage.stats.coverage_stats.compute_exon_level_statistics', mocked_compute_statistics)
    @patch('gelcoverage.stats.coverage_stats.compute_transcript_level_statistics', mocked_compute_statistics)
    @patch('gelcoverage.stats.coverage_stats.compute_coding_region_statistics', mocked_compute_coding_region_statistics)
    @patch('gelcoverage.stats.coverage_stats.find_gaps', mocked_find_gaps)
    def run(self):
        return super(GelCoverageMocker, self).run()
