import factory.fuzzy
import numpy.random
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

    avg = 0
    med = 0
    sd = 0
    gc = 0
    pct75 = 0
    pct25 = 0
    bases = factory.fuzzy.FuzzyInteger(400, 2000)
    gte50x = 0
    gte30x = 0
    gte15x = 0
    lt15x = 0
    rmsd = FuzzyFloatWithPrecision(0, 100)
    bases_gte_50x = 0
    bases_gte_30x = 0
    bases_gte_15x = 0
    bases_lt_15x = 0


class CoverageGapFactory(FactoryAvro):
    def __init__(self, *args, **kwargs):
        super(CoverageGapFactory, self).__init__(*args, **kwargs)

    class Meta:
        model = protocols.coverage_0_1_0.CoverageGap

    _version = '6.0.0'

    s = factory.fuzzy.FuzzyInteger(10000, 50000000)
    e = factory.fuzzy.FuzzyInteger(10000, 50000000)


def mocked_compute_statistics(*args, **kwargs):
    region_statistics_factory = GenericFactoryAvro.get_factory_avro(
        protocols.coverage_0_1_0.RegionStatistics,
        version='6.0.0'
    )
    mean = kwargs.get('mean', 50.0)
    sd = kwargs.get('sd', 20.0)
    region_statistics = region_statistics_factory.create()
    size = region_statistics.bases if region_statistics.bases is not None else 1
    randomised_mean = numpy.random.normal(loc=mean, scale=sd)
    coverages = numpy.random.normal(loc=randomised_mean, scale=sd, size=size)
    region_statistics.avg = numpy.mean(coverages)
    region_statistics.med = numpy.median(coverages)
    region_statistics.sd = numpy.std(coverages)
    region_statistics.pct75 = numpy.percentile(coverages, 75)
    region_statistics.pct25 = numpy.percentile(coverages, 25)
    region_statistics.bases_gte_50x = numpy.sum([1 for x in coverages if x >= 50])
    region_statistics.bases_gte_30x = numpy.sum([1 for x in coverages if x >= 30])
    region_statistics.bases_gte_15x = numpy.sum([1 for x in coverages if x >= 15])
    region_statistics.bases_lt_15x = numpy.sum([1 for x in coverages if x < 15])
    region_statistics.gte50x = float(region_statistics.bases_gte_50x) / region_statistics.bases
    region_statistics.gte30x = float(region_statistics.bases_gte_30x) / region_statistics.bases
    region_statistics.gte15x = float(region_statistics.bases_gte_15x) / region_statistics.bases
    region_statistics.gte15x = float(region_statistics.bases_lt_15x) / region_statistics.bases
    region_statistics.gc = float(numpy.random.normal(loc=50, scale=20))/100
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

    def __init__(self, config, mocking_config):

        # mock everything that needs mocking
        GelCoverageMocker.prepare_mockers()
        super(GelCoverageMocker, self).__init__(config)
        self.lower_coverage_genes = mocking_config['lower_coverage_genes']
        self.lower_coverage_mean = mocking_config['lower_coverage_mean']
        self.lower_coverage_sd = mocking_config['lower_coverage_sd']
        self.higher_coverage_genes = mocking_config['higher_coverage_genes']
        self.higher_coverage_mean = mocking_config['higher_coverage_mean']
        self.higher_coverage_sd = mocking_config['higher_coverage_sd']

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

    def set_differential_coverage_genes(self, results):
        if self.lower_coverage_genes is not None or self.higher_coverage_genes is not None:
            for gene in results["results"]["genes"]:
                if gene["name"] in self.lower_coverage_genes:
                    gene["union_tr"]["stats"] = mocked_compute_statistics(
                        mean=float(self.lower_coverage_mean), sd=float(self.lower_coverage_sd))
                if gene["name"] in self.higher_coverage_genes:
                    gene["union_tr"]["stats"] = mocked_compute_statistics(
                        mean=float(self.higher_coverage_mean), sd=float(self.higher_coverage_sd))
        return results

    @patch.object(gelcoverage.tools.bigwig_reader.BigWigReader, 'read_bigwig_coverages', lambda x, y, z, w: None)
    @patch('gelcoverage.stats.coverage_stats.compute_exon_level_statistics', mocked_compute_statistics)
    @patch('gelcoverage.stats.coverage_stats.compute_transcript_level_statistics', mocked_compute_statistics)
    @patch('gelcoverage.stats.coverage_stats.compute_coding_region_statistics', mocked_compute_coding_region_statistics)
    @patch('gelcoverage.stats.coverage_stats.find_gaps', mocked_find_gaps)
    def run(self):
        results, bed = super(GelCoverageMocker, self).run()
        results = self.set_differential_coverage_genes(results)
        return results, bed
