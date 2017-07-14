import pybedtools
from pycellbase.cbclient import CellBaseClient
from pycellbase.cbconfig import ConfigClient
import logging

import gelcoverage.stats.sequence_stats as sequence_stats
import gelcoverage.tools.backoff_retrier as backoff_retrier


chromosome_mapping = {
    "hsapiens": {
        "grch38": {
            "1": "chr1",
            "2": "chr2",
            "3": "chr3",
            "4": "chr4",
            "5": "chr5",
            "6": "chr6",
            "7": "chr7",
            "8": "chr8",
            "9": "chr9",
            "10": "chr10",
            "11": "chr11",
            "12": "chr12",
            "13": "chr13",
            "14": "chr14",
            "15": "chr15",
            "16": "chr16",
            "17": "chr17",
            "18": "chr18",
            "19": "chr19",
            "20": "chr20",
            "21": "chr21",
            "22": "chr22",
            "X": "chrX",
            "Y": "chrY",
            "MT": "chrM"
        },
        "grch37": {
            "1": "chr1",
            "2": "chr2",
            "3": "chr3",
            "4": "chr4",
            "5": "chr5",
            "6": "chr6",
            "7": "chr7",
            "8": "chr8",
            "9": "chr9",
            "10": "chr10",
            "11": "chr11",
            "12": "chr12",
            "13": "chr13",
            "14": "chr14",
            "15": "chr15",
            "16": "chr16",
            "17": "chr17",
            "18": "chr18",
            "19": "chr19",
            "20": "chr20",
            "21": "chr21",
            "22": "chr22",
            "X": "chrX",
            "Y": "chrY",
            "MT": "chrM"
        }
    }
}


class CellbaseHelper:

    def __init__(self, species, version, assembly, host, retries, filter_flags, filter_biotypes):
        """
        Initializes the CellBase helper.
        :param species: The species
        :param version: The API version
        :param assembly: The species reference genome assembly
        :param host: The webservices host
        :param filter_flags: list of transcript flags to consider in analysis
        :param filter_biotypes: list of biotypes to consider in analysis
        """
        # Saves transcript filtering configuration
        self.filter_flags = filter_flags
        self.filter_biotypes = filter_biotypes
        self.assembly = assembly
        self.species = species
        self.retries = retries
        # Builds JSON configuration
        json_config = {
            "species": species,
            "version": version,
            "rest": {
                "hosts": [
                    host
                ]
            }
        }
        # Initializes the CellBase client
        self.__secure_initialise_cellbase_client = backoff_retrier.wrapper(self.__initialise_cellbase_client,
                                                                           self.retries)
        self.cellbase_client = self.__secure_initialise_cellbase_client(json_config)
        self.cellbase_gene_client = self.cellbase_client.get_gene_client()
        # Wraps the CB search call into the truncated binary backoff implementation
        self.cellbase_search = backoff_retrier.wrapper(self.cellbase_gene_client.search, self.retries)
        # Loads chromosome mapping for this specific reference
        self.chromosome_mapping = chromosome_mapping[self.species.lower()][self.assembly.lower()]

    def __initialise_cellbase_client(self, json_config):
        config = ConfigClient(json_config)
        return CellBaseClient(config)

    def __is_any_flag_included(self, flags):
        """
        Returns a boolean indicating if the list of input flags contain anyone not to be filtered out
        :param flags: The input list of flags
        :return: boolean
        """
        return len(set(self.filter_flags).intersection(set(flags))) > 0

    @staticmethod
    def __get_all_flags_for_gene(transcripts):
        """
        For a list of transcripts, each with a set of flags gets a flattened list
        :param transcripts: the transcripts data structure
        :return: flattened list of transcripts
        """
        return sum([y["annotationFlags"] if "annotationFlags" in y else [] for y in transcripts], [])

    def get_all_gene_names(self, _filter=True):
        """
        Gets all existing HGNC gene names
        :param _filter: flag indicating if filtering should be applied
        :return: list of HGNC gene names
        """
        logging.debug("Getting gene list from CellBase...")
        # calls to CB
        cellbase_result = self.cellbase_search(
                    assembly=self.assembly,
                    include=",".join(["name", "transcripts.annotationFlags"]),
                    **{"transcripts.biotype": ",".join(self.filter_biotypes) if _filter else ""}
                )
        gene_list = [x["name"] for x in cellbase_result[0]["result"]
                     if self.__is_any_flag_included(CellbaseHelper.__get_all_flags_for_gene(x["transcripts"]))]
        logging.debug("Gene list obtained from CellBase of %s genes" % str(len(gene_list)))
        return gene_list

    def __get_gene_info(self, gene_list, _filter=True):
        """
        For a list of HGNC gene names queries CellBase for all information about transcripts and
        their corresponding exons, including the exonic sequence.
        :param gene_list: the list of HGNC gene names
        :param _filter: flag indicating whether to filter by biotypes
        :return: the data structure returned by CellBase
        """
        # calls to CB
        cellbase_genes = self.cellbase_search(
            name=gene_list,
            assembly=self.assembly,
            include=["name", "chromosome", "transcripts.exons.start", "transcripts.exons.exonNumber",
                     "transcripts.id,transcripts.strand", "transcripts.exons.end", "transcripts.exons.sequence",
                     "exonNumber", "transcripts.annotationFlags"],
            **{"transcripts.biotype": self.filter_biotypes if _filter else []}
        )
        # TODO: check for errors and empty results
        # TODO: unit test
        return cellbase_genes

    def make_exons_bed(self, gene_list, _filter=True, has_chr_prefix=False):
        """
        Gets all exons from cellbase and makes a bed - also calculates gc content, returns a valid bed with gc in the
        score column
        :param gene_list The list of genes to be analysed
        :param _filter: flag indicating if filtering should be applied
        :param has_chr_prefix: flag indicating if the chromosomes have the chr prefix or not
        :return: pybedtools object
        """
        logging.info("Building gene annotations bed file from CellBase...")
        if gene_list is None or len(gene_list) == 0:
            raise SystemError("Input gene list is empty")
        number_genes = len(gene_list)

        # TODO: Verify that all genes in the list are present in the reference

        # Iterates through genes in 1000 genes batches
        gene_batch = 1000
        gene_count = 0
        all_exons = list()
        while gene_count < number_genes:
            limit = min(len(gene_list) - gene_count, gene_batch)
            current_list = gene_list[gene_count: limit]
            gene_count += limit
            # Query CellBase for a subset of genes
            cellbase_genes = self.__get_gene_info(current_list)
            # Parse results into BED fields
            for gene in cellbase_genes[0]["result"]:
                gene_name = gene["name"]
                for transcript in gene["transcripts"]:
                    filtered_out = "annotationFlags" not in transcript or \
                                   not self.__is_any_flag_included(transcript["annotationFlags"])
                    if _filter and filtered_out:
                        # We ignore transcripts not flagged as any of a set of flags in the config file
                        continue
                    txid = transcript["id"]
                    strand = transcript["strand"]
                    chromosome = gene["chromosome"]
                    if has_chr_prefix:
                        if chromosome not in self.chromosome_mapping:
                            continue  # skips genes not in canonical chromosomes when requires transformation
                        chromosome = self.chromosome_mapping[chromosome]  # transforms chromosome
                    for exon in transcript["exons"]:
                        gc = sequence_stats.compute_gc_content(exon["sequence"]) if "sequence" in exon else None
                        exon_number = exon["exonNumber"]
                        row_id = gene_name + "|" + txid + "|exon" + str(exon_number)
                        all_exons.append(
                            (chromosome, exon["start"], exon["end"], row_id, str(gc), strand))
        # Build BED file
        bed = pybedtools.BedTool(all_exons)
        logging.info("Gene annotations bed file built!")
        return bed

    def get_gene_list_from_transcripts(self, transcripts_list):
        # TODO: implement it
        # gene_list = []
        # return gene_list
        raise NotImplemented
