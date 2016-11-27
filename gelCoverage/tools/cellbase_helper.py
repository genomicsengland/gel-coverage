import pybedtools
from pycellbase.cbclient import CellBaseClient
from pycellbase.cbconfig import ConfigClient

import gelCoverage.stats.sequence_stats as sequence_stats


class CellbaseHelper:

    def __init__(self, species, version, assembly, host, filter_flags, filter_biotypes):
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
        config = ConfigClient(json_config)
        self.cellbase_client = CellBaseClient(config)
        self.cellbase_gene_client = self.cellbase_client.get_gene_client()

    def is_any_flag_included(self, flags):
        """
        Returns a boolean indicating if the list of input flags contain anyone not to be filtered out
        :param flags: The input list of flags
        :return: boolean
        """
        return len(set(self.filter_flags).intersection(set(flags))) > 0

    def get_all_flags_for_gene(self, transcripts):
        """
        For a list of transcripts, each with a set of flags gets a flattened list
        :param transcripts: the transcripts data structure
        :return: flattened list of transcripts
        """
        return sum([y["annotationFlags"] if "annotationFlags" in y else [] for y in transcripts], [])

    def get_all_gene_names(self, filter=True):
        """
        Gets all existing HGNC gene names
        :param filter: flag indicating if filtering should be applied
        :return: list of HGNC gene names
        """
        cellbase_result = self.cellbase_gene_client.search(
            None,
            assembly=self.assembly,
            include=",".join(["name", "transcripts.annotationFlags"]),
            **{"transcripts.biotype": ",".join(self.filter_biotypes) if filter else ""}
        )
        gene_list = [x["name"] for x in cellbase_result[0]["result"]
                     if self.is_any_flag_included(self.get_all_flags_for_gene(x["transcripts"]))]

        return gene_list

    def get_gene_info(self, gene_list, filter = True):
        """
        For a list of HGNC gene names queries CellBase for all information about transcripts and
        their corresponding exons, including the exonic sequence.
        :param gene_list: the list of HGNC gene names
        :param filter: flag indicating whether to filter by biotypes
        :return: the data structure returned by CellBase
        """
        cellbase_genes = self.cellbase_gene_client.search(
            None,
            name=",".join(gene_list),
            assembly=self.assembly,
            include=",".join(["name", "chromosome", "transcripts.exons.start",
                              "transcripts.exons.exonNumber",
                              "transcripts.id,transcripts.strand",
                              "transcripts.exons.end", "transcripts.exons.sequence",
                              "exonNumber", "transcripts.annotationFlags"]),
            **{"transcripts.biotype": ",".join(self.filter_biotypes) if filter else ""}
        )
        # TODO: check for errors and empty results
        # TODO: unit test
        return cellbase_genes

    def make_exons_bed(self, gene_list, filter=True):
        """
        Gets all exons from cellbase and makes a bed - also calculates gc content, returns a valid bed with gc in the
        score column
        :param gene_list The list of genes to be analysed
        :param filter: flag indicating if filtering should be applied
        :return: pybedtools object
        """
        if gene_list is None or len(gene_list) == 0:
            raise SystemError("Input gene list is not correct")
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
            cellbase_genes = self.get_gene_info(current_list)
            # Parse results into BED fields
            for gene in cellbase_genes[0]["result"]:
                gene_name = gene["name"]
                for transcript in gene["transcripts"]:
                    filtered_out = "annotationFlags" not in transcript or \
                                   not self.is_any_flag_included(transcript["annotationFlags"])
                    if filter and filtered_out:
                        # We ignore transcripts not flagged as any of a set of flags in the config file
                        continue
                    txid = transcript["id"]
                    strand = transcript["strand"]
                    for exon in transcript["exons"]:
                        gc = sequence_stats.compute_gc_content(exon["sequence"]) if "sequence" in exon else None
                        exon_number = exon["exonNumber"]
                        row_id = gene_name + "|" + txid + "|exon" + str(exon_number)
                        all_exons.append(
                            (gene["chromosome"], exon["start"], exon["end"], row_id, str(gc), strand))
        # Build BED file
        bed = pybedtools.BedTool(all_exons)
        return bed

    def get_gene_list_from_transcripts (self, transcripts_list):
        # TODO:
        # gene_list = []
        # return gene_list
        raise NotImplemented
