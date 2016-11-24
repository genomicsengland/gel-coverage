import pybedtools
from pycellbase.cbconfig import *
from pycellbase.cbclient import *
import re

class CellbaseHelper():

    def __init__(self, species, version, assembly, host, filter_basic_flag, filter_biotypes):
        """
        Initializes the CellBase helper.
        :param species: The species
        :param version: The API version
        :param assembly: The species reference genome assembly
        :param host: The webservices host
        :param filter_basic_flag: boolean indicating wether to filter out transcripts without the basic flag
        :param filter_biotypes: list of biotypes to be filtered out
        """
        # Saves transcript filtering configuration
        self.filter_basic_flag = filter_basic_flag
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
        conf = ConfigClient("../../resources/cellbase_config.json")
        self.cellbase_client = CellBaseClient(conf)
        self.cellbase_gene_client = self.cellbase_client.get_gene_client()

    def gc_content(self, sequence):
        """
        Calculates gc content of a DNA sequence

        :param sequence: a string of ATGC's
        :return: the gc fraction of the region
        """
        gcCount = 0
        totalBaseCount = 0
        gcCount += len(re.findall("[GC]", sequence))
        totalBaseCount += len(re.findall("[GCTA]", sequence))
        gcFraction = float(gcCount) / totalBaseCount
        return (gcFraction)

    def get_all_genes (self):
        """
        Gets all existing genes ENSEMBL ids
        :param assembly: the assembly
        :return: list of ENSEMBL identifiers
        """
        # TODO: filter by biotype and basic flag
        cellbase_result = self.cellbase_client.get("feature", "gene", "list", assembly = self.assembly)
        gene_list = [x["id"] for x in cellbase_result [0]["result"]]
        return gene_list

    def get_gene_list_from_transcripts (self, transcripts_list):
        # TODO:
        #gene_list = []
        #return gene_list
        raise NotImplemented

    def make_exons_bed(self, gene_list):
        """
        Gets all exons from cellbase and makes a bed - also calculates gc content, returns a valid bed with gc in the
        score column
        :param gene_list The list of genes to be analysed
        :return: pybedtools object
        """
        # TODO: filter by biotype and basic flag
        if gene_list is None or len(gene_list) == 0:
            raise SystemError("Input gene list is not correct")
        number_genes = len(gene_list)

        # TODO: Verify that all genes in the list are present in the reference

        # Iterates through genes in 1000 genes batches
        GENE_BATCH = 1000
        gene_count = 0
        all_exons = list()
        while gene_count < number_genes:
            limit = min( len(gene_list) - gene_count, GENE_BATCH)
            current_list = gene_list[gene_count : limit]
            gene_count += limit
            # Query CellBase for a subset of genes
            # TODO: filter by biotype and basic flag
            cellbase_exons = self.cellbase_gene_client.search(
                None,
                name = ",".join(current_list),
                assembly = self.assembly,
                include = ",".join(["name","chromosome","transcripts.exons.start",
                                    "transcripts.exons.exonNumber", "transcripts.id,transcripts.strand",
                                    "transcripts.exons.end", "transcripts.exons.sequence", "exonNumber"]))
            # TODO: check for errors and empty results
            # Parse results into BED fields
            for gene in cellbase_exons[0]["result"]:
                gene_name = gene["name"]
                for transcript in gene["transcripts"]:
                    txid = transcript["id"]
                    strand = transcript["strand"]
                    for exon in transcript["exons"]:
                        gc = self.gc_content(exon["sequence"]) if "sequence" in exon else None
                        exon_number = exon["exonNumber"]
                        row_id = gene_name + "|" + txid + "|exon" + str(exon_number)
                        all_exons.append(
                            (gene["chromosome"], exon["start"], exon["end"], row_id, str(gc), strand))
        # Build BED file
        bed = pybedtools.BedTool(all_exons)
        return bed