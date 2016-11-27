import ConfigParser
import argparse
import json
import os

from gelCoverage.gel_coverage import GelCoverageEngine


def main():

    config_file = '../resources/exon_coverage_summary.config'

    parser = argparse.ArgumentParser(description = 'Coverage summary for exons')
    parser.add_argument('--bw', metavar='bw', help = 'This is the bigwig file')
    parser.add_argument('--panel', metavar = 'panel',
                        help = 'The panel name or identifier in PanelApp /'
                             '(see https://bioinfo.extge.co.uk/crowdsourcing/WebServices/list_panels)',
                        default = None)
    parser.add_argument('--panel-version', metavar = 'panel_version',
                        help='The panel version',
                        default = None)
    parser.add_argument('--gene_list', metavar = 'gene_list',
                        default = None,
                        help = 'Comma separated list of genes (HGNC gene symbols) to analyse. Will be masked by a panel')
    #parser.add_argument('--transcripts', metavar='transcripts',
    #                    help='Comma separated list of transcripts to analyse. Will be masked by a panel or a list of genes')
    #parser.add_argument('--canonical_transcript', metavar='canonical_transcript',
    #                    help='collapse transcripts to one cannonical transcript - only used when using panel or genes as input',
    #                    default=0)
    #parser.add_argument('--output', metavar='output',
    #                    help='Output format',
    #                    choices = ["json", "tabular"],
    #                    default = "json")
    parser.add_argument('--coverage-threshold', metavar='coverage_threshold',
                        help='The coverage threshold used to compute continuous gaps with low coverage (0 = disabled)',
                        default = 15)
    #TODO: add parameter to add flanking regions to genes
    #TODO: add parameter to return only information about canonical transcripts
    #parser.add_argument('--cnv', metavar='cnv', help='cnv vcf - so that losses can be indicated', default=0)

    args = parser.parse_args()

    # Reads configuration file
    config_parser = ConfigParser.ConfigParser()
    config_filepath = os.path.join(os.path.dirname(os.path.realpath(__file__)), config_file)
    config_parser.readfp(open(config_filepath))

    # Creates a data structure with all config parameters
    config = {
        # Sets parameters from CLI
        "bw" : args.bw,
        "panel" : args.panel,
        "panel_version": args.panel_version,
        "gene_list": args.gene_list,
        "coverage_threshold": args.coverage_threshold,
        # Sets parameters from config file
        "cellbase_species": config_parser.get('cellbase', 'species'),
        "cellbase_version": config_parser.get('cellbase', 'version'),
        "cellbase_assembly": config_parser.get('cellbase', 'assembly'),
        "cellbase_host": config_parser.get('cellbase', 'host'),
        "panelapp_host": config_parser.get('panelapp', 'host'),
        "panelapp_gene_confidence": config_parser.get('panelapp', 'gene_confidence'),
        "transcript_filtering_flags": config_parser.get('transcript_filtering', 'flags'),
        "transcript_filtering_biotypes": config_parser.get('transcript_filtering', 'biotypes')
    }
    # Calls the GEL coverage engine
    gel_coverage_engine = GelCoverageEngine(config)
    output = gel_coverage_engine.run()
    # Prints output to stdout
    # TODO: we may want to write it to a file. Check that
    print json.dumps(output, indent=4)

    # TODO: output results in different formats
    #if args.output == "json":
    #    print json.dumps(output, indent = 4)
    #elif args.output == "tabular":
        # TODO: flatten the dict into a table
    #    pass

if __name__ == '__main__':
    main()