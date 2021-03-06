# Coverage analysis development status (05/12/2016)

## Added features

1. JSON output (the downside is not having an easy way to explore results)
2. 3 supported running modes:
	* Panel from PanelApp (name and version)
	* Custom gene list
	* None of the above => Whole exome analysis
3. Filtering transcripts by:
	* Transcripts flags
	* Biotypes
2. Configuration from command line parameters and from a config file
4. Coverage metrics at different genomic levels
	* Exon level (include gaps of low coverage)
	* Transcript level
	* Union transcript level (exons + full transcript)
	* Panel / gene list level (aggregated from the union transcript level)
	* Whole genome level
5. Exon padding (this is taken into account for the union transcript)
6. Tests on reduced data


## Next steps

1. Run tests on whole genome data
	* Already found a bug in whole exome analysis mode
	* Retrieve appropriate testing data for cancer and rare diseases
	* Implement some verification code (specially for the union transcript)
	* Test GRCh38 (should be supported out-of-the-box)
2. Make use of BED with N regions for whole genome metrics (can we retrieve this from CellBase??? Otherwise add a config parameter)
4. % of bases in a BED region below whateverX
5. Aggregate stats at chromosome level (including autosomes aggregation, only mean)

## Nice to have

1. Transcript list analysis mode
2. Cache queries to CellBase
3. Create the coverage distribution tables and plots created by the old script
4. Median and percentiles at whole genome level (this will require some weighting)
5. Automate test data download using OpenCGA
6. Simple tool to inspect coverage analysis output and even transform into tables ingestable by R