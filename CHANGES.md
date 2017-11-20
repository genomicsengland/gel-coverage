version 1.4.0 (20 November 2017)
----------------------------

* Major changes
    - Support for assembly when querying PanelApp

version 1.3.1 (31 August 2017)
----------------------------

* Minor changes
    - Avoid retrieving gene list for whole genome analysis

version 1.3.0 (24 August 2017)
----------------------------

* Major changes
    - The coverage module now accepts a bed file defining the coding regions (--coding-regions) to avoid connection to CellBase.

version 1.2.3 (24 May 2017)
----------------------------

* Minor changes
    - The binary backoff policy has been added to the cellbase client initialisation (https://jira.extge.co.uk/browse/BERTHA-239)

version 1.2.2 (3 April 2017)
----------------------------

* Minor changes
    - Implemented a truncated binary backoff retry policy for CellBase and PanelApp connections.
    - Fixed a bug with some percentages truncated to only 1 decimal place.
    - Logs are now more controlled when missing regions are found. Also, log level is configurable when using the module as a Python library, default log level is INFO.
    - PEP8 fixes

version 1.2.1 (23 January 2017 a bit later...)
----------------------------

* Minor changes
    - Upgrading required numpy version to 1.10.4

version 1.2.0 (23 January 2017)
----------------------------

* Major changes:
    - Coverage standard deviation is added to set of coverage metrics
    - Coverage statistics aggregated for autosomes
    - Some fields are renamed and every field name is not stored as a constant in the code, so renaming is not so painful.

* Minor Changes:
    - Fixed bug when all transcripts in a panel have been filtered out
    - Longs list of genes (>99) are never printed to logs


version 1.1.2 (19 December 2016)
----------------------------

* Minor Changes:
    - Using setuptools find_package() to create a proper python package



version 1.1.1 (19 December 2016)
----------------------------

* Minor Changes:
    - Fixed bug when iterating whole genome in chunks
    - Fixed bug when trying to create an unexisting BED file for whole genome analysis
    - Add flag to avoid running analysis on the coding region
    - Added execution times to documentation
    - Ported shell script for biohpc


version 1.1.0 (14 December 2016)
----------------------------

* Major Changes:
    - Output in machine readable format, JSON.
    - All required (so far) coverage statistics merged in the same module.
    - Reading panels from PanelApp webservices.
    - Reference from CellBase is now read using the Python client, pycellbase.
    - Added percentages of bases over certain coverage thresholds and a measure of uneveness at chromosome and whole genome levels.
