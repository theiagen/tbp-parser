# tbp-parser

!!! warning "Not for Diagnostic Use"
    **CAUTION**: The information produced by this program should **not** be used for clinical reporting unless and until extensive validation has occured in ==_your_== laboratory on a stable version. Otherwise, the outputs of tbp-parser are for research use only.

!!! dna "TBProfiler Compatibility"
    **Versions of TBProfiler prior to v6.0.0 are not compatible with v2+ of tbp-parser.** Please ensure that you are using the correct version of tbp-parser for your version of TBProfiler.

## Overview

`tbp-parser` is a tool developed in partnership with the California Department of Health (CDPH) to parse the output of [Jody Phelan’s TBProfiler tool](https://github.com/jodyphelan/TBProfiler) into four additional files:

1. A *Laboratorian* report, which contains information about each mutation detected and its associated drug resistance profile in a CSV file.
2. A *LIMS* report, formatted specifically for CDPH’s STAR LIMS, which summarizes the highest severity mutations for each antimicrobial drug and the relevant mutations.
3. A *Looker* report, which condenses the information contained in the Laboratorian report into a format suitable for generating a dashboard in Google’s Looker Studio.
4. A *coverage* report, which contains the percent coverage of each gene relative to the H37Rv reference genome in addition to any warnings, such as any deletions identified in the gene that might have contributed to a reduced percent coverage

Please reach out to us at <support@theiagen.com> if you would like any custom file formats and/or changes to these output files that suit your individual needs.
