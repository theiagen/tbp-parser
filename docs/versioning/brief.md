---
title: Brief Description of Versions
---

!!! warning "Validate Before Use"
    **CAUTION**: The information produced by this program should **not** be used for clinical reporting unless and until extensive validation has occured in ==_**your**_== laboratory on a stable version. Otherwise, the outputs of tbp-parser are for research use only.

You may notice there are many releases; `tbp-parser` is in active development and each release is "use at your own risk." We highly recommend upgrading to the latest release as they include important bug fixes. In order to help track the different changes, we have included a brief description of each release:

- **v1.2.x _& below_** - the initial developmental stages of tbp-parser for WGS data
- **v1.3.x** - the addition of tNGS data parsing and includes some updates applicable to WGS parsing
- **v1.4.x** - reworks how QC is performed (changes in order of operations)
    - **v1.4.3+** - changes how tNGS lineage determination is performed
    - **v1.4.4+** - changes how nonsynonymous mutations are interpretted; major interpretation differences between earlier versions
- **v1.6.x** - only considers the genes included in the LIMS report to determine the drug output in the LIMS report
- **v1.5.x+ and v2.0.0** - major changes to code in due to using results from TBProfiler v6.2.0+; compatible with WHO v2
- **v2.1.0** - ==_v1.6.0 and earlier versions are no longer supported_==; v2.1+ changes are included on `main` branch moving foward.
- **v3.0.0+** - a complete refactor of the codebase to improve maintainability and add new features

For a more exhaustive list, please visit [the Exhaustive List of Versions](exhaustive.md).
