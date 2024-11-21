# Versioning and Releases

## Validated Versions

The California Department of Public Health has clinically validated the following versions:

- **v1.2.2 for WGS**, and
- **v1.4.4.8 for tNGS**

!!! warning "Validate Before Use"
    **CAUTION**: The information produced by this program should **not** be used for clinical reporting unless and until extensive validation has occured in ==_**your**_== laboratory on a stable version. Otherwise, the outputs of tbp-parser are for research use only.

## Interpretation Documents

Interpretation documents for v1.2.2 and v1.4.4.8 are available in the [root directory](https://www.github.com/theiagen/tbp-parser) of the `tbp-parser` repository.

Interpretation documents for other versions are available in the [interpretation_docs](https://github.com/theiagen/tbp-parser/tree/main/interpretation_docs) directory on GitHub.

## When Using TheiaProk

!!! dna "FUTURE DEPRECATION NOTICE"
    ==**At the time of the PHB v2.3.0 release:**==

    - **all** branches on Terra that have been mentioned in this documentation will be deleted. Please use the v2.3.0 version of TheiaProk moving forward.
    - the `main` branch of tbp-parser will host v2.1.0 and above; earlier versions of tbp-parser will no longer be supported

If you are running tbp-parser as part of the TheiaProk pipeline(s) with Terra, the following branches are recommended:

- To run v1.2.2 on Terra, please use the **smw-tb-2024-01-16-dev** branch.
- To run v1.4.4.8+ and v1.6.x+, please use the **smw-tb-2024-05-03-dev** branch.
- To run v1.5.x+ and v2.x+, please use the **smw-tb-2024-05-03-*who2*-dev** branch.
- To run v2.1.0+, please use the **smw-tbprofiler-updates-dev** branch _until the time of the v2.3.0 PHB release_.

## Version Differences

For more information on the differences between versions, you can see the [Brief Description of Versions](brief.md) or the [Exhaustive List of Versions](exhaustive.md).
