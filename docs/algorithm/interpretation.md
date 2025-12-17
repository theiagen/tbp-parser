---
title: "The Interpretation Document"
---

Resistance calls are made in either one of two ways. The first is using the WHO annotation, which is output directly from the TBProfiler. The WHO has a catalogue of mutations and how they may confer antimicrobial resistance. If this annotation is present, it will *always* be used.

In the case where the WHO annotation is missing, either due to novel mutations or mutations with unclear significance in the literature, `tbp-parser` will apply *expert rules*. These expert rules are additional conditions used to decide if a mutation is considered to confer resistance or not. **These expert rules come from the CDC** and can be found documented in the `tbp-parser` [GitHub repository](https://github.com/theiagen/tbp-parser) inside the interpretation logic PDFs.

When an expert rule is applied, the `rationale` field of the laboratorian report will indicate which expert rule was used (the number prefacing the rule directly correlates to the appropriate section in the interpretation logic PDF) and also indicate that there was no WHO annotation.

Interpretation documents for v1.2.2 and v1.4.4.8 are available in the [root directory](https://www.github.com/theiagen/tbp-parser) of the `tbp-parser` repository. Versions that correspond to different releases are available in the [interpretation_docs](https://github.com/theiagen/tbp-parser/tree/main/interpretation_docs) directory on GitHub.
