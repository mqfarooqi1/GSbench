## Update to version 0.2.0

This is a feature update. New functionality:

* `gwas()` — single-marker genome-wide association scan with optional
  principal-component structure correction, returning a `gs_gwas` object with
  `print()` and `plot()` (Manhattan / QQ) methods.
* `Dmatrix()` — dominance genomic relationship matrix (Vitezica et al. 2013),
  usable in `gblup()` as `K` to fit a dominance model.

No user-facing changes to existing functions; the update is backward compatible.

## R CMD check results

Local `R CMD check --as-cran` (Windows 11, R 4.5.2): 0 errors | 0 warnings |
2 notes.

* The only WARNING seen locally is `'qpdf' is needed for checks on size
  reduction of PDFs` — an artifact of qpdf not being installed on the local
  machine; it does not occur on systems with qpdf (e.g. CRAN).
* NOTE "unable to verify current time" is a local network/clock artifact.
* NOTE (incoming feasibility) reports the maintainer and days since the last
  update; the update is justified by the substantial new features above.

## Notes

* Any "possibly misspelled words" are domain acronyms (GBLUP, GWAS, QQ) and
  author surnames (VanRaden, Endelman, Meuwissen, Vitezica); these are listed in
  `inst/WORDLIST`.
* Machine-learning backends remain optional (Suggests) and thread-limited; tests
  using them skip on CRAN so the test run stays single-core.
* Functions write only to `tempdir()` and do not modify the global environment.
  Method references carry DOIs.

## Test environments

* Local: Windows 11, R 4.5.2
