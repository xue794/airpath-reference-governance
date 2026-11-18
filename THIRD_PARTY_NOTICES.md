# Third-party notices and model attribution

This document covers third-party material identified in the archived computational sources. Any repository-wide license for study-authored material does not replace the terms of the components identified below.

## Bundled software: mpmath 1.3.0

The following three frozen source trees contain mpmath 1.3.0:

- `source/contribution_probe_20261106/feedback/vendor/mpmath/`
- `source/revision_wp3_20261107/vendor/mpmath/`
- `source/revision_wp4_20261107/analysis/vendor/mpmath/`

Copyright (c) 2005-2021 Fredrik Johansson and mpmath contributors. mpmath is distributed under its three-clause BSD license. The complete, unmodified copyright notice, conditions and disclaimer are in [LICENSES/mpmath-1.3.0-LICENSE.txt](LICENSES/mpmath-1.3.0-LICENSE.txt). The identical license is also supplied in each adjacent `vendor/mpmath-1.3.0.dist-info/LICENSE` directory. Existing per-file author and copyright notices remain unchanged.

Each tree contains 87 Python source files. All 261 included file copies were compared byte-for-byte with the corresponding files in the official 1.3.0 source distribution; all matched. The separate frozen copies are retained because historical analyses refer to their own relative directories. No mpmath implementation was changed for this release.

Upstream sources:

- [mpmath 1.3.0 release on PyPI](https://pypi.org/project/mpmath/1.3.0/)
- [Upstream 1.3.0 license](https://github.com/mpmath/mpmath/blob/1.3.0/LICENSE)
- [Upstream project and contributor credits](https://github.com/mpmath/mpmath/tree/1.3.0)

The source archive `mpmath-1.3.0.tar.gz` has SHA-256 `7a28eb2a9774d00c7bc92411c19a89209d5da7c4c9a9e227be8330a23a25b91f`, as published by PyPI. The license was taken from that archive and independently matched against the upstream GitHub 1.3.0 tag. Its SHA-256 is `c26cae81da4508e5e249985777a33821f183223ebb74d7f8cfbf90fe7eef2fb7`.

## Aircraft-model source

The study's `t2_gna_params.m` and `t2_gna_coeff.m` implementations identify their model source as:

Grauer, J. A.; Morelli, E. A. *A Generic Nonlinear Aerodynamic Model for Aircraft*. AIAA Atmospheric Flight Mechanics Conference, 2014, AIAA Paper 2014-0542. [DOI: 10.2514/6.2014-0542](https://doi.org/10.2514/6.2014-0542). [NASA NTRS record 20140011902](https://ntrs.nasa.gov/citations/20140011902).

The existing source comments identify Table 1 for rigid-body and geometry values, Tables 3–4 for the 45 aerodynamic coefficients, Table 5 for the reference trim case, and Table 6 for modal reference values. These comments are preserved in the `core` snapshot and its historical WP3/WP4 copies.

The NASA NTRS record marks the cited paper as public and states “Work of the US Gov. Public Use Permitted.” This attribution concerns the published model and parameters. The MATLAB implementation and the governor study supplied here are study materials; they are not presented as an official NASA software distribution or as NASA-endorsed results. No upstream paper PDF or NASA software package is bundled.

## External software dependencies

MATLAB and its toolboxes are not included. The representative MATLAB entry requires a separately installed MATLAB environment; MATLAB R2018a was used for the recorded check.

Some historical analysis or document-generation scripts import NumPy, SciPy, pandas, Matplotlib, python-docx and lxml. These packages are external environment dependencies; their implementations are not bundled in this repository. Their own distributions provide their license terms. Historical report scripts can also refer to environment-specific document-rendering helpers. The portable saved-output reconstruction entries under `python/` use the Python standard library only.

## Audit scope

The release audit inspected the source inventory, vendor directories, existing copyright/license/source comments, and the nested S11 data archive. mpmath was the vendored software library identified. The S11 archive contains saved data and text records, with no embedded Python or MATLAB implementation. Model attribution above follows the existing source comments and the official NASA record.

License files and this notice were added without changing the archived scientific source files, their paths, or `SOURCE_MANIFEST.json`. The current `RELEASE_MANIFEST.json` also covers the added notices and release tooling; the original v6 manifest is preserved under `validation/original_v6/`. This notice does not enlarge the numerical validation scope stated in the README.
