# Changelog

## v1.1.0 — 2026-11-18

- Published the revised manuscript's reserve bounds, two wind-direction and recovery experiments, sampled finite-horizon predictor and QP comparison, and event and uncertainty diagnostics in `Extension_Data_and_Code.zip`.
- Added the revised figure and table crosswalk and instructions for checking the archived results.
- Provided the original `Computational_Materials.zip` with the archive layout required to replay the four accepted S11.2 interval trees.

## v1.0.0 — 2026-11-18

Initial repository packaging of the computational materials for *Air-Path Reference Governance for Fixed-Wing Tracking under Wind*.

- Built the release from the original 1,491-file computational archive, preserving the scientific files in `source/`, archived data and 594 frozen protocol dependencies byte for byte.
- Retained the original text documentation and saved the original archive manifest and MATLAB release wrapper under `validation/original_v6/`.
- Added a repository overview, citation metadata, a reproducibility guide and the `scripts/verify_release.py` entry for the existing Python check sequence.
- Updated `RELEASE_MANIFEST.json` to cover the current release files, including added documentation, tooling and licensing files. Added a guarded manifest-refresh entry for intentional release metadata updates.
- Fixed `matlab/run_representative_direct.m` to resolve a caller-relative output directory to an absolute path before changing the working directory. This changes release orchestration only; the copied scientific functions and parameters are unchanged.
- Preserved the 43-entry map between manuscript and supplementary figures or tables and their supporting records.
- Included the recorded successful reconstruction check comprising 1,101 numerical comparisons.
- Included the recorded MATLAB R2018a reintegration of the original C1 Direct case: 35 s, maximum step 0.005 s, 701 saved points and 21,030 compared numeric values.

Repository packaging does not introduce a new research case, change the scientific models or parameters, or establish an end-to-end rerun of all historical workflows. Manuscript publication status and licensing terms are not inferred from this version number; consult the citation metadata and licensing files.

- Fixed the release wrapper's MATLAB R2018a cleanup callback to restore the caller's path and working directory after the run. The frozen scientific source remains unchanged.

- Applied the author-confirmed MIT license to original code and CC BY 4.0 to original data and standalone documentation; retained all third-party terms. See LICENSING.md.
