# Revised manuscript evidence guide

Download the two named assets from the [v1.1.0 release](https://github.com/xue794/airpath-reference-governance/releases/tag/v1.1.0). File paths in the table are relative to the extracted extension.

| Revised manuscript item | Computational record |
| --- | --- |
| Figure 2 and S13.5, C1 inversion event and angle sensitivity | Original Figure 9 records in v1.0.0 and `audit/results/fixed_event_angles.csv` |
| Figure 3 and S13.1, reserve sizing and trajectory coverage | `reserve/output/reserve_design.csv`, `reserve/output/reserve_design.pdf`, `reserve/output/trajectory_coverage.csv` and `reserve/output/witnesses/` |
| Figure 9 and S13.2, wind direction and return to trim | `wind/results/wind_recovery.pdf`, `wind/results/all_cases.csv` and `wind/results/attribution.csv` |
| Table 10 and S13.3–S13.4, sampled predictor and QP comparison | `predictive/metrics.csv`, `predictive/failure_check.json` and `predictive/timing_blocks.csv` |
| Table 12, addition of C3 Direct to the original deceleration table | `wind/results/all_cases.csv` and original Table 11 records in v1.0.0 |
| S13.5, 500 paired outcomes and two re-integrated examples | `audit/results/input_groups.csv`, `audit/results/discordant_021.csv` and `audit/results/discordant_026.csv` |

The [original 43-entry evidence index](../MAIN_AND_SUPPLEMENT_MAP.csv) uses v1.0.0 numbering. Original main Figures 1–9 appear in the revised manuscript as Figures 1, 4, S7, 5, 6, 7, S8, 8 and 2, respectively. Original main Tables 10 and 11 appear as Tables 11 and 12.

## Verify the archived results

Extract `Extension_Data_and_Code.zip` and run the standard-library check in the extracted directory.

```sh
python check_extension.py
```

Use the dependencies listed in the archive README for the reserve calculation and the additional scientific checks. Supply the named `Computational_Materials.zip` asset as an extra input to replay all four original S11.2 interval trees.

```sh
python check_extension.py --science --design --legacy-archive "/path/to/Computational_Materials.zip"
```

The extension README lists the entries for re-integrating the MATLAB studies and describes what each Python check verifies. The original repository checks remain available through `python scripts/verify_release.py`.

| Release asset | SHA-256 |
| --- | --- |
| `Extension_Data_and_Code.zip` | `a72ac6b53aec3a0d8f81013b8d87ccdeb44be85d32c6cb9e5cc4d53db731af1f` |
| `Computational_Materials.zip` | `3844fcb36b1f623cacd6cfc87e392277ec030b2ed1f092891abe629a4d6611cd` |
