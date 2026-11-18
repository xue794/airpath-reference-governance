# Air-Path Reference Governance for Fixed-Wing Tracking under Wind

Repository: [https://github.com/xue794/airpath-reference-governance](https://github.com/xue794/airpath-reference-governance)

Computational materials accompanying the manuscript by **Tianming Ren, Ou Li, Guo Jiang and Xue Zhong**. The **v1.1.0** release adds the experiments and reserve calculations in the revised manuscript.

The study constructs reference constraints for coupled lift–thrust inversion, quantifies how downstream feedback consumes the reference force reserve, and evaluates reference selection and prescribed deceleration schedules. This repository connects the reported quantities to archived data, frozen scientific source and executable checks.

## Revised manuscript materials

Download `Extension_Data_and_Code.zip` from the [v1.1.0 release](https://github.com/xue794/airpath-reference-governance/releases/tag/v1.1.0), extract it, and run `python check_extension.py` from the extracted directory. The archive contains the reserve design, wind and recovery simulations, sampled predictor and QP comparisons, and event and uncertainty diagnostics. The [revised manuscript evidence guide](docs/REVISED_STUDY.md) maps the new figure and table numbers to the supporting records and gives the commands for additional verification.

The [v1.0.0 release](https://github.com/xue794/airpath-reference-governance/releases/tag/v1.0.0) identifies the original computational materials. The sections below describe how to check those results.

## Start with the saved-output checks

Use **Python 3.10 or later**. The verification entry uses only the Python standard library. From the repository root, run:

```sh
python scripts/verify_release.py
```

This checks the release manifest, reconstructs selected summaries from saved outputs, and compares them with archived results. It writes results under `generated/` and does not launch MATLAB. The recorded preparation check passed **1,101 numerical comparisons**; see [the validation record](validation/SAVED_OUTPUT_RECONSTRUCTION.json).

The reconstructed quantities include all 500 paired uncertainty outcomes, common-completion RMS differences, ordinary-case tracking and thrust diagnostics, endpoint-deletion counts, computation-time summaries, and C1/C2 comparisons between the static backup, moving backup and direct center following.

For the commands, interpretation and optional MATLAB run, read [Reproducibility and verification scope](docs/REPRODUCIBILITY.md).

## Find the evidence

| Location | Purpose |
| --- | --- |
| [MAIN_AND_SUPPLEMENT_MAP.csv](MAIN_AND_SUPPLEMENT_MAP.csv) | 43 entries linking the v1.0.0 figure and table numbers to their supporting records and verification scope. |
| [DATA_DICTIONARY.txt](DATA_DICTIONARY.txt) | Method names, units, time windows, missing values and comparison conventions. |
| [data/](data/) | Archived outcomes, saved trajectories, event records, numerical certificates and the focused S11 archive. |
| [source/](source/) | Original scientific modules, frozen inputs and source snapshots in their retained relative layout. |
| [python/](python/) | Saved-output reconstruction, comparison and integrity checks. |
| [matlab/](matlab/) | Isolated entry for the representative C1 Direct reintegration. |
| [validation/](validation/) | Records of the checks performed during package preparation. |
| [SOURCE_WORKFLOWS.txt](SOURCE_WORKFLOWS.txt) | Original study entry points and their execution requirements. |

This release is built from the original **1,491-file computational archive**. Its scientific files in `source/`, archived data and **594 frozen protocol dependencies** preserve their original bytes. Repository documentation and tooling have been added, and the MATLAB release wrapper has been fixed to resolve a relative output path before changing directories. The original wrapper and archive manifest are retained in [validation/original_v6/](validation/original_v6/). The original [README.txt](README.txt) also remains unchanged as an archive record; the current execution instructions are in this page and the reproducibility guide.

[SOURCE_MANIFEST.json](SOURCE_MANIFEST.json) records provenance and frozen source hashes. [RELEASE_MANIFEST.json](RELEASE_MANIFEST.json) covers the current v1.0.0 files, including added documentation, tooling and licensing files; it excludes itself, generated outputs and development caches. The original archive counts are not the total file count of this repository. Maintainers updating release metadata should follow [the release checklist](docs/RELEASE_CHECKLIST.md), including its manifest refresh and verification steps.

## What has been reproduced

Two activities have been checked:

1. **Reconstruction from saved outputs.** Selected counts and numerical summaries pass the 1,101 comparisons described above. These checks do not regenerate the underlying trajectories or interval trees.
2. **One MATLAB reintegration.** The original C1 Direct case was rerun with MATLAB R2018a, unchanged scientific source and parameters, a maximum integration step of 0.005 s and a 35 s horizon. Its 701 saved points contain 21,030 compared numeric values, which matched the archive exactly on the checked workstation. The comparison tolerance is `1e-8` in each column's stated units. The metrics also matched; see [MATLAB_REINTEGRATION.json](validation/MATLAB_REINTEGRATION.json).

The other historical workflows are provided for provenance and further execution with their stated dependencies. Some launchers contain workstation-specific paths. This is **not a one-command reproduction of the entire study**, and the complete physical histories of all 1,000 uncertainty integrations are not included. Exact agreement for the representative run is not a claim of bitwise identity on other systems.

## Reading the results

- The source label `Proposed` denotes Scalar; `MatchedQP` denotes Matched-QP.
- The 408 common completions are the subset used for paired continuous metrics; all 500 pairs remain in the outcome accounting.
- `RAW_COMMAND_CHART_EXIT` denotes an upstream coordinate-domain exit, not a QP solver failure.
- Reference thrust, raw feedback thrust and modeled actuator thrust are distinct variables.
- The 0.606 m/s scale is a research comparison scale, not an independent mission acceptance threshold.
- The endpoint-protection experiments without a backup row and the regional certificates with a backup have different guarantees. Their results must retain those distinctions.

## Citation and licensing

Author and version metadata are provided in [CITATION.cff](CITATION.cff). The accompanying work is identified here as a manuscript; no publication DOI is assigned by this repository.

Study-authored code is licensed under [MIT](LICENSE). Study-authored data and standalone documentation are licensed under [CC BY 4.0](LICENSES/CC-BY-4.0.txt). See [LICENSING.md](LICENSING.md) for the scope of each license. Third-party components retain their original terms, including the BSD 3-Clause license for mpmath; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). The original provenance records are retained with the scientific files.

## 中文导航

修订稿新增实验和力裕量计算见 [v1.1.0](https://github.com/xue794/airpath-reference-governance/releases/tag/v1.1.0)。解压 `Extension_Data_and_Code.zip` 后运行 `python check_extension.py`。新版图号、数据文件和可选复算命令见[修订稿复现指南](docs/REVISED_STUDY.md)。

这是论文的计算资料与核查入口。先在仓库根目录运行 `python scripts/verify_release.py`，再用 [43项图表索引](MAIN_AND_SUPPLEMENT_MAP.csv) 查找论文结果对应的文件。[复现说明](docs/REPRODUCIBILITY.md) 区分“从存档输出复算指标”和“重新运行MATLAB仿真”：本版已核查1,101项数值比较，并重跑一条C1 Direct案例；不声称全研究已经一键复现。原始源文件、冻结输入和旧说明均保留。
