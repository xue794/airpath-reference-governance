# Execution scope

The release exposes a portable check of selected saved numerical outputs and a
separate representative MATLAB integration. It does not provide a single command
that regenerates every flight, uncertainty trial, interval tree, or figure.

## Saved-output verification

Requirements: Python 3.10 or newer, using only the standard library. No package
installation, MATLAB, or network access is needed once the repository is present.

From the repository directory:

```text
python scripts/verify_release.py
```

From another directory, pass the full path to the same script:

```text
python "/absolute/path/to/airpath-reference-governance/scripts/verify_release.py"
```

The wrapper locates the release relative to its own `__file__`, uses the active
Python interpreter in isolated mode, and runs these entries in order:

| Entry | Operation | Writes |
| --- | --- | --- |
| `python/check_manifest.py` | Verify every SHA-256 listed in `RELEASE_MANIFEST.json` | None |
| `python/reproduce.py` | Reconstruct selected paired, ordinary, endpoint-ablation, timing, and direct-center summaries from saved outputs | Root `generated/` |
| `python/verify.py` | Compare those reconstructed results with archived summaries and manuscript counts | `generated/verification.json` |

Each entry resolves its inputs from its own file location. The third requires the
second to have completed. The wrapper stops on the first failure, so stale output
cannot turn a failed reconstruction into a successful overall check. Existing
files in `generated/` may be replaced; frozen files in `data/` and `source/` are
not rewritten by these Python entries.

The check verifies the saved-output identities and selected statistics specified
by the scripts. It does not rerun the aircraft model or independently regenerate
the interval certificates. `MAIN_AND_SUPPLEMENT_MAP.csv` and
`DATA_DICTIONARY.txt` describe how the archived outputs relate to the manuscript.

## Representative MATLAB integration

The retained integration environment is MATLAB R2018a. Add only the release's
`matlab/` directory to the path and call the wrapper; do not recursively add all
historical study folders to the MATLAB path.

```matlab
addpath('/absolute/path/to/airpath-reference-governance/matlab');
run_representative_direct;
```

The default output is `generated/matlab_C1_005` inside this release. It must not
already exist. An optional character-vector argument selects a different output
directory, including a path with spaces:

```matlab
run_representative_direct('C1 check output');
```

A relative argument is resolved against the caller's current directory before
the wrapper changes folders. The wrapper copies `source/` into that new output
directory and invokes the original `run_direct('C1', .005)` on the copy. It restores
the caller's MATLAB path and current directory when execution returns or errors.
It refuses an existing output directory rather than overwriting a prior run.

The underlying case integrates the Direct configuration to 35 s, with 701 saved
points. It uses the frozen C1 F-configuration input needed by the original
protocol. After completion, verify the trace and derived metrics separately:

```text
python "/absolute/path/to/airpath-reference-governance/python/verify_matlab.py" --run-folder "/absolute/path/to/C1 check output"
```

The default `--run-folder` is the default directory above. The checker compares
all 30 saved numeric fields at 701 points and the selected speed-error,
position-RMS, and feedback-gap metrics with the archived C1 Direct record.
This is a representative one-case reintegration, not a new validation sample or
an all-study repeat. The historical check is recorded in
`validation/MATLAB_REINTEGRATION.json`.

The release wrapper resolves its output path before its internal directory
change. Its cleanup callback captures the caller's path and folder as values and
restores them through a file-local function, including in MATLAB R2018a. These
are wrapper fixes only. All files under `source/` retain their archived bytes;
no controller, model parameter, integration setting, or scientific function was
changed for the GitHub entry point.

## Historical workflow sources

`SOURCE_WORKFLOWS.txt` indexes the original stages. These files preserve the
original workflow and provenance; they are not additional portable release
commands. Several launchers or reports still refer to the original machine or
to input/output archives that are not all bundled here. Examples include:

| Archived file | Retained environment assumption |
| --- | --- |
| `source/endpoint_probe_20261106/static/run_static_endpoint_probe.m` | Source snapshot under `D:/Research/paper1_Submission/...` |
| `source/decel_d0_20261108/analysis/run_replays.py` | MATLAB executable at `C:/Program Files/MATLAB/R2018a/bin/matlab.exe` |
| `source/decel_confirm_20261108/run_all.py` | `matlab` on PATH or the original Windows MATLAB installation |
| `source/decel_center_20261108/report.py` | Windows font `C:/Windows/Fonts/msyh.ttc` for report rendering |
| `source/revision_wp3_20261107/` and `source/revision_wp4_20261107/` | Original stage paths and their input/output archives |

Historical drivers can write into their own `runs/` directories. Some of those
directories and `.mat` files are frozen inputs required by later stages, so they
remain tracked and are deliberately not globally ignored. Adapt or investigate
historical launchers only in a separate working copy. Their dependencies and
full-study execution have not all been validated by the release entry points.

## Continuous integration and byte preservation

The prepared GitHub Actions workflow runs the saved-output command on Linux and
Windows with Python 3.10 and 3.14, from a working directory outside the checkout.
It does not launch MATLAB. A workflow definition alone is not evidence of a
successful remote run; inspect the repository's Actions run after publication.

The workflow pins the official releases checked on 18 November 2026:
[checkout v7.0.1](https://github.com/actions/checkout/releases/tag/v7.0.1) and
[setup-python v7.1.0](https://github.com/actions/setup-python/releases/tag/v7.1.0).

`.gitattributes` disables line-ending normalization for all files so Git preserves
the archived SHA-256 identities across platforms. `.gitignore` excludes only
root-level generated output and local caches/editor files; it does not ignore
scientific CSV, JSON, MAT, or historical `runs/` inputs. Keep scientific source
and data changes separate from generated verification output, and regenerate
the release manifest only when preparing a deliberately updated release.
