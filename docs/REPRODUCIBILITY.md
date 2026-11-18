# Reproducibility and verification scope

This guide separates reconstruction of reported statistics from reintegration of the aircraft model. The original [README.txt](../README.txt), [SOURCE_WORKFLOWS.txt](../SOURCE_WORKFLOWS.txt) and scientific files are retained without alteration. The current MATLAB release wrapper additionally resolves relative output paths before changing directories; its original version is archived in [validation/original_v6/](../validation/original_v6/). Use the instructions below for the current release wrapper.

## 1. Reconstruct selected summaries from saved outputs

Run the following from the repository root with Python 3.10 or later:

```sh
python scripts/verify_release.py
```

The entry runs the existing sequence:

```sh
python python/check_manifest.py
python python/reproduce.py
python python/verify.py
```

No third-party Python package or MATLAB installation is needed for this sequence. The manifest check tests the identity of the current release files, including repository documentation and tooling. Reconstruction reads the saved outputs and writes derived summaries under `generated/`; it does not read the expected summary values. Verification then compares those derived results with separately archived summaries and declared counts. Inspect the command output and `generated/verification.json` for the result of the current run.

The recorded package check passed **1,101 numerical comparisons**, with a maximum absolute difference of approximately `7.64e-13`. This number describes selected saved-output comparisons, not a number of independently simulated cases. The record is in [SAVED_OUTPUT_RECONSTRUCTION.json](../validation/SAVED_OUTPUT_RECONSTRUCTION.json).

The reconstructed results cover:

- the 500-pair outcome table and RMS contrasts on the 408 common completions;
- the ordinary-case saved-grid position metrics and reference, feedback and actuator thrust diagnostics;
- the 40,817-state endpoint-deletion audit and its 21 changed updates;
- the archived timing-block summaries;
- C1/C2 speed-error contrasts between O, F and Direct on the declared common window.

Other included records document numerical certificates, coverage and events. Reading or hash-checking those records does not independently regenerate their interval proofs. The [43-entry evidence map](../MAIN_AND_SUPPLEMENT_MAP.csv) states the verification scope for each figure or table. The [data dictionary](../DATA_DICTIONARY.txt) specifies units, windows and unavailable-value conventions. Reading the archived NPZ figure inputs requires NumPy; those arrays are outside the standard-library reconstruction entry.

## 2. Rerun the representative MATLAB case

The checked reintegration uses **MATLAB R2018a**. Start a fresh MATLAB session at the repository root and run:

```matlab
addpath(fullfile(pwd,'matlab'))
run_representative_direct
```

The wrapper copies the frozen `source/` tree into `generated/matlab_C1_005/source/`, then calls the original `run_direct('C1',.005)`. It leaves the scientific functions and parameters unchanged and refuses to reuse an existing output directory.

After MATLAB finishes, run this command from the repository root:

```sh
python python/verify_matlab.py
```

The checker compares the new 701-point trace, completion and event status, ground-speed MAE over 3–16 s, full-flight position RMS, and minimum upper feedback gap over 14–20 s with the archive. It writes `verification.json` in the run folder. The default trace tolerance is `1e-8` in each column's stated units.

For another run, provide a **new absolute or relative output path** as a character vector to `run_representative_direct(outputRoot)`. A relative path is resolved from MATLAB's working directory at the time of the call. Pass the same directory to `python python/verify_matlab.py --run-folder PATH`; an absolute path avoids ambiguity if the terminal uses a different working directory. Quote paths containing spaces. Inspect an existing run rather than overwriting it.

In the recorded preparation check, the case completed to 35 s and all 21,030 saved numeric values matched exactly, as did the three metrics. Integration took 52.09 s on that workstation. Runtime and bitwise agreement are environment-dependent; the recorded result does not guarantee either on another system. See [MATLAB_REINTEGRATION.json](../validation/MATLAB_REINTEGRATION.json).

## 3. Use the historical workflows with their original scope

The source archive preserves all 594 frozen protocol dependencies at their recorded relative paths. [FROZEN_INPUT_HASH_CHECK.json](../validation/FROZEN_INPUT_HASH_CHECK.json) documents identity with the archived D0, confirmation and Direct protocol inputs. It is an input-identity check, not an independent simulation result.

[SOURCE_WORKFLOWS.txt](../SOURCE_WORKFLOWS.txt) identifies the original uncertainty, endpoint, local-transfer and moving-region entries. Some historical Python launchers contain workstation-specific paths, and some MATLAB drivers write into their study directories. They are supplied as provenance and require their documented inputs and environment before further execution. The complete physical histories of all 1,000 uncertainty runs are not included.

Do not add the whole archive recursively with `addpath(genpath('source'))`: study directories intentionally retain different versions of functions sharing the same names. Use an isolated copy and the paths specified for the intended study. The representative wrapper follows that approach.

No other MATLAB case or interval tree was regenerated for the recorded package validation. Packaging the old uncertainty records, selected diagnostics, confirmation maneuvers and post hoc Direct comparison together does not change their experimental roles or denominators. Full-study portable reintegration remains outside the verified scope of this release.

## 4. Follow the provenance

Use these records for different questions:

| Record | Question answered |
| --- | --- |
| [MAIN_AND_SUPPLEMENT_MAP.csv](../MAIN_AND_SUPPLEMENT_MAP.csv) | Which saved records support a figure or table, and what has been recomputed? |
| [SOURCE_MANIFEST.json](../SOURCE_MANIFEST.json) | Where did a copied archive file originate, and what is its SHA256 hash? |
| [RELEASE_MANIFEST.json](../RELEASE_MANIFEST.json) | Does the current repository file match its recorded bytes, including added documentation, tooling and licensing files? |
| [Original archive manifest](../validation/original_v6/RELEASE_MANIFEST.json) | What file identities were recorded for the original v6 computational archive? |
| [validation/](../validation/) | Which identity, saved-output and representative reintegration checks were recorded? |

The 1,491-file archive is the source baseline, not the current repository file count. Scientific files in `source/`, archived data and the 594 frozen protocol dependencies retain their original bytes. The release layer adds documentation and tooling and fixes the MATLAB wrapper's path handling; the original wrapper and manifest remain under `validation/original_v6/`.

For an intentional release metadata change, follow [RELEASE_CHECKLIST.md](RELEASE_CHECKLIST.md). Its `scripts/update_manifest.py` step checks the frozen provenance hashes before refreshing the current release manifest, then `scripts/verify_release.py` checks the result. Refreshing a manifest does not validate a scientific change. Generated outputs, `.git/` and Python caches are excluded from the release manifest.
