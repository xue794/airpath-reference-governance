Computational materials for aerospace-4638229
Air-Path Reference Governance for Fixed-Wing Tracking under Wind
Prepared from the existing numerical record, 18 November 2026.

START HERE
1. Run: python python/check_manifest.py
2. Run: python python/reproduce.py
3. Run: python python/verify.py
These entries use Python 3.10 or later and its standard library only.
They write to generated/. No MATLAB simulation is launched.

WHAT IS RECONSTRUCTED
The script recomputes the 500-pair completion table, common-completion RMS
differences, ordinary saved-grid position metrics and thrust/membership counts,
the 40,817-state endpoint deletion counts, timing medians, and the C1/C2
Direct/O/F speed contrasts from the included fine-step trajectories.
Verification compares them with separately archived summaries or stated counts.
Expected summaries are not read by reproduce.py.

The remaining CSV/JSON records provide the existing local/whole-region
certificate bounds, coverage, events, and contrasts. Including these records
is not a new independent interval proof. MAIN_AND_SUPPLEMENT_MAP.csv maps
the numbered figures/tables and their supporting records. DATA_DICTIONARY.txt
states the units, averaging windows and missing-value conventions.

SOURCE AND FROZEN INPUTS
source/ retains the original relative layout of scientific MATLAB/Python
modules and frozen source snapshots, plus all 594 protocol dependency checks
for D0, C1-C3 and Direct. See SOURCE_WORKFLOWS.txt for original entry points.
SOURCE_MANIFEST.json identifies every copied file and its source and SHA-256.
RELEASE_MANIFEST.json checks this release, including the documentation and
standalone Python reconstruction code.

REPRESENTATIVE MATLAB REINTEGRATION
In a fresh MATLAB R2018a session opened at this release root, run:
  addpath(fullfile(pwd,'matlab'))
  run_representative_direct
Then, from the release root in a terminal, run:
  python python/verify_matlab.py
The entry copies the frozen source tree to generated/matlab_C1_005/source
and calls the original run_direct('C1',.005) without changing its scientific
source or parameters. It refuses an existing output directory. To repeat a
run, give a new absolute output path to run_representative_direct and pass
the same path with --run-folder to verify_matlab.py. Original inputs and
archives are never overwritten by this wrapper.

This representative 35 s case completed in 52.09 s of integration during
release checking. All 21,030 values in its 701-point saved trace matched the
archive exactly, as did speed MAE, position RMS and the 14-20 s upper feedback
gap. The comparison tolerance is 1e-8 in each column's stated units; exact
equality on this workstation is not a claim of bitwise identity elsewhere.
See validation/MATLAB_REINTEGRATION.json.

The other original MATLAB workflows require MATLAB R2018a and their declared
dependencies and were not reintegrated during this release preparation.
Some historical Python launchers contain workstation-specific paths
and are supplied as provenance, not portable execution commands. Do not use
addpath(genpath(source)): study folders intentionally contain different
versions of functions with the same names. Read SOURCE_WORKFLOWS.txt first.

IMPORTANT DISTINCTIONS
data/original_comparison/ preserves W3 tables and certificate records. Its
ordinary summaries are historical; final ordinary layered data are in
data/ordinary/. Proposed is the source name for Scalar. MatchedQP is QP.
The 408 common completions are a conditional population, not all 500 trials.
RAW_COMMAND_CHART_EXIT is an upstream coordinate-domain exit, not a QP failure.
Reference, raw feedback and applied actuator thrust are different variables.
The 0.606 m/s scale is a research scale, not a mission acceptance threshold.
No new research case, hardware result, stability theorem or public release is
implied by this package. The complete physical histories of all 1000
uncertainty runs are not included. No all-study portable reintegration claim
is made. The unchanged Supplement_S11_Data.zip is retained in data/ and also
expanded under data/S11/ for navigation.
