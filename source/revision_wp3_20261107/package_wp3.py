from pathlib import Path
import json,hashlib,zipfile,datetime,shutil,csv
from lxml import etree
H=Path(__file__).resolve().parent;B=H.parent/'revision_stage1_20261107'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
baseline=json.loads((B/'baseline_manifest.json').read_text(encoding='utf-8'))
for r in baseline:
    assert sha(r['source'])==r['sha256'],r['source']
    assert sha(B/r['frozen'])==r['sha256'],r['frozen']
independent=json.loads((H/'checks/independent_verification.json').read_text())
assert independent['all_passed'] and independent['interval_enclosure_violations']==0
summary=json.loads((H/'results/transfer_summary.json').read_text())
assert sum(x['samples'] for x in summary)==14407
assert all(x['counts']['certified_J']==x['samples'] for x in summary)
assert all(x['counts']['inclusion_G']==x['samples'] for x in summary)
assert all(x['counts']['attitude_S4_ok']==x['samples'] for x in summary)
doc=H/'Aerospace_WP3反馈裕度传递核查_v1_0.docx'
with zipfile.ZipFile(doc) as z:xml=etree.fromstring(z.read('word/document.xml'))
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','m':'http://schemas.openxmlformats.org/officeDocument/2006/math'}
assert len(xml.findall('.//m:oMath',ns))==10
assert len(xml.findall('.//w:tbl',ns))==8
maxdiff=max(max(s['reconstruction_max_errors']) for s in summary)
assert maxdiff<8e-13
runtime=json.loads((H/'inputs/runtime.json').read_text())
runtime.update(matlab='MATLAB R2018a',matlab_executable=r'C:\Program Files\MATLAB\R2018a\bin\matlab.exe',
    document_python=r'C:\Users\IMUT\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe',
    interval_algorithm='mpmath.iv 40 decimal digits; independent mp 70 digits',new_integrations=0)
write(H/'inputs/runtime.json',runtime)
shutil.copy2(H.parent/'major_revision_plan_20261106/export_word.ps1',H/'export_word.ps1')
status=dict(date='2026-11-07',baseline=str(B),WP0='completed in stage1; preserved',WP1='completed bounded comparison; novelty adequacy remains open',
    WP2='completed in stage1; WP3 adds local-transfer evidence',WP3='completed: seven retained original trajectories, quantitative conditional instantaneous command transfer',
    WP4='pending: Transparent hold/noisy gaps; any integration or V refinement remains subject to frozen bounded protocol',
    WP5='pending manuscript integration; standalone English and S4 additions now supplied',WP6='pending',WP7='pending',
    new_trajectory_integrations=0,reused_trajectory_cells=7,base_protocol_cells=9,maximum_including_conditional_options=15,
    scientific_boundaries=['finite saved reference points only','no invariant tracking-error tube','no full-KB uniform feedback radius','no endpoint-specific engineering advantage established'])
write(H/'execution_status.json',status)
notes='''WP3 computational notes and change record

The protocol was frozen before the new numerical exports and interval run.
protocol.json SHA256: %s
No feedback gain, reference constraint, wind profile, initial state, backup set,
inversion domain, or sample grid was altered. No integrations were executed.

Compatibility correction: the first MATLAB export assumed that the historical
Transparent Y record contained V and thetaC. Its optional fields differ from
the W2C record. The corrected exporter checks fields only when they exist and
records OptionalSavedFieldAvailable explicitly; absent checks are not evidence
of agreement. Saved er was available; V, Tcmd and thetaC were absent from this
historical Y record and are explicitly flagged unavailable. The independent wind-aware
Python reconstruction checks all seven runs from saved 36-state arrays.
The failed and successful export logs are retained. Four existing exports were
regenerated algebraically; no scientific trajectory was rerun.

Scientific numerical work uses the already available Python 3.14 scientific
environment and the prior mpmath 1.3.0 distribution, copied with its license.
Document authoring and rendering use the bundled document runtime. No runtime
packages were installed. A MathJax rendering attempt found that the module was
unavailable; the final Word report uses structured editable OMML instead.
The unused attempted helper is excluded from the delivery package.

The full-grid maximum reconstruction difference is %.17g (the largest entry
is a force/thrust quantity in SI units). Earlier interim commentary referred to
approximately 6e-13 before all trajectories were combined; the final complete
record and report use the bound 8e-13. This does not affect any classification.

The primary attitude condition is the S4 sum, not the direct attitude distance.
Both primary and fixed-bound direct-distance inclusion counts coincide here.
The predeclared posthoc test using the observed direct attitude angle gives
1587/2401 for Transparent varying, compared with the primary 1582/2401; every
other count is unchanged. It does not replace the primary result, nor establish
an invariant or future error bound.

The independent 70-digit calculation checks 2268 error states (8 corners and
the center, 12 equally spaced saved indices, 3 layers, 7 trajectories), 6804
point derivative evaluations and 756 roots. Finite checks verify implementation
consistency; interval arithmetic, not the corners, supplies the whole-box bound.

Joint radius zeros, if returned, would be uncertified and excluded from all
inclusion counts. All 14407 primary joint radii are positive in this dataset.
Geometry radii reach the fixed ceiling; no maximal-radius claim is made.
Raw thrust positivity uses a strict test. Searches provide sufficient numerical
half-widths, not optimal radii or exact physical parameter certification.

No conclusions are drawn about actual aircraft hardware certification, all
between-sample times, or actual-error invariance under active constraints.
'''%(sha(H/'inputs/protocol.json'),maxdiff)
(H/'checks/technical_notes.txt').write_text(notes,encoding='utf-8')
(H/'qa/visual_review.txt').write_text('All ten report pages inspected at full image resolution. Native equations, tables and the figure are legible with no clipping or overlap. After the final numeric-precision wording correction, only page 10 changed; the unchanged nine page image hashes match the earlier inspected render. Page 10 was reinspected. Three standalone scientific PNGs were inspected; thrust limits now display the full 0–60 N range.\n',encoding='utf-8')
readme='''Aerospace WP3 feedback-margin transfer audit
Completed 2026-11-07. No new trajectory integrations.

Read Aerospace_WP3反馈裕度传递核查_v1_0.docx first (10 pages).
The TXT companion supplies searchable prose; native math is in the Word report
and the complete LaTeX addition is in S4_quantitative_transfer_draft.tex.
manuscript_insert_EN.txt contains proposed main-text and response wording.
These additions have not been merged into the original submitted manuscript.

Results: seven original trajectories, 14407 saved points, three nested sufficient
error-box half-widths E/G/J. Use summary_table.json and *_transfer.csv for counts;
*_bounds.csv contains interval bounds. Positive certificate and sample inclusion
are distinct fields. All claims concern the specified nominal command interface.
No continuous-time actual-aircraft safety or error-invariance claim is supplied.

Independent verification passed. Original baseline 472-file identities, all
seven input trajectories and 315 frozen source files remain unchanged.
The result of WP3 and remaining WP4/5/6/7 tasks is in execution_status.json.

Reproduction in the original workspace:
1. setup_wp3.py freezes inputs and sources; do not use it to redefine protocol.
2. MATLAB R2018a: addpath('D:/Research/paper1_Submission/revision_wp3_20261107');
   export_wp3;  This is algebraic postprocessing, not simulation.
3. Scientific Python with numpy/scipy/matplotlib:
   python interval_wp3.py
   python check_interval_wp3.py
   python assemble_wp3.py
   A local licensed mpmath 1.3.0 copy is supplied in vendor.
4. Bundled document Python: python build_wp3_report.py
   Export through export_word.ps1 and rasterize using render_report.py for QA.
   Document rendering helpers expect the original Codex skill/runtime location.

The package contains the exported 36-state arrays, configuration, algebraic
diagnostics, frozen code, interval results, scripts and figures. The interval
analysis and assembly can run using these exports without the large original
MAT archives. Source-integrity verification and re-export require the original
input paths recorded in inputs; original large MAT files are referenced by
SHA256 and are not duplicated here. It is not a full submission-data distribution.

PNG/SVG figures are scientific outputs. Red shading refers only to sampled
joint-condition gaps, never to observed aircraft failure. QA page images/PDF
and the unused MathJax helper are excluded from the ZIP.
'''
(H/'README.txt').write_text(readme,encoding='utf-8')
validation=dict(date='2026-11-07',original_stage1_baseline_files_verified=len(baseline),original_trajectories_verified=7,
    frozen_source_files_verified=315,samples=14407,independent_all_passed=True,max_reconstruction_difference=maxdiff,
    report_pages=10,native_equations=10,tables=8,new_integrations=0,all_scientific_inputs_preserved=True)
write(H/'checks/delivery_validation.json',validation)
def allowed(p):
    rel=p.relative_to(H)
    return p.is_file() and '__pycache__' not in rel.parts and rel.parts[0]!='qa' and p.name not in ['artifact_manifest.json','render_equations.cjs']
paths=sorted(p for p in H.rglob('*') if allowed(p))
manifest=[dict(path=p.relative_to(H).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in paths]
write(H/'artifact_manifest.json',manifest)
target=H.parent/'Aerospace_WP3反馈裕度核查材料包_v1_0.zip'
with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in paths+[H/'artifact_manifest.json']:z.write(p,arcname=H.name+'/'+p.relative_to(H).as_posix())
with zipfile.ZipFile(target) as z:
    assert z.testzip() is None
    for r in manifest:assert hashlib.sha256(z.read(H.name+'/'+r['path'])).hexdigest()==r['sha256']
out=dict(path=str(target),bytes=target.stat().st_size,sha256=sha(target),files=len(paths)+1,verified=True)
write(H/'qa/package_verification.json',out)
print(json.dumps(out,ensure_ascii=False))
