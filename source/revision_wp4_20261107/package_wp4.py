from pathlib import Path
import json,hashlib,zipfile,shutil,datetime
H=Path(__file__).resolve().parent;W=H.parent;P=W/'revision_wp3_20261107';B=W/'revision_stage1_20261107'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
baseline=read(B/'baseline_manifest.json')
for r in baseline:assert sha(r['source'])==r['sha256']==sha(B/r['frozen'])
prior=read(P/'artifact_manifest.json')
for r in prior:assert sha(P/r['path'])==r['sha256'],r['path']
source=read(H/'inputs/source_manifest.json')
for r in source:assert sha(r['source'])==r['sha256']==sha(r['snapshot'])
frozen=read(H/'inputs/run_fingerprint.json')
for r in frozen['files']:assert sha(r['path'])==r['sha256'],r['path']
records=read(H/'inputs/all_trajectories.json');catalog=[];(H/'inputs/reused_raw').mkdir(exist_ok=True)
for r in records:
    assert sha(r['path'])==r['sha256']
    new=r.get('new_integration',False);rel=('runs/' if new else 'inputs/reused_raw/')+r['method']+'_'+r['case']+'.mat'
    if not new:shutil.copy2(r['path'],H/rel)
    assert sha(H/rel)==r['sha256']
    catalog.append(dict(case=r['case'],method=r['method'],relative_path=rel,original_path=r['path'],sha256=r['sha256'],new_integration=new,
        final_status='completed',expected_horizon_s=r['duration_s']))
dump(H/'inputs/portable_trajectory_catalog.json',catalog)
shutil.copy2(P/'results/summary_table.json',H/'analysis/prior_results/summary_table.json')
shutil.copy2(P/'results/Proposed_hold_export.mat',H/'analysis/prior_results/Proposed_hold_export.mat')
run_status=[read(H/'runs'/f'Transparent_{x}_status.json') for x in ['hold','noisy']]
assert all(s['Completed'] and not s['ErrorId'] and not s['WarningId'] for s in run_status)
ind=read(H/'analysis/checks/independent_verification.json');gap=read(H/'checks/negative_reference_verification.json');assem=read(H/'checks/assembly_check.json')
assert ind['all_passed'] and ind['interval_enclosure_violations']==0 and gap['negative_samples']==17
assert assem['samples']==18009 and assem['all_complete'] and assem['all_EG_included'] and assem['all_theta_condition_observed']
assert not assem['all_J_positive']
summary=read(H/'results/summary_table.json');newN=next(r for r in summary if r['method']=='Transparent' and r['case']=='noisy')
assert newN['J_positive']==2384 and newN['J_inclusion']==1567
assert all(r['feedback_margin_violations']==r['modeled_actuator_violations']==r['study_domain_violations']==0 for r in summary)
assert all(r['thrust_target_clipped']==r['thrust_rate_limited']==r['surface_rate_limited']==r['allocator_active']==0 for r in summary)
assert sha(H/'analysis/interval_wp3.py')==sha(P/'interval_wp3.py')
margin=read(H/'results/margin_summary.json')
for r in summary:
    n0=next(x['NegativeSamples'] for x in margin if x['Case']==r['case'] and x['Method']==r['method'] and x['Category']=='reference' and x['Margin']=='reference_T_0N')
    assert n0==r['reference_other_state_violations']
    # In this fixed dataset all non-thrust reference rows are nonnegative.
    assert all(x['NegativeSamples']==0 for x in margin if x['Case']==r['case'] and x['Method']==r['method'] and x['Category']=='reference' and x['Margin'] not in ['reference_T_0N','reference_T_3N','backup'])
doc=H/'Aerospace_WP4记录补齐与分层核查_v1_0.docx'
with zipfile.ZipFile(doc) as z:
    assert 'word/document.xml' in z.namelist()
render=read(H/'qa/render_summary.json');assert render['pages']==7
before={x['name']:x['sha256'].lower() for x in read(H/'qa/reviewed_page_hashes_before_wording_fix.json')}
after={f'page-{i}.png':sha(H/'qa/report'/f'page-{i}.png') for i in range(1,8)}
changed=[k for k in after if before[k]!=after[k]];assert set(changed)<=set(['page-3.png','page-4.png'])
visual=dict(docx_sha256=sha(doc),pages=7,all_pages_visually_inspected=True,
    final_reinspected_changed_pages=changed,other_page_images_match_inspected_hashes=True,page_sha256=after,
    standalone_figures_inspected=['noisy_thrust_layers','noisy_tracking_backup','transparent_transfer','noisy_negative_reference'],
    layout_findings='No clipping, overlap, missing glyphs, split tables or unintended blank pages in the final render.',
    backend='Native Word read-only PDF export and bundled render_docx.py rasterizer; bundled LibreOffice unavailable')
dump(H/'checks/report_visual_check.json',visual)
validation=dict(completed_at=datetime.datetime.now().astimezone().isoformat(),all_required_checks_passed=True,
    baseline_sources_and_snapshots_unchanged=len(baseline),prior_WP3_artifacts_unchanged=len(prior),source_files_unchanged=len(source),
    pre_integration_fingerprint_files_unchanged=len(frozen['files']),raw_trajectory_hashes_verified=len(records),
    original_interval_kernel_unchanged=True,independent_checks_passed=True,negative_reference_check_passed=True,
    original_seven_metrics_unchanged=True,retained_samples=18009,new_integrations=2,all_nine_cells_complete=True,
    scientific_negative_result_retained='17 Transparent noisy points have no positive joint certificate; no independent endpoint benefit established',
    document_pages=7,visual_review_complete=True)
dump(H/'checks/final_validation.json',validation)
status=dict(date='2026-11-07',WP0='completed and preserved',WP1='completed bounded literature comparison; novelty adequacy remains open',
    WP2='completed evidence inventory, updated by WP3 and WP4',WP3='completed original seven; definitions reused unchanged',
    WP4='completed original nine-cell matrix by filling Transparent hold/noisy only',WP5='pending manuscript and supplementary-material integration',WP6='pending',WP7='pending',
    newly_integrated_cells=2,reused_cells=7,total_cells=9,retained_samples=18009,optional_L_executed=False,optional_V_refinement_executed=False,
    additional_MC_executed=False,original_sources_or_submission_modified=False,
    contribution_boundary='Reference-admissibility and selected-set tradeoff documented; endpoint-specific advantage beyond thrust/AoA protection remains unestablished')
dump(H/'execution_status.json',status)
runtime=read(P/'inputs/runtime.json');runtime['new_integrations']=2;runtime['postprocessing_integrations']=0;dump(H/'inputs/runtime.json',runtime)
(H/'checks/technical_notes.txt').write_text('''WP4 computational and provenance notes

Exactly two scientific trajectory integrations were executed: Transparent hold and noisy.
Both completed without exceptions or warnings. No tuning, retry integration, model
extension, new scenario, random resampling, or optional refinement was performed.
Integration wall times in the logs are not governor performance benchmarks.

The 65-ZIP inventory inspects ZIP central directories in the stated root/FWUAV
search scope. It does not claim recursive inspection of arbitrary nested ZIP content.
The matching BASELINE_020 hold/noisy names belong to the W1C Scalar runner with
rhoB=0.20 and are rejected as Transparent data. Original source witnesses are retained.

Preflight checks use both original initial states and 12 equally spaced saved
states per scenario. All downstream RHS entries agree; reference-acceleration
differences equal the selector-input differences exactly. New C, K, I0 and x0
and the original signal source text are frozen before integration.

The first algebraic export succeeded for both new cells. The following all-nine
audit initially stopped before loading the first trajectory because heterogeneous
JSON objects decoded to a MATLAB cell array. The wrapper was corrected to handle
both cell and struct arrays, and the audit then completed. export_audit.log retains
the failed diagnostic wrapper attempt; layered_audit.log retains the completed run.
No trajectory was reintegrated. The integration driver and all frozen inputs match
their pre-run hashes. This is the only postprocessing compatibility correction.

The WP3 interval evaluator is byte-identical. Only exporter function/catalog names,
the acceptable partial-horizon assertion, and the independent checker's input-count
label were adapted. Partial-horizon support was unused because both runs completed.
The interval definition JSON retains WP3's historical creation timestamp and
analysis-only integration count of zero. The authoritative WP4 integration plan is
inputs/protocol_manifest.json. Some copied historical trajectory metadata says
WP3 pending; this is provenance, not current status. Current status is execution_status.json.

Independent 70-digit checks covered 648 error states, 1944 point derivative values,
216 roots, and the saved certificate flags over all 3602 new points. Whole-box
sufficiency comes from interval enclosures, not from the finite corner checks.
The new negative-reference finding prompted an additional independent check of
all 17 negative saved points and their two neighboring points, without integration.
The maximum reference-thrust difference is 1.0245276849119023e-14 N.

All figures connect saved values only. Shading is a guide to sampled condition
inclusion, not an event-localized interval. Zero J returns are unverified, omitted
from positive-radius curves and separately marked. Posthoc observed-attitude
diagnostics never replace the fixed 0.70-degree S4 condition.

Document authoring uses the bundled runtime; native Word read-only export supplies
PDF conversion because bundled LibreOffice is unavailable. The bundled renderer
rasterizes it for full-page inspection. QA PDFs and PNGs are excluded from the bundle.
''',encoding='utf-8')
(H/'evidence_update.txt').write_text('''WP4 证据更新 2026年11月7日

原矩阵缺口已关闭：Transparent hold与noisy各一次新增积分，加原七条，共九条、18,009点。
原七条的数据与WP3局部传递数值保持不变，旧包不覆盖。

新增H：三方法全程位置RMS同为0.027789 m；Transparent与Scalar完整状态相差不超过1.53e-13；联合包含1201/1201。
新增N：Transparent位置RMS为2.416064 m，Scalar 7.750469 m，Matched-QP 8.409029 m。
Transparent N：157点参考推力低于3 N，17点低于0 N（最低-0.277596 N），2339点在KB外。
其原始反馈推力仍在10.681875–27.817714 N内，反演端点、执行器及研究界限无保存点越界。
Transparent N联合包含1567/2401，17个负参考推力点无正J证书，另外817点超出充分盒。

解释更新：参考层可行性和选定集合保持的收益有具体记录；必要物理保护及端点相对已有推力/AoA限制的独立收益仍未证明。
不能把WP3“原七条均有正J半宽”改写成新九条均成立；也不能把原点无正证书写成实际反馈推力越界。

后续限于WP5整合正文、S4、S9及审稿回复。可选L/V加密未启用，不增加恢复备份研究。
''',encoding='utf-8')
(H/'README.txt').write_text('''Aerospace WP4 Transparent record-gap completion
Completed 2026-11-07. Read the seven-page Chinese Word report first.

Scope
Two new integrations only: Transparent hold 60 s / 1201 points and noisy 120 s /
2401 points. Seven existing trajectories reused. Total: 9 cells / 18,009 points.
No new controller, backup design, tuning, Monte Carlo, L scenario or V refinement.

Main artifacts
- Aerospace_WP4记录补齐与分层核查_v1_0.docx and TXT: findings and revision implications.
- manuscript_insert_EN.txt: results and suggested R1/R3 response wording for later integration.
- results/summary_table.json and CSV: all nine cells, including zero certificates.
- results/margin_summary.json and CSV: named reference, derivative, feedback,
  actuator and study-domain margins, units, minima and first negative saved times.
- results/*_layered.mat: all margins, actual applied u, nominal u0, QP KKT values,
  solver paths and actuator/allocation diagnostics for each saved point.
- results/*_layered_timeseries.csv: 28 clearly named columns, including thrust layers.
- runs/: the two new raw 36-state MAT records, statuses and integration logs.
- inputs/reused_raw/: byte-identical copies of the seven reused raw MAT files.
- inputs/portable_trajectory_catalog.json: relative file paths and original hashes.
- analysis/results/: new two-cell interval results and independent verification inputs.
- analysis/prior_results/: original seven-cell transfer series and summary evidence.
- figures/: four scientific figures as PNG and SVG; sampled curves only.
- checks/: preflight, archive search, successful/failed algebraic logs, negative-reference
  check, source/baseline validation and document-review record.

Read the negative finding carefully
Transparent noisy has 17 negative reference-thrust samples and 17 zero J returns.
These are not actual feedback/actuator violations. Its full-record position RMS is
smaller than Scalar and Matched-QP. This audit does not establish an endpoint-specific
advantage over a matched thrust or angle-of-attack construction.

Reproduce postprocessing in the original workspace
Use the executables and library versions in inputs/runtime.json.
1. MATLAB: addpath(the WP4 root); export_wp4; audit_wp4;
2. Copy the two *_export.mat files and export_summary.json from results to analysis/results.
3. Scientific Python: analysis/interval_wp3.py
4. Scientific Python: analysis/check_interval_wp3.py
5. Scientific Python: assemble_wp4.py and check_reference_gap.py
These steps do not integrate. Original all_trajectories.json uses historical absolute
paths; audit's hold comparison and assembly's prior-metric check also reference the
adjacent WP3 directory. The package preserves the original executable audit, not a
path-rewritten claim of portable reruns. For another host, use the relative catalog
and included prior_results copies to relocate those paths in a separate working copy,
retaining original provenance and hashes. Raw inputs and source snapshots are included.

Do not execute run_wp4_case in the delivered folder to reproduce the existing result:
it intentionally refuses an existing final MAT. An independent scientific rerun
must use a fresh copy/output directory and be labeled as a new integration.
prepare_wp4.py/preflight_wp4.m/prepare_analysis.py document the original preparation
sequence; do not rerun them over this frozen evidence archive.

Bundle integrity
python verify_bundle.py validates all packaged relative-path checksums using the
standard library. artifact_manifest.json is the package inventory. Original source
and baseline checks were performed separately on the task host and are recorded in
checks/final_validation.json. Delivery ZIP checksum is recorded outside the ZIP to
avoid a recursive checksum. Historical status fields are retained as provenance;
execution_status.json gives the current phase status.

The original submission, stage1 baseline and WP3 deliverables remain unchanged.
WP5 manuscript integration is the next stage and has not been performed here.
''',encoding='utf-8')
exclude={'artifact_manifest.json','checks/delivery_check.json','analysis/assemble_wp3.py'}
files=[]
for p in sorted(H.rglob('*')):
    if not p.is_file():continue
    rel=p.relative_to(H).as_posix()
    if rel in exclude or rel.startswith('qa/') or '__pycache__' in p.parts or p.suffix=='.pyc' or '_checkpoint' in p.name:continue
    files.append(p)
manifest=[dict(path=p.relative_to(H).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in files]
dump(H/'artifact_manifest.json',manifest)
bundle=W/'Aerospace_WP4_Transparent记录补齐材料包_v1_0.zip'
with zipfile.ZipFile(bundle,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in files+[H/'artifact_manifest.json']:z.write(p,H.name+'/'+p.relative_to(H).as_posix())
with zipfile.ZipFile(bundle) as z:
    assert z.testzip() is None
    for r in manifest:assert hashlib.sha256(z.read(H.name+'/'+r['path'])).hexdigest()==r['sha256']
delivery=dict(file=str(bundle),bytes=bundle.stat().st_size,sha256=sha(bundle),entries=len(manifest)+1,CRC_passed=True,all_entry_hashes_passed=True,
    final_docx=str(doc),new_integrations=2,total_cells=9,retained_samples=18009)
dump(H/'checks/delivery_check.json',delivery)
print(json.dumps(delivery,ensure_ascii=False))
