from pathlib import Path
import json
h=Path(__file__).resolve().parent
s=(h/'ep6_snapshot.m').read_text()
assert s.count('function S=ep6_snapshot(')==1 and s.count('ep6_governor(')==1
s=s.replace('function S=ep6_snapshot(', 'function S=ep6_snapshot_warm(').replace('ep6_governor(', 'ep6_governor_warm(')
(h/'ep6_snapshot_warm.m').write_text(s)
s=(h/'run_ep6.m').read_text()
s=s.replace('function run_ep6(Vlow,courseDegrees,tag,maxStep,turnRise)', 'function run_ep6_warm(Vlow,courseDegrees,tag,maxStep,turnRise,methods)')
s=s.replace("if nargin<5,turnRise=3;end", "if nargin<5,turnRise=3;end\nif nargin<6,methods={'T','N0','N10'};end\nclear ep6_governor_warm;")
s=s.replace("methods={'T','N0','N10'};Results=cell(3,1);\nfor j=1:3", "Results=cell(numel(methods),1);\nfor j=1:numel(methods)")
s=s.replace('ep6_snapshot(', 'ep6_snapshot_warm(')
s=s.replace("'method',method,'status',status", "'method',method,'numerical_start','previous support proposal; unchanged w2b_qp3 verifies KKT','status',status")
(h/'run_ep6_warm.m').write_text(s)
(h/'WARM_START_NOTE.txt').write_text('The original solver and all objective/constraint tolerances are unchanged.\nOnly FeasibleStart proposes a current projection onto previously supporting rows.\nThe original solver independently checks feasibility and KKT and can enumerate.\nThis is a numerical performance change, not a control-law change.\nOriginal baseline runs are retained for same-trajectory verification.\n',encoding='utf-8')
print('Generated warm-start runner using unchanged historical QP solver.')
