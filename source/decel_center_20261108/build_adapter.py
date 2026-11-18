"""Build the direct-reference adapter from the archived production equations."""
from pathlib import Path
H=Path(__file__).resolve().parent
def single(src,old,new):
    assert src.count(old)==1,old
    return src.replace(old,new)
def build():
    root=H.parent
    ref=(root/'decel_region_20261108/moving_backup_affine.m').read_text(encoding='utf8')
    schedule=ref[ref.index('if t < 0'):ref.index('center=[v;')]
    blend=ref[ref.index('function [v,n,a]=blend'):]
    (H/'direct_reference.m').write_text(
        'function [xf,accel]=direct_reference(t,R)\n'+schedule+
        'xf=[v;0;0;n;0;0];accel=[a;0;0];\nend\n'+blend,encoding='utf8')
    s=(root/'decel_d0_20261108/d0_snapshot.m').read_text(encoding='utf8')
    s=single(s,'function S=d0_snapshot(t,x,mode,C,method,O)','function S=direct_snapshot(t,x,mode,C,O)')
    s=single(s,'if nargin<6||isempty(O),O=w2d_options();end',"assert(nargin==5 && isfield(O,'WindFun'));")
    s=single(s,'[cc,chart]=rg_raw_command_chart(Gd.vc,C.R);u0=rg_nominal_control(xf,cc,C.R);Gov=d0_governor(t,xf,u0,C.R,method,O);','% The analytic six-component reference was inserted in x upstream of this function.')
    s=single(s,'xfDot=[xf(4:6);Gov.u];','')
    s=single(s,'S.cc=cc;S.chart=chart;S.u0=u0;S.Gov=Gov;S.xfDot=xfDot;','')
    (H/'direct_snapshot.m').write_text(s,encoding='utf8')
if __name__=='__main__':build()
