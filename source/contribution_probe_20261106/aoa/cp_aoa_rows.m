function A=cp_aoa_rows(x,R,mode,alpha)
%CP_AOA_ROWS Exact coupled-angle CBF rows on reference acceleration.
% A=cp_aoa_rows(x,R,'A0') uses fixed [alphaLo,alphaHi], k=R.kN.
% A=cp_aoa_rows(x,R,'Am') uses force-reserve-matched MOVING boundaries:
%   M(aL;V,X)=M(alphaLo;V,X)+DeltaN,
%   M(aU;V,X)=M(alphaHi;V,X)-DeltaN.
% Output A.c+A.G*u>=0 contains LOWER then UPPER angle-CBF conditions.
% Append these to frozen rows [1:13,16:17] for the matched comparison.
% Optional alpha reuses the already-computed matching coupled root; residual
% validation protects against an unrelated root. It may be an original-endpoint
% EVENT-LOCATION trial root in [alphaPhysicalMin,alphaPhysicalMax]; the caller
% must terminate at the ORIGINAL endpoint and reject all post-event trajectory.
% This function does not authorize an enlarged operating interval.
% This is a nominal reference kernel, not a physical-aircraft certificate.
V=x(1);gam=x(2);nuV=x(4);nuG=x(5);nuChi=x(6);m=R.P.m;
X=m*(nuV+R.g*sin(gam));
ag=V*nuG+R.g*cos(gam);ac=V*cos(gam)*nuChi;rr=hypot(ag,ac);
if V<=0||rr<=0,error('CPAoA:Domain','Invalid reference speed/normal geometry.');end
N=m*rr;
Mlo=rg_N_map(R.alphaLo,V,X,R);Mhi=rg_N_map(R.alphaHi,V,X,R);
supplied=(nargin>=4&&~isempty(alpha));
if ~supplied
 [alpha,info]=rg_alpha_root(V,X,N,R);
 if ~info.bracket,error('CPAoA:Root','Current coupled root not bracketed.');end
end
Ma=rg_N_map(alpha,V,X,R);
checkLo=R.alphaLo;checkHi=R.alphaHi;
if supplied&&isfield(R,'alphaPhysicalMin')&&isfield(R,'alphaPhysicalMax')
 checkLo=min(checkLo,R.alphaPhysicalMin);checkHi=max(checkHi,R.alphaPhysicalMax);
end
if ~isfinite(alpha)||alpha<checkLo-1e-12||alpha>checkHi+1e-12|| ...
   abs(Ma.value-N)>1e-9*(1+abs(N))||Ma.alpha<=0
 error('CPAoA:SuppliedRoot','Invalid, mismatched or nonregular coupled root.');
end

dag_c=nuV*nuG-R.g*sin(gam)*nuG;dag_G=[0 V 0];
dac_c=nuV*cos(gam)*nuChi-V*sin(gam)*nuG*nuChi;
dac_G=[0 0 V*cos(gam)];
Ndot_c=m*(ag*dag_c+ac*dac_c)/rr;
Ndot_G=m*(ag*dag_G+ac*dac_G)/rr;
Xdot_c=m*R.g*cos(gam)*nuG;Xdot_G=[m 0 0];
adot_c=(Ndot_c-Ma.V*nuV-Ma.X*Xdot_c)/Ma.alpha;
adot_G=(Ndot_G-Ma.X*Xdot_G)/Ma.alpha;

switch upper(mode)
 case 'A0'
  delta=0;al=R.alphaLo;au=R.alphaHi;
  al_c=0;al_G=zeros(1,3);au_c=0;au_G=zeros(1,3);
  Ml=Mlo;Mu=Mhi;
 case 'AM'
  delta=R.DeltaN;
  if delta<0,error('CPAoA:Reserve','DeltaN must be nonnegative.');end
  if Mhi.value-Mlo.value<2*delta
   error('CPAoA:EmptyMatchedInterval','No interval has both requested force reserves.');
  end
  if delta==0
   al=R.alphaLo;au=R.alphaHi;Ml=Mlo;Mu=Mhi;
  else
   [al,il]=rg_alpha_root(V,X,Mlo.value+delta,R);
   [au,iu]=rg_alpha_root(V,X,Mhi.value-delta,R);
   if ~il.bracket||~iu.bracket||al>au
    error('CPAoA:MatchedBoundary','Force-reserve boundary inversion failed.');
   end
   Ml=rg_N_map(al,V,X,R);Mu=rg_N_map(au,V,X,R);
  end
  if Ml.alpha<=0||Mu.alpha<=0,error('CPAoA:BoundaryRegularity','Nonpositive boundary derivative.');end
  al_c=((Mlo.V-Ml.V)*nuV+(Mlo.X-Ml.X)*Xdot_c)/Ml.alpha;
  al_G=(Mlo.X-Ml.X)*Xdot_G/Ml.alpha;
  au_c=((Mhi.V-Mu.V)*nuV+(Mhi.X-Mu.X)*Xdot_c)/Mu.alpha;
  au_G=(Mhi.X-Mu.X)*Xdot_G/Mu.alpha;
 otherwise
  error('CPAoA:Mode','Use mode A0 or Am.');
end
h=[alpha-al;au-alpha];
A=struct('c',[adot_c-al_c;au_c-adot_c]+R.kN*h, ...
 'G',[adot_G-al_G;au_G-adot_G], ...
 'names',{{'alpha_lower';'alpha_upper'}}, ...
 'h',h,'alpha',alpha,'alphaLower',al,'alphaUpper',au, ...
 'alphaDot_c',adot_c,'alphaDot_G',adot_G, ...
 'alphaLowerDot_c',al_c,'alphaLowerDot_G',al_G, ...
 'alphaUpperDot_c',au_c,'alphaUpperDot_G',au_G, ...
 'forceMargins',[N-Mlo.value-delta;Mhi.value-N-delta], ...
 'forceReserve',delta,'k',R.kN,'mode',mode, ...
 'Malpha',Ma.alpha,'MalphaLower',Ml.alpha,'MalphaUpper',Mu.alpha, ...
 'suppliedRootOutsideOriginal',alpha<R.alphaLo||alpha>R.alphaHi);
end
