function Results = test_open_loop()
%TEST_OPEN_LOOP  P1--P3 validation for the T-2 GTM--GNA open-loop plant.
P=t2_gna_params();
fprintf('\nT-2 GTM--GNA open-loop validation\n');
fprintf('================================\n');

%% P1: GNA coefficient regression at source Table-5 condition
alpha=P.src.alpha_rad; de=P.src.deltae_rad; V=P.src.V_mps;
C=t2_gna_coeff(alpha,0,[0;0;0],[0;de;0],V,P);
Cexp=[0.039859948035330;0;0.477110161888924;0;-0.028956967816721;0];
err=max(abs(C-Cexp));
fprintf('P1 coefficient regression: max abs error = %.3e\n',err);
fprintf('   CD=%.9f, CL=%.9f, Cm=%.9f\n',C(1),C(3),C(5));
assert(err<1e-12,'P1 failed');

%% P2: trim of the study plant (body-axis point thrust, no propulsive moment)
TR=t2_trim_level(P.src.V_mps,P.src.h_m,P);
fprintf('P2 reduced-plant level trim:\n');
fprintf('   alpha = %.6f deg\n',rad2deg(TR.alpha));
fprintf('   delta_e = %.6f deg\n',rad2deg(TR.deltae));
fprintf('   T = %.6f N\n',TR.T);
fprintf('   residual norm = %.3e\n',norm(TR.residual));
assert(norm(TR.residual)<1e-6,'P2 failed');

%% P3: finite-difference local modes around the reduced-plant trim
u0=TR.V*cos(TR.alpha); w0=TR.V*sin(TR.alpha);
x0=[u0;0;w0;0;0;0;0;TR.alpha;0];
U=struct('delta',[0;TR.deltae;0],'T',TR.T,'wI',[0;0;0]);
fun=@(x)t2_reduced9_dynamics(x,U,TR.h,P);
A=t2_linearize_fd(fun,x0);
lambda=eig(A);

% Identify modes using robust nearest-pattern logic.
realEig=lambda(abs(imag(lambda))<1e-6 & abs(real(lambda))>1e-8);
[~,ir]=min(abs(abs(real(realEig))-P.src.mode.roll)); roll=realEig(ir);
real2=realEig; real2(ir)=[]; [~,is]=min(abs(abs(real(real2))-P.src.mode.spiral)); spiral=real2(is);
cmp=lambda(imag(lambda)>1e-6);
wn=abs(cmp); zeta=-real(cmp)./wn;
[~,ip]=min(abs(wn-P.src.mode.phugoid_wn)); ph=cmp(ip);
[~,id]=min(abs(wn-P.src.mode.dutch_wn)); dut=cmp(id);
[~,isp]=min(abs(wn-P.src.mode.short_wn)); sho=cmp(isp);

fprintf('P3 local modes (study plant vs Grauer GNA Table 6):\n');
fprintf('   spiral: |lambda| %.5f vs %.5f rad/s\n',abs(spiral),P.src.mode.spiral);
fprintf('   roll:   |lambda| %.5f vs %.5f rad/s\n',abs(roll),P.src.mode.roll);
fprintf('   phugoid: wn %.5f vs %.5f, zeta %.5f vs %.5f\n',abs(ph),P.src.mode.phugoid_wn,-real(ph)/abs(ph),P.src.mode.phugoid_zeta);
fprintf('   dutch:   wn %.5f vs %.5f, zeta %.5f vs %.5f\n',abs(dut),P.src.mode.dutch_wn,-real(dut)/abs(dut),P.src.mode.dutch_zeta);
fprintf('   short:   wn %.5f vs %.5f, zeta %.5f vs %.5f\n',abs(sho),P.src.mode.short_wn,-real(sho)/abs(sho),P.src.mode.short_zeta);

Results=struct('P',P,'coeff',C,'trim',TR,'A',A,'eig',lambda, ...
 'modes',struct('spiral',spiral,'roll',roll,'phugoid',ph,'dutch',dut,'short',sho));
end
