function [c,info] = rg_raw_command_chart(vc,R)
%RG_RAW_COMMAND_CHART v_c^I -> [V_c; gamma_c; chi_c] on a local S1 chart.
V = norm(vc);
Vh = hypot(vc(1),vc(2));
info.regular = (V>R.Vc_min) && (Vh>R.Vhc_min);
if ~info.regular
    error('rg_raw_command_chart:singular','Raw velocity command leaves certified chart.');
end
gamma = atan2(-vc(3),Vh);
chi = atan2(vc(2),vc(1));
c=[V;gamma;chi];
info.V=V; info.Vh=Vh;
end
