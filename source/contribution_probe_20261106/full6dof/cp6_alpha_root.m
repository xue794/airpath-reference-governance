function [alpha,info]=cp6_alpha_root(V,X,N,R,channel)
% Original root always attempted first, so the regular-domain RHS is exact.
% Extension ONLY lets ode45 locate terminal events at the original 0/9deg
% endpoints. No post-event trajectory is accepted as part of the experiment.
[alpha,info]=rg_alpha_root(V,X,N,R);info.originalBracket=info.bracket;info.extended=false;
if ~info.bracket
    global CP6_EXTENSION CP6_EVAL_TIME;
    E=R;E.alphaLo=R.alphaPhysicalMin;E.alphaHi=R.alphaPhysicalMax;
    [alpha,J]=rg_alpha_root(V,X,N,E);
    info=J;info.originalBracket=false;info.extended=true;
    if ~isempty(CP6_EXTENSION)
        CP6_EXTENSION.(channel)=CP6_EXTENSION.(channel)+1;
        CP6_EXTENSION.last_time_s=CP6_EVAL_TIME;
    end
end
end
