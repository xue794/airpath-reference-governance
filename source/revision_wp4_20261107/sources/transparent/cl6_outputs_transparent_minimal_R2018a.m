function Y = cl6_outputs_transparent_minimal_R2018a(t,x,mode,C)
%CL6_OUTPUTS_TRANSPARENT_MINIMAL_R2018A
% Diagnostics required for the controlled-ablation comparison figure.

n = numel(t);

Y.er = zeros(n,3);
Y.rd = zeros(n,3);
Y.r = x(:,1:3);
Y.lambda = ones(n,1);
Y.minPhi = zeros(n,1);
Y.xf = x(:,20:25);

for k = 1:n
    S = cl6_control_snapshot_transparent_R2018a(t(k),x(k,:).',mode,C);

    Y.er(k,:) = S.Gd.er.';
    Y.rd(k,:) = S.Gd.rd.';
    Y.lambda(k) = S.Gov.lambda;
    Y.minPhi(k) = min(S.Gov.phi);
end
end
