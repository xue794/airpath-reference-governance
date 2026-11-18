function J=viii_collect_macro_jacobians(Data,C,P8)
%VIII_COLLECT_MACRO_JACOBIANS Feasible-center normalized Jacobians.

Pts=viii_collect_witness_points(Data,C,P8);

A=cell(numel(Pts),1);
f=cell(numel(Pts),1);
y=cell(numel(Pts),1);

for k=1:numel(Pts)
    L=viii_local_jacobian(Pts(k).t,Pts(k).X,Pts(k).mode,C,P8);
    A{k}=L.A;
    f{k}=L.f;
    y{k}=L.y;
end

J.points=Pts;
J.A=A;
J.f=f;
J.y=y;
end
