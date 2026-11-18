function J=viii_local_jacobian(t,X,mode,C,P8)
%VIII_LOCAL_JACOBIAN Local Jacobian in normalized reduced coordinates.
%
% A central difference is attempted first.  Near a piecewise-smooth
% feasibility boundary the perturbation step is deterministically halved,
% up to eight times.  The center itself must remain feasible.

sc=P8.scaleVector(:);
y0=X(:)./sc;
n=numel(y0);
A=zeros(n,n);
usedStep=zeros(n,1);

% The center must be evaluable.
f0=viii_reduced_vector_field(t,X,mode,C,P8)./sc;

for k=1:n
    hk=P8.jacStep;
    ok=false;

    for it=1:8
        yp=y0; ym=y0;
        yp(k)=yp(k)+hk;
        ym(k)=ym(k)-hk;

        Xp=sc.*yp;
        Xm=sc.*ym;

        try
            fp=viii_reduced_vector_field(t,Xp,mode,C,P8);
            fm=viii_reduced_vector_field(t,Xm,mode,C,P8);

            A(:,k)=(fp./sc-fm./sc)/(2*hk);
            usedStep(k)=hk;
            ok=true;
            break;
        catch
            hk=hk/2;
        end
    end

    if ~ok
        error('viii_local_jacobian:Feasibility', ...
            ['Unable to obtain a feasible central finite difference at ' ...
             't=%g, mode=%s, coordinate=%d.'],t,mode,k);
    end
end

J.A=A;
J.f=f0;
J.X=X;
J.y=y0;
J.t=t;
J.mode=mode;
J.usedStep=usedStep;
end
