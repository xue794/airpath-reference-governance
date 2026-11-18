function A=viii_instance_block_gain_audit_rev2(Data,C,P8,P2)
%VIII_INSTANCE_BLOCK_GAIN_AUDIT_REV2 Sampled local block-interface gains.
%
% This computes blockwise derivatives ONLY for interfaces explicitly present
% in the repaired proof. It does not infer the theorem graph from a full
% frozen-state Jacobian.

B=viii_blocks();
Pts=viii_instance_points_rev2(Data,C,P8,P2);

% Receiving rows with nonlinear interfaces.
recv=[4 5 6 7 8];
field={'p','G','F','R','A'};
outdim=[3 6 3 3 4];

% Allowed source blocks, including self only for absorption diagnostics.
src=cell(8,1);
src{4}=[4 6 7 8];
src{5}=[1 2 3 4 5];
src{6}=[1 2 3 4 5 6 7 8];
src{7}=[1 2 3 4 5 6 7 8];
src{8}=[1 2 3 4 5 6 7 8];

L=zeros(8,8);
usedMin=inf;
roundTripMax=0;

for kp=1:numel(Pts)
    pt=Pts(kp);

    % Center must be exactly representable by the reduced lifting.
    [x0,~]=viii_lift_reduced(pt.t,pt.X,pt.mode,C,P8);
    Xback=viii_reduced_from_full(pt.t,x0,pt.mode,C,P8);
    roundTripMax=max(roundTripMax,norm((Xback-pt.X)./P8.scaleVector,inf));

    % Jacobian columns are computed once and reused for every receiving map.
    J=cell(5,1);
    for ir=1:5
        J{ir}=zeros(outdim(ir),B.n);
    end

    for k=1:B.n
        h=P2.relStep*P8.scaleVector(k);
        if h<=0, h=P2.relStep; end
        ok=false;

        for ih=1:P2.maxHalves
            Xp=pt.X; Xm=pt.X;
            Xp(k)=Xp(k)+h;
            Xm(k)=Xm(k)-h;

            try
                Op=viii_instance_interface_values_rev2(pt.t,Xp,pt.mode,C,P8,P2);
                Om=viii_instance_interface_values_rev2(pt.t,Xm,pt.mode,C,P8,P2);

                for ir=1:5
                    fp=Op.(field{ir});
                    fm=Om.(field{ir});
                    J{ir}(:,k)=(fp-fm)/(2*h);
                end
                usedMin=min(usedMin,h);
                ok=true;
                break;
            catch
                h=h/2;
            end
        end

        if ~ok
            error('viii_instance_block_gain_audit_rev2:Feasibility', ...
                'Could not obtain feasible block perturbation at t=%g, mode=%s, coordinate=%d.', ...
                pt.t,pt.mode,k);
        end
    end

    for ir=1:5
        i=recv(ir);
        for jj=src{i}
            val=norm(J{ir}(:,B.idx{jj}),2);
            L(i,jj)=max(L(i,jj),val);
        end
    end
end

% Exact analytic rows z and s.
L=max(L,P2.Lexact);

% Numerical validation margin, but never inflate exact z/s values twice.
Lexact=P2.Lexact;
maskExact=Lexact>0;
Lnum=L;
Lnum(~maskExact)=P2.interfaceInflation*Lnum(~maskExact);
Lnum(maskExact)=Lexact(maskExact);

A.Lraw=L;
A.L=Lnum;
A.points=Pts;
A.nPoints=numel(Pts);
A.roundTripMax=roundTripMax;
A.minStepUsed=usedMin;
A.labels=P2.labels;
end
