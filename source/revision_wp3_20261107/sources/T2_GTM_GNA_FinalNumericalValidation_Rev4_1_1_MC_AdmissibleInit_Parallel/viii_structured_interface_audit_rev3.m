function A=viii_structured_interface_audit_rev3(Data,C,P8,P3)
%VIII_STRUCTURED_INTERFACE_AUDIT_REV3 Structured theorem-interface audit.
%
% The dependency graph is fixed by the repaired proof.  Numerical
% differentiation is used only to instantiate the explicit receiving
% interfaces p,G,F,R,A at feasible accepted T-2 centers.

B=viii_blocks();
P2=viii_instance_gain_params_rev2(C,P8);
P2.relStep=P3.relStep;
P2.maxHalves=P3.maxHalves;
P2.sample=P3.sample;

Pts=viii_instance_points_rev2(Data,C,P8,P2);

recv=[4 5 6 7 8];
field={'p','G','F','R','A'};
outdim=[3 6 3 3 4];

% The self block is included ONLY for structured absorption.
src=cell(8,1);
src{4}=[4 6 7 8];
src{5}=[1 2 3 4 5];
src{6}=[1 2 3 4 5 6 7 8];
src{7}=[1 2 3 4 5 6 7 8];
src{8}=[1 2 3 4 5 6 7 8];

alpha=inf(8,1);
alpha(1:3)=P3.alphaExact(1:3);
cij=zeros(8,8);
cij=max(cij,P3.cExact);

rawInterfaceNorm=zeros(8,8);
minNominal=inf(8,1);
minNet=inf(8,1);
worstPoint=zeros(8,1);
roundTripMax=0;
minStep=inf;

% Store center-by-center structured self margins for diagnostics.
selfMargins=nan(numel(Pts),8);

for kp=1:numel(Pts)
    pt=Pts(kp);

    [x0,~]=viii_lift_reduced(pt.t,pt.X,pt.mode,C,P8);
    Xback=viii_reduced_from_full(pt.t,x0,pt.mode,C,P8);
    roundTripMax=max(roundTripMax,norm((Xback-pt.X)./P8.scaleVector,inf));

    O0=viii_instance_interface_values_rev2(pt.t,pt.X,pt.mode,C,P8,P2);
    SS=viii_structured_storage_coeff_rev3(pt.X,O0,C,P8,P3);

    H=cell(5,8);
    for ir=1:5
        i=recv(ir);
        for j=src{i}
            H{ir,j}=zeros(outdim(ir),B.dim(j));
        end
    end

    % Differentiate by SOURCE BLOCK rather than by a monolithic Jacobian.
    for j=1:8
        if ~any(cellfun(@(v) any(v==j),src(recv)))
            continue;
        end

        inds=B.idx{j};
        for kc=1:numel(inds)
            k=inds(kc);
            h=P3.relStep*P8.scaleVector(k);
            if h<=0, h=P3.relStep; end
            ok=false;

            for ih=1:P3.maxHalves
                Xp=pt.X; Xm=pt.X;
                Xp(k)=Xp(k)+h;
                Xm(k)=Xm(k)-h;

                try
                    Op=viii_instance_interface_values_rev2(pt.t,Xp,pt.mode,C,P8,P2);
                    Om=viii_instance_interface_values_rev2(pt.t,Xm,pt.mode,C,P8,P2);

                    for ir=1:5
                        i=recv(ir);
                        if any(src{i}==j)
                            fp=Op.(field{ir});
                            fm=Om.(field{ir});
                            H{ir,j}(:,kc)=(fp-fm)/(2*h);
                        end
                    end

                    minStep=min(minStep,h);
                    ok=true;
                    break;
                catch
                    h=h/2;
                end
            end

            if ~ok
                error('viii_structured_interface_audit_rev3:Feasibility', ...
                    'No feasible perturbation at t=%g, mode=%s, block=%s, col=%d.', ...
                    pt.t,pt.mode,B.labels{j},kc);
            end
        end
    end

    % Structured self absorption and structured cross coefficients.
    for ir=1:5
        i=recv(ir);
        Ri=SS.R{i};
        Ni=SS.N{i};
        minNominal(i)=min(minNominal(i),min(eig(0.5*(Ri+Ri.'))));

        Hii=H{ir,i};
        Qnet=Ri-0.5*(Ni*Hii+(Ni*Hii).');
        aHere=min(eig(0.5*(Qnet+Qnet.')));
        selfMargins(kp,i)=aHere;

        if aHere<minNet(i)
            minNet(i)=aHere;
            worstPoint(i)=kp;
        end

        for j=src{i}
            rawInterfaceNorm(i,j)=max(rawInterfaceNorm(i,j),norm(H{ir,j},2));

            if j~=i
                cij(i,j)=max(cij(i,j),norm(Ni*H{ir,j},2));
            end
        end
    end
end

% Apply numerical instance margin to NON-EXACT cross coefficients.
for i=4:8
    for j=1:8
        if j~=i && cij(i,j)>0
            cij(i,j)=P3.interfaceInflation*cij(i,j);
        end
    end
end

% Self margins are DEFLATED, not cross-inflated.
for i=4:8
    minNet(i)=minNet(i)/P3.interfaceInflation;
end

alpha(4:8)=minNet(4:8);

A.alpha=alpha;
A.c=cij;
A.minNominal=minNominal;
A.minNet=minNet;
A.rawInterfaceNorm=rawInterfaceNorm;
A.selfMargins=selfMargins;
A.worstPoint=worstPoint;
A.points=Pts;
A.nPoints=numel(Pts);
A.roundTripMax=roundTripMax;
A.minStepUsed=minStep;
A.labels=P3.labels;
end
