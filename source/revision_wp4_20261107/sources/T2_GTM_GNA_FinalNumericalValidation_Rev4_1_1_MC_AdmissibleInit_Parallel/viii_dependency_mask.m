function M=viii_dependency_mask()
%VIII_DEPENDENCY_MASK Theory dependency pattern in order (z,s,w,p,G,F,R,A).
%
% M(i,j)=true means subsystem i may depend internally on subsystem j.

M=false(8,8);

% I_z={s,G}
M(1,[2 5])=true;

% I_s={w,p}
M(2,[3 4])=true;

% I_w={}
% none

% I_p={F,R,A}
M(4,[6 7 8])=true;

% I_G={z,s,w,p}
M(5,[1 2 3 4])=true;

% I_F={z,s,w,p,G,R,A}
M(6,[1 2 3 4 5 7 8])=true;

% I_R={z,s,w,p,G,F,A}
M(7,[1 2 3 4 5 6 8])=true;

% I_A={z,s,w,p,G,F,R}
M(8,[1 2 3 4 5 6 7])=true;
end
