function phi = rg_constraint_direct(x,u,R)
%RG_CONSTRAINT_DIRECT Evaluate the affine CBF/HOCBF constraints.
C=rg_constraint_affine(x,R);
phi=C.c+C.G*u;
end
