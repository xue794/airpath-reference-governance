function U=cl6_plant_input(delta,T,mode,C)
%CL6_PLANT_INPUT Build the exact input structure for t2_gna_6dof_quat.

K=cl6_case_info(mode,C);

U.delta=delta(:);
U.T=T;
U.wind_fun=@(tt)cl6_wind_profile(tt,K.windMode,C);
end
