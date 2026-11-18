function [cls,domainExit,numericalFailure]=w2e_classify_error(id)
%W2E_CLASSIFY_ERROR Frozen error taxonomy; all failures remain in the ledger.
if nargin<1||isempty(id),id='';end
domainExit=false;numericalFailure=false;
if strcmp(id,'rg_raw_command_chart:singular')
 cls='RAW_COMMAND_CHART_EXIT';domainExit=true;
elseif strcmp(id,'rg_scalar_governor:backup')
 cls='PROPOSED_BACKUP_FEASIBILITY_EXIT';domainExit=true;
elseif strcmp(id,'rg_constraint_affine:root') || strcmp(id,'rg_feedforward_rate:root') || strncmp(id,'rg_alpha_root:',14)
 cls='REFERENCE_INVERSION_DOMAIN_EXIT';domainExit=true;
elseif strcmp(id,'viA_feedback_command:root')
 cls='FEEDBACK_INVERSION_DOMAIN_EXIT';domainExit=true;
elseif strcmp(id,'W2B:QPSolver')
 cls='MATCHED_QP_SOLVER_FAILURE';numericalFailure=true;
elseif strcmp(id,'W2E:kappa') || strcmp(id,'cl6_airpath_state:V') || strcmp(id,'t2_gna_aero:V') || strcmp(id,'t2_gna_coeff:V') || strcmp(id,'vaB_to_airdata:ZeroSpeed') || strcmp(id,'quat_normalize:ZeroNorm') || strcmp(id,'dcm_to_euler321:Singularity')
 cls='PHYSICAL_GEOMETRY_DOMAIN_EXIT';domainExit=true;
else
 cls='NUMERICAL_OR_OTHER_FAILURE';numericalFailure=true;
end
end
