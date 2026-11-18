try
 addpath(fileparts(mfilename('fullpath')));preflight_comparison;
 run_cp6_warm(31,60,'V31_C60',.01,3,{'Am'});
 run_cp6_warm(31,60,'V31_C60_fine',.005,3,{'Am'});
catch ME
 disp(getReport(ME,'extended'));
end
exit;
