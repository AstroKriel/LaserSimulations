%% PREPARE WORKSPACE
clc, clear, close
format compact
format shortg

%% Compile C-Programs (if changes have been made):
% mex fast_sim_laser_prof.c
% mex fast_sim_laser_prof_noise.c
% mex fast_sim_laser_prpcf.c
% mex fast_sim_laser_prpcf_noise.c
% mex fast_sim_laser_prpcfuf_noise.c
% mex fast_sim_laser_pcf.c

%% CHOOSE SYSTEM
% DIM     = 5; BIF_EQN = 1;
% name_sys    = 'PCF';

% DIM     = 4; BIF_EQN = 2;
% % name_sys    = 'PROFN';
% name_sys    = 'PRPCFUFN';

DIM     = 7; BIF_EQN = 3;
name_sys    = 'PRPCFN';

%% EVALUATE and CALL SCRIPTS
evalParamSweep()
plotParamSweep()
plotTime()

%% END SCRIPT