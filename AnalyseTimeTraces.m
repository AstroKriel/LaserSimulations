%% Script beginning
% Script created to integrate and analyse a single system (to test if
% saving issues could be fixed).
format compact
format short
clc, clear, close all

tic
disp('Script Started.')

%% Compile C code:
% mex fast_sim_laser_prof.c
% mex fast_sim_laser_prof_noise.c
DIM = 4;
BIF_EQN = 2;
% name = 'PROF';
name = 'PROFN';

% mex fast_sim_laser_prpcf.c
% mex fast_sim_laser_prpcf_noise.c
% DIM = 7;
% BIF_EQN = 3;
% name = 'PRPCF';
% name = 'PRPCFN';

% mex fast_sim_laser_pcf.c
% mex fast_sim_laser_pcf_noise.c
% DIM = 5;
% BIF_EQN = 1;
% name = 'PCF';
% name = 'PCFN';

%% INITIALISE PARAMETERS
% Initialise System's Parameter Values
P       = 0.6;
T       = 250;
theta   = 7e3;
eta     = 0.08;
ka      = 0.96;
beta    = (1-ka)/(2*ka);
alpha   = 2;
omega   = 0;
tau_R   = 0;
tau_P   = 1.4e-3;
R       = 1e-12; % {0, 1e-12}

% Initialise Analysis Parameters
h         = 1;
horizon   = 0.2e6;
transient = 1e6;
delay     = floor(theta/h);

%% Initialise Analysis:
param_set = [P, T, theta, eta, beta, ka, alpha, tau_R, omega, R];

% Initialise System's Past
sim_past = [];

%% Integrate
if ~isequal(transient, 0)
    sim_past = integ_sim_laser(param_set, [h transient], sim_past, DIM,...
        name, BIF_EQN);
end
sim_past = integ_sim_laser(param_set, [h horizon], sim_past, DIM,...
    name, BIF_EQN);

hor = (sim_past(:, 1)).^2;
ver = (sim_past(:, BIF_EQN)).^2;
[itter_start, itter_end] = findIndices(hor, ver);

x = (1:length(hor(itter_start:(itter_start + 2*(itter_end - itter_start)))))./delay;

% Plot E_1 and E_2
subplot(2, 1, 1), hold on
plot(x, hor(itter_start:(itter_start + 2*(itter_end - itter_start))))
title('Horizontal')
subplot(2, 1, 2), hold on
plot(x, ver(itter_start:(itter_start + 2*(itter_end - itter_start))))
title('Vertical')

toc
disp('Script finished.')

%% Define Function
function [itter_start, itter_end_2] = findIndices(hor, ver)
    points_found = false;
    itter_start = 1;
    while (itter_start < length(hor) && ~points_found)
        cur_hor = hor(itter_start);
        cur_ver = ver(itter_start);
        nex_hor = hor(itter_start + 1);
        nex_ver = ver(itter_start + 1);
        if (cur_hor <= cur_ver && nex_hor >= nex_ver && ~points_found)
            % hor grad > 0 & overlapping -> start period found
            itter_end_1 = itter_start + 1;
            while (itter_end_1 < length(hor) && ~points_found)
                cur_hor = hor(itter_end_1);
                cur_ver = ver(itter_end_1);
                nex_hor = hor(itter_end_1 + 1);
                nex_ver = ver(itter_end_1 + 1);
                if (cur_hor >= cur_ver && nex_hor <= nex_ver && ~points_found)
                    % hor grad < 0 & overlapping -> end period found
                    itter_end_2 = itter_end_1 + 1;
                    while (itter_end_2 < length(hor) && ~points_found)
                        cur_hor = hor(itter_end_2);
                        cur_ver = ver(itter_end_2);
                        nex_hor = hor(itter_end_2 + 1);
                        nex_ver = ver(itter_end_2 + 1);
                        if (cur_hor <= cur_ver && nex_hor >= nex_ver && ~points_found)
                            % hor grad > 0 & overlapping -> end period found
                            points_found = true;
                        end

                        % not overlapping
                        itter_end_2 = itter_end_2 + 1;
                    end
                end

                % not overlapping
                itter_end_1 = itter_end_1 + 1;
            end
        end

        % not overlapping
        itter_start = itter_start + 1;
    end
end
