%% START
format compact
clc, clear, close all

tic
%% COMPILE C-FILE
mex fast_sim_laser_SYS4.c
    % system 1: DIM = 4, Bif_num = 2
    % system 3: DIM = 7, Bif_num = 3
    % system 4: DIM = 5, Bif_num = 1
mex fast_detect_extrema.c

%% INITIALISE PARAMETERS
% System Parameters (that change)
DIM     = 5;        % <- CHANGES: for each system
bif_eqn = 1;        % <- CHANGES: for each system

% System Parameters
% PCF Paper Parameters Values
eta     = 0;
omega   = 0;
alpha   = 2;
ka      = 0.96;
beta    = (1-ka)/(2*ka);
T       = 1200;
P       = 0.6016;
theta   = 1143;
tau_R   = 20;

% Typical Parameter Values
% eta     = 0;                          % sweep [0, 0.06]
% omega   = 0;
% alpha   = 2;                          % in [2, 3, 5]
% ka      = 0.96;                       % in [0.96, 1]
% beta    = (1-ka)/(2*ka);
% T       = 1200;                       % in [1000, 500, 100]
% P       = 0.6016;                     % in [0.6, 1, 2]
% theta   = 1143;                       % in [1000, 2500, 5000, 7000]
% tau_R   = 20;                         % in {20, 50}

% Create Folder
Name = ['test_SYS4_tau=', num2str(tau_R)]; % <- CHANGES: system, and bif var
folder = [Name, '/'];
mkdir(folder);

% Parameters: Evaluation
params_laser = [eta, omega, alpha, beta, ka, T, P, theta, tau_R];
num_itter    = 100;     % <- CHANGES: resolution
param_num    = 1;       % <- (can) CHANGE: bifurcation parameter
param_vals   = linspace(0, 0.06, num_itter);
h            = 1;
horizon      = 0.2e6;
transient    = 1e6;
direction    = sign(param_vals(end) - param_vals(1));

% Bifurcation Parameters
num_extrema = 105;
bif_epsilon = 1e-3;
num_val     = num_extrema*num_itter;
bif_gamma   = zeros(num_val, 1);
bif_extrema = zeros(num_val, 1);
index       = 1;

sim_past = [];

%% Simulation
disp('Simulation Started.')
tic

for itter = 1:num_itter
    params_laser(param_num) = param_vals(itter);
    
    if ~isequal(transient, 0)
        sim_past = sim_laser_PCF(params_laser, [h transient], sim_past, DIM);
    end
    
    sim_past = sim_laser_PCF(params_laser, [h horizon], sim_past, DIM);

    bif_sim_past = abs(sim_past(1:end, bif_eqn));
    temp_output = detect_extrema(bif_sim_past, num_extrema, bif_epsilon);
    temp_len = length(temp_output);
    bif_gamma(index:index+temp_len-1) = param_vals(itter);
    bif_extrema(index:index+temp_len-1) = temp_output;
    index = index + temp_len;
end

%% Bifurcation Diagram
bif_gamma = bif_gamma(1:index-1);
bif_extrema = bif_extrema(1:index-1);

h2 = figure('Renderer', 'painters', 'Position', [10 10 900 500]);
plot(bif_gamma, bif_extrema, '.k', 'MarkerSize', 1)
grid on, axis tight
title(['Bifurcation Diagram: SYS1, \tau_R =', num2str(tau_R)], 'fontsize', 16)
xlabel('\eta', 'fontsize', 16)
ylabel('Output power (a. u.)', 'fontsize', 16)

saveas(h, [folder, 'bif_diag.png'], 'png');
save([folder, 'bif_extrema.mat'], 'EXTREMA');
save([folder, 'bif_gamma.mat'], 'GAMMA');

close all

%% END
disp('Simulation Completed')
toc
