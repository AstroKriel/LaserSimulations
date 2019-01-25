%% Script Begining
%{ 
    Authors:    Neco Kriel  (2019)
                Chi Hak     (2017)

    Purpose:    
    Analyse and save the bifurcation and frequency spectrum data for a 
    chosen system and sarameter space, over a chosen domain for a specified 
    parameter.

    Function Dependencies:
    - IntegSimLaser
    - findIndices

    Sub-Function Dependencies:
    - fast_sim_laser_pcf
    - fast_sim_laser_pcf_noise
    - fast_sim_laser_prof
    - fast_sim_laser_prof_noise
    - fast_sim_laser_prpcf
    - fast_sim_laser_prpcf_noise
    
    Special Actions:
    - creates a folder if the particular one does not already exist
    - save data (bif_extrema, bif_eta, fft_display)
%}

format compact
format short
clc, clear, close all

tic
disp('Program Started: AnalyseTimeTraces.')

%% Compile C code:
% mex fast_sim_laser_prof.c
% mex fast_sim_laser_prof_noise.c
DIM     = 4;
BIF_EQN = 2;
name    = 'PROF';
% name    = 'PROFN';

% mex fast_sim_laser_prpcf.c
% mex fast_sim_laser_prpcf_noise.c
% DIM     = 7;
% BIF_EQN = 3;
% name    = 'PRPCF';
% name    = 'PRPCFN';

% mex fast_sim_laser_pcf.c
% mex fast_sim_laser_pcf_noise.c
% DIM     = 5;
% BIF_EQN = 1;
% name    = 'PCF';
% name    = 'PCFN';

%% INITIALISE PARAMETERS
num_cycles = 2; % number of cycles of hor. and ver. polar.

% Initialise System's Parameter Values
P         = 0.6;
T         = 250;
theta     = 7e3;
eta       = 0.08;
ka        = 0.96;
beta      = (1-ka)/(2*ka);
alpha     = 2;
omega     = 0;
tau_R     = 0;
tau_P     = 1.4e-3;
R         = 0;
param_set = [P, T, theta, eta, beta, ka, alpha, tau_R, omega, R];

% Initialise Analysis Parameters
h         = 1;
horizon   = 0.2e6;
transient = 1e6;
delay     = floor(theta/h);

% Initialise System's Past
sim_past = [];

%% Integrate and Analyse
if ~isequal(transient, 0)
    sim_past = IntegSimLaser(param_set, [h transient], sim_past, DIM,...
        name, BIF_EQN);
end
sim_past = IntegSimLaser(param_set, [h horizon], sim_past, DIM,...
    name, BIF_EQN);

% calculate indices of one cycle of hor. and ver. polarizations
hor = (sim_past(:, 1)).^2;
ver = (sim_past(:, BIF_EQN)).^2;
[itter_start, itter_end] = findIndices(hor, ver);

% normalise time wrt. delay
x = (1:length(hor(itter_start:(itter_start + num_cycles*(itter_end - itter_start)))))./delay;

%% Plot 
figure('Renderer', 'painters', 'Position', [50 50 800 500])
subplot(2, 1, 1), hold on
plot(x, hor(itter_start:(itter_start + 2*(itter_end - itter_start))))
title('Horizontal')

subplot(2, 1, 2), hold on
plot(x, ver(itter_start:(itter_start + 2*(itter_end - itter_start))))
title('Vertical')

%% End Script
disp('Program Finished: AnalyseTimeTraces.')
toc
