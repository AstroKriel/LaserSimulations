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

    Sub-Function Dependencies:
    - fast_sim_laser_pcf
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
disp('Program Started: AnalysingParamSweep.')

%% Compile C-Programs (if changes have been made):
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
% DIM     = 5;
% BIF_EQN = 1;
% name    = 'PCF';

%% INITIALISE PARAMETERS
% Logical Variables
log_bif_data_save = true;
log_fft_save      = true;

% Sweeping Parameters
param_itter   = 100; % bif. resolution
param_elem    = 4;   % bif. w/respect to param_set
param_perturb = 1e-5;
param_start   = 0;
param_end     = 0.99;
param_vals    = linspace(param_start, param_end, param_itter) + param_perturb;

% Initialise System's Parameter Values
paper_name = 'Typ';
P       = 0.6;
T       = 1000;
theta   = 7000;
eta     = 0;
ka      = 0.96;
beta    = (1-ka)/(2*ka);
alpha   = 3;
omega   = 0;
tau_R   = 50;

% PCF
% paper_name = 'HO_ECM';
% P         = 0.6016;           % <-
% T         = 1200;             % <-
% theta     = 1143;             % <-
% eta       = 0;                % <-
% ka        = 0;
% beta      = 0;
% alpha     = 2;                % <-
% omega     = 0;
% tau_R     = 50;               % <-

% paper_name = 'M_ECM';
% P         = 1;                % <-
% T         = 500;              % <-
% theta     = 1000;             % <-
% eta       = 0;                % <-
% ka        = 0;
% beta      = 0;
% alpha     = 3;                % <-
% omega     = 0;
% tau_R     = 50;               % <-


% PROF
% paper_name = 'SOAPS';
% P         = 0.6;              % <-
% T         = 250;              % <-
% theta     = 7e3;              % <-
% eta       = 0;                % <-
% ka        = 0.96;             % <-
% beta      = (1-ka)/(2*ka);    % <-
% alpha     = 2;                % <-
% omega     = 0;                % <-
% tau_R     = 50;

% paper_name = 'R_SW_O_SL_1';
% P         = 0.5;              % <-
% T         = 150;              % <-
% theta     = 1e3;              % <-
% eta       = 0.125;            % <-
% ka        = 0.9;              % <-
% beta      = (1-ka)/(2*ka);    % <-
% alpha     = 2;                % <-
% omega     = -0.15;            % <-
% tau_R     = 50;

% paper_name = 'R_SW_O_SL_2';
% P         = 0.5;              % <-
% T         = 1500;             % <-
% theta     = 2e3;              % <-
% eta       = 0.4;              % <-
% ka        = 0.9;              % <-
% beta      = (1-ka)/(2*ka);    % <-
% alpha     = 2;                % <-
% omega     = -0.15;            % <-
% tau_R     = 50;

% paper_name = 'R_SW_O_SL_3';
% P         = 0.5;            % <-
% T         = 1500;           % <-
% theta     = 2e3;            % <-
% eta       = 0;              % <-
% ka        = 0.33;           % <-
% beta      = (1-ka)/(2*ka);  % <-
% alpha     = 2;              % <-
% omega     = -0.15;          % <-
% tau_R     = 50;

tau_P     = 1.4e-3;
R         = 1e-6;
param_set = [P, T, theta, eta, beta, ka, alpha, tau_R, omega, R];

% Initialise Analysis Parameters
h         = 1;
horizon   = 0.2e6;
transient = 1e6;
delay     = floor(theta/h);

% Initialise Bifurcation Parameters
num_extrema = 105;
bif_epsilon = 1e-3;
num_val     = num_extrema*param_itter;
bif_eta     = zeros(num_val, 1);
bif_extrema = zeros(num_val, 1);
bif_index   = 1;

% Initialise Frequency Analysis:
tau_p       = 1.4e-3;
L           = floor(horizon/h) + 1;
NFFT        = 2^nextpow2(L);
f           = 1/(2*h*tau_p)*linspace(0, 1/16, NFFT/32+1);
DISPLAY     = zeros(NFFT/32+1, param_itter);
SPECTRA     = zeros(NFFT,1);

% Initialise System's Past
sim_past = [];

% Create Folder to Save Data
mkdir(['Param_', paper_name])
folder = ['Param_', paper_name, '/SYS_', name];
mkdir(folder);

%% Simulation
for itter = 1:param_itter
    %% Integrate
    if (mod(itter, 5) == 0)
        disp(['Simulation: ', num2str(100*itter/param_itter), '% Complete'])
    end
    
    param_set(param_elem) = param_vals(itter);
    
    if ~isequal(transient, 0)
        sim_past = IntegSimLaser(param_set, [h transient], sim_past,...
            DIM, name, BIF_EQN);
    end
    
    sim_past = IntegSimLaser(param_set, [h horizon], sim_past, DIM,...
        name, BIF_EQN);
    
    if any(any(isnan(sim_past)))
        disp(['System ', name,...
            ' has NaN values at itter = ', num2str(itter)])
        break
    end
    % input: param_set, transient, h, sim_past, DIM, name, BIF_EQN
    % output: sim_past
    
    %% Bifurcation Analysis
    bif_sim_past = abs(sim_past(1:end, BIF_EQN));
    temp_output = detect_extrema(bif_sim_past, num_extrema, bif_epsilon);
    temp_len = length(temp_output);
    bif_eta(bif_index:bif_index+temp_len-1) = param_vals(itter);
    bif_extrema(bif_index:bif_index+temp_len-1) = temp_output;
    bif_index = bif_index + temp_len;
    % input: sim_past, BIF_EQN, num_extrema, bif_epsilon, bif_index
    % output: bif_eta, bif_extrema, bif_index
    
    %% Frequency Analysis
    SPECTRA = fft(sim_past(:,BIF_EQN).^2, NFFT)/L;
    DISPLAY(:,itter) = SPECTRA(1:NFFT/32+1);
    % input: sim_past, BIF_EQN, NFFT, L
    % output: SPECTRA, DISPLAY
end

%% Save Bifurcation
if (log_bif_data_save && ~any(any(isnan(sim_past))))
    bif_eta = bif_eta(1:bif_index-1);
    bif_extrema = bif_extrema(1:bif_index-1);
    
    csvwrite([folder, '/', 'bif_extrema.txt'], bif_extrema);
    csvwrite([folder, '/', 'bif_eta.txt'], bif_eta);
end
% output: bif_eta, bif_extrema

%% Save Frequency
if (log_fft_save)
    csvwrite([folder, '/', 'fft_display.txt'], DISPLAY);
end
% output: DISPLAY

%% End Script
disp('Program Finished: AnalysingParamSweep.')
toc
