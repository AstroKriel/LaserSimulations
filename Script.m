%% Script beginning
% Script created to integrate and analyse a single system (to test if
% saving issues could be fixed).
format compact
format short
clc, clear, close all

tic
disp('Script Started.')

%% Compile C-Programs (if changes have been made):
% mex fast_detect_extrema.c

% mex fast_sim_laser_prof.c
% mex fast_sim_laser_prof_noise.c
DIM = 4; 
BIF_EQN = 2;  
name = 'PROF';
% name = 'PROFN';

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

%% Variables that change
tau_R   = 0;
theta   = 7e3;
% bifurcation Domain
bif_start   = 0;
bif_end     = 0.2;

%% INITIALISE PARAMETERS
% Sweeping Parameters:
num_itter   = 10;  % bif. resolution
% Tuning Parameters
bif_elem    = 4;    % bif. w/respect to param_set
bif_perturb = 1e-5;

% Logical Variables
log_bif_data_save = true;
log_fft_save      = true;

% Initialise System's Parameter Values
P       = 0.6;
T       = 250;
% theta   = 7e3;
eta     = 0.08;
ka      = 0.96;
beta    = (1-ka)/(2*ka);
alpha   = 2;
omega   = 0;
% tau_R   = 0;
tau_P   = 1.4e-3;
R       = 0;

% Initialise Analysis Parameters
h         = 1;
horizon   = 0.2e6;
transient = 1e6;
delay     = floor(theta/h);

%% Initialise Analysis:
param_vals = linspace(bif_start, bif_end, num_itter) + bif_perturb;
param_set  = [P, T, theta, eta, beta, ka, alpha, tau_R, omega, R];

% Create Folder to Save Data
folder = ['SYS_', name, '_tau=', num2str(tau_R), '_theta=', num2str(theta)];
mkdir(folder);

% Initialise Bifurcation Parameters
num_extrema     = 105;
bif_epsilon     = 1e-3;
num_val         = num_extrema*num_itter;
bif_eta         = zeros(num_val, 1);
bif_extrema     = zeros(num_val, 1);
bif_index       = 1;

% Initialise System's Past
sim_past = [];

% Initialise Frequency Analysis:
tau_p       = 1.4e-3;
L           = floor(horizon/h) + 1;
NFFT        = 2^nextpow2(L);
f           = 1/(2*h*tau_p)*linspace(0, 1/16, NFFT/32+1);
DISPLAY     = zeros(NFFT/32+1, num_itter);
SPECTRA     = zeros(NFFT,1);

%% Simulation
for itter = 1:num_itter
    %% Show Progress
    if (mod(itter, 1) == 0)
        disp(itter)
    end
    
    %% Integrate
    param_set(bif_elem) = param_vals(itter);
    
    if ~isequal(transient, 0)
        sim_past = integ_sim_laser(param_set, [h transient], sim_past,...
            DIM, name, BIF_EQN);
    end
    
    sim_past = integ_sim_laser(param_set, [h horizon], sim_past, DIM,...
        name, BIF_EQN);
    
    if any(any(isnan(sim_past)))
        disp(['System ', name,...
            ' has NaN values at itter = ', num2str(itter)])
        break
    end
    % output: sim_past
    
    %% Bifurcation Analysis
    bif_sim_past = abs(sim_past(1:end, BIF_EQN));
    temp_output = detect_extrema(bif_sim_past, num_extrema, bif_epsilon);
    temp_len = length(temp_output);
    bif_eta(bif_index:bif_index+temp_len-1) = param_vals(itter);
    bif_extrema(bif_index:bif_index+temp_len-1) = temp_output;
    bif_index = bif_index + temp_len;
    % ouput: bif_eta, bif_extrema, bif_index
    
    %% Frequency Analysis
    SPECTRA = fft(sim_past(:,BIF_EQN).^2, NFFT)/L;
    DISPLAY(:,itter) = SPECTRA(1:NFFT/32+1);
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

toc
disp('Script finished.')

%% Spript end
