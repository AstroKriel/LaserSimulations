%% Script beginning
% Script created to integrate and analyse a single system (to test if
% saving issues could be fixed).
format compact
format short
clc, clear, close all

tic
disp('Script Started.')

%% Compile C code:
% Uncomment System of Interest
mex fast_detect_extrema.c

% Phase-Conjugate Feedback
% mex fast_sim_laser_pcf.c
% DIM = 5; BIF_EQN = 1; name = 'PCF';

% Polarization-Rotated Optical Feedback
mex fast_sim_laser_prof.c
DIM = 4; BIF_EQN = 2;  name = 'PROF';

% % Polarization-Rotated Phase-Conjugate Feedback
% mex fast_sim_laser_prpcf.c
% DIM = 7; BIF_EQN = 3; name = 'PRPCF';

tau_R = 50;

% bif_length      = 'Long';
% bif_start_1     = 0;
% bif_end         = 0.95;

bif_length      = 'Shorter';
bif_start_1     = 0.06;
bif_end         = 0.07;

%% INITIALISE PARAMETERS
% Sweeping Parameters:
num_itter       = 100;  % resolution
param_set_num   = 2;    % in {1, 2}
% Tuning Parameters
bif_elem        = 1;    % bif. w/respect to param.
dir             = '_F';
bif_start_2     = 1e-5;


% Logical Variables
% Bifurcation Vars:
log_bif_data_save   = true;
% Frequency Analysis Vars:
log_fft_save        = true;

%% System Parameter Sets:
switch param_set_num
    case 1
        % Values as seen in paper
        omega   = 0;
        alpha   = 2;
        ka      = 0.96;
        T       = 1200;
        P       = 0.6016;
        theta   = 1143;
    case 2
        % Typical Parameter Values
        omega   = 0;
        alpha   = 3;        % in [2, 3, 5]
        ka      = 0.96;     % in [0.96, 1]
        T       = 1000;     % in [1000, 500, 100]
        P       = 0.6;      % in [0.6, 1, 2]
        theta   = 7000;     % in [1000, 2500, 5000, 7000]
end

% Integrating Forward
param_vals = linspace(bif_start_1, bif_end, num_itter) + bif_start_2;

% Create Folder to Save Data
if log_bif_data_save || log_fft_save
    folder = ['SYS_', name, '_Param=', num2str(param_set_num),...
        '_tau=', num2str(tau_R), dir];
    mkdir(folder);
end

eta       = 0; % initialise, but sweep parameter
tau_P     = 1.4e-3;
beta      = (1-ka)/(2*ka);
param_set = [eta, omega, alpha, beta, ka, T, P, theta, tau_R];

% Analysis Parameters
h         = 1;
horizon   = 0.2e6;
transient = 1e6;
delay     = floor(theta/h);

% Initialise Bifurcation Parameters
num_extrema     = 105;
bif_epsilon     = 1e-3;
num_val         = num_extrema*num_itter;
bif_eta         = zeros(num_val, 1);
bif_extrema     = zeros(num_val, 1);
bif_index       = 1;

% Initialise Frequency Analysis Parameters
Chaos_bw            = zeros(num_itter, 1);
Chaus_eff_bw        = zeros(num_itter, 1);
Freq_ext_cav_node   = zeros(num_itter, 1);

sim_past = [];

% Frequency Analysis:
tau_p       = 1.4e-3;
L           = floor(horizon/h) + 1;
NFFT        = 2^nextpow2(L);
f           = 1/(2*h*tau_p)*linspace(0, 1/16, NFFT/32+1);
DISPLAY     = zeros(NFFT/32+1, num_itter);
SPECTRA     = zeros(NFFT,1);

%% Simulation
for itter = 1:num_itter
    %% Show Progress
    if (mod(itter, 10) == 0)
        disp(itter)
    end
    
    %% Integrate
    param_set(bif_elem) = param_vals(itter);
    
    if ~isequal(transient, 0)
        sim_past = integ_sim_laser(param_set, [h transient], sim_past, DIM);
    end
    
    sim_past = integ_sim_laser(param_set, [h horizon], sim_past, DIM);
    
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
    
    csvwrite([folder, '/', 'bif_extrema_', bif_length,'.txt'], bif_extrema);
    csvwrite([folder, '/', 'bif_eta_', bif_length,'.txt'], bif_eta);
end
% output: bif_eta, bif_extrema

%% Save Frequency
if (log_fft_save)
    csvwrite([folder, '/', 'fft_display_', bif_length, '.txt'], DISPLAY);
end

toc
disp('Script finished.')

%% Spript end
