%% Script Begining
%{
    Authors:    Neco Kriel  (2019)

    Purpose:
    The purpose of this script is to investigate a laser as the bifurcation
    parameter is swept. As the system is integrated, for every point in the
    bifurcation parameter's domain, the bifurcation structure of the system
    is analysed, frequency analyses is performed, and the length of the
    period for the vertical oscilations are performed. At the end of the
    script, these variables are stored in a folder (folder) created unique
    to the parameters chosen.

    Function Dependencies:
    - IntegSimLaser
    - findIndices

    Sub-Function Dependencies:
    - fast_sim_laser_pcf.c
    - fast_sim_laser_prof.c
    - fast_sim_laser_prof_noise.c
    - fast_sim_laser_prpcf.c
    - fast_sim_laser_prpcf_noise.c
    
    Special Actions:
    - creates a folder if the particular one does not already exist
    - save data:
        > ParamSet
        > bif_extrema
        > bif_eta
        > fft_display
        > max_freq
        > period_length
%}

%% PREPARE WORKSPACE
format compact
format short
clc, clear, close all

tic
disp('Started: ParamSweep')

%% Compile C-Programs (if changes have been made):
% mex fast_sim_laser_prof.c
% mex fast_sim_laser_prof_noise.c
% mex fast_sim_laser_prpcf.c
% mex fast_sim_laser_prpcf_noise.c
% mex fast_sim_laser_pcf.c

%% CHOOSE SYSTEM
% DIM     = 5; BIF_EQN = 1;
% name_sys    = 'PCF';
% name_sys    = 'PCFN';

DIM     = 4; BIF_EQN = 2;
% name_sys    = 'PROF';
name_sys    = 'PROFN';

% DIM     = 7; BIF_EQN = 3;
% % name_sys    = 'PRPCF';
% name_sys    = 'PRPCFN';

%% INITIALISE PARAMETERS
% Param Sets:
% PRPCFN and PROFN
% k = 0.96, a = 2, theta = 7000
% PRPCFN
% 1 . k = 0.96, a = {[1, 3], 5}, theta = 7000
% 2 . k = 0.96, a = 2, theta = {[800, 1e4], 2e4}
% 3 . k = [0.8, 1], a = 2, theta = 7000
% PRPCFN and PROFN
% 4 . k = 0.96, a = 3, theta = {[800, 1e4], 2e4}

% Sweeping Parameters
param_iter   = 100;   % bif. resolution
param_elem    = 4;    % bif. w/respect to param_set
param_perturb = 1e-5; % perturb ini. val. (for when param_start = 0)
param_start   = 0;    % first param. val.
param_end     = 0.25; % final param. val.
param_vals    = linspace(param_start, param_end, param_iter) +...
    param_perturb; % create array of param. vals.

% Initialise System Parameters
for theta = [800, 1e3:1e3:1e4, 2e4]
    param_set_name = 'SOAPS_4'; % name_sys of parameter set (file-storage purposes)
    P         = 0.6;    % pump parameter (above threshold)
    T         = 250;    % ratio of carrier to cavity lifetime
%     theta     = 7000;   % normalised (time) delay
    eta       = 0;      % feedback rate
    ka        = 0.96;   % gain coefficient ratio between the TM and TE modes
    beta      = (1-ka)/(2*ka); % TM mode additional losses,
    alpha     = 3;      % linewidth enhancement factor
    omega     = 0;      %
    tau_R     = 50;     %
    tau_P     = 1.4e-3; % photon lifetime (scale so freq-spectrum is in GHz)
    R         = 1e-12;  % variance of white Gaussian noise
    param_set = [P, T, theta, eta, beta, ka, alpha, tau_R, omega, R]; % store params
    
    % Initialise Analysis Parameters
    h         = 1;
    horizon   = 0.2e6;
    transient = 1e6;
    delay     = floor(theta/h);
    
    % Initialise Bifurcation Parameters
    num_extrema = 105;
    bif_epsilon = 1e-3;
    num_val     = num_extrema*param_iter;
    bif_eta     = zeros(num_val, 1);
    bif_extrema = zeros(num_val, 1);
    bif_index   = 1;
    
    % Initialise Frequency Analysis:
    L           = floor(horizon/h) + 1;
    NFFT        = 2^nextpow2(L);
    DISPLAY     = zeros(NFFT/32+1, param_iter);
    SPECTRA     = zeros(NFFT,1);
    
    % Initialise Data Variables
    max_freq = zeros(param_iter, 1);
    period_length = zeros(param_iter, 3);
    
    % Initialise System's Past
    sim_past = [];
    
    % Create Folder (to save data)
    mkdir(['Param_', param_set_name])
    folder = ['Param_', param_set_name, '/SYS_', name_sys, '_theta=', num2str(theta)];
%     folder = ['Param_', param_set_name, '/SYS_', name_sys, '_alpha=', num2str(alpha)];
%     folder = ['Param_', param_set_name, '/SYS_', name_sys, '_kappa=', num2str(ka)];
    mkdir(folder);
    
    %% SIMULATION
    for iter = 1:param_iter
        %% INTEGRATE
        if (mod(iter, 5) == 0)
            disp(['Simulation: ', num2str(100*iter/param_iter), '% Complete'])
        end
        param_set(param_elem) = param_vals(iter);
        if ~isequal(transient, 0)
            sim_past = IntegSimLaser(param_set, [h transient], sim_past,...
                DIM, name_sys, BIF_EQN);
        end
        sim_past = IntegSimLaser(param_set, [h horizon], sim_past, DIM,...
            name_sys, BIF_EQN);
        if any(any(isnan(sim_past)))
            disp(['System ', name_sys,...
                ' has NaN values at itter = ', num2str(iter)])
            break
        end
        
        %% ANALYSE INTEGRAL
        % Analyse Bifurcation Structure
        bif_sim_past = abs(sim_past(1:end, BIF_EQN));
        temp_output = detect_extrema(bif_sim_past, num_extrema, bif_epsilon);
        temp_len = length(temp_output);
        bif_eta(bif_index:bif_index+temp_len-1) = param_vals(iter);
        bif_extrema(bif_index:bif_index+temp_len-1) = temp_output;
        bif_index = bif_index + temp_len;
        
        % Analyse Power Spectrum for Vertical Polarisation
        SPECTRA = fft(sim_past(:, BIF_EQN).^2, NFFT)/L;
        DISPLAY(:, iter) = SPECTRA(1:NFFT/32+1);
        
        % Analyse Maximum Frequency of Vertical Polarisation
        FFT_array_ver = abs(fft(sim_past(:, BIF_EQN).^2));
        FFT_array_ver = FFT_array_ver(1:floor(length(FFT_array_ver)/2));
        FFT_array_ver(1:5) = NaN;
        [~, val_itter] = max(FFT_array_ver);
        x_ver = (1:length(FFT_array_ver))/(2*h*tau_P*length(FFT_array_ver));
        max_freq(iter) = x_ver(val_itter);
        
        % Analyse Period of Polarisation Oscilations
        hor = (sim_past(:, 1)).^2;
        ver = (sim_past(:, BIF_EQN)).^2;
        if (any(abs(hor-ver) < 0.1))
            % calculate indices of one cycle of hor. and ver. polarizations
            [itter_start, itter_mid, itter_end] = findIndices(hor, ver);
            period_length(iter, :) = [(itter_mid - itter_start),...
                (itter_end - itter_mid), (itter_end - itter_mid - theta)];
        else
            % no period detected, so set to default
            period_length(iter, :) = [-1, -1, -1];
        end
    end
    
    %% SAVE DATA
    % Save Copy of Parameter Set
    save(['Param_', param_set_name, '/', 'ParamSet.mat'], 'P', 'T', 'theta', 'beta', 'ka',...
        'alpha', 'tau_R', 'omega', 'R', 'h', 'tau_P', 'horizon', 'param_start', 'param_end')
    
    % Save Bifurcation
    bif_eta = bif_eta(1:bif_index-1);
    bif_extrema = bif_extrema(1:bif_index-1);
    csvwrite([folder, '/', 'bif_extrema.txt'], bif_extrema); % y-var (maxima)
    csvwrite([folder, '/', 'bif_eta.txt'], bif_eta); % x-var (eta)
    
    % Save Frequency
    csvwrite([folder, '/', 'fft_display.txt'], DISPLAY); % frequency vals
    csvwrite([folder, '/', 'max_freq.txt'], max_freq); % maximum freq per eta
    
    % Save Periods: length of vert. polar. osc. period
    csvwrite([folder, '/', 'period_length.txt'], period_length);
end

%% End Script
disp('Finished: ParamSweep')
disp(['System: ', name_sys]), disp(['Params: ', param_set_name])
toc
