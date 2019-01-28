%% Script Begining
%{
    Authors:    Neco Kriel  (2019)
                Chi Hak     (2017)

    Purpose:

    Function Dependencies:
    - IntegSimLaser
    - findIndices

    Sub-Function Dependencies:
    - fast_sim_laser_pcf
    - fast_sim_laser_prof
    - fast_sim_laser_prof_noise
    - fast_sim_laser_prpcf
    - fast_sim_laser_prpcf_noise

    TODO:
    - Impliment PCF
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
% DIM     = 5;
% BIF_EQN = 1;
% name    = 'PCF';

%% INITIALISE PARAMETERS
num_cycles = 2; % number of cycles of hor. and ver. polar.

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

% Sweeping Parameters
param_itter   = 100; % bif. resolution
param_elem    = 4;   % bif. w/respect to param_set
param_perturb = 1e-5;
param_start   = 0;
param_end     = 0.99;
param_vals    = linspace(param_start, param_end, param_itter) + param_perturb;

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

folder = ['Param_', paper_name, '/SYS_', name];
eta    = importdata([folder, '/bif_eta.txt']);
maxima = importdata([folder, '/bif_extrema.txt']);

% Initialise Analysis Parameters
h         = 1;
horizon   = 0.2e6;
transient = 1e6;
delay     = floor(theta/h);

% Initialise System's Past
sim_past = [];

%% Integrate, Analyse, and Plot
figure('Renderer', 'painters', 'Position', [10 10 1200 700]);
vidObj = VideoWriter([folder, '.avi']);
open(vidObj); % start video

for itter = 1:param_itter
    param_set(param_elem) = param_vals(itter);
    
    if ~isequal(transient, 0)
        sim_past = IntegSimLaser(param_set, [h transient], sim_past, DIM,...
            name, BIF_EQN);
    end
    sim_past = IntegSimLaser(param_set, [h horizon], sim_past, DIM,...
        name, BIF_EQN);
    
    if any(any(isnan(sim_past)))
        disp(['System ', name, ' has NaN values.'])
    else
        hor = (sim_past(:, 1)).^2;
        ver = (sim_past(:, BIF_EQN)).^2;
        
        if (any(abs(hor-ver) < 0.1))
            % calculate indices of one cycle of hor. and ver. polarizations
            [itter_start, itter_mid, itter_end] = findIndices(hor, ver);
            
            % normalise time wrt. delay
            x = (1:length(hor(itter_start:(itter_start + num_cycles*(itter_end - itter_start)))))./delay;
            
            % Plot Time
            subplot(3, 3, 4:9)
            plot(x, hor(itter_start:(itter_start + 2*(itter_end - itter_start))), 'b')
            hold on
            plot(x, ver(itter_start:(itter_start + 2*(itter_end - itter_start))), 'r')
            title(name)
            legend({'Horizontal', 'Vertical'}, 'FontSize', 11)
            ylabel('Power (10 log-dB)', 'fontsize',14)
        else            
            % normalise time wrt. delay
            x = (1:length(hor(5*delay:10*delay)))./delay;
            
            % Plot Time
            subplot(3, 3, 4:9)
            plot(x, hor(5*delay:10*delay), 'b')
            hold on
            plot(x, ver(5*delay:10*delay), 'r')
            title(name)
            legend({'Horizontal', 'Vertical'}, 'FontSize', 11)
            ylabel('Power (10 log-dB)', 'fontsize',14)
        end
        
        
        % Plot FFT
        FFT_array_hor = abs(fft(sim_past(:, 1).^2));
        FFT_array_hor = FFT_array_hor(1:floor(length(FFT_array_hor)/2));
        x_hor = (1:length(FFT_array_hor))/(2*h*tau_P*length(FFT_array_hor));
        subplot(3, 3, 2), hold on
        plot(x_hor, FFT_array_hor, 'b-', 'MarkerSize', 1)
        title('Frequency: Horizontal', 'FontSize', 16)
        ylabel('Magnitude')
        xlabel('Frequencies (GHz)')
        xlim([1, 50])
        
        FFT_array_ver = abs(fft(sim_past(:, BIF_EQN).^2));
        FFT_array_ver = FFT_array_ver(1:floor(length(FFT_array_ver)/2));
        x_ver = (1:length(FFT_array_ver))/(2*h*tau_P*length(FFT_array_ver));
        max_ver = max(FFT_array_ver);
        subplot(3, 3, 3), hold on
        plot(x_ver, FFT_array_ver, 'r-', 'MarkerSize', 1)
        title('Frequency: Vertical', 'FontSize', 16)
        ylabel('Magnitude')
        xlabel('Frequencies (GHz)')
        xlim([1, 50])
        
        % Plot Bifurcation
        subplot(3, 3, 1)
        plot(eta, maxima, '.k', 'MarkerSize', 5)
        hold on
        line([param_set(param_elem), param_set(param_elem)], [0, max(maxima)])
        title('Bifurcation Diagram', 'FontSize', 16)
        xlabel('\eta', 'FontSize', 16)
        ylabel('Power Output', 'FontSize', 16)
        grid on
        axis tight
        
        % Save Current Frame
        writeVideo(vidObj, getframe(gcf));
        clf
    end
end

close(vidObj); % end video
close all

%% End Script
disp('Program Finished: AnalyseTimeTraces.')
toc
