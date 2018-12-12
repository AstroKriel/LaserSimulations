%% START MATLAB
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

% Logical Variables
log_plot_sim        = true;
log_fig_sim_save    = true;
log_fig_bif_save    = true;
log_fig_close       = true;
log_data_save       = true;

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
Name = ['test_SYS4_tau=', num2str(tau_R)]; % <- CHANGES: system and bif var
folder = [Name, '/'];
mkdir(folder);

% Parameters: Evaluation
param_list      = [eta, omega, alpha, beta, ka, T, P, theta, tau_R];
% Analysis Parameters
test_num        = [0.0006, 0.001212, 0.0067, 0.014, 0.017, 0.03, 0.052,...
    0.055, 0.058];      % analysis points
tol             = 1e-2;
num_itter       = 1000;  % <- CHANGES: resolution
param_num       = 1;    % <- CHANGEs: bif parameter
param_vals      = linspace(0, 0.15, num_itter);
h               = 1;
horizon         = 0.2e6;
transient       = 1e6;
direction       = sign(param_vals(end) - param_vals(1));
delay           = floor(theta/h);

% Bifurcation Parameters
num_extrema     = 105;
bif_epsilon     = 1e-3;
num_val         = num_extrema*num_itter;
bif_gamma       = zeros(num_val, 1);
bif_extrema     = zeros(num_val, 1);
index           = 1;

sim_past = [];

%% Simulation
disp('Simulation Started.')
tic

for itter = 1:num_itter
    %% Integrate
    param_list(param_num) = param_vals(itter);
    
    if ~isequal(transient, 0)
        sim_past = sim_laser_PCF(param_list, [h transient], sim_past, DIM);
    end
    
    sim_past = sim_laser_PCF(param_list, [h horizon], sim_past, DIM);
    
    if any(any(isnan(sim_past)))
        disp('System has NaN values')
        break
    end
    
    %% Bifurcation Analysis
    bif_sim_past = abs(sim_past(1:end, bif_eqn));
    temp_output = detect_extrema(bif_sim_past, num_extrema, bif_epsilon);
    temp_len = length(temp_output);
    bif_gamma(index:index+temp_len-1) = param_vals(itter);
    bif_extrema(index:index+temp_len-1) = temp_output;
    index = index + temp_len;
    
    %% Plot
    if (log_plot_sim && (any((1-tol)*test_num <= param_list(param_num) & param_list(param_num) <= (1+tol)*test_num)))
        % Plot: Time Traces
        h1 = figure('Renderer', 'painters', 'Position', [10 10 900 600]);
        subplot(2, 2, 1)
        plot(sim_past(:, 1).^2, 'b')
        ylabel('Power (10 log-dB)', 'fontsize',14)
        title('Horizontal Polarisation', 'fontsize',18)
        axis([15*delay, 20*delay, 0, 1])
        
        subplot(2, 2, 3)
        plot(sim_past(:, bif_eqn).^2, 'r')
        ylabel('Power (10 log-dB)', 'fontsize',14)
        title('Vertical Polarisation', 'fontsize',18)
        axis([15*delay, 20*delay, 0, 1])
        
        % Plot radio frequency/RF spectrum FFT in dB
        subplot(2, 2, [2, 4])
        plot(log(abs(fft(sim_past(:, bif_eqn).^2))))
        xlim([5, 6e3])
        title(['RF spectrum FFT in dB: eta=',...
            num2str(param_list(param_num))], 'fontsize',18)
        
        % Save: Fig Time Trace
        if log_fig_sim_save
            saveas(h1,[folder, 'Plot_', num2str(itter),'.png'],'png');
        end
        
        % Save: Time Trace Data
        if log_data_save
            csvwrite([folder, 'Data_', num2str(itter),'_hor_polar.txt'],...
                sim_past(:, 1).^2)
            csvwrite([folder, 'Data_', num2str(itter),'_ver_polar.txt'],...
                sim_past(:, bif_eqn).^2)
            csvwrite([folder, 'Data_', num2str(itter),'_FFT.txt'],...
                log(abs(fft(sim_past(:, bif_eqn).^2))))
        end
        
        % Close figs in large parameter sweeps
        if log_fig_close
            close all
        end
    end
end

%% Bifurcation Diagram
bif_gamma = bif_gamma(1:index-1);
bif_extrema = bif_extrema(1:index-1);

h2 = figure('Renderer', 'painters', 'Position', [10 10 900 500]);
plot(bif_gamma, bif_extrema, '.k', 'MarkerSize', 1)
grid on, axis tight
title(['Bifurcation Diagram: SYS1 ', Name], 'fontsize', 16)
xlabel('\eta', 'fontsize', 16)
ylabel('Output power (a. u.)', 'fontsize', 16)

% Save: Bifurcation Plot
if log_fig_bif_save
    saveas(h, [folder, 'bif_diag.png'], 'png');
end

% Save: Bifurcation Data
if log_data_save
    csvwrite([folder, 'bif_extrema.txt'], bif_extrema);
    csvwrite([folder, 'bif_gamma.txt'], bif_gamma);
end

if log_fig_close
    close all
end

%% END MATLAB
disp('Simulation Completed')
toc

