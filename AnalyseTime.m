%% START SCRIPT
format compact
format shortg
clc, clear, close all

tic
disp('Program Starting')

%% Compile C code:
% Uncomment System of Interest
mex fast_detect_extrema.c

% mex fast_sim_laser_prof.c
% DIM = 4; BIF_EQN = 2;  name = 'PROF';
% range = [0.011051, 0.03026, 0.042105, 0.10132]; % PROF

mex fast_sim_laser_prpcf.c
DIM = 7; BIF_EQN = 3; name = 'PRPCF';
% range = [0.019005, 0.028534, 0.031461, 0.057124, 0.0599]; % PRPCF, tau = 50
range = [0.013377, 0.031686, 0.054872, 0.062451]; % PRPCF, tau = 1

% mex fast_sim_laser_pcf.c
% DIM = 5; BIF_EQN = 1; name = 'PCF';

%% INITIALISE PARAMETERS
% Sweeping Parameters:
num_itter       = 100;  % resolution
tau_R           = 1;    % in {20, 50}
param_set_num   = 2;

% Tuning Parameters
bif_elem        = 1;    % bif. w/respect to param.
bif_start       = 1e-5; % bif. param. starting val.

% System Parameter Sets:
omega   = 0;
alpha   = 3;
ka      = 0.96;
T       = 1000;
P       = 0.6;
theta   = 7000;
eta     = 0;            % initialise, but sweep parameter
beta    = (1-ka)/(2*ka);

% Integrating Forward
param_vals  = linspace(0, 0.15, num_itter) + bif_start;
param_set   = [eta, omega, alpha, beta, ka, T, P, theta, tau_R];

% Analysis Parameters
h           = 1;
horizon     = 0.2e6;
transient   = 1e6;
delay       = floor(theta/h);
sim_past    = [];

% Initialise Bifurcation Parameters
num_extrema     = 105;
bif_epsilon     = 1e-3;
num_val         = num_extrema*num_itter;
bif_eta         = zeros(num_val, 1);
bif_extrema     = zeros(num_val, 1);
bif_index       = 1;

%% Bifurcation
filename = ['SYS_',name,'_Param=',num2str(param_set_num),'_tau=',num2str(tau_R),'_F'];

x = importdata([filename, '/bif_eta.txt']);
x = x(1:20:end);
y = importdata([filename, '/bif_extrema.txt']);
y = y(1:20:end);

%% Simulation
figure('Renderer', 'painters', 'Position', [10 10 1200 700]);
vidObj = VideoWriter(['SYS_', name, '_Param=', num2str(param_set_num),...
    '_tau=', num2str(tau_R), '.avi']);
open(vidObj); % start video

for itter = 1:num_itter
    param_set(bif_elem) = param_vals(itter);
    
    % Integrate
    if ~isequal(transient, 0)
        sim_past = integ_sim_laser(param_set, [h transient], sim_past, DIM);
    end
    
    sim_past = integ_sim_laser(param_set, [h horizon], sim_past, DIM);
    
    if any(any(isnan(sim_past)))
        disp(['System ', name,...
            ' has NaN values at itter = ', num2str(itter)])
        break
    end
    
    % Plot Time Series
    subplot(3, 3, 4:9), hold on
    plot(sim_past(:, 1).^2, 'b.', 'MarkerSize', 5)
    plot(sim_past(:, BIF_EQN).^2, 'r.', 'MarkerSize', 1)
    
    title([name, ' Time Series: eta = ', num2str(param_set(bif_elem))], 'fontsize',18)
    legend({'Horizontal', 'Vertical'}, 'FontSize', 11)
    ylabel('Power (10 log-dB)', 'fontsize',14)
    xlim([15*delay, 20*delay])
    
    % Plot Bifurcation
    subplot(3, 3, 1), hold on
    plot(x(x > range(end)), y(x > range(end)), '.g', 'MarkerSize', 5)
    for i = 0:length(range)-1
        if (mod(i, 2) == 0)
            plot(x(x < range(end-i)), y(x < range(end-i)), '.k', 'MarkerSize', 5)
        else
            plot(x(x < range(end-i)), y(x < range(end-i)), '.g', 'MarkerSize', 5)
        end
    end
    line([param_set(bif_elem), param_set(bif_elem)], [0, 2])
    title('Bifurcation Diagram', 'FontSize', 16)
    xlabel('\eta', 'FontSize', 16)
    ylabel('Power Output', 'FontSize', 16)
    grid on
    xlim([param_vals(1), param_vals(end)])
    
    % Plot FFT
    FFT_array = abs(fft(sim_past(:, BIF_EQN).^2));
    FFT_array = FFT_array(1:floor(length(FFT_array)/2));
    subplot(3, 3, 3), hold on
    plot(1:length(FFT_array), FFT_array, 'r.', 'MarkerSize', 1)
    title('Frequency', 'FontSize', 16)
    ylabel('Magnitude')
    xlim([5, 6e3])
    
    % Save Current Frame
    writeVideo(vidObj, getframe(gcf));
    clf
end

% Save Video
close(vidObj); % end video
close all

disp('Program Finished')
toc

%% END SCRIPT