%% Script Begining
%{
    Authors:    Neco Kriel  (2019)

    Purpose:
    The purpose of this script is to create and save plots for the data
    created by ParamSweep.m.

    Dependencies:
    - The data and folders exist
%}

%% PREPARE WORKSPACE
format compact
format short
clc, clear, close all

tic
disp('Started: PlotParamSweep')

%% INITIALISE AND DEFINE PARAMETERS
% System and Parameter set name
name_sys       = 'PROFN';   % system to integrate
% name_sys       = 'PRPCFN'; 
name_param_set = 'SOAPS_4'; % foldername to save data/figures

load(['Param_', name_param_set, '/ParamSet.mat']);

% FFT variables
L    = floor(horizon/h) + 1;
NFFT = 2^nextpow2(L);
f    = 1/(2*h*tau_P)*linspace(0, 1/16, NFFT/32+1);

for theta = [800, 1e3:1e3:1e4, 2e4]
    % file save path
    filename = ['Param_', name_param_set, '/SYS_', name_sys, '_theta=',...
        num2str(theta)];
%     filename = ['Param_', name_param_set, '/SYS_', name_sys, '_alpha=',...
%         num2str(alpha)];
%     filename = ['Param_', name_param_set, '/SYS_', name_sys, '_kappa=',...
%         num2str(ka)];
    
    % Import data
    fft = importdata([filename, '/', 'fft_display', '.txt']);
    xF  = importdata([filename, '/', 'bif_eta', '.txt']);
    yF  = importdata([filename, '/', 'bif_extrema', '.txt']);
    % Truncate data
    log = (xF <= param_end & xF >= param_start);
    xF  = xF(log);
    yF  = yF(log);
    
    %% PLOTTING
    % Plot Bifurcation Diagram and Power Spectrum
    param_res   = 1e3; % frequency resolution
    param_vals  = linspace(param_start, param_end, param_res);
    h1 = figure('Renderer', 'painters', 'Position', [10 10 1200 700]);
    subplot(2, 1, 1)
    imagesc(param_vals, f, log10(2*abs(fft)));
    colormap(hot)
    caxis([-4 0])
    cb = colorbar;
    set(cb,'Position',[0.92 0.6 .01 0.3])
    set(gca,'YDir','normal')
    axis tight
    title('Power Spectrum', 'FontSize', 16)
    subplot(2, 1, 2)
    plot(xF, yF, 'k.', 'MarkerSize', 3)
    axis tight
    title(['Bifurcation Diagram: ',...
                             name_sys,...
                           ': \{\theta=', num2str(theta),...
                             ', \kappa=', num2str(ka),...
                             ', \alpha=', num2str(alpha),...
                             ', \tau_R=', num2str(tau_R),...
                                  ', R=', num2str(R),' \}'],...
                             'FontSize', 16)
    saveas(h1, ['Param_', name_param_set, '/', 'SYS_', name_sys,...
        '_theta=', num2str(theta), '.png'], 'png');
%     saveas(h1, ['Param_', name_param_set, '/', 'SYS_', name_sys,...
%         '_alpha=', num2str(alpha), '.png'], 'png');
%     saveas(h1, ['Param_', name_param_set, '/', 'SYS_', name_sys,...
%         '_kappa=', num2str(ka), '.png'], 'png');
    close all
    
    %     % prepare bifurcation values for FFT and Period plot
    %     param_res     = 100;
    %     param_perturb = 1e-5;
    %     param_vals    = linspace(param_start, param_end, param_res) + param_perturb;
    %
    %     % Plot Maximum Frequency
    %     max_freq = 0;
    %     h1 = figure('Renderer', 'painters', 'Position', [10 10 1200 700]);
    %     hold on
    %     filename = ['Param_', name_param_set, '/SYS_', name_sys, '_theta=',...
    %         num2str(theta)];
    %     freq = importdata([filename, '/', 'max_freq', '.txt']);
    %     freq = freq*10;
    %     plot(param_vals, freq, 'o', 'MarkerSize', 3)
    %     max_freq = max([max_freq, freq']);
    %     title(name_sys)
    %     ylim([0, 40])
    %     xlabel('\eta')
    %     ylabel('max frequency (GHz)')
    %     % legend(['\theta = ', num2str(100)],...
    %     %     ['\theta = ', num2str(500)],...
    %     %     ['\theta = ', num2str(1000)],...
    %     %     ['\theta = ', num2str(3000)],...
    %     %     ['\theta = ', num2str(5000)],...
    %     %     ['\theta = ', num2str(7000)], 'Location','southeast')
    %     saveas(h1, ['Param_', name_param_set, '/', 'SYS_', name_sys,...
    %         '_maxfreq'], 'png');
    %     close all
    %
    %     % Plot Length of Vertical Polarisation Period
    %     h1 = figure('Renderer', 'painters', 'Position', [10 10 1200 700]);
    %     hold on
    %     filename = ['Param_', name_param_set, '/SYS_', name_sys, '_theta=',...
    %         num2str(theta)];
    %     time = importdata([filename, '/', 'period_length', '.txt']);
    %     time = time(:, 2);
    %     plot(param_vals, time, 'o', 'MarkerSize', 3)
    %     title(name_sys)
    %     xlabel('\eta')
    %     ylabel('period')
    %     % legend(['\theta = ', num2str(100)],...
    %     %     ['\theta = ', num2str(500)],...
    %     %     ['\theta = ', num2str(1000)],...
    %     %     ['\theta = ', num2str(3000)],...
    %     %     ['\theta = ', num2str(5000)],...
    %     %     ['\theta = ', num2str(7000)], 'Location','southeast')
    %
    %     saveas(h1, ['Param_', name_param_set, '/', 'SYS_', name_sys,...
    %         '_period'], 'png');
    %     close all
end

%% End Script
disp('Finished: PlotParamSweep')
disp(['System: ', name_sys]), disp(['Params: ', name_param_set])
toc
