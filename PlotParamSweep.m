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
clc, clear, close all

tic
disp('Started: PlotParamSweep')

%% INITIALISE AND DEFINE PARAMETERS
% System and Parameter set name
% name_sys       = 'PRPCFN';
name_sys    = 'PRPCFUFN';

name_param_set = 'alpha'; % foldername to save data/figures

load(['Param_', name_param_set, '/ParamSet.mat']);

% FFT variables
L    = floor(horizon/h) + 1;
NFFT = 2^nextpow2(L);
f    = 1/(2*h*tau_P)*linspace(0, 1/16, NFFT/32+1);

for alpha = [2:5]
    % file save path
    filename = ['Param_', name_param_set, '/SYS_', name_sys, '_',...
        name_param_set, '=',...
        num2str(alpha)];
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
        ', \alpha=', num2str(ka),...
        ', \tau_R=', num2str(tau_R),...
        ', R=', num2str(R),' \}'],...
        'FontSize', 16)
    saveas(h1, ['Param_', name_param_set, '/', 'SYS_', name_sys, '_',...
        name_param_set,...
        num2str(alpha), '.png'])
    savefig(['Param_', name_param_set, '/', 'SYS_', name_sys, '_',...
        name_param_set,...
        num2str(alpha), '.fig'])
    close all
end

%% End Script
disp('Finished: PlotParamSweep')
disp(['System: ', name_sys]), disp(['Params: ', name_param_set])
toc
