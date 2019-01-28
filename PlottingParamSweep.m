%% Script Begining
%{
    Authors:    Neco Kriel  (2019)
                Chi Hak     (2017)

    Purpose:    

%}

format compact
format short
clc, clear, close all

tic
disp('Program Started: PlottingParamSweep.')

%% Define Variables
name = 'PROF';
% name = 'PROFN';
% name = 'PRPCF';
% name = 'PRPCFN';
% name = 'PCF';

paper_name = 'Typ';
% paper_name = 'HO_ECM';
% paper_name = 'M_ECM';
% paper_name = 'SOAPS';
% paper_name = 'R_SW_O_SL_1';
% paper_name = 'R_SW_O_SL_2';
% paper_name = 'R_SW_O_SL_3';

pre_name = '';
filename = ['Param_', paper_name, '/SYS_', name];

param_start = 0;
param_end   = 0.99;
param_res   = 1e3;

param_vals  = linspace(param_start, param_end, param_res);
h           = 0.5;
horizon     = 0.2e6;
tau_p       = 1.4e-3;
L           = floor(horizon/h) + 1;
NFFT        = 2^nextpow2(L);
f           = 1/(2*h*tau_p)*linspace(0, 1/16, NFFT/32+1);

% Import data
fft = importdata([filename, '/', 'fft_display', '.txt']);
xF  = importdata([filename, '/', 'bif_eta', '.txt']);
yF  = importdata([filename, '/', 'bif_extrema', '.txt']);
% Truncate data
log = (xF <= param_end & xF >= param_start);
xF  = xF(log);
yF  = yF(log);

h1 = figure('Renderer', 'painters', 'Position', [10 10 1200 700]);
subplot(2, 1, 1)
imagesc(param_vals, f, log10(2*abs(fft)));
colormap(hot)
caxis([-4 0])
cb = colorbar;
set(cb,'Position',[0.92 0.6 .01 0.3])
set(gca,'YDir','normal')
axis tight
title(['Power Spectrum: ', name, ' System'], 'FontSize', 16)

subplot(2, 1, 2)
plot(xF, yF, 'k.', 'MarkerSize', 3)
axis tight
title('Bifurcation Diagram', 'FontSize', 16)

saveas(h1, filename, 'png');
close all

%% End Script
disp('Program Finished: PlottingParamSweep.')
toc 
