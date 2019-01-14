%% Bifurcation and Frequency Plots
format compact 
clc, clear, close all

% name = 'PCF';
name = 'PROF';
% name = 'PRPCF';

param = 2;
tau = 50;
filename = ['SYS_',name,'_Param=',num2str(param),'_tau=',num2str(tau),'_F'];

% bif_length = 'Long';
% param_start = 0;
% param_end = 0.95;

bif_length = 'Shorter';
param_start = 0.06;
param_end = 0.07;

h           = 0.5;
horizon     = 0.2e6;
tau_p       = 1.4e-3;
L           = floor(horizon/h) + 1;
NFFT        = 2^nextpow2(L);
f = 1/(2*h*tau_p)*linspace(0, 1/16, NFFT/32+1);
param_vals = linspace(param_start, param_end, 1000);

fft = importdata([filename, '/', 'fft_display_', bif_length, '.txt']);

xF = importdata([filename, '/', 'bif_eta_', bif_length,'.txt']);
yF = importdata([filename, '/', 'bif_extrema_', bif_length,'.txt']);
log = (xF >= param_start & xF <= param_end);
xF = xF(log);
yF = yF(log);

h1 = figure('Renderer', 'painters', 'Position', [10 10 1200 700]);

subplot(2, 1, 1) 
imagesc(param_vals, f, log10(2*abs(fft)));
colormap(hot)
caxis([-4 0])
cb = colorbar;
set(cb,'Position',[0.92 0.6 .01 0.3])
set(gca,'YDir','normal')
axis tight

title(['SYS ',name,', Param=',num2str(param),', tau=',num2str(tau), ', ',...
    bif_length], 'fontsize',18)

subplot(2, 1, 2)
plot(xF, yF, 'k.', 'MarkerSize', 3)
axis tight

saveas(h1, ['Plots/FFT_', filename, '_', bif_length], 'png');
% close all