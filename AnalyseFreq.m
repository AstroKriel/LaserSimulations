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

% mex fast_sim_laser_prpcf.c
% DIM = 7; BIF_EQN = 3; name = 'PRPCF';
%
mex fast_sim_laser_pcf.c
DIM = 5; BIF_EQN = 1; name = 'PCF';

%% INITIALISE PARAMETERS
% Sweeping Parameters:
num_itter       = 50;  % resolution
tau_R           = 50;    % in {20, 50}
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

% omega   = 0;
% alpha   = 2;
% ka      = 0.96;
% T       = 250;
% P       = 0.6;
% theta   = 7000;

% omega   = 0;
% alpha   = 5;
% ka      = 0.98;
% T       = 250;
% P       = 0.6;
% theta   = 7000;

eta     = 0;            % initialise, but sweep parameter
beta    = (1-ka)/(2*ka);

% Integrating Forward
param_vals  = linspace(0, 0.95, num_itter) + bif_start;
param_set   = [eta, omega, alpha, beta, ka, T, P, theta, tau_R];

% Analysis Parameters
h           = 0.5;
horizon     = 0.2e6;
transient   = 1e6;
delay       = floor(theta/h);
sim_past    = [];

tau_p       = 1.4e-3;
L           = floor(horizon/h) + 1;
NFFT        = 2^nextpow2(L);
f = 1/(2*h*tau_p)*linspace(0, 1/16, NFFT/32+1);
DISPLAY = zeros(NFFT/32+1, num_itter);
SPECTRA = zeros(NFFT,1);

%% Simulation
for itter = 1:num_itter
    if (mod(itter, 50) == 0)
        disp(itter)
    end
    
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
    
    SPECTRA = fft(sim_past(:,BIF_EQN).^2, NFFT)/L;
    DISPLAY(:,itter) = SPECTRA(1:NFFT/32+1);
end

h1 = figure('Renderer', 'painters', 'Position', [10 10 1200 700]);
imagesc(param_vals, f, log10(2*abs(DISPLAY)));
colormap(hot)
caxis([-4 0])
colorbar
set(gca,'YDir','normal')
axis tight

% folder = ['SYS_', name, '_Param=', num2str(param_set_num),...
%     '_tau=', num2str(tau_R), '_F'];
% mkdir(folder);
% csvwrite([folder, '/fft_display.txt'], DISPLAY);

% saveas(h1, [folder, '_FFT'], 'png');

disp('Program Finished')
toc

%% END SCRIPT