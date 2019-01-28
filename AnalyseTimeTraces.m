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

% Initialise Analysis Parameters
h         = 1;
horizon   = 0.2e6;
transient = 1e6;
delay     = floor(theta/h);

% Initialise System's Past
sim_past = [];

%% Integrate, Analyse, and Plot
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
    
    % calculate indices of one cycle of hor. and ver. polarizations
    [itter_start, itter_mid, itter_end] = findIndices(hor, ver);
    
    % normalise time wrt. delay
    x = (1:length(hor(itter_start:(itter_start + num_cycles*(itter_end - itter_start)))))./delay;
    
    % Plot
    h1 = figure('Renderer', 'painters', 'Position', [100 150 800 500]);
    subplot(2, 1, 1), hold on
    plot(x, hor(itter_start:(itter_start + 2*(itter_end - itter_start))), 'b')
    title([name, 'Horizontal'])
    
    subplot(2, 1, 2), hold on
    plot(x, ver(itter_start:(itter_start + 2*(itter_end - itter_start))), 'r')
    title('Vertical')
    
    folder = ['SYS_', name, '_tau=', num2str(tau_R), '_theta=', num2str(theta)];
    mkdir(folder);
    csvwrite([folder, '/periodLength.txt'],...
        [(itter_mid - itter_start), (itter_end - itter_mid)]);
    saveas(h1, [folder, '/TimeTrace'], 'png');
    close all
end

%% End Script
disp('Program Finished: AnalyseTimeTraces.')
toc
