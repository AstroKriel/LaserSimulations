format compact
format short
clc, clear, close all

tic
disp('Program Started: PlottingParamSweep.')

%% Define Variables
% % name = 'PROF';
name = 'PROFN';
% % name = 'PRPCF';
% name = 'PRPCFN';
% % name = 'PCF';

paper_name = 'SOAPS_2';

% Sweeping Parameters
param_itter   = 100; % bif. resolution
param_perturb = 1e-5;
param_start   = 0;
param_end     = 0.2;
param_vals    = linspace(param_start, param_end, param_itter) + param_perturb;

max_freq = 0;
h1 = figure('Renderer', 'painters', 'Position', [10 10 1200 700]);
hold on
for theta = [5000, 7000]
    pre_name = '';
    filename = ['Param_', paper_name, '/SYS_', name, '_theta=', num2str(theta)];
    time = importdata([filename, '/', 'period_length', '.txt']);
    time = time(:, 2);
    
    plot(param_vals, time, 'o', 'MarkerSize', 3)
end
title(name)
xlabel('\eta')
ylabel('period')
legend(['\theta = ', num2str(100)],...
    ['\theta = ', num2str(500)],...
    ['\theta = ', num2str(1000)],...
    ['\theta = ', num2str(3000)],...
    ['\theta = ', num2str(5000)],...
    ['\theta = ', num2str(7000)], 'Location','southeast')

saveas(h1, ['Param_', paper_name, '/', 'SYS_', name, '_period'], 'png');
close all

h1 = figure('Renderer', 'painters', 'Position', [10 10 1200 700]);
hold on
for theta = [5000, 7000]
    pre_name = '';
    filename = ['Param_', paper_name, '/SYS_', name, '_theta=', num2str(theta)];
    time = importdata([filename, '/', 'period_length', '.txt']);
    time = time(:, 3);
    
    plot(param_vals, time, 'o', 'MarkerSize', 3)
end
title(name)
xlabel('\eta')
ylabel('epsilon (a.u.)')
legend(['\theta = ', num2str(100)],...
    ['\theta = ', num2str(500)],...
    ['\theta = ', num2str(1000)],...
    ['\theta = ', num2str(3000)],...
    ['\theta = ', num2str(5000)],...
    ['\theta = ', num2str(7000)], 'Location','southeast')

saveas(h1, ['Param_', paper_name, '/', 'SYS_', name, '_epsilon'], 'png');
close all

%% End Script
disp('Program Finished: PlottingParamSweep.')
toc
