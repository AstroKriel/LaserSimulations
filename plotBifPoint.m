clc, clear, close all

prpcf_1   = [0.037889; 0.035364; 0.037889; 0.035364; 0.035364; 0.035364; 0.035364]; % Tau = 1
prpcf_2   = [0.040414; 0.040414; 0.040414; 0.040414; 0.040414; 0.040414; 0.040414]; % Tau = 20
prpcf_3   = [0.045465; 0.042939; 0.037889; 0.037889; 0.037889; 0.037889; 0.037889]; % Tau = 50
prpcf_4   = [0.042939; 0.037889; 0.035364; 0.032838; 0.032838; 0.032838; 0.035364]; % Tau = 100
prpcf_5   = [0.032838; 0.032838; 0.032838; 0.032838; 0.032838; 0.032838; 0.032838]; % Tau = 500
prpcf_6   = [0.032838; 0.032838; 0.032838; 0.032838; 0.032838; 0.032838; 0.032838]; % Tau = 1000

prof    = [0.037889; 0.015162; 0.015162; 0.015162; 0.015162; 0.015162; 0.015162];  % Tau neg.

pcf_1 = [0.0050605; 0.0050605; 0.0025353; 0.0025353; 0.0025353; 0.0025353; 0.0025353];  % Tau = 1
pcf_2 = [0.0075858; 0.0075858; 0.0050605; 0.0050605; 0.0050605; 0.0050605; 0.0050605];  % Tau = 50
pcf_3 = [0.0075858; 0.0075858; 0.0075858; 0.0075858; 0.0075858; 0.0075858; 0.0075858];  % Tau = 1000

theta       = [100, 500, 1e3:2e3:1e4];

figure('Renderer', 'painters', 'Position', [10 10 1000 500]), hold on

% PCF, Tau = 1
plot(pcf_1, theta,  'bo') % 1
plot(pcf_2, theta,  'b^') % 50
plot(pcf_3, theta,  'bs') % 1000

% PROF
plot(prof, theta,  'ko', 'MarkerSize',15)

% PRPCF
plot(prpcf_1, theta,  'ro', 'MarkerSize', 5) % 1
plot(prpcf_2, theta,  'co', 'MarkerSize', 10) % 20
plot(prpcf_3, theta,  'r^', 'MarkerSize', 5) % 50
plot(prpcf_4, theta,  'c^', 'MarkerSize', 10) % 100
plot(prpcf_5, theta,  'cs', 'MarkerSize', 10) % 500
plot(prpcf_6, theta,  'rs', 'MarkerSize', 5) % 1000

legend('PCF: Tau_R = 1', 'PCF: Tau_R = 50', 'PCF: Tau_R = 1000',...
    'PROF',...
    'PRPCF: Tau_R = 1', 'PRPCF: Tau_R = 20', 'PRPCF: Tau_R = 50',...
    'PRPCF: Tau_R = 100', 'PRPCF: Tau_R = 500', 'PRPCF: Tau_R = 1000')
xlabel('Bif Point')
ylabel('Theta')

