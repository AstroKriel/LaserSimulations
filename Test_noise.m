%% Start Script
format compact
format shortg
clc, clear, close all

% Compile C-Programs (if changes have been made)
% mex fast_sim_laser_prof.c
% mex fast_sim_laser_prof_noise.c

tic
disp('Program Started')

%% Initialise System
P       = 0.6;
T       = 250;
theta   = 7e3;
eta     = 0.08;
ka      = 0.96;
beta    = (1-ka)/(2*ka);
alpha   = 2;
omega   = 0;
tau_R   = 0;
tau_P   = 1.4e-3;
R       = 0;

params_laser = [P, T, theta, eta, beta, ka, alpha, tau_R, omega, R];
params_temp  = [1, 1e7];

delay = 7000;
steps = 1e7;

%% PROF Without Noise
disp('Plot Without Noise')

PAST    = 0.1 + zeros(delay+1, 4) + 1e-5;
RESULTS = fast_sim_laser_prof(params_laser, params_temp, PAST);
RESULTS = RESULTS(4*(delay+1)-3:end);
OUTPUT  = abs(reshape(RESULTS, 4, steps+1).');
hor     = OUTPUT(2*delay:end, 1);
ver     = OUTPUT(2*delay:end, 2);

[itter_start, itter_end] = findIndices(hor, ver);

x = (1:length(hor(itter_start:(itter_start + 2*(itter_end - itter_start)))))./delay;

% Plot E_1 and E_2
figure('Renderer', 'painters', 'Position', [50 50 800 500])
subplot(2, 1, 1)
plot(x, hor(itter_start:(itter_start + 2*(itter_end - itter_start))))
title('Horizontal')
subplot(2, 1, 2)
plot(x, ver(itter_start:(itter_start + 2*(itter_end - itter_start))))
title('Vertical')

%% PROF With noise
R = [1e-6, 1e-12, 1e-20];

for i = 1:length(R)
    disp(['Plot With Noise, R = ', num2str(R(i))])
    
    params_laser(9) = R(i);     % add noise
    NOISE = randn(steps+1, 2);  % initialise the noise
    
    PAST    = 0.1 + zeros(delay+1, 4) + 1e-5;
    RESULTS = fast_sim_laser_prof_noise(params_laser, params_temp, PAST, NOISE);
    RESULTS = RESULTS(4*(delay+1)-3:end);
    OUTPUT  = abs(reshape(RESULTS, 4, steps+1).');
    hor     = OUTPUT(2*delay:end, 1);
    ver     = OUTPUT(2*delay:end, 2);
    
    [itter_start, itter_end] = findIndices(hor, ver);
    
    x = (1:length(hor(itter_start:(itter_start + 2*(itter_end - itter_start)))))./delay;
    
    % Plot E_1 and E_2
    subplot(2, 1, 1), hold on
    plot(x, hor(itter_start:(itter_start + 2*(itter_end - itter_start))))
    subplot(2, 1, 2), hold on
    plot(x, ver(itter_start:(itter_start + 2*(itter_end - itter_start))))
end

subplot(2, 1, 1), legend('No Noise',...
    ['R = ', num2str(R(1))],...
    ['R = ', num2str(R(2))],...
    ['R = ', num2str(R(3))], 'Location', 'northeast')
subplot(2, 1, 2), legend('No Noise',...
    ['R = ', num2str(R(1))],...
    ['R = ', num2str(R(2))],...
    ['R = ', num2str(R(3))], 'Location', 'southeast')

%% End Script
disp('Program Finished')
toc

%% Define Function
function [itter_start, itter_end_2] = findIndices(hor, ver)
    points_found = false;
    itter_start = 1;
    while (itter_start < length(hor) && ~points_found)
        cur_hor = hor(itter_start);
        cur_ver = ver(itter_start);
        nex_hor = hor(itter_start + 1);
        nex_ver = ver(itter_start + 1);
        if (cur_hor <= cur_ver && nex_hor >= nex_ver && ~points_found)
            % hor grad > 0 & overlapping -> start period found
            itter_end_1 = itter_start + 1;
            while (itter_end_1 < length(hor) && ~points_found)
                cur_hor = hor(itter_end_1);
                cur_ver = ver(itter_end_1);
                nex_hor = hor(itter_end_1 + 1);
                nex_ver = ver(itter_end_1 + 1);
                if (cur_hor >= cur_ver && nex_hor <= nex_ver && ~points_found)
                    % hor grad < 0 & overlapping -> end period found
                    itter_end_2 = itter_end_1 + 1;
                    while (itter_end_2 < length(hor) && ~points_found)
                        cur_hor = hor(itter_end_2);
                        cur_ver = ver(itter_end_2);
                        nex_hor = hor(itter_end_2 + 1);
                        nex_ver = ver(itter_end_2 + 1);
                        if (cur_hor <= cur_ver && nex_hor >= nex_ver && ~points_found)
                            % hor grad > 0 & overlapping -> end period found
                            points_found = true;
                        end
                        
                        % not overlapping
                        itter_end_2 = itter_end_2 + 1;
                    end
                end
                
                % not overlapping
                itter_end_1 = itter_end_1 + 1;
            end
        end
        
        % not overlapping
        itter_start = itter_start + 1;
    end
end
