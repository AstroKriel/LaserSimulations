%% Script Begining
%{
    Authors:    Neco Kriel  (2019)

    Purpose:
    We care to see how the time traces change with respect to the changing
    bifurcation parameter. Here we integrate our system, and plot two
    periods of the solution. We also plot the frequency information for the
    horizontal and vertical polarisations.

    Function Dependencies:
    - IntegSimLaser
    - findIndices

    Sub-Function Dependencies:
    - fast_sim_laser_pcf.c
    - fast_sim_laser_prof.c
    - fast_sim_laser_prof_noise.c
    - fast_sim_laser_prpcf.c
    - fast_sim_laser_prpcf_noise.c
%}

%% PREPARE WORKSPACE
clc, clear, close all

tic
disp('Started: PlotTime')

%% CHOOSE SYSTEM
% DIM     = 5; BIF_EQN = 1;
% name_sys    = 'PCF';

% DIM     = 4; BIF_EQN = 2;
% name_sys    = 'PROFN';


DIM     = 7; BIF_EQN = 3;
name_sys    = 'PRPCFN';

%% INITIALISE PARAMETERS
param_set_name = 'alpha'; % name_sys of parameter set (file-storage purposes)

% Initialise System's Parameter Values
load(['Param_', param_set_name, '/ParamSet.mat']);

num_cycles = 2; % number of cycles of hor. and ver. polar. plotted
eta_save   = [5]; % save figure for these vals. of eta
eta_tol    = 1e-2; % the tol. of eta for which fig. will be saved

% Sweeping Parameters
param_iter    = 30; % bif. resolution
param_elem    = 4;   % bif. w/respect to param_set
param_perturb = 1e-5;
param_vals    = linspace(param_start, param_end, param_iter) + param_perturb;

for alpha = [2:5]
    beta      = (1-ka)/(2*ka);
    param_set = [P, T, theta, 0, beta, ka, alpha, tau_R, omega, R]; % store params
    
    % Create Folder (load data from/save video to)
    folder = ['Param_', param_set_name, '/SYS_', name_sys, '_',...
        param_set_name, '=',...
        num2str(alpha)];

    eta     = importdata([folder, '/bif_eta.txt']);
    maxima  = importdata([folder, '/bif_extrema.txt']);
    log     = (eta <= param_end & eta >= param_start);
    eta     = eta(log);
    maxima  = maxima(log);
    
    % Initialise Analysis Parameters
    transient = 1e6;
    delay     = floor(theta/h);
    
    % Initialise System's Past
    sim_past = [];
    
    % Open figure and prepare video
    h1 = figure('Renderer', 'painters', 'Position', [10 10 1200 700]);
    vidObj = VideoWriter([folder, '.avi']); % create the video obj.
    open(vidObj); % start video
    
    %% INTEGRATE, ANALYSE and PLOT
    for iter = 1:param_iter
        param_set(param_elem) = param_vals(iter);
        if ~isequal(transient, 0)
            sim_past = integSimLaser(param_set, [h transient], sim_past, DIM,...
                name_sys, BIF_EQN);
        end
        sim_past = integSimLaser(param_set, [h horizon], sim_past, DIM,...
            name_sys, BIF_EQN);
        if any(any(isnan(sim_past)))
            disp(['System ', name_sys, ' has NaN values.'])
        else
            hor = (sim_past(:, 1)).^2;
            ver = (sim_past(:, BIF_EQN)).^2;
            if (any(abs(hor-ver) < 0.1))
                % calculate indices of one cycle of hor. and ver. polarizations
                [itter_start, itter_mid, itter_end] = findIndices(hor, ver);
                % normalise time wrt. delay
                x = (1:length(hor(itter_start:(itter_start +...
                    num_cycles*(itter_end - itter_start)))))./delay;
                % Plot Time
                subplot(3, 3, 4:9)
                plot(x, hor(itter_start:(itter_start + 2*(itter_end -...
                    itter_start))), 'b')
                hold on
                plot(x, ver(itter_start:(itter_start + 2*(itter_end -...
                    itter_start))), 'r')
            else
                % normalise time wrt. delay
                x = (1:length(hor(5*delay:10*delay)))./delay;
                % Plot Time
                subplot(3, 3, 4:9)
                plot(x, hor(5*delay:10*delay), 'b')
                hold on
                plot(x, ver(5*delay:10*delay), 'r')
            end
            title([name_sys, ': \{\eta=', num2str(param_set(param_elem)),...
                             ', \theta=', num2str(theta),...
                             ', \kappa=', num2str(ka),...
                             ', \alpha=', num2str(alpha),...
                             ', \tau_R=', num2str(tau_R),...
                             ', T=', num2str(T),...
                             ', R=', num2str(R),' \}'])
            legend({'Horizontal', 'Vertical'}, 'FontSize', 11)
            ylabel('Power (10 log-dB)', 'fontsize',14)
            ylim([0, 2])
            
            
            % Plot FFT - horizontal
            FFT_array_hor = abs(fft(sim_past(:, 1).^2));
            FFT_array_hor = FFT_array_hor(1:floor(length(FFT_array_hor)/2));
            x_hor = (1:length(FFT_array_hor))/(2*h*tau_P*length(FFT_array_hor));
            subplot(3, 3, 2), hold on
            plot(x_hor, FFT_array_hor, 'b-', 'MarkerSize', 1)
            title('Frequency: Horizontal', 'FontSize', 16)
            ylabel('Magnitude')
            xlabel('Frequencies (GHz)')
            xlim([1, 50])
            
            % Plot FFT - vertical
            FFT_array_ver = abs(fft(sim_past(:, BIF_EQN).^2));
            FFT_array_ver = FFT_array_ver(1:floor(length(FFT_array_ver)/2));
            x_ver = (1:length(FFT_array_ver))/(2*h*tau_P*length(FFT_array_ver));
            subplot(3, 3, 3), hold on
            plot(x_ver, FFT_array_ver, 'r-', 'MarkerSize', 1)
            title('Frequency: Vertical', 'FontSize', 16)
            ylabel('Magnitude')
            xlabel('Frequencies (GHz)')
            xlim([1, 50])
            
            % Plot Bifurcation
            subplot(3, 3, 1)
            plot(eta, maxima, '.k', 'MarkerSize', 5)
            hold on
            line([param_set(param_elem), param_set(param_elem)], [0, max(maxima)])
            title('Bifurcation Diagram', 'FontSize', 16)
            xlabel('\eta', 'FontSize', 16)
            ylabel('Power Output', 'FontSize', 16)
            grid on
            axis tight
            
            if (any(abs(param_set(param_elem) - eta_save) <= eta_tol))
                saveas(h1, [folder, '/', 'SYS_', name_sys, '_timetraces_eta=',...
                    num2str(param_set(param_elem))], 'png');
                savefig([folder, '/', 'SYS_', name_sys, '_timetraces_eta=',...
                    num2str(param_set(param_elem)), '.fig'])
            end
            
            % Save Current Frame
            writeVideo(vidObj, getframe(gcf));
            clf
        end
    end
    
    close(vidObj); % end video
    close all
end

%% End Script
disp('Finished: PlotTime')
disp(['System: ', name_sys]), disp(['Params: ', param_set_name])
toc
