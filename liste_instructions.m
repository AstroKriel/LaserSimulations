%% START MATLAB
format compact
format short
clc, clear, close all

tic
%% COMPILE C-FILE
mex fast_sim_laser_SYS1.c % DIM = 4, Bif_num = 2
mex fast_sim_laser_SYS3.c % DIM = 7, Bif_num = 3
mex fast_sim_laser_SYS4.c % DIM = 5, Bif_num = 1
mex fast_detect_extrema.c

%% INITIALISE PARAMETERS
% Logical Variables
log_plot_sim        = false;    % simulation logical variables
log_fig_sim_save    = false;
log_sim_data_save   = false;
log_plot_bif        = false;
log_fig_bif_save    = false;    % bifurcation logical variables
log_bif_data_save   = false;
log_fig_close       = false;    % close fig
log_freq_analysis   = false;    % analyse frequency

% System Parameters
% % Param Set 1: PCF Paper Parameters Values
% var_param = 1;
% eta     = 0;
% omega   = 0;
% alpha   = 2;
% ka      = 0.96;
% beta    = (1-ka)/(2*ka);
% T       = 1200;
% P       = 0.6016;
% theta   = 1143;
% tau_R   = 20;

% % Param Set 2: Typical Parameter Values
var_param = 2;
eta     = 0;                          % sweep [0, 0.06]
omega   = 0;
alpha   = 3;                          % in [2, 3, 5]
ka      = 0.96;                       % in [0.96, 1]
beta    = (1-ka)/(2*ka);
T       = 1000;                       % in [1000, 500, 100]
P       = 0.6;                        % in [0.6, 1, 2]
theta   = 7000;                       % in [1000, 2500, 5000, 7000]
tau_R   = 50;                         % in {20, 50}
tau_P   = 1.4e-12;

for SYS = 3
    switch SYS
        case 1
            % System 1
            VAR = 1;
            DIM = 4;
            bif_eqn = 2;
        case 2
            % System 2
            VAR = 3;
            DIM = 7;
            bif_eqn = 3;
        case 3
            % System 3
            VAR = 4;
            DIM = 5;
            bif_eqn = 1;
        otherwise
            % Eror: no system defined
            disp('System number is wrong!')
    end
    
    % Create Folder
    Name = ['Param', num2str(var_param),'_SYS', num2str(VAR),'_tau=', num2str(tau_R)];
    folder = [Name, '/'];
    mkdir(folder);
    
    % Evaluation Parameters
    param_list      = [eta, omega, alpha, beta, ka, T, P, theta, tau_R];
    
    % Analysis Parameters
    test_num        = 0.09:0.01:0.15;   % analysis points
    tol             = 1e-3;
    num_itter       = 100;                   % <- CHANGES: resolution
    param_num       = 1;                    % <- CHANGES: sys bif param
    param_vals      = linspace(0, 0.15, num_itter);
    h               = 1;
    horizon         = 0.2e6;
    transient       = 1e6;
    direction       = sign(param_vals(end) - param_vals(1));
    delay           = floor(theta/h);
    
    % Bifurcation Parameters
    num_extrema     = 105;
    bif_epsilon     = 1e-3;
    num_val         = num_extrema*num_itter;
    bif_gamma       = zeros(num_val, 1);
    bif_extrema     = zeros(num_val, 1);
    bif_index       = 1;
    
    % Frequency Analysis Parameters
    Chaos_bw        = zeros(num_itter, 1);
    Chaus_eff_bw    = zeros(num_itter, 1);
    Freq_ext_cav_node = zeros(num_itter, 1);
    
    sim_past = [];
    
    %% Simulation
    disp('Simulation Started.')
    tic
    
    for itter = 1:num_itter
        %% Integrate
        % dependencies: param_list, param_num, itter, param_vals, h,
        % transient, sim_past, DIM, horizon, VAR
        param_list(param_num) = param_vals(itter);
        
        if ~isequal(transient, 0)
            sim_past = sim_laser_PCF(param_list, [h transient], sim_past, DIM);
        end
        
        sim_past = sim_laser_PCF(param_list, [h horizon], sim_past, DIM);
        
        if any(any(isnan(sim_past)))
            disp(['System ', num2str(VAR), ' has NaN values'])
            break
        end
        % output: sim_past, param_list
        
        %% Bifurcation Analysis
        % dependencies: bif_eqn, sim_past, num_extrema, bif_epsilon,
        % bif_index, param_vals, itter
        bif_sim_past = abs(sim_past(1:end, bif_eqn));
        temp_output = detect_extrema(bif_sim_past, num_extrema, bif_epsilon);
        temp_len = length(temp_output);
        bif_gamma(bif_index:bif_index+temp_len-1) = param_vals(itter);
        bif_extrema(bif_index:bif_index+temp_len-1) = temp_output;
        bif_index = bif_index + temp_len;
        % ouput: bif_gamma, bif_extrema, bif_index
        
        %% Plot
        % dependencies: log_plot_sim, tol, [param_list, param_num], delay,
        % bif_eqn, sim_past, log_fig_sim_save, itter, log_sim_data_save,
        % folder, log_fig_close
        if (log_plot_sim && (any((1-tol)*test_num <= param_list(param_num)...
                & param_list(param_num) <= (1+tol)*test_num)))
            % Plot: Time Traces
            h1 = figure('Renderer', 'painters', 'Position', [10 10 900 600]);
            subplot(2, 2, 1)
            plot(sim_past(:, 1).^2, 'b')
            ylabel('Power (10 log-dB)', 'fontsize',14)
            title('Horizontal Polarisation', 'fontsize',18)
            axis([15*delay, 20*delay, 0, 1])
            
            subplot(2, 2, 3)
            plot(sim_past(:, bif_eqn).^2, 'r')
            ylabel('Power (10 log-dB)', 'fontsize',14)
            title('Vertical Polarisation', 'fontsize',18)
            axis([15*delay, 20*delay, 0, 1])
            
            % Plot radio frequency/RF spectrum FFT in dB
            subplot(2, 2, [2, 4])
            plot(log(abs(fft(sim_past(:, bif_eqn).^2))))
            xlim([5, 6e3])
            title(['RF spectrum FFT in dB: eta=',...
                num2str(param_list(param_num))], 'fontsize',18)
            
            % Save: Fig Time Trace
            if log_fig_sim_save
                saveas(h1,[folder, 'Plot_', num2str(itter),'.png'],'png');
            end
            
            % Save: Time Trace Data
            if log_sim_data_save
                csvwrite([folder, 'Data_', num2str(itter),'_hor_polar.txt'],...
                    sim_past(:, 1).^2)
                csvwrite([folder, 'Data_', num2str(itter),'_ver_polar.txt'],...
                    sim_past(:, bif_eqn).^2)
                csvwrite([folder, 'Data_', num2str(itter),'_FFT.txt'],...
                    log(abs(fft(sim_past(:, bif_eqn).^2))))
            end
            
            % Close figs in large parameter sweeps
            if log_fig_close
                close all
            end
        end
        
        %% Evaluate/Analyse Frequency
        % dependencies: log_freq_analysis, sim_past, bif_eqn, tau_P,
        % Chaos_bw, eff_chaosBD, itter, Freq_ext_cav_node
        FFT_array = abs(fft(sim_past(:, bif_eqn)));
        P_array = FFT_array(1:ceil(length(FFT_array)/2));
        F_array = (1:length(P_array))/(tau_P*length(sim_past(:, bif_eqn)));
                
        if log_freq_analysis
            [tab_periode, ~] = detectionPERIODE(sim_past(:, bif_eqn));
            
            if min(tab_periode) == -10 % stable
                Chaos_bw(itter) = 0;
                Chaus_eff_bw(itter) = 0;
            else
                if min(tab_periode) == 0 % ap�riodique
                    Chaos_bw(itter) = chaosBD(P_array, F_array); % La bande passante de chaos
                    Chaus_eff_bw(itter) = eff_chaosBD(P_array, F_array);
                                        
                    if Chaos_bw(itter) < 10e6
                        % pour r�gler le probl�me des transitoires infinis de la restabilisation
                        Chaos_bw(itter)     = 0;
                        Chaus_eff_bw(itter) = 0;
                    end
                else % p�riodique
                    unevaleur = 1/(mean(tab_periode)*tau_P);
                    Freq_ext_cav_node(itter) = unevaleur;
                end
            end
        end
        % output: Freq_ext_cav_node, Chaos_bw, eff_chaosBD
    end
    
    %% Bifurcation Diagram
    % dependencies: log_plot_bif, sim_past, bif_gamma, bif_index,
    % bif_extrema, log_fig_bif_save, folder, log_fig_close
    if (log_plot_bif && ~any(any(isnan(sim_past))))
        bif_gamma = bif_gamma(1:bif_index-1);
        bif_extrema = bif_extrema(1:bif_index-1);
        
        h2 = figure('Renderer', 'painters', 'Position', [10 10 900 500]);
        plot(bif_gamma, bif_extrema, '.k', 'MarkerSize', 1)
        grid on, axis tight
        title(['Bifurcation Diagram: ', Name], 'fontsize', 16)
        xlabel('\eta', 'fontsize', 16)
        ylabel('Output power (a. u.)', 'fontsize', 16)
        
        % Save: Bifurcation Plot
        if log_fig_bif_save
            saveas(h, [folder, 'bif_diag.png'], 'png');
        end
        
        % Save: Bifurcation Data
        if log_bif_data_save
            csvwrite([folder, 'bif_extrema.txt'], bif_extrema);
            csvwrite([folder, 'bif_gamma.txt'], bif_gamma);
        end
        
        if log_fig_close
            close all
        end
    end
    
    %% Plot Frequency
    figure, hold on
    plot(param_vals, Chaos_bw)
    plot(param_vals, Freq_ext_cav_node, 'o')
end

%% END MATLAB
disp('Simulation Completed')
toc
