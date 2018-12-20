%% Script beginning
% Script created to integrate and analyse a single system (to test if
% saving issues could be fixed).
tic
disp('Script Started.')

format compact
format short
clc, clear, close all

%% Compile C code:
% Uncomment System of Interest
mex fast_detect_extrema.c

% mex fast_sim_laser_prof.c
% DIM = 4; BIF_EQN = 2;  name = 'PROF';

mex fast_sim_laser_prpcf.c
DIM = 7; BIF_EQN = 3; name = 'PRPCF';

% mex fast_sim_laser_pcf.c
% DIM = 5; BIF_EQN = 1; name = 'PCF';

%% INITIALISE PARAMETERS
% Sweeping Parameters:
num_itter              = 500;       % resolution
parfor param_set_num   = 1:2        % in {1, 2}
    for dir            = 1:2        % in {1, 2}
        for tau_R      = [20, 50]   % in {20, 50}
            % Tuning Parameters
            sim_tol         = 1e-2;
            sim_test_vals   = []; % simu. analysis points
            bif_elem        = 1;  % bif. w/respect to param.
            bif_start       = 1e-5;
            
            % Logical Variables
            % Simulation Vars:
            log_plot_sim        = true;
            log_fig_sim_save    = true;
            log_sim_data_save   = true;
            log_plot_bif        = true;
            % Bifurcation Vars:
            log_fig_bif_save    = true;
            log_bif_data_save   = true;
            % Frequency Analysis Vars:
            log_freq_analysis   = false;
            log_plot_freq       = false;
            log_fig_freq_save   = false;
            % Close all figures:
            log_fig_close       = false;
            
            % System Parameter Sets:
            switch param_set_num
                case 1
                    % Values as seen in paper
                    omega   = 0;
                    alpha   = 2;
                    ka      = 0.96;
                    T       = 1200;
                    P       = 0.6016;
                    theta   = 1143;
                case 2
                    % Typical Parameter Values
                    omega   = 0;
                    alpha   = 3;        % in [2, 3, 5]
                    ka      = 0.96;     % in [0.96, 1]
                    T       = 1000;     % in [1000, 500, 100]
                    P       = 0.6;      % in [0.6, 1, 2]
                    theta   = 7000;     % in [1000, 2500, 5000, 7000]
            end
            switch dir
                case 1
                    % Integrating Forward
                    param_vals = linspace(0, 0.15, num_itter) + bif_start;
                    
                    % Create Folder to Save Data and Plots
                    if log_fig_sim_save || log_sim_data_save || log_fig_bif_save ||...
                            log_bif_data_save
                        folder = ['Param', num2str(param_set_num),'_SYS_', name,...
                            '_tau=', num2str(tau_R), '_F/'];
                        mkdir(folder);
                    end
                case 2
                    % Integrating Backward
                    param_vals = linspace(0.15, 0, num_itter) + bif_start;
                    
                    % Create Folder to Save Data and Plots
                    if log_fig_sim_save || log_sim_data_save || log_fig_bif_save ||...
                            log_bif_data_save
                        folder = ['Param', num2str(param_set_num),'_SYS_', name,...
                            '_tau=', num2str(tau_R), '_B/'];
                        mkdir(folder);
                    end
            end
            eta       = 0; % initialise, but sweep parameter
            tau_P     = 1.4e-12;
            beta      = (1-ka)/(2*ka);
            param_set = [eta, omega, alpha, beta, ka, T, P, theta, tau_R];
            
            % Analysis Parameters
            h         = 1;
            horizon   = 0.2e6;
            transient = 1e6;
            delay     = floor(theta/h);
            
            % Initialise Bifurcation Parameters
            num_extrema     = 105;
            bif_epsilon     = 1e-3;
            num_val         = num_extrema*num_itter;
            bif_eta       = zeros(num_val, 1);
            bif_extrema     = zeros(num_val, 1);
            bif_index       = 1;
            
            % Initialise Frequency Analysis Parameters
            Chaos_bw            = zeros(num_itter, 1);
            Chaus_eff_bw        = zeros(num_itter, 1);
            Freq_ext_cav_node   = zeros(num_itter, 1);
            
            sim_past = [];
            
            %% Simulation
            for itter = 1:num_itter
                %% Integrate
                % dependencies: param_vals, bif_elem, param_vals, itter, transient,
                % sim_past, h, DIM
                param_set(bif_elem) = param_vals(itter);
                
                if ~isequal(transient, 0)
                    sim_past = integ_sim_laser(param_set, [h transient], sim_past, DIM);
                end
                
                sim_past = integ_sim_laser(param_set, [h horizon], sim_past, DIM);
                
                if any(any(isnan(sim_past)))
                    disp(['2. System ', name,...
                        ' has NaN values at itter = ', num2str(itter)])
                    break
                end
                % output: sim_past
                
                %% Bifurcation Analysis
                % dependencies: sim_past, BIF_EQN, num_extrema, bif_epsilon, bif_index,
                % param_vals, itter, bif_eta, bif_extrema
                bif_sim_past = abs(sim_past(1:end, BIF_EQN));
                temp_output = detect_extrema(bif_sim_past, num_extrema, bif_epsilon);
                temp_len = length(temp_output);
                bif_eta(bif_index:bif_index+temp_len-1) = param_vals(itter);
                bif_extrema(bif_index:bif_index+temp_len-1) = temp_output;
                bif_index = bif_index + temp_len;
                % ouput: bif_eta, bif_extrema, bif_index
                
                %% Plot
                % dependencies: log_plot_sim, tol, test_num, param_set, bif_elem,
                % delay, sim_past, log_fig_sim_save, folder, itter, log_sim_data_save,
                % log_fig_close
                if (log_plot_sim && (any((1-sim_tol)*sim_test_vals <= param_set(bif_elem)...
                        & param_set(bif_elem) <= (1+sim_tol)*sim_test_vals)))
                    % Plot: Time Traces
                    h1 = figure('Renderer', 'painters', 'Position', [10 10 900 600]);
                    subplot(2, 2, 1)
                    plot(sim_past(:, 1).^2, 'b')
                    ylabel('Power (10 log-dB)', 'fontsize',14)
                    title('Horizontal Polarisation', 'fontsize',18)
                    axis([15*delay, 20*delay, 0, 1])
                    
                    subplot(2, 2, 3)
                    plot(sim_past(:, BIF_EQN).^2, 'r')
                    ylabel('Power (10 log-dB)', 'fontsize',14)
                    title('Vertical Polarisation', 'fontsize',18)
                    axis([15*delay, 20*delay, 0, 1])
                    
                    % Plot radio frequency/RF spectrum FFT in dB
                    subplot(2, 2, [2, 4])
                    plot(log(abs(fft(sim_past(:, BIF_EQN).^2))))
                    xlim([5, 6e3])
                    title(['RF spectrum FFT in dB: eta=',...
                        num2str(param_set(bif_elem))], 'fontsize',18)
                    
                    % Save: Fig Time Trace
                    if log_fig_sim_save
                        saveas(h1,[folder, 'Plot_', num2str(itter),'.png'],'png');
                    end
                    
                    % Save: Time Trace Data
                    if log_sim_data_save
                        csvwrite([folder, 'Data_', num2str(itter),'_hor_polar.txt'],...
                            sim_past(:, 1).^2)
                        csvwrite([folder, 'Data_', num2str(itter),'_ver_polar.txt'],...
                            sim_past(:, BIF_EQN).^2)
                        csvwrite([folder, 'Data_', num2str(itter),'_FFT.txt'],...
                            log(abs(fft(sim_past(:, BIF_EQN).^2))))
                    end
                    
                    % Close figs in large parameter sweeps
                    if log_fig_close
                        close all
                    end
                end
                % output: none
                
                %% Evaluate/Analyse Frequency
                % dependencies: log_freq_analysis, sim_past, BIF_EQN, tau_P, Chaos_bw,
                % Chaus_eff_bw, itter, Freq_ext_cav_node
                if log_freq_analysis
                    FFT_array = abs(fft(sim_past(:, BIF_EQN)));
                    P_array = FFT_array(1:ceil(length(FFT_array)/2));
                    F_array = (1:length(P_array))/(tau_P*length(sim_past(:, BIF_EQN)));
                    
                    [tab_periode, ~] = detectionPERIODE(sim_past(:, BIF_EQN));
                    
                    if min(tab_periode) == -10
                        Chaos_bw(itter) = 0;
                        Chaus_eff_bw(itter) = 0;
                    else
                        if min(tab_periode) == 0
                            Chaos_bw(itter) = chaosBD(P_array, F_array);
                            Chaus_eff_bw(itter) = eff_chaosBD(P_array, F_array);
                            
                            if Chaos_bw(itter) < 10e6
                                Chaos_bw(itter)     = 0;
                                Chaus_eff_bw(itter) = 0;
                            end
                        else
                            unevaleur = 1/(mean(tab_periode)*tau_P);
                            Freq_ext_cav_node(itter) = unevaleur;
                        end
                    end
                end
                % output: Chaos_bw, Chaus_eff_bw, Freq_ext_cav_node
            end
            
            %% Bifurcation Diagram
            % dependencies: log_plot_bif, sim_past, bif_eta, bif_extrema, bif_index,
            % log_fig_bif_save, folder, log_bif_data_save
            if (log_plot_bif && ~any(any(isnan(sim_past))))
                bif_eta = bif_eta(1:bif_index-1);
                bif_extrema = bif_extrema(1:bif_index-1);
                
                h2 = figure('Renderer', 'painters', 'Position', [10 10 900 500]);
                plot(bif_eta, bif_extrema, '.k', 'MarkerSize', 1)
                grid on, axis tight
                title(['Bifurcation Diagram: System ', name], 'fontsize', 16)
                xlabel('\eta', 'fontsize', 16)
                ylabel('Output power (a. u.)', 'fontsize', 16)
                
                % Save: Bifurcation Plot
                if log_fig_bif_save
                    saveas(h2, [folder, 'bif_diag_', name,'.png'], 'png');
                end
                
                % Save: Bifurcation Data
                if log_bif_data_save
                    csvwrite([folder, 'bif_extrema.txt'], bif_extrema);
                    csvwrite([folder, 'bif_eta.txt'], bif_eta);
                end
            end
            % output: bif_eta, bif_extrema
            
            %% Plot Frequency
            % dependencies: log_plot_freq, param_vals, Chaos_bw, Freq_ext_cav_node,
            % log_fig_freq_save, folder, log_fig_close
            if log_plot_freq
                h3 = figure;
                hold on
                plot(param_vals, Chaos_bw)
                plot(param_vals, Freq_ext_cav_node, 'o')
            end
            
            if log_fig_freq_save
                saveas(h3, [folder, 'freq_analysis_system_', name,'.png'], 'png');
            end
            
            if log_fig_close
                close all
            end
            % output: none
        end
    end
end

toc
disp('Script finished.')

%% Spript end