function PAST = sim_laser_loop(params_laser, variable_param, params_temp, folder)
    tic

    % Get parameter to scan and range values
    parameter = variable_param(1);
    start_value = variable_param(2);
    end_value = variable_param(3);
    num_values = variable_param(4);

    % Temporal parameters
    h = params_temp(1);			%iteration step for RK4
    horizon = params_temp(2);	%total time of simulation
    transient = params_temp(3); %transient time (discarded)

    % Values of the varying parameter to scan
    values = linspace(start_value, end_value, num_values);
    direction = sign(end_value-start_value);

    PAST = [];

    if nargin < 4 || isempty(folder)
        folder = [uigetdir('', 'Select saving directory') '/'];
    end

    mkdir(folder);	% Create the folder
    progressbar;

    for inc=1:num_values
        params_laser(parameter) = values(inc);

        if ~isequal(params_temp(3),0)
            PAST = sim_laser_PCF(params_laser, [h transient], PAST);
        end

        PAST = sim_laser_PCF(params_laser, [h horizon], PAST);

        if ~isequal(direction, -1)
            save([folder,'sim_laser_pcf_',num2str(inc),'.mat'], 'PAST');
        else
            save([folder,'sim_laser_pcf_',num2str(num_values-inc+1),'.mat'], 'PAST');
        end

        progressbar(inc/num_values);
    end

    toc	% if you are curious about execution time
end
