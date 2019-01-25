function OUTPUT = integ_sim_laser(params_laser, params_temp, INPUT, DIM, name, BIF)
    theta   = params_laser(3);
    h       = params_temp(1);
    horizon = params_temp(2);
    delay   = floor(theta/h);
    steps   = floor(horizon/h);
    NOISE   = randn(steps+1, BIF);

    if (nargin < 3) || isempty(INPUT)
        PAST = 0.1+zeros(delay+1, DIM) + 1e-5;
    elseif length(INPUT) < delay+2
        PAST = vertcat(zeros(delay+1-size(INPUT, 1), DIM) + 1e-5, INPUT);
    else
        PAST = INPUT(end-delay:end,:);
    end

    PAST = reshape(PAST.', DIM*(delay+1), 1);

    switch name
        % PCF
        case 'PCF'
            RESULTS = fast_sim_laser_pcf(params_laser, params_temp, PAST);
        case 'PCFN'
            RESULTS = fast_sim_laser_pcf_noise(params_laser, params_temp,...
                PAST, NOISE);
        % PROF
        case 'PROF'
            RESULTS = fast_sim_laser_prof(params_laser, params_temp, PAST);
        case 'PROFN'
            RESULTS = fast_sim_laser_prof_noise(params_laser, params_temp,...
                PAST, NOISE);
        % PRPCF
        case 'PRPCF'
            RESULTS = fast_sim_laser_prpcf(params_laser, params_temp, PAST);
        case 'PRPCFN'
            RESULTS = fast_sim_laser_prpcf_noise(params_laser, params_temp,...
                PAST, NOISE);
    end

    RESULTS = RESULTS(DIM*delay + 1:end);
    OUTPUT = reshape(RESULTS, DIM, steps+1).';
end
