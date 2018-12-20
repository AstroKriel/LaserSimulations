function OUTPUT = integ_sim_laser(params_laser, params_temp, INPUT, DIM)
    theta = params_laser(8);
    h = params_temp(1);
    horizon = params_temp(2);
    delay = floor(theta/h);
    steps = floor(horizon/h);

    if (nargin<3) || isempty(INPUT)
        PAST = 0.1+zeros(delay+1, DIM) + 1e-5;
    elseif length(INPUT) < delay+2
        PAST = vertcat(zeros(delay+1-size(INPUT,1),DIM) + 1e-5, INPUT);
    else
        PAST = INPUT(end-delay:end,:);
    end

    PAST = reshape(PAST.', DIM*(delay+1), 1);
    
    switch DIM
        case 4
            RESULTS = fast_sim_laser_prof(params_laser, params_temp, PAST);
        
        case 7
            RESULTS = fast_sim_laser_prpcf(params_laser, params_temp, PAST);
        case 5
            RESULTS = fast_sim_laser_pcf(params_laser, params_temp, PAST);
    end
        
    RESULTS = RESULTS(DIM*delay + 1:end);
    OUTPUT = reshape(RESULTS, DIM, steps+1).';
end
