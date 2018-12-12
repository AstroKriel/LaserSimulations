function bif_diag(num_extrema, epsilon, folder, file)
    tic

    if (nargin < 4) || isempty(file)
        file = 'parameters.mat';
    end

    if (nargin < 3) || isempty(folder)
        [file, folder] = uigetfile;
    end

    load([folder, file]);

    start_value = PARAMS.variable_param(2);
    end_value = PARAMS.variable_param(3);
    num_values = PARAMS.variable_param(4);

    if end_value < start_value
        temp = start_value;
        start_value = end_value;
        end_value = temp;
        direction = 'backward scan';
    else
        direction = 'forward scan';
    end

    X = linspace(start_value, end_value, num_values);
    L = num_values*num_extrema;
    GAMMA = zeros(L,1);
    EXTREMA = zeros(L,1);
    index = 1;

    for inc=1:num_values
        load([folder,'sim_laser_pcf_',num2str(inc),'.mat']);
        PAST = abs(PAST(1:end,2));
        TEMP = detect_extrema(PAST,num_extrema, epsilon);
        l = length(TEMP);
        GAMMA(index:index+l-1) = X(inc);
        EXTREMA(index:index+l-1) = TEMP;
        index = index+l;
    end

    GAMMA = GAMMA(1:index-1);
    EXTREMA = EXTREMA(1:index-1);

    h = figure;
    hold on
    plot(GAMMA, EXTREMA,'.k','MarkerSize',1)

    grid on
    axis tight
    ylabel('Output power (a. u.)')
    title(direction)

    saveas(h,[folder,'bif_diag.fig'],'fig');
    saveas(h,[folder,'bif_diag.png'],'png');

    toc
end
