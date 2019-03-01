format compact
clc, clear, close all

name_sys = 'PCF';
name_sys = 'PROFN';
% name_sys    = 'PRPCFN';

param_set_name = 'bifpoint_analysis';

for tau_R = [1]
    for theta = [7000]
%         folder = ['Param_', param_set_name, '/SYS_', name_sys, '_theta=', num2str(theta)];
%         folder = ['Param_', param_set_name, '/SYS_', name_sys, '_alpha=', num2str(alpha)];
%         folder = ['Param_', param_set_name, '/SYS_', name_sys, '_kappa=', num2str(ka)];

        folder = ['Param_', param_set_name, '/SYS_', name_sys, '_theta=', num2str(theta), '_tau=', num2str(tau_R)];
        
        xF = importdata([folder, '/bif_eta.txt']);
        
        for i = 2:length(xF)
            if (xF(i) == xF(i-1))
                disp([num2str(xF(i))])
                break
            end
        end
    end
    disp(' ')
end

