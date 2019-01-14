format compact
clc, clear, close all

% name = 'PCF';


% name = 'PRPCF';

name = 'PROF';

param = 2;
tau = 50;
filename = ['SYS_',name,'_Param=',num2str(param),'_tau=',num2str(tau),'_F'];

xF = importdata([filename, '/bif_eta.txt']);
yF = importdata([filename, '/bif_extrema.txt']);

h1 = figure('Renderer', 'painters', 'Position', [10 10 1200 700]);
plot(xF, yF, 'k.', 'MarkerSize', 3)

% plot(xF(xF > range(end)), yF(xF > range(end)), '.r', 'MarkerSize', 5)
% hold on 
% for i = 0:length(range)-1
%     if (mod(i, 2) == 0)
%         plot(xF(xF < range(end-i)), yF(xF < range(end-i)), '.b', 'MarkerSize', 5)
%     else
%         plot(xF(xF < range(end-i)), yF(xF < range(end-i)), '.r', 'MarkerSize', 5)
%     end
% end


title(['Bifurcation Diagram: System \{', name, ', \tau = ', num2str(tau), ...
    ', param = ', num2str(param), '\}'], 'FontSize', 16)
xlabel('\eta', 'FontSize', 16)
ylabel('Power Output', 'FontSize', 16)
grid on
axis tight

savefig(h1, filename)
saveas(h1, filename, 'png');

close all