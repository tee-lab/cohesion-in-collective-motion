close all
clear
clc

%% Plotting CC and avg.conn on the same plot

load('n_1_np.mat', 'avg_conn', 'K', 'no_it', 'zor', 'rad_rep', ...
    'Time','n','time_scale', 'r_align', 'conncomp_size')
load('cluster_ana.mat', 'coh_coeff')

figure(1)

% for t_s = 1:numel(time_scale)
%     
%     avg_conn_temp = squeeze(avg_conn(t_s,:,:));
%     mean_conn = mean(avg_conn_temp,1);
%     std_conn = std(avg_conn_temp,0,1);
%     err_conn = std_conn./sqrt(size(avg_conn_temp,1));
%     
%     errorbar(K.', mean_conn/(n-1), err_conn/(n-1), 'Marker', 'o', 'MarkerSize', 6, ...
%         'MarkerFaceColor', [0.8 0.8 0.8], 'LineStyle', ':', 'LineWidth', 2, 'CapSize', 4)
%     hold on
%     
% end

for t_s = 1:numel(time_scale)
    
    conncomp_size_temp = squeeze(conncomp_size(t_s,:,:));
    mean_conncomp_size = mean(conncomp_size_temp,1);
    std_conncomp_size = std(conncomp_size_temp,0,1);
    err_conncomp_size = std_conncomp_size./sqrt(size(conncomp_size_temp,1));
    
    errorbar(K.', mean_conncomp_size/(n), err_conncomp_size/(n),...
        'Marker', 'o', 'MarkerSize', 6, ...
        'MarkerFaceColor', [0.8 0.8 0.8], 'LineStyle', ':', 'LineWidth',...
        2, 'CapSize', 4)
    hold on
    
end

axis([0 n 0 1])
xlabel('K', 'FontName', 'Serif', 'FontSize', 10, 'FontWeight', 'bold')
ylabel('Avg agents in fully conn net', 'FontName', 'Serif', 'FontSize', 10, 'FontWeight', 'bold')

yyaxis right


mean_coh_coeff = mean(coh_coeff,1);
std_coh_coeff = std(coh_coeff,0,1);
err_coh_coeff = std_coh_coeff./sqrt(size(coh_coeff, 1));
errorbar(K, mean_coh_coeff, err_coh_coeff, 'Marker', 'o', 'MarkerSize', 6, ...
    'MarkerFaceColor', [0.8 0.8 0.8], 'LineStyle', ':', 'LineWidth', 0.8, 'CapSize', 4);
hold on
    
axis([0 n 0 1])
ylabel('Cohesion Coefficient', 'FontName', 'TimesRoman', 'FontSize', 10, 'FontWeight', 'bold')

ax = gca;
ax.LineWidth = 1.5;

legend({'ts = 1', 'ts = 5', 'ts = 10', 'ts = 15',  'ts = 50',  'ts = 100',...
    'CC'}, 'FontName', 'Serif', 'FontSize', 8, ...
    'FontWeight', 'bold', 'Location', 'best')

% title('No.of Comps (n = 30)', 'FontName', 'Serif', 'FontSize', 14, 'FontWeight', 'bold')

% print -depsc2 path_conn_n2.eps
print -depsc2 np_n1.eps
       
toc