close all
clear
clc

%% %% Plotting NP, CC vs K

load('n_1_np.mat', 'K', 'conncomp_size', 'n')
load('cluster_ana.mat', 'coh_coeff')

conncomp_size_temp = squeeze(conncomp_size(3,:,:));
mean_conncomp_size = mean(conncomp_size_temp,1);
std_conncomp_size = std(conncomp_size_temp,0,1);
err_conncomp_size = std_conncomp_size./sqrt(size(conncomp_size_temp,1));

errorbar(K.', mean_conncomp_size/(n), err_conncomp_size/(n),...
    'Marker', 'o', 'MarkerSize', 6, 'MarkerFaceColor', [0.8 0.8 0.8],...
    'LineStyle', '--', 'LineWidth', 2, 'CapSize', 4)
hold on
mean_coh_coeff = mean(coh_coeff,1);
std_coh_coeff = std(coh_coeff,0,1);
err_coh_coeff = std_coh_coeff./sqrt(size(coh_coeff, 1));
errorbar(K, mean_coh_coeff, err_coh_coeff, 'Marker', 'o',...
    'MarkerSize', 6, 'MarkerFaceColor', [0.8 0.8 0.8], 'LineStyle', ...
    '-.', 'LineWidth', 0.8, 'CapSize', 4);
hold off

axis([0 n 0 1])
xlabel('K', 'FontName', 'Serif', 'FontSize', 10, 'FontWeight', 'bold')
ylabel('', 'FontName', 'Serif', 'FontSize', 10, 'FontWeight', 'bold')

legend({'NP', 'CC'}, 'Location', 'southeast', 'FontName', 'serif', ...
    'FontSize',9, 'TextColor', 'k')
legend('boxoff')