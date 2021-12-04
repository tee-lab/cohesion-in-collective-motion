clc
clear
close all
%%

all_marks = {'o','+','*','x','s','d','^','v','>','<','p','h'};
col=[0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840];

newfolder = sprintf('/Users/vivekj_19/OneDrive - Indian Institute of Science/IISc/MS Thesis/MS/commented_codes/figure_2');
cd(newfolder)
load('cluster_ana_pw.mat', 'coh_coeff', 'K')
leg = zeros(1,2); % For legend

% Your data...
x = K;
y = ((mean(coh_coeff,1)));
% Now make the 'line' (actually a surface)...
z = zeros(size(x));
S = surface([x;x],[y;y],[z;z],...
    'facecol','no',...
    'edgecol','interp',...
    'linew',3,...
    'edgealpha',.2,...
    'edgecolor',col(1,:));
hold all

leg(1) = scatter(K,mean(coh_coeff,1),'filled','SizeData',100);
hold all
alpha(.6)

newfolder = sprintf('/Users/vivekj_19/OneDrive - Indian Institute of Science/IISc/MS Thesis/MS/commented_codes/figure_4');
cd(newfolder)
load('cluster_ana_at.mat', 'coh_coeff', 'K')

% Your data...
x = K;
y = ((mean(coh_coeff,1)));
% Now make the 'line' (actually a surface)...
z = zeros(size(x));
S = surface([x;x],[y;y],[z;z],...
    'facecol','no',...
    'edgecol','interp',...
    'linew',3,...
    'edgealpha',.2,...
    'edgecolor',col(2,:));
hold all

leg(2) = scatter(K,mean(coh_coeff,1),'filled','SizeData',100);
hold all
alpha(.6)

xlabel('$K$','Interpreter','latex')
ylabel('$\mathcal{C}$','Interpreter','latex')
yticks([0 0.25 0.5 0.75 1])
xticks([0 10 20 30 40 50])
lgd = legend(leg, 'Stoc Pairwise Model', 'Averaging Type Model');
legend('boxoff')
lgd.FontSize = 8;
lgd.Location = 'southeast';