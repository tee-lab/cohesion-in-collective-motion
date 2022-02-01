clc
clear
close all
%%
load('n_pw.mat', 'conncomp_size', 'K', 'n') % Load data files
load('cluster_ana_pw.mat')

all_marks = {'o','+','*','x','s','d','^','v','>','<','p','h'};
col=[0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840];


% Load datas saved using the variable names used in n_pw.mat and
% cluster_ana.mat
N = n; % Group Size
cc = coh_coeff;
cs = avg_clus_size;
nc = num_clus;
np = conncomp_size;

% Plot for K vs Cohesion parameter
figure(1)
% Your data...
x = K;
y = ((mean(cc,1)));
% Now make the 'line' (actually a surface)...
z = zeros(size(x));
S = surface([x;x],[y;y],[z;z],...
    'facecol','no',...
    'edgecol','interp',...
    'linew',3,...
    'edgealpha',.2,...
    'edgecolor',col(1,:));
hold all

scatter(K,mean(cc,1),'filled','SizeData',100);
hold all
alpha(.6)

% Plot for K/N vs Cohesion parameter
figure(2)
x = K/N;
y = ((mean(cc,1)));
% Now make the 'line' (actually a surface)...
z = zeros(size(x));
S = surface([x;x],[y;y],[z;z],...
    'facecol','no',...
    'edgecol','interp',...
    'linew',2,...
    'edgealpha',.2,...
    'edgecolor',col(1,:));
hold all

scatter(K/N,mean(cc,1),'filled','SizeData',50)
hold all
alpha(.6)

% Plot for number of clusters
figure(3)
x = K;
y = mean(nc,1);
% Now make the 'line' (actually a surface)...
z = zeros(size(x));
S = surface([x;x],[y;y],[z;z],...
    'facecol','no',...
    'edgecol','interp',...
    'linew',3,...
    'edgealpha',.2,...
    'edgecolor',col(1,:));
hold all

scatter(K,mean(nc,1),'filled','SizeData',100)
hold all
alpha(.6)

figure(4)
% Your data...
x = K;
y = mean(cc,1);
% Now make the 'line' (actually a surface)...
z = zeros(size(x));
S1 = surface([x;x],[y;y],[z;z],...
    'facecol','no',...
    'edgecol','interp',...
    'linew',3,...
    'edgealpha',.2,...
    'edgecolor',col(1,:));
hold all

sc(1) = scatter(K,mean(cc,1),'filled','SizeData',100);
hold all
alpha(.6)


% Plot for network parameter and cohesion parameter vs K
x = K;
y = mean(np,1)/N;
% Now make the 'line' (actually a surface)...
z = zeros(size(x));
S2 = surface([x;x],[y;y],[z;z],...
    'facecol','no',...
    'edgecol','interp',...
    'linew',3,...
    'edgealpha',.2,...
    'edgecolor',col(2,:));
hold all

sc(2) = scatter(K,mean(np,1)/N,'filled','SizeData',100);
hold all
alpha(.6)


figure(1)
xlabel('$\textbf{K}$','Interpreter','latex', 'FontSize', 14)
ylabel('$\mathcal{C}$','Interpreter','latex', 'FontSize', 14)
axis([0 n 0 1])
yticks([0 0.25 0.5 0.75 1])
xticks([0 10 20 30 40 50])

figure(2)
xlabel('$\frac{K}{N}$','Interpreter','latex')
ylabel('$\mathcal{C}$','Interpreter','latex')
yticks([0 0.5 1])
xticks([0 0.5 1])

figure(3)
xlabel('$K$','Interpreter','latex')
ylabel('Number of clusters')
xticks([0 10 20 30 40 50])

figure(4)
xlabel('$\textbf{K}$','Interpreter','latex', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('$\mathcal{C},~ \mathcal{N}_p$','Interpreter','latex', 'FontSize', 14, 'FontWeight', 'bold')
axis([0 N 0 1])
yticks([0 0.25 0.5 0.75 1])
xticks([0 10 20 30 40 50])
lgd = legend(sc, '$\mathcal{C}$', ' $\mathcal{N}_p$');
legend('boxoff')
lgd.FontSize = 8;
lgd.Location = 'northeast';
lgd.InterpreterMode = 'manual';
lgd.Interpreter = 'latex';