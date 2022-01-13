close all
clear
clc

%% Alternate measure of cohesion

tic;

load('n_pw.mat'); % Load the file that has all the simulation data

n_iter = round((Time - st_t)/0.1); % No.of data points in a simulation after removing data till st_t
srt_p = st_t/0.1;

no_K = numel(K); % Length of K values explored

% As defined in long_sim_data.m
theta_t = theta_t(:,srt_p+1:end,:,:);
pos_t = pos_t(:,:,srt_p+1:end,:,:);


%% Mean distance between agents and Expanse

mean_neigh_dist = zeros(n_iter, no_it, no_K); % mean distance between agents
expanse = zeros(n_iter, no_it, no_K); % Expanse, as defined by Huth (1992)

for k = 1:no_K

    parfor i = 1:no_it

        pos_temp = pos_t(:,:,:,i,k); % Position vector for given realisation and K

        for t = 1:n_iter
            
            grp_cent = mean(pos_temp(:,:,t)); % Calculate centroid
            dist_gc = pos_temp(:,:,t) - grp_cent; % Relative position to centroid
            dist_gc = sqrt(dist_gc(:,1).^2 + dist_gc(:,2).^2); % Distance of each agent from centroid
            expanse(t,i,k) = mean(dist_gc); % Mean distance of agents from centroid
            dist_temp = squareform(pdist(pos_temp(:,:,t))); % Distance between agents
            mean_neigh_dist(t,i,k) = mean(sum(dist_temp, 2)/(n-1)); % Mean distance btw agents
            
        end

    end

end

% Merge all cohesion parameter values for all realisations for given K
mean_neigh_dist = reshape(mean_neigh_dist, n_iter*no_it, no_K);
expanse = reshape(expanse, n_iter*no_it, no_K);

neigh_dist = struct('mean_neigh_dist', mean_neigh_dist, 'expanse', expanse);
save('alt_coh_msr.mat', '-struct', 'neigh_dist', '-v7.3')

toc;