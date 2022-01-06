close all
clear
clc

%%

tic;

load('n_at.mat'); % Load the file that has all the simulation data

n_iter = round((Time - st_t)/0.1); % No.of data points in a simulation after removing data till st_t
srt_p = st_t/0.1;

no_K = numel(K); % Length of K values explored

% As defined in long_sim_data.m
theta_t = theta_t(:,srt_p+1:end,:,:);
pos_t = pos_t(:,:,srt_p+1:end,:,:);
vel_t = vel_t(:,:,srt_p+1:end,:,:);

%% Cluster Analysis

num_clus = zeros(n_iter, no_it, no_K); % Stores no.of clusters
coh_coeff = zeros(n_iter, no_it, no_K); % Cohesion parameter
noise_clus = zeros(n_iter, no_it, no_K); % No.of lone agents (refer SI)
pol = zeros(n_iter, no_it, no_K); % Group Polarisation 
avg_clus_size = zeros(n_iter, no_it, no_K); % Average cluster size

for i = 1:no_K % Run over all K's

    parfor k = 1:no_it % Run over all realisations

        pos_temp = pos_t(:,:,:,k,i); % Position vector for given realisation and K
        theta_temp = theta_t(:,:,k,i); % Heading angle for given realisation and K

        for j = 1:n_iter % Run over a given simulation

            x = pos_temp(:,1,j); % x coordinates of all agents
            y = pos_temp(:,2,j); % y coordinates of all agents

            ids = DBSCAN([x y], 2*zor, 1); % IDs of clusters. Refer main text for the definition of a cluster
            num_clus(j,k,i) = max(ids); % No.of clusters formed
            size_large_clus = 0; % Size of the largest cluster
            avg_clus_size_temp = 0; % Average size of the cluster

            % This loop calculates largest cluster and average cluster size
            for m1 = 1:max(ids) 
                size_large_clus_temp = numel(ids(ids == m1));
                avg_clus_size_temp = avg_clus_size_temp + size_large_clus_temp;
                if size_large_clus_temp > size_large_clus
                    size_large_clus = size_large_clus_temp;
                end
            end

            avg_clus_size(j,k,i) = avg_clus_size_temp/max(ids);
            coh_coeff(j,k,i) = size_large_clus/n; % Cohesion parameter

            noise_clus(j,k,i) = numel(ids(ids == 0)); % no.of lone agents
            
            % Calculates group polarisation (Refer main text for definition)
            pol_temp_x = mean(cos(theta_temp(:,j)));
            pol_temp_y = mean(sin(theta_temp(:,j)));
            pol_temp = sqrt(pol_temp_x^2 + pol_temp_y^2);
            pol(j,k,i) = pol_temp;
        end
    end
end

% Merge all cohesion parameter values for all realisations for given K
coh_coeff = reshape(coh_coeff, n_iter*no_it, no_K);
num_clus = reshape(num_clus, n_iter*no_it, no_K);
noise_clus = reshape(noise_clus, n_iter*no_it, no_K);
avg_clus_size = reshape(avg_clus_size, n_iter*no_it, no_K);
pol = reshape(pol, n_iter*no_it, no_K);

clusters = struct('num_clus', num_clus, 'coh_coeff', coh_coeff, 'pol', pol, ...
    'noise_clus', noise_clus, 'avg_clus_size', avg_clus_size, 'K', K, 'n', n);
save('cluster_ana_at.mat', '-struct', 'clusters', '-v7.3')

toc;