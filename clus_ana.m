close all
clear
clc

%%

tic;
var = struct();
load('n_1_np.mat');

st_t = 0;
n_iter = round((Time - st_t)/0.1);
srt_p = st_t/0.1;

no_K = numel(K);
no_ra = numel(r_atr);

theta_t = theta_t(:, srt_p+1:end,:,:,:);
pos_t = pos_t(:,:,srt_p+1:end,:,:,:);
vel_t = vel_t(:,:,srt_p+1:end,:,:,:);

%% Cluster Analysis

num_clus = zeros(n_iter, no_it, no_K, no_ra);
coh_coeff = zeros(n_iter, no_it, no_K, no_ra);
num_in_clus = zeros(n_iter, no_it, no_K, no_ra);
avg_num_in_clus = zeros(no_it, no_K, no_ra);
noise_clus = zeros(n_iter, no_it, no_K, no_ra);
noise_coeff = zeros(no_it, no_K, no_ra);
m1_size_large_clus = zeros(n_iter, no_it, no_K, no_ra);
m2_size_other_clus = zeros(n_iter, no_it, no_K, no_ra);
pol = zeros(n_iter, no_it, no_K, no_ra);
avg_clus_size = zeros(n_iter, no_it, no_K, no_ra);

for nro = 1:no_ra
    
    for i = 1:no_K
        
        parfor k = 1:no_it
            
            pos_temp = pos_t(:,:,:,k,i,nro);
            theta_temp = theta_t(:,:,k,i,nro);
            
            for j = 1:n_iter
                
                x = pos_temp(:,1,j);
                y = pos_temp(:,2,j);
                
                ids = DBSCAN([x y], 2*zor, 1);
                num_clus(j,k,i,nro) = max(ids);
                size_large_clus = 0;
                ids_lc = 0;
                avg_clus_size_temp = 0;
                
                for m1 = 1:max(ids)
                    size_large_clus_temp = numel(ids(ids == m1));
                    avg_clus_size_temp = avg_clus_size_temp + size_large_clus_temp;
                    if size_large_clus_temp > size_large_clus
                        size_large_clus = size_large_clus_temp;
                        ids_lc = m1;
                    end
                end
                
                m1_size_large_clus(j,k,i,nro) = size_large_clus;
                avg_clus_size(j,k,i,nro) = avg_clus_size_temp/max(ids);
                coh_coeff(j,k,i,nro) = size_large_clus/n;
                
                if max(ids) == 1
                    num_in_clus(j,k,i,nro) = numel(ids(ids == 1));
                else
                    num_in_clus(j,k,i,nro) = NaN;
                end
                
                noise_clus(j,k,i,nro) = numel(ids(ids == 0));
                m2_size_other_clus(j,k,i,nro) = n - m1_size_large_clus(j,k,i,nro) - noise_clus(j,k,i,nro);
                
                pol_temp_x = mean(cos(theta_temp(:,j)));
                pol_temp_y = mean(sin(theta_temp(:,j)));
                pol_temp = sqrt(pol_temp_x^2 + pol_temp_y^2);
                pol(j,k,i,nro) = pol_temp;
            end
        end
    end
    
end

for nro = 1:no_ra
    for i = 1:no_K
        for j = 1:no_it
%             coh_coeff(j,i,nro) = (sum(num_clus(:,j,i,nro) == 1))/n_iter;
            noise_coeff(j,i,nro) = (sum(noise_clus(:,j,i,nro)))/n_iter;
            num_in_clus_temp = num_in_clus(:,j,i,nro);
            num_in_clus_temp = num_in_clus_temp(~isnan(num_in_clus_temp));
            avg_num_in_clus(j,i,nro) = mean(num_in_clus_temp);
        end
    end
end

coh_coeff = reshape(coh_coeff, n_iter*no_it, no_K, no_ra);
num_clus = reshape(num_clus, n_iter*no_it, no_K, no_ra);
noise_clus = reshape(noise_clus, n_iter*no_it, no_K, no_ra);
num_in_clus = reshape(num_in_clus, n_iter*no_it, no_K, no_ra);
avg_clus_size = reshape(avg_clus_size, n_iter*no_it, no_K, no_ra);
m1_size_large_clus = reshape(m1_size_large_clus, n_iter*no_it, no_K, no_ra);
m2_size_other_clus = reshape(m2_size_other_clus, n_iter*no_it, no_K, no_ra);
pol = reshape(pol, n_iter*no_it, no_K, no_ra);

clusters = struct('num_clus', num_clus, 'coh_coeff', coh_coeff, 'pol', pol, 'noise_clus', noise_clus,...
    'noise_coeff', noise_coeff, 'num_in_clus', num_in_clus, 'avg_num_in_clus', avg_num_in_clus, ...
    'm1_size_large_clus', m1_size_large_clus, 'm2_size_other_clus', m2_size_other_clus,...
    'avg_clus_size', avg_clus_size, 'K', K, 'r_align', r_align, 'r_atr', r_atr, 'n', n);
save('cluster_ana.mat', '-struct', 'clusters', '-v7.3')

toc;