clear
close all
clc

%%

tic;

n = 10; % NUMBER OF FISH
K = 1:2:9;
no_K = numel(K);
no_it = 2;

sight = 270; % In degrees

rad_rep = 0.2;  
zor = 1.2;
omega = 2*pi/3;

S0=0.2;
k_r = -1*1e-3;
beta = 3; % power of the distance dependence

attr_c = 3*1e-2;
gamma = 3; latr=zor/2; Smax=5*0.2;

Time = 100;
dt = 0.05; n_iter = ceil(Time/dt); tau = 0.2; theta_tau = 0.5;

for i = 1:no_K
    
    r_spon = 1.0; sigma_t = (sight/2)*pi/180; % sigma is the width of the distribution
    r_align = 1.5; K_alg = K(i); k_alg = 1;
    r_atr = 1.0; K_atr = K(i); k_atr = 1;
    conn_time = [1 5 10 15 50 100];
    
    parfor j = 1:no_it
        
        [t_t, theta_t, pos_t, v, vel_t, v_mag, avg_conn_t, conncomp_size_t] = n_particles(n, r_spon, ...
            r_align, r_atr, dt, n_iter, zor, rad_rep, tau, theta_tau, K_alg, k_alg, K_atr, k_atr, k_r,...
            beta, sight, gamma, attr_c, latr, S0, Smax, conn_time)
        
        pos(:,:,:,j,i) = pos_t;
        theta(:,:,j,i) = theta_t;
        vel(:,:,:,j,i) = vel_t;
        avg_conn(:,j,i) = avg_conn_t;
        conncomp_size(:,j,i) = conncomp_size_t;
        
    end
    
end

n_n = struct('pos_t', pos, 'theta_t', theta, 'vel_t', vel, 'sight', sight, 'S0', S0,...
    'K', K, 'no_it', no_it, 'zor', zor, 'rad_rep', rad_rep, ...
    'Time', Time, 'n', n, 'dt', dt, 'time_scale', conn_time, 'avg_conn', avg_conn,...
    'r_align', r_align, 'r_atr', r_atr, 'r_spon', r_spon, 'conncomp_size', conncomp_size);
save('n_1_np.mat','-struct', 'n_n', '-v7.3')

disp('Simulation complete')
toc;