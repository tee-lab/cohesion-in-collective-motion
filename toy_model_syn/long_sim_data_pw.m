clear
close all
clc

%%

tic;

n = 10; % NUMBER OF FISH
K = [2 4]; % Values of K to be explored
no_K = numel(K);
no_it = 2; % No.of realisations

sigma_t = pi;
omega = pi/2; % angular speed
S0=0.2; % Cruise speed
sight = 360; % Visual angle in degrees

latr = 1;
gamma = 3;

Time = 1500; % Simulation time
dt = 0.05; % Integration time. Take dt such that T/dt is an interger because of approximations.
n_iter = floor(Time/dt); 
theta_tau = 0.5; % Relaxation time for angular speed
st_t = 10; % Starting time to ignore to remove the effect of initial conditions


for i = 1:no_K
    
    r_int = 1.2; % Interaction rate
    K_alg = K(i); % K as defined in the main text 
    k_alg = 1; % k as defined in the main text 
    K_atr = K(i); % K as defined in the main text
    k_atr = 1; % k as defined in the main text
%     conn_time = 10; % t_{w} as defined in the main text
    conn_time = 8;
    parfor j = 1:no_it

        [t_t, theta_t, pos_t, conncomp_size_t, avg_uni_neigh_t] = n_particles(n, r_int, ...
            dt, n_iter, sigma_t, omega, sight, K_alg, k_alg, K_atr, k_atr, ...
            S0, conn_time, st_t)
        
        pos(:,:,:,j,i) = pos_t; % Position vector
        theta(:,:,j,i) = theta_t; % Direction vector
        conncomp_size(j,i) = conncomp_size_t; % Mean network parameter
        avg_uni_neigh(j,i) = avg_uni_neigh_t; % Average no.of unique neighbours
        
    end
    
end

% Store all data in .mat file as structure
n_n = struct('pos_t', pos, 'theta_t', theta, 'sigma_t', sigma_t, 'omega', omega, 'S0', S0,...
    'K', K, 'no_it', no_it, 'st_t', st_t, 'Time', Time, 'n', n, 'dt', dt, ...
    'time_scale', conn_time, 'r_int', r_int, 'sight', sight, 'conncomp_size', conncomp_size,...
    'avg_uni_neigh', avg_uni_neigh);
save('n_pw.mat','-struct', 'n_n', '-v7.3')

disp('Simulation complete')
toc;