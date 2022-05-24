clear
close all
clc

%%

tic;

n = 10; % NUMBER OF FISH
K = [2 4]; % Values of K to be explored
no_K = numel(K);
no_it = 2; % No.of realisations

sight = 270; % Sight as defined in the supplementary text. Units are mentioned that text.

rad_rep = 0.2; % Agent size 
zor = 1.2; % Zone of repulsion

S0=0.2; % Cruise speed
k_r = -1*1e-3;
beta = 3; % power of the distance dependence

attr_c = 3*1e-2;
gamma = 3; latr=zor/2; 
Smax=5*0.2; % Maximum speed an agent can travel at. 

Time = 1000; % Simulation time
dt = 0.05; % Integration time. Take dt such that T/dt is an interger because of approximations. 
n_iter = floor(Time/dt); 
tau = 0.2; % Relaxation time for speed
theta_tau = 0.5; % Relaxation time for angular speed
st_t = 10; % Starting time to ignore to remove the effect of initial conditions


for i = 1:no_K
    
    r_spon = 1.0; % Spontaneous interaction rate
    sigma_t = (sight/2)*pi/180; % sigma is the width of the distribution
    r_align = 1.5; % Rate of alignment
    K_alg = K(i); % K as defined in the main text 
    k_alg = 1; % k as defined in the main text
    r_atr = 1; % Rate of attraction
    K_atr = K(i); % K as defined in the main text
    k_atr = 1; % k as defined in the main text

    parfor j = 1:no_it

        [t_t, theta_t, pos_t, v, vel_t, v_mag] = n_particles(n, r_spon, ...
            r_align, r_atr, dt, n_iter, zor, rad_rep, tau, theta_tau, K_alg, k_alg, K_atr, k_atr, k_r,...
            beta, sight, gamma, attr_c, latr, S0, Smax)
        
        pos(:,:,:,j,i) = pos_t; % Position vector
        theta(:,:,j,i) = theta_t; % Direction vector
        vel(:,:,:,j,i) = vel_t; % velocity vector
        
    end
    
end

% Store all data in .mat file as structure
n_n = struct('pos_t', pos, 'theta_t', theta, 'vel_t', vel, 'sight', sight, 'S0', S0,...
    'K', K, 'no_it', no_it, 'zor', zor, 'rad_rep', rad_rep, 'st_t', st_t, ...
    'Time', Time, 'n', n, 'dt', dt, ...
    'r_align', r_align, 'r_atr', r_atr, 'r_spon', r_spon);
save('n_pw.mat','-struct', 'n_n', '-v7.3')

disp('Simulation complete')
toc;