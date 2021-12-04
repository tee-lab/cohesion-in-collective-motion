clear;
clc;
close all;
%% to load and read variables from saved mat files.
tic;

% Load data
load('n_at.mat', 'pos_t', 'theta_t', 'vel_t', 'rad_rep', 'zor', 'dt', 'n');
l = 2; % iteration 
m = 1; % K
st_time = 100; % start time

toc;

%% simulation
tic 
figure(4)
pos_t = pos_t(:,:,:,l,m); % first ith iteration and kth K. 
vel_t = vel_t(:,:,:,l,m);

pos_t = pos_t(:,:,st_time:end);
vel_t = vel_t(:,:,st_time:end);

run_time = length(vel_t);

for i = 1:20:run_time
    vel_x = vel_t(:,1,i);
    vel_y = vel_t(:,2,i);
    pos_x = pos_t(:,1,i);
    pos_y = pos_t(:,2,i);
    quiver(pos_x, pos_y, vel_x, vel_y,0.5);
    hold all
    viscircles([pos_x, pos_y], rad_rep*ones(n,1));
    hold off
    
    axis('equal');
    axis([min(min(pos_t(:,1,:)))-zor, max(max(pos_t(:, 1, :)))+zor,...
        min(min(pos_t(:,2,:)))-zor, max(max(pos_t(:,2,:)))+zor])
    drawnow('limitrate');

end

toc