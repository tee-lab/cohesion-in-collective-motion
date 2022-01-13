clear;
clc;
close all;
%% to load and read variables from saved mat files.
tic;

% Load data
load('n_pw.mat', 'pos_t', 'theta_t', 'dt', 'n');
l = 1; % iteration 
m = 1; % K
st_time = 100; % start time

toc;

%% simulation
tic 
figure(4)
pos_t = pos_t(:,:,:,l,m); % first ith iteration and kth K.
theta_t = theta_t(:,:,l,m);

pos_t = pos_t(:,:,st_time:end);
theta_t = theta_t(:,st_time:end);

run_time = length(pos_t);

for i = 1:40:run_time
    vel_x = cos(theta_t(:,i));
    vel_y = sin(theta_t(:,i));
    pos_x = pos_t(:,1,i);
    pos_y = pos_t(:,2,i);
    quiver(pos_x, pos_y, vel_x, vel_y,0.5);
    hold all
    viscircles([pos_x, pos_y], 0.2*ones(n,1));
    hold off
    
    axis('equal');
    axis([min(min(pos_t(:,1,:))), max(max(pos_t(:, 1, :))),...
        min(min(pos_t(:,2,:))), max(max(pos_t(:,2,:)))])
    drawnow('limitrate');

end

toc