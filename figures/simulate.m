clear;
clc;
close all;
%% to load and read variables from saved mat files.
tic;

% Load data
load('n_pw.mat', 'pos_t', 'theta_t', 'vel_t', 'rad_rep', 'zor', 'dt', 'n');
l = 1; % iteration 
m = 1; % K
st_time = 500; % start time
end_time = 8000;

toc;

%% simulation
tic 
figure(1)
pos_t = pos_t(:,:,:,l,m); % first ith iteration and kth K. 
vel_t = vel_t(:,:,:,l,m);
theta_t = theta_t(:,:,l,m);

pos_t = pos_t(:,:,st_time:end_time);
vel_t = vel_t(:,:,st_time:end_time);
theta_t = theta_t(:,st_time:end_time);

run_time = length(vel_t);

for i = 1:40:run_time
    vel_x = cos(theta_t(:,i));
    vel_y = sin(theta_t(:,i));
    pos_x = pos_t(:,1,i);
    pos_y = pos_t(:,2,i);
    quiver(pos_x, pos_y, -vel_x, -vel_y, 0.35, 'LineWidth', 2.5, 'ShowArrowHead','off',...
        'Color', '#B0E0E6');
    hold all
    plot(pos_x, pos_y, '.', 'Color', '#00A693', 'MarkerSize', 25);
    hold off
    
    axis('equal');
    axis([min(min(pos_t(:,1,:)))-zor, max(max(pos_t(:, 1, :)))+zor,...
        min(min(pos_t(:,2,:)))-zor, max(max(pos_t(:,2,:)))+zor])
    drawnow('limitrate');            

end

toc