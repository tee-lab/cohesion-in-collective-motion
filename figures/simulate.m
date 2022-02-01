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
end_time = 2400;

toc;

%% simulation
tic 
figure(1)
% figure('units','pixels','position',[0 0 1920 1080])
% set(gcf,'color','w');
pos_t = pos_t(:,:,:,l,m); % first ith iteration and kth K. 
vel_t = vel_t(:,:,:,l,m);
theta_t = theta_t(:,:,l,m);

pos_t = pos_t(:,:,st_time:end_time);
vel_t = vel_t(:,:,st_time:end_time);
theta_t = theta_t(:,st_time:end_time);

run_time = length(vel_t);

mo = VideoWriter('n_30_at_k', 'MPEG-4');
mo.FrameRate = 10;
mo.Quality = 100;
open(mo);

for i = 1:50:run_time
    vel_x = cos(theta_t(:,i));
    vel_y = sin(theta_t(:,i));
    pos_x = pos_t(:,1,i);
    pos_y = pos_t(:,2,i);
    quiver(pos_x, pos_y, -vel_x, -vel_y, 0.55, 'LineWidth', 2.5, 'ShowArrowHead','off',...
        'Color', '#B0E0E6');
%     alpha(0.1)
    hold all
%     viscircles([pos_x, pos_y], rad_rep*ones(n,1), 'LineWidth', 2,...
%         'Color', [0.9290 0.6940 0.1250]);
    plot(pos_x, pos_y, '.', 'Color', '#00A693', 'MarkerSize', 25);
%     alpha(1)
    hold off
    
    axis('equal');
    axis([min(min(pos_t(:,1,:)))-zor, max(max(pos_t(:, 1, :)))+zor,...
        min(min(pos_t(:,2,:)))-zor, max(max(pos_t(:,2,:)))+zor])
    drawnow('limitrate');            

    image = getframe(figure(1));
    writeVideo(mo, image);

end

close(mo)

toc