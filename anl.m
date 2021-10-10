clear;
clc;
close all;
%% to load and read variables from saved mat files.
tic;

load('n_4.mat', 'pos_t', 'theta_t', 'vel_t', 'rad_rep', 'zor', 'dt', 'n');
l = 1; % iteration 
m = 1; % K
st_time = 100;

toc;

%% simulation
tic 
figure(4)
pos_t = pos_t(:,:,:,l,m); % first ith iteration and kth K. 
vel_t = vel_t(:,:,:,l,m);

pos_t = pos_t(:,:,st_time:end);
vel_t = vel_t(:,:,st_time:end);

run_time = length(vel_t); % respective time
% mo = VideoWriter('n30_K7.avi');
% mo.FrameRate = 10;
% open(mo);

for i = 1:10:run_time
    vel_x = vel_t(:,1,i);
    vel_y = vel_t(:,2,i);
    pos_x = pos_t(:,1,i);
    pos_x_c = mean(pos_x);
    pos_y = pos_t(:,2,i);
    pos_y_c = mean(pos_y);
    quiver(pos_x, pos_y, vel_x, vel_y,0.5);
%     quiver(pos_x - pos_x_c, pos_y - pos_y_c, vel_x, vel_y, 1,'ShowArrowHead', 'on', 'LineWidth',2);
    hold all
%     viscircles([pos_x - pos_x_c , pos_y - pos_y_c], rad_rep*ones(n,1));
    viscircles([pos_x, pos_y], rad_rep*ones(n,1));
    hold off
    
    axis('equal');
    axis([min(min(pos_t(:,1,:)))-zor, max(max(pos_t(:, 1, :)))+zor,...
        min(min(pos_t(:,2,:)))-zor, max(max(pos_t(:,2,:)))+zor])
    drawnow('limitrate');
    
%     image = getframe(figure(4));
%     writeVideo(mo, image);

end

% close(mo)
toc

%% INTER-DISTANCE ANALYSIS

% load('n60.mat', 'pos_t');
% pos_t = pos_t(:,:,:,l,m); % first ith iteration and kth K.
% 
% dij_min=zeros(length(pos_t),1);
% dij_mean=zeros(length(pos_t),1);
% 
% for t=1:length(pos_t)
%     dij=zeros(n);
%     for j=1:n
%         dij(:,j)=sqrt((pos_t(j,1,t)-pos_t(:,1,t)).^2+(pos_t(j,2,t)-pos_t(:,2,t)).^2);
%     end
%     dij_min(t)=min(dij(dij>0))-2*rad_rep;
%     dij_mean(t)=mean(dij(dij>0));
% end
% figure(13)
% plot(dij_min)
% hold all
% %plot(dij_mean)
% %hold all
% 
% plot(zor*ones(size(dij_min)),'--')
% plot(rad_rep*ones(size(dij_min)),'--')
% disp(min(dij_min))