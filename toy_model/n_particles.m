% Date: 26-10-2021
% Edited code with comments
% Written by Vivek J and Danny Raj M

function [t_t, theta_t, pos_t, conncomp_size_t, avg_uni_neigh_t] = n_particles(n, r_spon, ...
    r_align, r_atr, dt, n_iter, latr, gamma, sigma_t, K_alg, k_alg, K_atr, k_atr, ...
    S0, conn_time, st_t)
 
% INITIAL CONDITIONS - Putting agents in a lattice with orientation around
% a random angle.

m=ceil(sqrt(n));
pos(:,1)=(rem((1:n)-1,m)); % x coordinates 
pos(:,2)= floor(((1:n)-1)/m); % y coordinates
theta = 2*pi*rand(1) + pi*randn(n,1); % Distribution of heading angles around a 
% randomly chosen angle

% Checking if all the heading angles are between [0-2pi] 
for j=1:n
    if theta(j)>2*pi
        theta(j)=theta(j)-2*pi;
    elseif theta(j)<0
        theta(j)=2*pi+theta(j);
    end
end

v0 = S0; % speed

% Truncated Distributions of angles for spontaneous reaction

nor_dis_ang = makedist('Normal', 'sigma', sigma_t);
trunc_dis_ang = truncate(nor_dis_ang, -pi, pi); 

% Defining speed, velocity for n agents

vel = [v0*cos(theta), v0*sin(theta)];

d_t = vel; % Desired heading direction 
theta_d = theta; % Desired heading angle 

sk_t = 0.1/dt; % Store data after sk_t time steps

% position, velocity, speed and direction at time t-1

pos_t_1 = pos;
theta_t_1 = theta;

pos_t = zeros(n, 2, ceil(n_iter/sk_t)); % Stores positions of agents over the simulation
theta_t = zeros(n, ceil(n_iter/sk_t)); % Stores heading angles of agents over the simulation

t_t = 0.1:0.1:n_iter*dt; % Time points at which data is stored 

% Reaction rates
te_spon = (1/r_spon) * log(1./rand(n,1)); % time at which spontaneous reaction happens
te_align = (1/r_align) * log(1./rand(n,1)); % time at which alignment interaction happens
te_atr = (1/r_atr) * log(1./rand(n,1)); % time at which attraction interaction happens

T = n_iter*dt; 

% For time window t_{w} what is the size of the array to store max
% size of connected cluster. Connection times 
% (time interval over which networks are constructed) are defined in
% long_sim_data.m file. 
conncomp_size_1 = zeros(floor(T/conn_time),1);

% Average no.of unique neighbour an agent interacts with in the time
% interval t_{w}
avg_uni_neigh_1 = zeros(floor(T/conn_time),1);
 
% Arrays to store which agent interacts with which agents with conn_time
% interval 
connect_agents_1 = [];
agents_1 = [];

tic;

% Start of the simulation 

for t = 2:n_iter
    
    % Displays time elapsed
    if rem(t,1e4)==0
        disp(t*dt)
    end
    
    % Constructs networks over conn_time time interval only time elapsed is
    % greater than st_t (to account for initial conditions)
    if rem(t,round(conn_time*(1/dt))) == 0 
        
        con_graph = avg_connections(agents_1, connect_agents_1, n);
%         ids = dbscan(pos, 2*zor, 1);
        
        connect_agents_1 = [];
        agents_1 = [];
        [~, clussize] = conncomp(con_graph);
        if isempty(clussize) == 0
            conncomp_size_1((t/(conn_time*(1/dt))),1) = max(clussize);
            avg_uni_neigh_1((t/(conn_time*(1/dt))),1) = mean(outdegree(con_graph));
        end

        figure(1)
        
%         plot(con_graph, 'XData', pos(:,1), 'YData', pos(:,2))
%         hold all
%         sc = gscatter(pos(:,1),  pos(:,2), net_bin);
%         hold all
%         xlabel('X', 'FontName', 'Serif', 'FontSize', 10, 'FontWeight', 'bold')
%         ylabel('Y', 'FontName', 'Serif', 'FontSize', 10, 'FontWeight', 'bold')
%         legend(sc)

    end
    
    for i = 1:n
        
        e_spon = te_spon(i) <= dt; % Check if spontaneous reaction happens at this time step
        e_align = te_align(i) <= dt; % Check if alignment interaction happens at this time step
        e_atr = te_atr(i) <= dt; % Check if attraction interaction happens at this time step
        te_spon(i) = te_spon(i) - dt;
        te_align(i) = te_align(i) - dt;
        te_atr(i) = te_atr(i) - dt;
        
        % If the interaction happened at this time step, chose a next time
        % point where the interaction will happen.
        if te_spon(i) < 0
            te_spon(i) = (1/r_spon) * log(1/rand());
        end
        
        if te_align(i) < 0
            te_align(i) = (1/r_align) * log(1/rand());
        end
        
        if te_atr(i) < 0
            te_atr(i) = (1/r_atr) * log(1/rand());
        end
        
        dis_vec = pos_t_1 - repmat(pos_t_1(i,:), size(pos_t_1,1), 1); % rij vector (as defined in the main text)
        mag_vec = sqrt(dis_vec(:,1).^2 + dis_vec(:,2).^2); % magnitude of rij
        mag_vec(i) = Inf;
        
        % Check if any event is happening in this time interval

        %if no event, continue to move in the desired direction of previous step

        if (e_spon+e_align+e_atr) == 0
            
            theta_d(i) = theta_d(i);
            
        else % in case of an event 
            
            % SPONTANEOUS REACTION

            % Chose a random angle from truncated normal dist 
            delta = random(trunc_dis_ang);
            
            theta_d_s = theta_t_1(i) + delta; % Desired direction due to spontaneous rotation

            % Ensure that the angles are with [0, 2pi]
            if theta_d_s >= 2*pi
                theta_d_s = theta_d_s - 2 * pi;
            elseif theta_d_s < 0
                theta_d_s = theta_d_s + 2 * pi;
            end
            
            % ALIGNMENT 
            
            % sorting agents in the ascending order of their
            % distance from focal individual
            [~, s_ind_st] = sort(mag_vec, 'ascend');

            neighbours_kalg=s_ind_st(randperm(K_alg, k_alg)); %Picking a random neighbor to align
            d_align_t = mean([cos(theta_t_1(neighbours_kalg,1)) sin(theta_t_1(neighbours_kalg,1))],1); % Desired direction
            theta_d_a = atan2(d_align_t(1,2),d_align_t(1,1)); % Desired direction
            
            % ATTRACTION

            % sorting agents in the ascending order of their
            % distance from focal individual
            [~, s_ind_atr] = sort(mag_vec, 'ascend');

            %Picking a random neighbor to attract
            neighbours_katr = s_ind_atr(randperm(K_atr, k_atr));
            % Desired direction
            d_atr_t = mean([(((mag_vec(neighbours_katr,1))/latr).^gamma)...
                    .*dis_vec(neighbours_katr,:)./(norm(dis_vec(neighbours_katr,:))+eps) ; [cos(theta_t_1(i)) sin(theta_t_1(i))]],1);
%             d_atr_t = dis_vec(neighbours_katr,:)./(norm(dis_vec(neighbours_katr,:))+eps);
            theta_d_atr = atan2(d_atr_t(1,2),d_atr_t(1,1));

            %  If attraction interaction has happened and time is
            %  greater than st_t then record which agent interacts with
            %  which agent
            if e_atr == 1 && dt*t > st_t
                nc = ones(1,length(neighbours_katr));
                agents_1 = [agents_1, i*nc];
                connect_agents_1 = [connect_agents_1, neighbours_katr.'];
            end
            
            % Take the average of all the interactions that happended in
            % this time interval dt
            theta_d(i) = atan2((sin(theta_d_atr) * e_atr + sin(theta_d_a) * e_align + sin(theta_d_s) * e_spon),...
                (e_atr * cos(theta_d_atr) + cos(theta_d_s) * e_spon + cos(theta_d_a) * e_align)); % Desired direction due to all interactions
            
            if theta_d(i) < 0
                theta_d(i) = theta_d(i) + 2 * pi;
            end
            
        end
        
        % Turning rate depends on the difference between the current direction & desired direction
%         omega = abs(theta_d(i) - theta(i)); 
%         
%         % Selecting the smallest angle to turn towards
%         if omega > pi
%             omega = 2*pi - omega;
%         end
%         omega = omega/theta_tau;
        
        omega = pi/3;

        % Turn towards the side with shortest angular distance. 
        if theta(i) < theta_d(i)
            if theta_d(i) - theta(i) >= pi
                theta(i) = theta(i) - omega * dt;
                if theta(i) < 0
                    theta(i) = theta(i) + 2 * pi;
                    if theta(i) < theta_d(i)
                        theta(i) = theta_d(i);
                    end
                end
            else
                theta(i) = theta(i) + omega * dt;
                if theta(i) > theta_d(i)
                    theta(i) = theta_d(i);
                end
            end
            
        elseif theta(i) > theta_d(i)
            if theta(i) - theta_d(i) >= pi
                theta(i) = theta(i) + omega * dt;
                if theta(i) >= 2 * pi
                    theta(i) = theta(i) - 2 * pi;
                    if theta(i) > theta_d(i)
                        theta(i) = theta_d(i);
                    end
                end
            else
                theta(i) = theta(i) - omega * dt;
                if theta(i) < theta_d(i)
                    theta(i) = theta_d(i);
                end
            end
        end
        
        d_t(i,:) = v0 * [cos(theta(i)), sin(theta(i))];
        
        vel(i,:) = d_t(i,:); % Update velocity
        pos(i,:) = pos(i,:) + vel(i,:)*dt; % Update position
        
    end

    theta_t_1 = theta;
    pos_t_1 = pos;
    
    % Store data after sk_t time steps
    if rem(t,sk_t) == 0
        pos_t(:,:,t/sk_t) = pos;
        theta_t(:,t/sk_t) = theta;
    end
    
end

% Mean size of largest network cluster for given conn_time
conncomp_size_t = mean(conncomp_size_1(floor(st_t/(conn_time))+1:end));
% Mean no.of unique near neighbours for given conn_time
avg_uni_neigh_t = mean(avg_uni_neigh_1(floor(st_t/(conn_time))+1:end));