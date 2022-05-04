% Date: 14-1-2022
% Edited code with comments
% Written by Vivek J and Danny Raj M

function [t_t, theta_t, pos_t, conncomp_size_t, avg_uni_neigh_t] = n_particles(n, r_int, ...
    dt, n_iter, sigma_t, omega, sight, K_alg, k_alg, K_atr, k_atr, ...
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
    if theta(j)>=2*pi
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

sk_t = ceil(0.1/dt); % Store data after sk_t time steps with apprx every 0.1s

% position, velocity, speed and direction at time t-1

pos_t_1 = pos;
theta_t_1 = theta;
vel_t_1 = vel;

pos_t = zeros(n, 2, floor(n_iter/sk_t)); % Stores positions of agents over the simulation
theta_t = zeros(n, floor(n_iter/sk_t)); % Stores heading angles of agents over the simulation

t_t = (sk_t*dt):(sk_t*dt):n_iter*dt; % Time points at which data is stored  

% Reaction rates
te_int = (1/r_int) * log(1/rand()); % time at which interaction happens

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
    
    % Constructs networks over conn_time time interval only if time elapsed is
    % greater than st_t (to account for initial conditions)
    if rem(t,floor(conn_time*(1/dt))) == 0 
        
        con_graph = avg_connections(agents_1, connect_agents_1, n);
%         ids = dbscan(pos, 2*zor, 1);
        
        connect_agents_1 = [];
        agents_1 = [];
        [~, clussize] = conncomp(con_graph);
        if isempty(clussize) == 0
            conncomp_size_1((t/floor(conn_time*(1/dt))),1) = max(clussize);
            avg_uni_neigh_1((t/floor(conn_time*(1/dt))),1) = mean(outdegree(con_graph));
        end

%         figure(1)
        
%         plot(con_graph, 'XData', pos(:,1), 'YData', pos(:,2))
%         hold all
%         sc = gscatter(pos(:,1),  pos(:,2), net_bin);
%         hold all
%         xlabel('X', 'FontName', 'Serif', 'FontSize', 10, 'FontWeight', 'bold')
%         ylabel('Y', 'FontName', 'Serif', 'FontSize', 10, 'FontWeight', 'bold')
%         legend(sc)

    end

    e_int = te_int <= dt; % Check if interaction happens at this time step
    te_int = te_int - dt;

    % If the interaction happened at this time step, chose a next time
    % point where the interaction will happen.
    if te_int < 0
        te_int = (1/r_int) * log(1/rand());
    end

    if e_int == 0

        for i = 1:n

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

    else

        for i = 1:n

            dis_vec = pos_t_1 - repmat(pos_t_1(i,:), size(pos_t_1,1), 1); % rij vector (as defined in the main text)
            mag_vec = sqrt(dis_vec(:,1).^2 + dis_vec(:,2).^2); % magnitude of rij
            %         mag_vec(i) = Inf;

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

            % Check if the neighbors are in the visual zone
            dth_mag = (vel_t_1(i,1) * dis_vec(:,1) + vel_t_1(i,2) * dis_vec(:,2))./(mag_vec*v0+eps);

            th_jk = acosd(dth_mag);
            th_jk(i)=1e3; % Not considering the focal individual for the interaction
            n_array = find(th_jk<sight/2);

            % sorting those in visual zone in the ascending order of their
            % distance from focal individual
            [~, s_ind_st] = sort(mag_vec(n_array), 'ascend');
            % Take the minimum of K (as defined in the main text) and agents in the visual field
            K_alg_st = min(K_alg, numel(s_ind_st));
            neighbours_alg = n_array(s_ind_st(1:K_alg_st));

            if isempty(n_array)==0 % If there are neighbors in visual field to interact with
                neighbours_kalg=neighbours_alg(randperm(K_alg_st, min(k_alg, K_alg_st))); %Picking a random neighbor to align
                d_align_t = mean([cos(theta_t_1(neighbours_kalg,1)) sin(theta_t_1(neighbours_kalg,1))],1); % Desired direction
                theta_d_a = atan2(d_align_t(1,2),d_align_t(1,1)); % Desired direction
            else % Else continue to move towards the desired direction and speed from previous time
                theta_d_a=theta_d(i);
            end

            % ATTRACTION

            % Check if the neighbors are in the visual zone
            th_at = abs(acosd(dth_mag)); th_at(i)=1e3;
            n_array_atr = find(th_at<sight/2);

            % sorting those in visual zone in the ascending order of their
            % distance from focal individual
            [~, s_ind_atr] = sort(mag_vec(n_array_atr), 'ascend');
            % Take the minimum of K (as defined in the main text) and agents in the visual field
            K_atr_st = min(K_atr, numel(s_ind_atr));
            neighbours_atr = n_array_atr(s_ind_atr(1:K_atr_st));

            if isempty(n_array_atr)==0
                %Picking a random neighbor to align
                neighbours_katr = neighbours_atr(randperm(K_atr_st, min(k_atr, K_atr_st)));
                % Desired direction
                %             d_atr_t = mean([(((mag_vec(neighbours_katr,1))/latr).^gamma)...
                %                     .*dis_vec(neighbours_katr,:)./(norm(dis_vec(neighbours_katr,:))+eps) ; [cos(theta_t_1(i)) sin(theta_t_1(i))]],1);
                d_atr_t = mean(dis_vec(neighbours_katr,:)./(mag_vec(neighbours_katr,1)+eps),1);
                theta_d_atr = atan2(d_atr_t(1,2),d_atr_t(1,1));

                %  If attraction interaction has happened and time is
                %  greater than st_t then record which agent interacts with
                %  which agent
                if e_int == 1 && floor(dt*t) > st_t
                    nc = ones(1,length(neighbours_katr));
                    agents_1 = [agents_1, i*nc];
                    connect_agents_1 = [connect_agents_1, neighbours_katr.'];
                end
            else
                theta_d_atr=theta_d(i);
            end

            % Take the average of all the interactions that happended in
            % this time interval dt
            theta_d(i) = atan2((sin(theta_d_atr) + sin(theta_d_a) + sin(theta_d_s)),...
                (cos(theta_d_atr) + cos(theta_d_s) + cos(theta_d_a))); % Desired direction due to all interactions

            if theta_d(i) < 0
                theta_d(i) = theta_d(i) + 2 * pi;
            end

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

    end

    theta_t_1 = theta;
    pos_t_1 = pos;
    vel_t_1 = vel;

    % Store data after sk_t time steps
    if rem(t,sk_t) == 0
        pos_t(:,:,t/sk_t) = pos;
        theta_t(:,t/sk_t) = theta;
    end

end

% Mean size of largest network cluster for given conn_time by considering
% datas only after st_t
conncomp_size_t = mean(conncomp_size_1(ceil(st_t/(conn_time))+1:end));
% Mean no.of unique near neighbours for given conn_time by considering
% datas only after st_t
avg_uni_neigh_t = mean(avg_uni_neigh_1(ceil(st_t/(conn_time))+1:end));