% Date: 14-1-2022
% Edited code with comments
% Written by Vivek J and Danny Raj M 

function [t_t, theta_t, pos_t, v, vel_t, v_mag, conncomp_size_t, avg_uni_neigh_t] = n_particles(n, r_spon, ...
    r_align, r_atr, dt, n_iter, zor, rad_rep, tau, theta_tau, K_alg, k_alg, K_atr, k_atr, k_r,...
    beta, sight, gamma, attr_c, latr, S0, Smax, conn_time, st_t)
 
% INITIAL CONDITIONS - Putting agents in a lattice with orientation around
% a random angle.

m=ceil(sqrt(n));
pos(:,1)=(rem((1:n)-1,m))*(zor*0.6); % x coordinates 
pos(:,2)= floor(((1:n)-1)/m)*(zor*0.6); % y coordinates
theta = 2*pi*rand(1) + (sight/2)*(pi/180)*randn(n,1); % Distribution of heading angles around a 
% randomly chosen angle

% Checking if all the heading angles are between [0-2pi] 
for j=1:n
    if theta(j)>=2*pi
        theta(j)=theta(j)-2*pi;
    elseif theta(j)<0
        theta(j)=theta(j)+2*pi;
    end
end

v0 = S0; % Cruise speed
vmax=Smax; % Maximum speed an agent can move with 

% Truncated Distributions of angles and speeds for spontaneous reaction

nor_dis_spd = makedist('Normal', 'sigma', v0); % v0 is the mean velocity of the fish also used as stdiv
trunc_dis_spd = truncate(nor_dis_spd, -v0, vmax-v0);

nor_dis_ang = makedist('Normal', 'sigma', (sight/6)*(pi/180));
trunc_dis_ang = truncate(nor_dis_ang, (-sight/2)*(pi/180), (sight/2)*(pi/180)); 

% Defining speed, velocity for n agents
v = v0 +  0.1*rand(n,1);
vel = [v.*cos(theta), v.*sin(theta)];

d_t = vel; % Desired heading direction 
s_d = zeros(n,1); % Desired chanege in speed
theta_d = theta; % Desired heading angle 

sk_t = ceil(0.1/dt); % Store data after sk_t time steps with apprx every 0.1s

% position, velocity, speed and direction at time t-1

pos_t_1 = pos;
vel_t_1 = vel;
v_t_1 = v;
theta_t_1 = theta;

pos_t = zeros(n, 2, floor(n_iter/sk_t)); % Stores positions of agents over the simulation
theta_t = zeros(n, floor(n_iter/sk_t)); % Stores heading angles of agents over the simulation
v_mag = zeros(n, floor(n_iter/sk_t)); % Stores speeds of agents over the simulation
vel_t = zeros(n, 2, floor(n_iter/sk_t)); % Stores velocities of agents over the simulation

t_t = (sk_t*dt):(sk_t*dt):n_iter*dt; % Time points at which data is stored 

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
        mag_vec_r = sqrt(dis_vec(:,1).^2 + dis_vec(:,2).^2) - 2 * rad_rep; %magnitude of rij minus space oppucied by agents
        
        % Check if any event is happening in this time interval

        %if no event, continue to move in the desired direction of previous step

        if (e_spon+e_align+e_atr) == 0
            
            theta_d(i) = theta_d(i);
            s_d(i) = s_d(i);
            
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
            
            s_d_s = random(trunc_dis_spd); % Change in speed due to spontaneous rotation 
            
            % ALIGNMENT 
            
            % Check if the neighbors are in the visual zone
            dth_mag = (vel_t_1(i,1) * dis_vec(:,1) + vel_t_1(i,2) * dis_vec(:,2))./(mag_vec*v_t_1(i)+eps); 
                      
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
                s_d_a = mean(v_t_1(neighbours_kalg,1),1) - v0; % Desired change in speed due to alignment interaction
                d_align_t = mean([cos(theta_t_1(neighbours_kalg,1)) sin(theta_t_1(neighbours_kalg,1))],1); % Desired direction
                theta_d_a = atan2(d_align_t(1,2),d_align_t(1,1)); % Desired direction
            else % Else continue to move towards the desired direction and speed from previous time
                theta_d_a=theta_d(i);
                s_d_a=s_d(i);
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
                %Picking a random neighbor to attract
                neighbours_katr = neighbours_atr(randperm(K_atr_st, min(k_atr, K_atr_st)));
                % Desired change in speed due to attraction interaction
                s_d_attr = min(attr_c * (mean(((mag_vec(neighbours_katr,1) - 2 *rad_rep)/latr).^gamma,1)), vmax-v0);
                % Desired direction
                d_atr_t = mean([attr_c * (((mag_vec(neighbours_katr,1) - 2 *rad_rep)/latr).^gamma)...
                    .*(dis_vec(neighbours_katr,:)./(mag_vec(neighbours_katr,1)+eps)) ; [cos(theta_t_1(i)) sin(theta_t_1(i))]],1);
                theta_d_atr = atan2(d_atr_t(1,2),d_atr_t(1,1));

                %  If attraction interaction has happened and time is
                %  greater than st_t then record which agent interacts with
                %  which agent
                if e_atr == 1 && floor(dt*t) > st_t
                    nc = ones(1,length(neighbours_katr));
                    agents_1 = [agents_1, i*nc];
                    connect_agents_1 = [connect_agents_1, neighbours_katr.'];
                end
            else
                theta_d_atr=theta_d(i);
                s_d_attr=s_d(i);
            end
            
            % Take the average of all the interactions that happended in
            % this time interval dt
            theta_d(i) = atan2((sin(theta_d_atr) * e_atr + sin(theta_d_a) * e_align + sin(theta_d_s) * e_spon),...
                (e_atr * cos(theta_d_atr) + cos(theta_d_s) * e_spon + cos(theta_d_a) * e_align)); % Desired direction due to all interactions
            s_d(i) = ((e_spon * s_d_s) + (e_align * (s_d_a)) + (e_atr * s_d_attr))/(e_spon + e_align + e_atr); % Desired direction due to all interactions
            
            if theta_d(i) < 0
                theta_d(i) = theta_d(i) + 2 * pi;
            end
            
        end
        
        % REPULSION
        
        % Check if the distance between agents is decreasing (Refer main
        % text for details).
        fapp = max((dot(repmat(vel_t_1(i,:), size(vel_t_1,1), 1) ,dis_vec, 2) ...
            + dot(vel_t_1, -dis_vec, 2)), zeros(size(pos_t_1,1),1)); 
        
        % Check if the focal individual is moving towards the neighbor
        dth_mag = (vel_t_1(i,1) * dis_vec(:,1) + vel_t_1(i,2) * dis_vec(:,2))./(mag_vec*v_t_1(i)+eps);
        dth_mag(i) = -1; % This check if two are approaching each other
        
        % Select only those agents which are within colision zone and the
        % distances between them is decreasing
        z = (mag_vec_r < zor) & (fapp > 0);
        z(i) = 0;
        z1 = sum(z);
        
        if z1 ~= 0 && k_r~=0
            
            % variation of speed because of repulsion
            fprox = k_r./(mag_vec_r(z,1).^beta + eps); 
            
            % Only if agent is moving towards the neighbor this term is
            % non-zero.
            s_r = mean(fapp(z) .* fprox .*(dth_mag(z)>0)); % ASYMMETRIC VERSION: using the direction of movement of agent
            
            vec_r_1 = [-(dis_vec(z,2)), dis_vec(z,1)]; % Vector-1 perpendiculat to rij
            vec_r_1 = vec_r_1 ./ (sqrt(vec_r_1(:,1).^2 + vec_r_1(:,2).^2)+eps);
            vec_r_2 = [(dis_vec(z,2)), -dis_vec(z,1)]; % Vector-2 perpendiculat to rij
            vec_r_2 = vec_r_2 ./ (sqrt(vec_r_2(:,1).^2 + vec_r_2(:,2).^2)+eps);
            
            v_norm = [cos(theta_t_1(i)) sin(theta_t_1(i))]; % Heading direction to find the nearest perpendicular
            
            % Find the perpendicular vector with lesser angular distance. 
            theta_r_1_c = acos(dot(repmat(v_norm, z1, 1), vec_r_1, 2));
            theta_r_2_c = acos(dot(repmat(v_norm, z1, 1), vec_r_2, 2));
            
            vec_r = zeros(z1,2); % Perpendicular vector it choses to move
            for jr = 1:z1
                [~,ind] = min([theta_r_1_c(jr), theta_r_2_c(jr)], [], 2);
                if ind == 1
                    vec_r(jr,:) = vec_r_1(jr,:);
                elseif ind == 2
                    vec_r(jr,:) = vec_r_2(jr,:);
                end
            end
            
            %desired direction after including repulsion, s_r depends on 
            % the distance between the particles.
            s_d(i) = s_d(i) + s_r; 
            wts=max([zeros(z1,1) fapp(z) .* fprox .*-sign(dth_mag(z))], [], 2);
            vec_r = mean([wts wts].*vec_r, 1);
            vec_r = vec_r/(vecnorm(vec_r)+eps);
            theta_r_n = [cos(theta_d(i))+vec_r(1) sin(theta_d(i))+vec_r(2)];
            theta_r_n = theta_r_n/(vecnorm(theta_r_n)+eps);
            theta_d(i) = atan2(theta_r_n(1,2), theta_r_n(1,1));
            if theta_d(i) < 0
                theta_d(i) = theta_d(i) + 2 * pi;
            end
        end
        
        % Calculated speed due to all interactions
        v(i) = v(i) * (1 - (dt/tau)) + s_d(i) * dt/tau + v0 * dt/tau; % Calculated speed due to all interactions

        % Ensuring that the speeds are within [0 Smax]
        if v(i) < 0
            v(i) = 0;
        end
        if v(i) > Smax
            v(i) = Smax;
        end
        
        % Turning rate depends on the difference between the current direction & desired direction
        omega = abs(theta_d(i) - theta(i)); 
        
        % Selecting the smallest angle to turn towards
        if omega > pi
            omega = 2*pi - omega;
        end
        omega = omega/theta_tau;

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
        
        d_t(i,:) = v(i) * [cos(theta(i)), sin(theta(i))];
        
        vel(i,:) = d_t(i,:); % Update velocity
        pos(i,:) = pos(i,:) + vel(i,:)*dt; % Update position
        
    end

    theta_t_1 = theta;
    vel_t_1 = vel;
    v_t_1 = v;
    pos_t_1 = pos;
    
    % Store data after sk_t time steps
    if rem(t,sk_t) == 0
        vel_t(:,:,t/sk_t) = vel;
        pos_t(:,:,t/sk_t) = pos;
        theta_t(:,t/sk_t) = theta;
        v_mag(:,t/sk_t) = v;
    end
    
end

% Mean size of largest network cluster for given conn_time by considering
% datas only after st_t
conncomp_size_t = mean(conncomp_size_1(ceil(st_t/(conn_time))+1:end));
% Mean no.of unique near neighbours for given conn_time by considering
% datas only after st_t
avg_uni_neigh_t = mean(avg_uni_neigh_1(ceil(st_t/(conn_time))+1:end));