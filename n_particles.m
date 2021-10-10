% Date: 22-7-21
% Corrected Code

function [t_t, theta_t, pos_t, v, vel_t, v_mag, avg_conn_t, conncomp_size_t] = n_particles(n, r_spon, ...
    r_align, r_atr, dt, n_iter, zor, rad_rep, tau, theta_tau, K_alg, k_alg, K_atr, k_atr, k_r,...
    beta, sight, gamma, attr_c, latr, S0, Smax, conn_time)
 
% INITIAL CONDITIONS - Putting agents in a lattice with orientation around
% a random angle.

m=ceil(sqrt(n));
pos(:,1)=(rem((1:n)-1,m))*(zor*0.6);
pos(:,2)= floor(((1:n)-1)/m)*(zor*0.6);
theta = 2*pi*rand(1) + (sight/2)*(pi/180)*randn(n,1);

for j=1:n
    if theta(j)>2*pi
        theta(j)=theta(j)-2*pi;
    elseif theta(j)<0
        theta(j)=2*pi+theta(j);
    end
end

v0 = S0;
vmax=Smax;

% Truncated Distribution

nor_dis_spd = makedist('Normal', 'sigma', v0); % v0 is the mean velocity of the fish also used as stdiv
trunc_dis_spd = truncate(nor_dis_spd, -v0, vmax-v0);

nor_dis_ang = makedist('Normal', 'sigma', (sight/6)*(pi/180));
trunc_dis_ang = truncate(nor_dis_ang, (-sight/2)*(pi/180), (sight/2)*(pi/180)); 

v = v0 +  0.1*rand(n,1);
vel = [v.*cos(theta), v.*sin(theta)];
d_t = vel;
s_d = zeros(n,1);
theta_d = theta;

sk_t = 0.1/dt;

pos_t = zeros(n, 2, ceil(n_iter/sk_t));
theta_t = zeros(n, ceil(n_iter/sk_t));
v_mag = zeros(n, ceil(n_iter/sk_t));
vel_t = zeros(n, 2, ceil(n_iter/sk_t));

t_t = 0.1:0.1:n_iter*dt;

% Reaction rates
te_spon = (1/r_spon) * log(1./rand(n,1));
te_align = (1/r_align) * log(1./rand(n,1));
te_atr = (1/r_atr) * log(1./rand(n,1));

T = n_iter*dt;
avg_conn_t_1 = zeros(floor(T/conn_time(1)),1);
avg_conn_t_2 = zeros(floor(T/conn_time(2)),1);
avg_conn_t_3 = zeros(floor(T/conn_time(3)),1);
avg_conn_t_4 = zeros(floor(T/conn_time(4)),1);
avg_conn_t_5 = zeros(floor(T/conn_time(5)),1);
avg_conn_t_6 = zeros(floor(T/conn_time(6)),1);

conncomp_size_1 = zeros(floor(T/conn_time(1)),1);
conncomp_size_2 = zeros(floor(T/conn_time(2)),1);
conncomp_size_3 = zeros(floor(T/conn_time(3)),1);
conncomp_size_4 = zeros(floor(T/conn_time(4)),1);
conncomp_size_5 = zeros(floor(T/conn_time(5)),1);
conncomp_size_6 = zeros(floor(T/conn_time(6)),1);
 
connect_agents_1 = [];
agents_1 = [];

connect_agents_2 = [];
agents_2 = [];

connect_agents_3 = [];
agents_3 = [];

connect_agents_4 = [];
agents_4 = [];

connect_agents_5 = [];
agents_5 = [];

connect_agents_6 = [];
agents_6 = [];

tic;

for t = 2:n_iter
    
    if rem(t,1e4)==0
        disp(t*dt)
    end
    
    % Ignore this part for now
    if rem(t,round(conn_time(1)*(1/dt))) == 0
        
        [avg_conn_temp, con_graph] = avg_connections(agents_1, connect_agents_1, n);
        avg_conn_t_1((t/(conn_time(1)*(1/dt))),1) = avg_conn_temp;
        
        connect_agents_1 = [];
        agents_1 = [];
        [~, clussize] = conncomp(con_graph);
        if isempty(clussize) == 0
            conncomp_size_1((t/(conn_time(1)*(1/dt))),1) = max(clussize);
        end
        
    end
    
    if rem(t,round(conn_time(2)*(1/dt))) == 0
        
        [avg_conn_temp, con_graph] = avg_connections(agents_2, connect_agents_2, n);
        avg_conn_t_2((t/(conn_time(2)*(1/dt))),1) = avg_conn_temp;
        
        connect_agents_2 = [];
        agents_2 = [];
        [~, clussize] = conncomp(con_graph);
        if isempty(clussize) == 0
            conncomp_size_2((t/(conn_time(2)*(1/dt))),1) = max(clussize);
        end
        
    end
    
    if rem(t,round(conn_time(3)*(1/dt))) == 0
        
        [avg_conn_temp, con_graph] = avg_connections(agents_3, connect_agents_3, n);
        avg_conn_t_3((t/(conn_time(3)*(1/dt))),1) = avg_conn_temp;
        
        connect_agents_3 = [];
        agents_3 = [];
        [~, clussize] = conncomp(con_graph);
        if isempty(clussize) == 0
            conncomp_size_3((t/(conn_time(3)*(1/dt))),1) = max(clussize);
        end
        
    end
    
    if rem(t,round(conn_time(4)*(1/dt))) == 0
        
        [avg_conn_temp, con_graph] = avg_connections(agents_4, connect_agents_4, n);
        avg_conn_t_4((t/(conn_time(4)*(1/dt))),1) = avg_conn_temp;
        
        connect_agents_4 = [];
        agents_4 = [];
        [~, clussize] = conncomp(con_graph);
        if isempty(clussize) == 0
            conncomp_size_4((t/(conn_time(4)*(1/dt))),1) = max(clussize);
        end
        
    end
    
    if rem(t,round(conn_time(5)*(1/dt))) == 0
        
        [avg_conn_temp, con_graph] = avg_connections(agents_5, connect_agents_5, n);
        avg_conn_t_5((t/(conn_time(5)*(1/dt))),1) = avg_conn_temp;
        
        connect_agents_5 = [];
        agents_5 = [];
        [~, clussize] = conncomp(con_graph);
        if isempty(clussize) == 0
            conncomp_size_5((t/(conn_time(5)*(1/dt))),1) = max(clussize);
        end
        
    end
    
    if rem(t,round(conn_time(6)*(1/dt))) == 0
        
        [avg_conn_temp, con_graph] = avg_connections(agents_6, connect_agents_6, n);
        avg_conn_t_6((t/(conn_time(6)*(1/dt))),1) = avg_conn_temp;
%         ids = dbscan(pos, 2*zor, 1);
        
        connect_agents_6 = [];
        agents_6 = [];
        [~, clussize] = conncomp(con_graph);
        if isempty(clussize) == 0
            conncomp_size_6((t/(conn_time(6)*(1/dt))),1) = max(clussize);
        end
        
%         figure(1)
%         plot(con_graph, 'XData', pos(:,1), 'YData', pos(:,2))
%         hold all
%         gscatter(pos(:,1),  pos(:,2), ids)
%         hold off
        
    end
    
    for i = 1:n
        
        e_spon = te_spon(i) <= dt;
        e_align = te_align(i) <= dt;
        e_atr = te_atr(i) <= dt;
        te_spon(i) = te_spon(i) - dt;
        te_align(i) = te_align(i) - dt;
        te_atr(i) = te_atr(i) - dt;
        
        if te_spon(i) < 0
            te_spon(i) = (1/r_spon) * log(1/rand());
        end
        
        if te_align(i) < 0
            te_align(i) = (1/r_align) * log(1/rand());
        end
        
        if te_atr(i) < 0
            te_atr(i) = (1/r_atr) * log(1/rand());
        end
        
        dis_vec = pos - repmat(pos(i,:), size(pos,1), 1); % rij vector
        mag_vec = sqrt(dis_vec(:,1).^2 + dis_vec(:,2).^2); % magnitude of rij
        mag_vec_r = sqrt(dis_vec(:,1).^2 + dis_vec(:,2).^2) - 2 * rad_rep; %magnitude of rij - space oppucied by agents
        
        if (e_spon+e_align+e_atr) == 0
            
            theta_d(i) = theta_d(i);
            s_d(i) = s_d(i);
            
        else
            
            % SPONTANEOUS REACTION
            delta = random(trunc_dis_ang);
            
            theta_d_s = theta(i) + delta; % Desired direction due to spontaneous rotation
            if theta_d_s >= 2*pi
                theta_d_s = theta_d_s - 2 * pi;
            elseif theta_d_s < 0
                theta_d_s = theta_d_s + 2 * pi;
            end
            
            s_d_s = random(trunc_dis_spd); % Change in speed due to spontaneous rotation 
            
            % ALIGNMENT 
            
            % Check if the neighbors are in the visual zone
            dth_mag = (vel(i,1) * dis_vec(:,1) + vel(i,2) * dis_vec(:,2))./(mag_vec*v(i)+eps); 
                      
            th_jk = acosd(dth_mag); th_jk(i)=1e3;
            n_array = find(th_jk<sight/2);
            
            % sorting those in visual zone in the ascending order of their
            % distance from focal individual
            [~, s_ind_st] = sort(mag_vec(n_array), 'ascend');
            K_alg_st = min(K_alg, numel(s_ind_st));
            neighbours_alg = n_array(s_ind_st(1:K_alg_st));
            
            if isempty(n_array)==0
                neighbours_kalg=neighbours_alg(randperm(K_alg_st, min(k_alg, K_alg_st))); %Picking a random neighbor to align
                s_d_a = mean(v(neighbours_kalg,1),1) - v0; % Desired change in speed due to alignment interaction
                d_align_t = mean([cos(theta(neighbours_kalg,1)) sin(theta(neighbours_kalg,1))],1); % Desired direction
                theta_d_a = atan2(d_align_t(1,2),d_align_t(1,1)); % Desired direction
            else
                theta_d_a=theta(i);
                s_d_a=v(i)-v0;
            end
            
            % ATTRACTION
            
            % Check if the neighbors are in the visual zone
            th_at = abs(acosd(dth_mag)); th_at(i)=1e3;
            n_array_atr = find(th_at<sight/2);
            
            % sorting those in visual zone in the ascending order of their
            % distance from focal individual
            [~, s_ind_atr] = sort(mag_vec(n_array_atr), 'ascend');
            K_atr_st = min(K_atr, numel(s_ind_atr));
            neighbours_atr = n_array_atr(s_ind_atr(1:K_atr_st));
            
            if isempty(n_array_atr)==0
                neighbours_katr = neighbours_atr(randperm(K_atr_st, min(k_atr, K_atr_st)));
                s_d_attr = min(attr_c * (mean(((mag_vec(neighbours_katr,1) - 2 *rad_rep)/latr).^gamma,1)), vmax-v0);
                d_atr_t = mean([attr_c *(((mag_vec(neighbours_katr,1) - 2 *rad_rep)/latr).^gamma)...
                    .*dis_vec(neighbours_katr,:)./(norm(dis_vec(neighbours_katr,:))+eps) ; [cos(theta(i)) sin(theta(i))]],1);
                theta_d_atr = atan2(d_atr_t(1,2),d_atr_t(1,1));
                if e_atr == 1
                    nc = ones(1,length(neighbours_katr));
                    agents_1 = [agents_1, i*nc];
                    connect_agents_1 = [connect_agents_1, neighbours_katr.'];
                    agents_2 = [agents_2, i*nc];
                    connect_agents_2 = [connect_agents_2, neighbours_katr.'];
                    agents_3 = [agents_3, i*nc];
                    connect_agents_3 = [connect_agents_3, neighbours_katr.'];
                    agents_4 = [agents_4, i*nc];
                    connect_agents_4 = [connect_agents_4, neighbours_katr.'];
                    agents_5 = [agents_5, i*nc];
                    connect_agents_5 = [connect_agents_5, neighbours_katr.'];
                    agents_6 = [agents_6, i*nc];
                    connect_agents_6 = [connect_agents_6, neighbours_katr.'];
                end
            else
                theta_d_atr=theta(i);
                s_d_attr=0;
            end
            
            theta_d(i) = atan2((sin(theta_d_atr) * e_atr + sin(theta_d_a) * e_align + sin(theta_d_s) * e_spon),...
                (e_atr * cos(theta_d_atr) + cos(theta_d_s) * e_spon + cos(theta_d_a) * e_align)); % Desired direction due to all interactions
            s_d(i) = ((e_spon * s_d_s) + (e_align * (s_d_a)) + (e_atr * s_d_attr))/(e_spon + e_align + e_atr); % Desired direction due to all interactions
            
            if theta_d(i) < 0
                theta_d(i) = theta_d(i) + 2 * pi;
            end
            
        end
        
        % REPULSION
        
        %check if the distance between agents is decreasing
        fapp = max((dot(repmat(vel(i,:), size(vel,1), 1) ,dis_vec, 2) ...
            + dot(vel, -dis_vec, 2)), zeros(size(pos,1),1)); 
        
        % Check if the focal individual is moving towards the neighbor
        dth_mag = (vel(i,1) * dis_vec(:,1) + vel(i,2) * dis_vec(:,2))./(mag_vec*v(i)+eps);
        dth_mag(i) = -1; % This check if two are approaching each other
        
        % Select only those agents which are within colision zone and the
        % distances between them is decreasing
        z = (mag_vec_r < zor) & (fapp > 0);
        z(i) = 0;
        z1 = sum(z);
        
        if z1 ~= 0 && k_r~=0
            
            % variation of speed because of repulsion
            fprox = k_r./(mag_vec_r(z,1).^beta + eps); 
%             c1=10;
%             c2=10;
%             fc = 1 - 1./(exp(-c1.*mag_vec_r(z,1) + c2) + 1);
            
            s_r = mean(fapp(z) .* fprox .*(dth_mag(z)>0)); % ASYMMETRIC VERSION: using the direction of movement of agent
            
            vec_r_1 = [-(dis_vec(z,2)), dis_vec(z,1)]; % Vector-1 perpendiculat to rij
            vec_r_1 = vec_r_1 ./ (sqrt(vec_r_1(:,1).^2 + vec_r_1(:,2).^2)+eps);
            vec_r_2 = [(dis_vec(z,2)), -dis_vec(z,1)]; % Vector-2 perpendiculat to rij
            vec_r_2 = vec_r_2 ./ (sqrt(vec_r_2(:,1).^2 + vec_r_2(:,2).^2)+eps);
            
            v_norm = vel(i,:)/norm(vel(i,:));
            
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
            
            s_d(i) = s_d(i) + s_r; %desired direction after including repulsion, s_r depends on the distance between the particles.
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
        
        v(i) = v(i) * (1 - (dt/tau)) + s_d(i) * dt/tau + v0 * dt/tau; % Calculated speed due to all interactions
        
        if v(i) < 0
            v(i) = 0;
        end
        if v(i) > Smax
            v(i) = Smax;
        end
        
        omega = abs(theta_d(i) - theta(i)); % Turning rate depends on the difference between the current direction & desired direction
        
        % Selecting the smallest angle 
        if omega > pi
            omega = 2*pi - omega;
        end
        omega = omega/theta_tau;

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
        
        vel(i,:) = d_t(i,:);
        pos(i,:) = pos(i,:) + vel(i,:)*dt;
        
    end
    
    if rem(t,sk_t) == 0
        vel_t(:,:,t/sk_t) = vel;
        pos_t(:,:,t/sk_t) = pos;
        theta_t(:,t/sk_t) = theta;
        v_mag(:,t/sk_t) = v;
    end
    
end

avg_conn_t = [mean(avg_conn_t_1); mean(avg_conn_t_2); mean(avg_conn_t_3); mean(avg_conn_t_4); mean(avg_conn_t_5); mean(avg_conn_t_6)];
conncomp_size_t = [mean(conncomp_size_1); mean(conncomp_size_2); mean(conncomp_size_3); mean(conncomp_size_4);...
    mean(conncomp_size_5); mean(conncomp_size_6)];