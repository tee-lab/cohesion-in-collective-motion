function [avg_conn_temp, con_graph] = avg_connections(agents, connect_agents, n)


if isempty(agents) == 0 && isempty(connect_agents) == 0
    
    graph_cons = [agents.', connect_agents.'];
    graph_cons = unique(graph_cons, 'rows');
    
    con_graph = digraph(graph_cons(:,1), graph_cons(:,2)); % Creates directional graph
    p_d = distances(con_graph); % Distance matrix
    
    for p1 = 1:size(p_d,1)
        
        for p2 = 1:size(p_d,1)
            
            if ~isinf(p_d(p1,p2)) && p_d(p1,p2) > 0
                p_d(p1,p2) = 1;
            else
                p_d(p1,p2) = 0;
            end
            
        end
        
    end
    
    avg_conn_temp = sum(sum(p_d,2))/n;
    
else
    
    avg_conn_temp = 0;
    con_graph = digraph(agents, connect_agents);
    
end