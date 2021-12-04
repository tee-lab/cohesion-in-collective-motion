%  Takes in agents and the agents they form edges with and gives the
%  connected graph. Also, refer main text - network analysis

% Outputs directed graph
function con_graph = avg_connections(agents, connect_agents, n)

% Check if the input is a non-empty set
if isempty(agents) == 0 && isempty(connect_agents) == 0
    
    graph_cons = [agents.', connect_agents.'];
    graph_cons = unique(graph_cons, 'rows');% Remove the same repeated interaction 
    % within the time interval of conn_time
    
    weights = ones(1, size(graph_cons,1));
    con_graph = digraph(graph_cons(:,1), graph_cons(:,2), weights, n); % Creates directional graph
    
else % Else null graph 

    weights = ones(1, size(agents,1));
    con_graph = digraph(agents, connect_agents, weights, n);
    
end