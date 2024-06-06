%%INITIALIZE ^~

addpath('SupportingFiles/')

%Read data
Nodes = load('SupportingFiles/Nodes.txt');
Links = load('SupportingFiles/Links.txt');
L = load('SupportingFiles/L.txt');

%Get sizes of input data
nNodes = size(Nodes,1);
nLinks = size(Links, 1);

%Create graph
G = graph(L);

%Show graph
figure(1)
selected = [];
plotTopology(Nodes, Links,selected)

n = 8;
t = 30;

[nodes,result] = Grasp(G,n,t);
plotTopology(Nodes, Links, nodes);
disp("GRASP BEST RESULT: " + result)

function [nodes_best, value_best] = Grasp(G,n,time)
t = tic;
[nodes_best] = GreedyRandomized(G,n);
[nodes_best, value_best] = Adaptative_Search(G, nodes_best,toc(t));


while toc(t) < time 
    [nodes] = GreedyRandomized(G,n);
    [nodes, value] = Adaptative_Search(G,nodes,toc(t));

    if value < value_best
        nodes_best = nodes;
        value_best = value;
    end

end

end

function s = GreedyRandomized(G,n)

r = 3;
E = 1:numnodes(G);
s = [];

for i=1:n
    R = [];
    for j=E
        R = [R; j ConnectedNP(G,[s j])];
    end
    R = sortrows(R,2);
    e = R(randi(r),1);
    s = [s e];
    E = setdiff(E,e);
end

end

function [result_nodes, result] = Adaptative_Search(G, initial,time)
    all_nodes = size(G.Nodes,1);
    best_result = ConnectedNP(G,initial);
    best_nodes = initial;
    improved = 1;
    t = tic;
    while improved && toc(t) < time
        other_nodes = setdiff(1:all_nodes,best_nodes);
        best_neigh_nodes = best_nodes;
        best_neigh_result = best_result;
        for i1 = best_nodes
            % other_nodes = setdiff(neighbors(G,i1)',best_nodes);
            for i2 = other_nodes
                neigh_nodes = [setdiff(best_nodes,i1) i2];
                neigh_result = ConnectedNP(G,neigh_nodes);
                if neigh_result < best_neigh_result
                    best_neigh_nodes = neigh_nodes;
                    best_neigh_result = neigh_result;
                end
            end
        end
        if best_neigh_result <= best_result && ~isequal(best_nodes,best_neigh_nodes)
            best_nodes = best_neigh_nodes;
            best_result = best_neigh_result;
        else
            improved = 0;
        end
    end
    result = best_result;
    result_nodes = best_nodes;
end



