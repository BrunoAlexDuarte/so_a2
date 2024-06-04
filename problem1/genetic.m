
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
P = 30;
q = 0.5;
m = 3;
t = 30;

[nodes,result] = Genetic(G,n,P,q,m,t);
plotTopology(Nodes, Links, nodes);
disp("BEST RESULT: " + result)

%%%

function [s_best,result_best] = Genetic(G,n,P,q,m,time)

%Create variables to use
Pop = zeros(n,P); %Population
all_nodes = size(G.Nodes,1); %Size of all nodes

%Initialize initial population with random nodes
for i=1:P
    Pop(:,i) = randperm(all_nodes,n);
end

t = tic;
count = 0;
while toc(t) < time % count < 100%

    Pop2 = zeros(n,P);
    for i=1:P
        s = Crossover(Pop,G);
        if rand() < q
            s = Mutation(s,all_nodes);
        end
        Pop2(:,i) = s;
    end
    Pop = Selection(Pop,Pop2,m, G);
    count = count + 1;
end

[index_best,result_best] = Best(Pop, G);
s_best = Pop(:, index_best);

end

%%%

function s = Crossover(Pop,G)

n = size(Pop,1);
P = size(Pop,2);
Nparents = 2; %Number of parents in the tourney
Parents = zeros(n,Nparents);

%Torney for the 2 parents
for i=1:2
    %Create candidate parents
    Parents_aux = randperm(P, Nparents);
    best_parents = 1;
    best_value = ConnectedNP(G, Pop(:, Parents_aux(best_parents)));
    %Get the best parent of the candidates
    for j=2:Nparents
        value = ConnectedNP(G, Pop(:, Parents_aux(j)));
        if value <= best_value
            best_value = value;
            best_parents = j;
        end
    end
    Parents(:, i) = Pop(:, Parents_aux(best_parents));
end

all_parents = unique([Parents(:, 1); Parents(:, 2)]);
tam = size(all_parents,1);
s = all_parents(randperm(tam,n));

end

%%%

function s = Mutation(S, all_nodes)

n = size(S, 1);
aux= setdiff(1:all_nodes,S);
s= [S(randperm(n,n-1));aux(randperm(all_nodes-n,1))];

end

%%%

function Pop = Selection(Pop1, Pop2, m, G)

n = size(Pop1,1);
P = size(Pop1,2);
Pop = zeros(n, P);
% disp("SELECTION WILL START")
for i=1:P
    % disp("POP 2")
    [index_2, s2] = Best(Pop2,G);
    if m > 0 %In case i can still be elitist
        % disp("POP 1")
        [index_1, s1] = Best(Pop1,G);
        % Check if there is a solution in the older population better 
        % than in the new population 
        if s1 <= s2 
            Pop(:, i) = Pop1(:, index_1); 
            Pop1(:, index_1) = [];
            m = m - 1;
        else
            Pop(:, i) = Pop2(:, index_2);
            Pop2(:, index_2) = [];
        end
    else %If we can't be elitist, we just take the better solution in the newer population
        Pop(:, i) = Pop2(:, index_2);
        Pop2(:, index_2) = [];    
    end
end

end


function [index, value] = Best(Pop,G)
%Gets the best solution for a problem with graph G in population Pop

P = size(Pop,2);
index = 1;
value = ConnectedNP(G, Pop(:, index)); %In this problem we only use this
%Otherwise we would need a parameter to pass a function to evaluate
for i = 2:P
    new_value = ConnectedNP(G, Pop(:, i));
    if new_value <= value
        value = new_value;
        index = i;
    end
end
% disp("BEST VALUE: " + value + " AND INDEX:" + index)

end
