clc; clear all; close all;
%% Example 1 Evaluation
% nodes are ordered begining from source in anticlock wise direction in
% Fig. 2

source=7; target=10; 
directed=true; 
nodes_neglected=true;
terminals_excluded=true;
 
source_nodes=[7 8 8 9 7 9 ];
target_nodes=[8 10 9 8 9 10];
net=reliability_net(source, target, directed, nodes_neglected, ...
        terminals_excluded, source_nodes, target_nodes);
 
%_______________ Entering Maximal Capacity for each element _______________
CP=cell(1,6);           
CP{1}=[0.05 0.1 0.25 0.6]; 
CP{2}=[0.1 0.3 0.6];       
CP{3}=[0.1 0.9];           
CP{4}=CP{3};
CP{5}=CP{3};
CP{6}=[0.05 0.25 0.7];
demand=randi([1 10], 1, max(cellfun(@(h) length(h), CP)));
demand=demand/sum(demand);
take_capacity (net, CP, demand);
shortest_paths (net, true);
net.flow_constraints={'flow=demand', 'less than Lj', 'maximal capacity constraint'};
evaluate_reliability(net, true, [], [], false);
%% Displaying network data
disp(net);
%% Plotting the Reliability curve
x = net.reliability.demand;
y = net.reliability.rel;
plot(x,y, "LineWidth",3, "Marker","x","MarkerSize",10)
grid on;
xlabel("Demand");
ylabel("Reliability");
