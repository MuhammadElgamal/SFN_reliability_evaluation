clc; clear all; close all;
%% Example 1 Evaluation
% This example runs on the paper of 
% Lin, Yi-Kuei, Cheng-Ta Yeh, and Cheng-Fu Huang. "Reliability evaluation of a 
% stochastic-flow distribution network  with delivery spoilage." 
% Computers & Industrial Engineering 66.2 (2013): 352-359.
% We are concerned about the flow from Kaosiung to Hangzhou city 
% through arcs of 1 , 2, 9, 10 which we reorder as 1, 2, 3, 4 respectively and 
% rename the source as 5 and Shanghai as 6, Wenzhou as 7 and Hangzhou as 8
% hi babe


source=5; target=8; 
directed=true; 
nodes_neglected=true;
terminals_excluded=true;
 
source_nodes=[5 5 6 7];
target_nodes=[6 7 8 8];
net=reliability_net(source, target, directed, nodes_neglected, ...
        terminals_excluded, source_nodes, target_nodes);
 
%_______________ Entering Maximal Capacity for each element _______________
CP=cell(1,4);           
CP{1}=[0.002 0.004 0.005 0.006 0.007 0.01 0.013 0.033 0.92]; 
CP{2}=[0.006 0.007 0.007 0.011 0.013 0.016 0.016 0.024 0.9];       
CP{3}=[0.01 0.01 0.05 0.013 0.013 0.018 0.886];           
CP{4}=[0.008 0.008 0.01 0.011 0.013 0.03 0.039 0.881];

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
