function [n_nodes, n_edges, mean_degree, degree_distribution] = basic_graph_properties(A)
%BASIC_GRAPH_PROPERTIES Compute basic graph properties
% ARGUMENTS:
% - A: adjacency matrix of graph to analyse

n_nodes = size(A,1);  % number of nodes
degree_distribution = sum(A);  % degree distribution of undirected network
n_edges = sum(degree_distribution)/2;
mean_degree = mean(degree_distribution); % mean degree of network
end

