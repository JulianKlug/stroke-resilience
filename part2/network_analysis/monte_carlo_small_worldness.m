function [S_MC]= monte_carlo_small_worldness(A, Num_ER_repeats, FLAG)
%monte_carlo_small_worldness: Output Small-Worldness stat coefficient according to
%Monte Carlo estimates of path-length and clustering of an Erdos-Renyi random graph
% 
% ARGUMENTS:
% - A: adjacency matrix of graph to analyse
% - Num_ER_repeats: number of random graphs generated (to estimate C and L
% numerically for E-R random graph)
% - FLAG is a number indicating which clustering coefficient to compute:
%   1 - Cws 
%   2 - transitivity C (no. of triangles)


% get its basic properties
n = size(A,1);  % number of nodes
k = sum(A);  % degree distribution of undirected network
m = sum(k)/2;

%% computing small-world-ness by estimating L_rand and C_rand from an ensemble of random graphs
[Lrand,Crand] = NullModel_L_C(n,m,Num_ER_repeats,FLAG);
Lrand_mean = mean(Lrand);

% compute small-world-ness using mean value over Monte-Carlo realisations
[S_MC,~,~] = small_world_ness(A,Lrand_mean,mean(Crand),FLAG);  
