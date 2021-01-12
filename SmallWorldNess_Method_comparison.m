function [small_worldness_results] = SmallWorldNess_Method_comparison(A, Num_ER_repeats, Num_S_repeats, I, compute_significance)
%SMALLWORLDNESS_METHOD_COMPARISON Output Small-Worldness stats according to
%multiple methods
%   Methods used: 
%       - analytical approximation if network was Erdos-Renyi random graph
%       - Monte Carlo estimates of path-length and clustering of an Erdos-Renyi random graph
% ARGUMENTS:
% - A: adjacency matrix of graph to analyse

% get its basic properties
n = size(A,1);  % number of nodes
k = sum(A);  % degree distribution of undirected network
m = sum(k)/2;
K = mean(k); % mean degree of network

FLAG_Cws = 1;
FLAG_Ctransitive = 2;

%% computing small-world-ness using the analytical approximations for the E-R graph

[expectedC,expectedL] = ER_Expected_L_C(K,n);  % L_rand and C_rand

[S_ws,C_ws,L_ws,Cs_ws, Ls_ws] = custom_small_world_ness(A,expectedL,expectedC,FLAG_Cws);  % Using WS clustering coefficient
[S_trans,C_trans,L_trans, Cs_trans, Ls_trans] = custom_small_world_ness(A,expectedL,expectedC,FLAG_Ctransitive);  %  Using transitive clustering coefficient

%% computing small-world-ness by estimating L_rand and C_rand from an ensemble of random graphs
% check when using small networks...

[Lrand,CrandWS] = NullModel_L_C(n,m,Num_ER_repeats,FLAG_Cws);
[Lrand,CrandTrans] = NullModel_L_C(n,m,Num_ER_repeats,FLAG_Ctransitive);

% Note: if using a different random graph null model, e.g. the
% configuration model, then use this form

% compute small-world-ness using mean value over Monte-Carlo realisations

% NB: some path lengths in L will be INF if the ER network was not fully
% connected: we disregard these here if the network is fully
% connected.
% 
%% TODO CHECK If graph is fully connected
Lrand_mean = mean(Lrand(Lrand < inf));
%% 
[S_ws_MC,C_ws_MC,L_ws_MC, Cs_ws_MC, Ls_ws_MC] = custom_small_world_ness(A,Lrand_mean,mean(CrandWS),FLAG_Cws);  % Using WS clustering coefficient
[S_trans_MC,C_trans_MC,L_trans_MC, Cs_trans_MC, Ls_trans_MC] = custom_small_world_ness(A,Lrand_mean,mean(CrandTrans),FLAG_Ctransitive);  %  Using transitive clustering coefficient

headers = {'graph_method', 'clustering_coefficient_type', ... 
'mean_path_length', 'mean_clustering_coefficient', 'small_worldness_coefficient', ...
'normalized_mean_path_length', 'mean_clustering_coefficient', ...
'random_graph_path_length', 'random_graph_clustering_coefficient'};
results = {
                'analytical approximation','WS clustering coefficient',L_ws,C_ws,S_ws, Ls_ws, Cs_ws, expectedL, expectedC;
                'analytical approximation','transitive clustering coefficient',L_trans,C_trans,S_trans, Ls_trans, Cs_trans, expectedL, expectedC;
                'Monte Carlo estimate of random graph','WS clustering coefficient',L_ws_MC,C_ws_MC,S_ws_MC, Ls_ws_MC, Cs_ws_MC, Lrand,CrandWS;
                'Monte Carlo estimate of random graph','transitive clustering coefficient',L_trans_MC,C_trans_MC,S_trans_MC, Ls_trans, Cs_trans, Lrand,CrandTrans;
                };
small_worldness_results = cell2table(results);
small_worldness_results.Properties.VariableNames = headers;

%% do Monte-Carlo test of distribution of S in ER random graph
if compute_significance

[I_ws,P_ws,Sb_ws] = SsampleER(n,K,m,I,Num_S_repeats,S_ws,FLAG_Cws);  % for Sws here
[I_trans,P_trans,Sb_trans] = SsampleER(n,K,m,I,Num_S_repeats,S_trans,FLAG_Cws);  % for Strans here
[I_ws_MC,P_ws_MC,Sb_ws_MC] = SsampleER(n,K,m,I,Num_S_repeats,S_ws_MC,FLAG_Cws);  % for Sws MC here
[I_trans_MC,P_trans_MC,Sb_trans_MC] = SsampleER(n,K,m,I,Num_S_repeats,S_trans_MC,FLAG_Cws);  % for Strans MC here

small_worldness_results.confidence_interval = {I_ws, I_trans, I_ws_MC, I_trans_MC};
small_worldness_results.p_value = {P_ws, P_trans, P_ws_MC, P_trans_MC};

end

            
% NOT RETURNED FOR NOW
% check how many samples ended up being used to calculate P-value
% Nsamps = numel(Sb);
% Pmax = 1 / Nsamps;
end

