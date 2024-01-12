from bct import modularity_und
from tqdm import tqdm
import networkx as nx
import os
import bct
from sklearn.metrics import auc
import math
import pandas as pd
import numpy as np

from part2.network_analysis.network_construction_tools import to_unweighted_graph

def thresholded_auc(min_threshold, thresholds, metric_dict):
    return auc([threshold for threshold in thresholds if threshold >= min_threshold],
               [metric_dict[threshold] for threshold in thresholds if
                threshold >= min_threshold])

class ModularityCalculator:
    def __init__(self, index_to_region_mapping):
        # Modules are based on predefined Yeo 17 networks (which have been mapped to Brainnetome)
        bn_to_yeo_mapping_path = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'utils', 'BN_to_Yeo_mapping.csv')
        bn_to_yeo_mapping = pd.read_csv(bn_to_yeo_mapping_path)
        # Regions of brainnetome without a clear Yeo network assignment are assigned to network 18 (subcortical)
        bn_to_yeo_mapping['atlas2_region'] = bn_to_yeo_mapping['atlas2_region'].fillna(18)
        # if fractional region overlap < 0.5, set to 18
        bn_to_yeo_mapping.loc[bn_to_yeo_mapping['fractional_region_overlap'] < 0.5, 'atlas2_region'] = 18
        self.index_to_module_mapping = [
            bn_to_yeo_mapping.loc[
                bn_to_yeo_mapping['atlas1_region'] == index_to_region_mapping[i], 'atlas2_region'].values[
                0] for i in range(len(index_to_region_mapping))]

    def newmann_q(self, graph):
        # Compute modularity using Newman's method
        # Q=  1/2m * ∑_ij〖(A_ij-(k_i k_j)/2m) * δ_(c_i c_j )〗
        # Where Q is the modularity of the network, with m edges and an adjacency matrix A. The expected number of edges falling between two vertices i and j in the configuration model is equal to (k_i k_j)/2m, where ki is the degree of node i. ci and cj are the cluster assignments of nodes I and j respectively. δ_(c_i c_j ) is the Kronecker delta function, which equals 1 if  ci equals c , and 0 otherwise.
        _, q = modularity_und(graph, kci=self.index_to_module_mapping)
        return q

    def participation_coefficients_per_node(self, graph):
        # Compute participation coefficient per node
        pc_per_node = bct.participation_coef(graph, self.index_to_module_mapping)
        return pc_per_node

    def participation_coefficients_per_module(self, graph):
        # Compute median participation coefficient per module

        pc_per_node = self.participation_coefficients_per_node(graph)
        modules = np.unique(self.index_to_module_mapping)
        pc_per_module = {}
        for module in modules:
            pc_per_module[module] = np.median(pc_per_node[self.index_to_module_mapping == module])

        return pc_per_module


def global_efficiency(graph):
    """Compute global efficiency of a graph.

    Args:
        graph (np.ndarray): Graph.

    Returns:
        float: Global efficiency.
    """
    # compute global efficiency
    binary_reachability_matrix, distance_matrix = bct.reachdist(graph)
    characteristic_path_length, global_efficiency, eccentricity, radius, diameter = bct.charpath(distance_matrix)

    return global_efficiency


def clusttriang(A):
    """
    C = CLUSTTRIANG(A) compute the clustering coefficient C of the
    adjacency matrix A, based on the ratio (number of triangles)/(number of triples).

    Notes:
    (1) This is a subtly different definition than that of Watts & Strogatz (1998).
    In the present form, it does give the fraction of neighboring nodes
    that are also neighbors of each other. (Also note: for directed
    graphs, this means *any* direction of connection, as along as all three
    nodes are connected).

    (2) This function uses the method of Keeling (1999) to compute C.

    Original matlab function by Mark Humphries (22/8/2006)
    """

    A2 = np.dot(A, A)
    A3 = np.dot(A2, A)
    sumA2 = np.sum(A2)

    C = np.trace(A3) / (sumA2 - np.trace(A2))

    return C

def ER_Expected_L_C(k, n):
    """
    The function computes the expected clustering coefficient (expectedC) and expected path-length (expectedL) of an Erdos-Renyi random graph with n nodes and mean degree k.
    Original matlab function by Mark Humphries (3/2/2017)
    :param k: mean degree
    :param n: number of nodes
    :return: expectedC, expectedL
    """
    expectedC = k / n
    z1 = k
    z2 = k ** 2
    expectedL = (math.log((n - 1) * (z2 - z1) + z1 ** 2) - math.log(z1 ** 2)) / math.log(z2 / z1)
    return expectedC, expectedL

def small_world_ness(A, LR, CR, FLAG):
    # Compute mean shortest path length L
    _, D = bct.reachdist(A)
    L = np.mean(D.flatten())

    # Calculate required form of C
    if FLAG == 1:
        c = bct.clustering_coef_bu(A)
        C = np.mean(c)
    elif FLAG == 2:
        C = clusttriang(A)

    # Compute small-world-ness score S
    Ls = L / LR
    Cs = C / CR
    S = Cs / Ls

    return S

def small_worldness_sigma(adjacency_matrix_graph: np.ndarray, method:str='analytical', niter: int = 100, nrand: int = 5):
    """
    Returns small-worldness sigma coefficient given the adjacency matrix of a graph
    Speed: approximately 1min per iteration (1 iter / 1 nrand)
    :param adjacency_matrix_graph: adjacency matrix of given graph
    :param niter Approximate number of rewiring per edge to compute the equivalent random graph.
    :param nrand Number of random graphs generated to compute the average clustering coefficient (Cr) and average shortest path length (Lr).
    :return: sigma "small-worldness" coefficient
    """
    if method == 'analytical':
        md = np.mean(bct.degrees_und(adjacency_matrix_graph))
        n_nodes = adjacency_matrix_graph.shape[0]
        # compute clustering coefficient and path length of random graph
        Cr, Lr = ER_Expected_L_C(md, n_nodes)
        return small_world_ness(adjacency_matrix_graph, Lr, Cr, FLAG=1)
    elif method == 'monte_carlo':
        networkX_graph = nx.from_numpy_array(adjacency_matrix_graph)
        return nx.sigma(networkX_graph, niter=niter, nrand=nrand)
    else:
        raise ValueError('method must be either "analytical" or "monte_carlo"')


def order_parameter(adjacency_matrix_graph: np.ndarray) -> float:
    """
    Returns order parameter of a graph (probability of a random node being in the giant component)
    Ref: http://networksciencebook.com/chapter/8#percolation-theory (8.2.1)
    :param adjacency_matrix_graph: adjacency matrix of given graph
    :return: order parameter
    """
    networkX_graph = nx.from_numpy_array(adjacency_matrix_graph)
    num_nodes = networkX_graph.number_of_nodes()
    # find giant component
    Gcc = sorted(nx.connected_components(networkX_graph), key=len, reverse=True)
    G0 = networkX_graph.subgraph(Gcc[0])
    G0_num_nodes = G0.number_of_nodes()
    # probability of a node being in the giant component
    p = G0_num_nodes / num_nodes
    # order parameter
    return p


def analyze_connectivity_graph(connectivity_matrix: np.ndarray, minimum_connectivity_threshold: float = 0.3,
                                index_to_region_mapping:np.ndarray=None,
                               binned_thresholding: bool = False,
                               compute_smallwordness: bool = 0, sigma_niter: int = 100,
                               sigma_nrand: int = 5,
                               compute_participation_coefficients: bool = 0,
                               ) -> list:
    """Analyze connectivity matrix, transform into a graph at multiple thresholds and analyze each graph. For each metric, the AUC over all thresholds is returned.

    Steps:
    - Take absolute value of connectivity matrix
    - Transform into graph at multiple thresholds (0 - 1)
    - For each graph, compute the following metrics:
        - Degree
        - Clustering
        - Global efficiency
    - Compute AUC for each metric over all thresholds

    Args:
        connectivity_matrix (np.ndarray): Connectivity matrix.
        minimum_connectivity_threshold (float, optional): Minimum threshold to include for AUC computation (default: 0.3, i.e. [0.3-1]).
        index_to_region_mapping (np.ndarray, optional): Mapping from index to region name (default: None).
        binned_thresholding
        compute_smallwordness: bool = 0
            - 0: do not compute small-worldness
            - 1: compute small-worldness using analytical formula
            - 2: compute small-worldness using Monte Carlo simulation (slow)

    """
    # compute overall functional connectivity (FC): mean of all positive values across all elements of the  matrix
    temp_CM = connectivity_matrix.copy()
    temp_CM[temp_CM < 0] = 0
    overall_functional_connectivity = np.nanmean(temp_CM)

    # transform into graph at multiple thresholds [0.1, 0.2, ..., 1]
    thresholds = np.arange(0.1, 1.1, 0.1)
    # only keep one decimal place
    thresholds = np.around(thresholds, decimals=1)
    graphs = [to_unweighted_graph(connectivity_matrix, threshold, binned_thresholding=binned_thresholding) for threshold
              in thresholds]
    # build dictionary with graph and threshold
    graphs = dict(zip(thresholds, graphs))

    # compute metrics for each graph
    # compute mean & median degree with degrees_und()
    mean_degrees = {threshold: np.mean(bct.degrees_und(graph)) for threshold, graph in graphs.items()}
    median_degrees = {threshold: np.median(bct.degrees_und(graph)) for threshold, graph in graphs.items()}

    # compute mean & median clustering coefficient with clustering_coef_bu()
    mean_clustering_coefficients = {threshold: np.mean(bct.clustering_coef_bu(graph)) for threshold, graph in
                                    graphs.items()}
    median_clustering_coefficients = {threshold: np.median(bct.clustering_coef_bu(graph)) for threshold, graph in
                                      graphs.items()}

    # compute betweenness centrality with betweenness_bin()
    mean_betweenness_centrality = {threshold: np.mean(bct.betweenness_bin(graph)) for threshold, graph in
                                    graphs.items()}
    median_betweenness_centrality = {threshold: np.median(bct.betweenness_bin(graph)) for threshold, graph in
                                        graphs.items()}

    # compute global efficiency
    global_efficiencies = {threshold: global_efficiency(graph) for threshold, graph in graphs.items()}

    # compute order parameter
    order_parameters = {threshold: order_parameter(graph) for threshold, graph in graphs.items()}

    # compute modularity
    modularity_calculator = ModularityCalculator(index_to_region_mapping)
    modularities = {threshold: modularity_calculator.newmann_q(graph) for threshold, graph in graphs.items()}

    if compute_participation_coefficients:
        # compute participation coefficients
        participation_coefficients_per_node = {threshold: modularity_calculator.participation_coefficients_per_node(graph) for threshold, graph in graphs.items()}
        participation_coefficients_per_module = {threshold: modularity_calculator.participation_coefficients_per_module(graph) for threshold, graph in graphs.items()}

    # compute small-worldness sigma
    if compute_smallwordness:
        small_worldness_sigmas = {}
        for threshold in tqdm(thresholds, "Calculating small worldness sigma"):
            if threshold < minimum_connectivity_threshold:
                continue

            if compute_smallwordness == 1:
                # use analytical formula
                small_worldness_sigmas[threshold] = small_worldness_sigma(graphs[threshold], niter=sigma_niter,
                                                                          nrand=sigma_nrand)
            elif compute_smallwordness == 2:
                # use Monte Carlo simulation
                small_worldness_sigmas[threshold] = small_worldness_sigma(graphs[threshold], method='monte_carlo',
                                                                          niter=sigma_niter, nrand=sigma_nrand)

    # compute AUC for each metric over all thresholds above minimum_connectivity_threshold
    mean_degree_auc = thresholded_auc(minimum_connectivity_threshold, thresholds, mean_degrees)
    median_degree_auc = thresholded_auc(minimum_connectivity_threshold, thresholds, median_degrees)
    mean_betweenness_centrality_auc = thresholded_auc(minimum_connectivity_threshold, thresholds,
                                                      mean_betweenness_centrality)
    median_betweenness_centrality_auc = thresholded_auc(minimum_connectivity_threshold, thresholds,
                                                        median_betweenness_centrality)

    mean_clustering_coefficient_auc = thresholded_auc(minimum_connectivity_threshold, thresholds,
                                                        mean_clustering_coefficients)
    median_clustering_coefficient_auc = thresholded_auc(minimum_connectivity_threshold, thresholds,
                                                            median_clustering_coefficients)
    global_efficiency_auc = thresholded_auc(minimum_connectivity_threshold, thresholds, global_efficiencies)
    order_parameter_auc = thresholded_auc(minimum_connectivity_threshold, thresholds, order_parameters)
    modularity_auc = thresholded_auc(minimum_connectivity_threshold, thresholds, modularities)

    if compute_smallwordness:
        small_worldness_sigma_auc = thresholded_auc(minimum_connectivity_threshold, thresholds, small_worldness_sigmas)
    else:
        small_worldness_sigma_auc = np.nan

    if compute_participation_coefficients:
        # compute aucs over thresholds for all nodes
        participation_coefficients_auc_per_node = (pd.DataFrame(participation_coefficients_per_node)
                                                   .apply(lambda x: thresholded_auc(minimum_connectivity_threshold,
                                                                                    thresholds, x.to_dict()), axis=1))
        # compute aucs over thresholds for all modules
        participation_coefficients_auc_per_module = (pd.DataFrame(participation_coefficients_per_module)
                                                     .apply(lambda x: thresholded_auc(minimum_connectivity_threshold,
                                                                                      thresholds, x.to_dict()), axis=1))
    else:
        participation_coefficients_auc_per_node = np.nan
        participation_coefficients_auc_per_module = np.nan

    return graphs, mean_degree_auc, median_degree_auc, \
        mean_betweenness_centrality_auc, median_betweenness_centrality_auc, \
          mean_clustering_coefficient_auc, median_clustering_coefficient_auc, \
        global_efficiency_auc, global_efficiencies, \
        order_parameter_auc, \
        modularity_auc, \
        overall_functional_connectivity, small_worldness_sigma_auc, \
        participation_coefficients_auc_per_node, participation_coefficients_auc_per_module
