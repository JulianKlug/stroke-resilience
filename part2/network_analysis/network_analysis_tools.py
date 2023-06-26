from tqdm import tqdm
import networkx as nx
import numpy as np
import bct
from sklearn.metrics import auc

from part2.network_analysis.network_construction_tools import to_unweighted_graph


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


def small_worldness_sigma(adjacency_matrix_graph: np.ndarray, niter: int = 100, nrand: int = 5):
    """
    Returns small-worldness sigma coefficient given the adjacency matrix of a graph
    Speed: approximately 1min per iteration (1 iter / 1 nrand)
    :param adjacency_matrix_graph: adjacency matrix of given graph
    :param niter Approximate number of rewiring per edge to compute the equivalent random graph.
    :param nrand Number of random graphs generated to compute the average clustering coefficient (Cr) and average shortest path length (Lr).
    :return: sigma "small-worldness" coefficient
    """
    networkX_graph = nx.from_numpy_array(adjacency_matrix_graph)
    return nx.sigma(networkX_graph, niter=niter, nrand=nrand)

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
                               binned_thresholding: bool = False,
                               compute_smallwordness: bool = False, sigma_niter: int = 100,
                               sigma_nrand: int = 5) -> list:
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
        binned_thresholding

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

    # compute small-worldness sigma
    if compute_smallwordness:
        small_worldness_sigmas = {}
        for threshold in tqdm(thresholds, "Calculating small worldness sigma"):
            if threshold < minimum_connectivity_threshold:
                continue
            small_worldness_sigmas[threshold] = small_worldness_sigma(graphs[threshold], niter=sigma_niter,
                                                                      nrand=sigma_nrand)

    # compute AUC for each metric over all thresholds above minimum_connectivity_threshold
    def thresholded_auc(min_threshold, thresholds, metric_dict):
        return auc([threshold for threshold in thresholds if threshold >= min_threshold],
                          [metric_dict[threshold] for threshold in thresholds if
                           threshold >= min_threshold])

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

    if compute_smallwordness:
        small_worldness_sigma_auc = thresholded_auc(minimum_connectivity_threshold, thresholds, small_worldness_sigmas)
    else:
        small_worldness_sigma_auc = np.nan

    return mean_degree_auc, median_degree_auc, \
        mean_betweenness_centrality_auc, median_betweenness_centrality_auc, \
          mean_clustering_coefficient_auc, median_clustering_coefficient_auc, \
        global_efficiency_auc, global_efficiencies, \
        order_parameter_auc, \
        overall_functional_connectivity, small_worldness_sigma_auc
