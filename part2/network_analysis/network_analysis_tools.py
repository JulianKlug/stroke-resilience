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


def analyze_connectivity_graph(connectivity_matrix: np.ndarray, minimum_connectivity_threshold: float = 0.3,
                               binned_thresholding: bool = False) -> None:
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
    graphs = [to_unweighted_graph(connectivity_matrix, threshold, binned_thresholding=binned_thresholding) for threshold in thresholds]
    # build dictionary with graph and threshold
    graphs = dict(zip(thresholds, graphs))

    # compute metrics for each graph
    # compute mean & median degree with degrees_und()
    mean_degrees = {threshold: np.mean(bct.degrees_und(graph)) for threshold, graph in graphs.items()}
    median_degrees = {threshold: np.median(bct.degrees_und(graph)) for threshold, graph in graphs.items()}

    # compute mean & median clustering coefficient with clustering_coef_bu()
    mean_clustering_coefficients = {threshold: np.mean(bct.clustering_coef_bu(graph)) for threshold, graph in graphs.items()}
    median_clustering_coefficients = {threshold: np.median(bct.clustering_coef_bu(graph)) for threshold, graph in graphs.items()}

    # compute global efficiency
    global_efficiencies = {threshold: global_efficiency(graph) for threshold, graph in graphs.items()}

    # compute AUC for each metric over all thresholds above minimum_connectivity_threshold
    mean_degree_auc = auc([threshold for threshold in thresholds if threshold >= minimum_connectivity_threshold],
                            [mean_degrees[threshold] for threshold in thresholds if threshold >= minimum_connectivity_threshold])
    median_degree_auc = auc([threshold for threshold in thresholds if threshold >= minimum_connectivity_threshold],
                            [median_degrees[threshold] for threshold in thresholds if threshold >= minimum_connectivity_threshold])
    mean_clustering_coefficient_auc = auc([threshold for threshold in thresholds if threshold >= minimum_connectivity_threshold],
                            [mean_clustering_coefficients[threshold] for threshold in thresholds if threshold >= minimum_connectivity_threshold])
    median_clustering_coefficient_auc = auc([threshold for threshold in thresholds if threshold >= minimum_connectivity_threshold],
                            [median_clustering_coefficients[threshold] for threshold in thresholds if threshold >= minimum_connectivity_threshold])
    global_efficiency_auc = auc([threshold for threshold in thresholds if threshold >= minimum_connectivity_threshold],
                                [global_efficiencies[threshold] for threshold in thresholds if threshold >= minimum_connectivity_threshold])

    return mean_degree_auc, median_degree_auc, mean_clustering_coefficient_auc, median_clustering_coefficient_auc, global_efficiency_auc, overall_functional_connectivity



