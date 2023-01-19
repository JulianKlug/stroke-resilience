import numpy as np
import bct
from sklearn.metrics import auc

"""
Matrix distributions
CM3D: min: -0.73, max: 1.00
CM3D_z: min: -0.93, max: 7.25
CM3D_z_norm: min: -0.00, max: 0.00
CM_mu: min: -0.73, max: 1.00
CM_mu_z: min: -0.93, max: 7.25
CM_std: min: 0.00, max: 0.00
CM_std_z: min: 0.00, max: 0.00
"""


def to_unweighted_graph(connectivity_matrix, threshold):
    """Transform connectivity matrix into an unweighted graph at a given threshold of proportional thresholding.

    Steps:
    - Threshold connectivity matrix at given threshold, through proportional thresholding (as implemented in BCT)
        - Gist: retain only the top t% of the weights in the matrix, where t is the threshold
        - Note: For issues with regard to proportional thresholding, see Martijn P. van den Heuvel, Siemon C. de Lange, Andrew Zalesky, Caio Seguin, B.T. Thomas Yeo, Ruben Schmidt, Proportional thresholding in resting-state fMRI functional connectivity networks and consequences for patient-control connectome studies: Issues and recommendations, NeuroImage, Volume 152, 2017, 437-449, ISSN 1053-8119, https://doi.org/10.1016/j.neuroimage.2017.02.005.
    - Transform into unweighted graph by binarizing graph (as implemented in BCT)
    - Fix common problems (as implemented in BCT): remove Inf and NaN, ensure exact binariness and symmetry (i.e. remove floating point instability), and zero diagonal.

    Args:
        connectivity_matrix (np.ndarray): Connectivity matrix.
        threshold (float): Threshold for proportional thresholding - from t=0 (no connection) to t=1 (all connections retained)

    Returns:
        np.ndarray: Graph.
    """

    # set nan to 0
    connectivity_matrix[np.isnan(connectivity_matrix)] = 0
    # Rubinov M, Sporns O (2010) NeuroImage 52:1059-69: "all self-connections or negative connections (such as functional anticorrelations) must currently be removed from the networks prior to analysis"
    pos_connectivity_matrix = bct.threshold_absolute(connectivity_matrix, 0, copy=True)
    thresholded_graph = bct.threshold_proportional(pos_connectivity_matrix, threshold, copy=True)
    binarized_graph = bct.weight_conversion(thresholded_graph, 'binarize', copy=True)
    autofixed_graph = autofix(binarized_graph, copy=True)

    return autofixed_graph


def autofix(W, copy=True):
    '''
    This is a local version of the BCT function autofix. It is a copy of the original function with the following changes:
    - corrected setting to zero nans & infs

    Fix a bunch of common problems. More specifically, remove Inf and NaN,
    ensure exact binariness and symmetry (i.e. remove floating point
    instability), and zero diagonal.


    Parameters
    ----------
    W : np.ndarray
        weighted connectivity matrix
    copy : bool
        if True, returns a copy of the matrix. Otherwise, modifies the matrix
        in place. Default value=True.

    Returns
    -------
    W : np.ndarray
        connectivity matrix with fixes applied
    '''
    if copy:
        W = W.copy()

    # zero diagonal
    np.fill_diagonal(W, 0)

    # remove np.inf and np.nan
    W[np.where(np.isinf(W))] = 0
    W[np.where(np.isnan(W))] = 0

    # ensure exact binarity
    u = np.unique(W)
    if np.all(np.logical_or(np.abs(u) < 1e-8, np.abs(u - 1) < 1e-8)):
        W = np.around(W, decimals=5)

    # ensure exact symmetry
    if np.allclose(W, W.T):
        W = np.around(W, decimals=5)

    return W


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


def analyze_connectivity_graph(connectivity_matrix: np.ndarray, minimum_connectivity_threshold: float = 0.3) -> None:
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

    """
    # compute overall functional connectivity (FC): mean of all positive values across all elements of the  matrix
    temp_CM = connectivity_matrix.copy()
    temp_CM[temp_CM < 0] = 0
    overall_functional_connectivity = np.nanmean(temp_CM)

    # transform into graph at multiple thresholds [0.1, 0.2, ..., 1]
    thresholds = np.arange(0.1, 1.1, 0.1)
    # only keep one decimal place
    thresholds = np.around(thresholds, decimals=1)
    graphs = [to_unweighted_graph(connectivity_matrix, threshold) for threshold in thresholds]
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



