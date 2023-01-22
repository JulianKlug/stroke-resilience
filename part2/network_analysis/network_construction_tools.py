import numpy as np
import bct

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


def to_unweighted_graph(connectivity_matrix, threshold, binned_thresholding=False):
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
        binned_thresholding (bool)

    Returns:
        np.ndarray: Graph.
    """

    # set nan to 0
    connectivity_matrix[np.isnan(connectivity_matrix)] = 0
    # Rubinov M, Sporns O (2010) NeuroImage 52:1059-69: "all self-connections or negative connections (such as functional anticorrelations) must currently be removed from the networks prior to analysis"
    pos_connectivity_matrix = bct.threshold_absolute(connectivity_matrix, 0, copy=True)

    if binned_thresholding:
        thresholded_graph = threshold_by_proportional_bin(pos_connectivity_matrix, threshold, n_total_bins=10, copy=True)
    else:
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

def threshold_by_proportional_bin(W, p, n_total_bins=10, copy=True):
    '''
    This function "thresholds" the connectivity matrix by preserving a
    bin corresponding to a proportion p (0<p<1) to p + 1/n_total_bins of the strongest weights. All other weights, and
    all weights on the main diagonal (self-self connections) are set to 0.

    If copy is not set, this function will *modify W in place.*

    Parameters
    ----------
    W : np.ndarray
        weighted connectivity matrix
    p : float
        proportional weight threshold (0<p<1)
    n_total_bins : int
        number of bins to use for thresholding
    copy : bool
        if True, returns a copy of the matrix. Otherwise, modifies the matrix
        in place. Default value=True.

    Returns
    -------
    W : np.ndarray
        bin-thresholded connectivity matrix
    '''
    from bct.utils.miscellaneous_utilities import teachers_round as round

    if p == 0.9:
        print()

    if p > 1 or p < 0:
        raise ValueError('Threshold must be in range [0,1]')
    if copy:
        W = W.copy()
    n = len(W)						# number of nodes
    np.fill_diagonal(W, 0)			# clear diagonal

    if np.allclose(W, W.T):				# if symmetric matrix
        W[np.tril_indices(n)] = 0		# ensure symmetry is preserved
        ud = 2						# halve number of removed links
    else:
        ud = 1

    ind = np.where(W)
    # of all non zero elements in W, retain only those between percentile p and p+1/n_total_bins
    weights_with_initial_indices = np.vstack((W[ind], ind)).T
    sorted_weights_with_initial_indices = weights_with_initial_indices[weights_with_initial_indices[:, 0].argsort()]
    number_of_links_per_bin = round(len(sorted_weights_with_initial_indices) / n_total_bins)
    # lower bound of the bin is indicated by the percentile p
    lower_bound = int(p*n_total_bins*number_of_links_per_bin)
    # upper bound of the bin is indicated by the percentile p+1/n_total_bins
    upper_bound = int(lower_bound + number_of_links_per_bin)
    # set elements outside of bounds to zero
    sorted_weights_with_initial_indices[0:lower_bound, 0] = 0
    sorted_weights_with_initial_indices[upper_bound:, 0] = 0

    # replaces values in original form of matrix
    for idx in range(len(sorted_weights_with_initial_indices)):
        W[int(sorted_weights_with_initial_indices[idx, 1]), int(sorted_weights_with_initial_indices[idx, 2])] = sorted_weights_with_initial_indices[idx, 0]

    if ud == 2:						# if symmetric matrix
        W[:, :] = W + W.T						# reconstruct symmetry

    return W


def test_threshold_by_proportional_bin():
    W = np.array(
        [[0, 54, 50, 73, 0],
         [96, 0, 17, 81, 84],
         [1, 65, 0, 62, 25],
         [21, 36, 39, 0, 71],
         [82, 44, 59, 89, 0]]
    )

    W = threshold_by_proportional_bin(W, 0.5)
    assert np.allclose(W, np.array([[ 0,  0,  0,  0,  0],
                                     [ 0,  0,  0,  0,  0],
                                     [ 0, 65,  0, 62,  0],
                                     [ 0,  0,  0,  0,  0],
                                     [ 0,  0,  0,  0,  0]]))
