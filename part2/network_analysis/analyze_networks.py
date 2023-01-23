import argparse
import os
import warnings

import pandas as pd
import scipy.io as sio

from part2.network_analysis.network_analysis_tools import analyze_connectivity_graph


def analyze_networks(data_dir:str, matrix_name: str = 'CM3D_z_norm', minimum_connectivity_threshold: float = 0.3,
                    binned_thresholding:bool = False,
                    compute_smallwordness:bool = False, sigma_niter:int = 100, sigma_nrand: int = 5,
                     connectivity_file_prefix:str = 'filtered_masked_transfer_') -> pd.DataFrame:
    """Analyze networks for all subjects in data_dir.

    Requires that the data_dir contains a subdirectory for each subject, and that each subject directory contains a .mat file with the connectivity matrix.

    Args:
        data_dir (str): Path to directory containing subject data.
        matrix_name (str, optional): Name of matrix to analyze (default: 'CM3D_z_norm').
        minimum_connectivity_threshold (float, optional): Minimum threshold to include for AUC computation (default: 0.3, i.e. [0.3-1]).
        connectivity_file_prefix (str, optional): Prefix of connectivity file (default: 'filtered_masked_transfer_').
        Smallworldness (sigma) parameters:
                niter Approximate number of rewiring per edge to compute the equivalent random graph.
                nrand Number of random graphs generated to compute the average clustering coefficient (Cr) and average shortest path length (Lr).

    Returns:
        pd.DataFrame: Dataframe containing all metrics for all subjects.
    """
    # get list of subjects
    subjects = [subject for subject in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, subject))]

    # initialize output_df
    output_df = pd.DataFrame(columns=['subject',
                                   'mean_degree_auc', 'median_degree_auc', 'mean_clustering_coefficient_auc',
                                  'median_clustering_coefficient_auc', 'global_efficiency_auc', 'overall_functional_connectivity'])

    # loop over subjects
    for subject in subjects:
        # find the connectivity matrix path
        connectivity_matrix_path_possibilities = [os.path.join(data_dir, subject, file) for file in os.listdir(os.path.join(data_dir, subject)) if file.startswith(connectivity_file_prefix) and file.endswith('.mat')]
        if len(connectivity_matrix_path_possibilities) == 0:
            warnings.warn(f'Could not find connectivity matrix for subject {subject}. Skipping subject.')
            continue
        elif len(connectivity_matrix_path_possibilities) > 1:
            warnings.warn(f'Multiple connectivity matrices found for subject {subject}, retaining only {os.path.basename(connectivity_matrix_path_possibilities[0])}')
        connectivity_matrix_path = connectivity_matrix_path_possibilities[0]

        connectivity_matrix = sio.loadmat(connectivity_matrix_path)[matrix_name]
        mean_degree_auc, median_degree_auc, mean_clustering_coefficient_auc, median_clustering_coefficient_auc, global_efficiency_auc, overall_functional_connectivity, small_worldness_sigma_auc = analyze_connectivity_graph(connectivity_matrix, minimum_connectivity_threshold,
                                                                                                                                                                                                    binned_thresholding=binned_thresholding,
                                                                                                                                                                                                   compute_smallwordness=compute_smallwordness,
                                                                                                                                                                                                    sigma_niter=sigma_niter, sigma_nrand=sigma_nrand)

        output_df = pd.concat([output_df, pd.DataFrame({'subject': subject,
                                'mean_degree_auc': mean_degree_auc, 'median_degree_auc': median_degree_auc,
                                'mean_clustering_coefficient_auc': mean_clustering_coefficient_auc,
                                'median_clustering_coefficient_auc': median_clustering_coefficient_auc,
                                'global_efficiency_auc': global_efficiency_auc,
                                'small_worldness_sigma_auc': small_worldness_sigma_auc,
                                'overall_functional_connectivity': overall_functional_connectivity},
                                 index=[0])], ignore_index=True)

    # add columns 'matrix_name', 'minimum_connectivity_threshold', 'connectivity_file_prefix' to output_df
    output_df['matrix_name'] = matrix_name
    output_df['sigma_parameters (niter, nrand)'] = f'({sigma_niter}, {sigma_nrand})'
    output_df['minimum_connectivity_threshold'] = minimum_connectivity_threshold
    output_df['binned_thesholding'] = binned_thresholding
    output_df['connectivity_file_prefix'] = connectivity_file_prefix
    output_df['connectivity_file_name'] = os.path.basename(connectivity_matrix_path)

    return output_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze networks for all subjects in data_dir.')
    parser.add_argument('-i', '--input_data_dir', type=str, help='Path to directory containing subject data.', required=True)
    parser.add_argument('-m', '--matrix_name', type=str, help='Name of matrix to analyze (default: CM3D_z_norm).', default='CM3D_z_norm')
    parser.add_argument('-t', '--minimum_connectivity_threshold', type=float, help='Minimum inclusive threshold to include for AUC computation (default: 0.3, i.e. [0.3-1]).', default=0.3)
    parser.add_argument('-b', '--binned_thresholding', required=False, action='store_true', help='use binned proportionnal threshold where each bin contains the links between pX and pX + 0.1')
    parser.add_argument('-p', '--connectivity_file_prefix', type=str, help='Prefix of connectivity file (default: filtered_masked_transfer_).', default='filtered_masked_transfer_')
    parser.add_argument('-o', '--output_dir', type=str, help='Path to output directory.', default='')
    # smallwordness arguments
    parser.add_argument('-s', '--smallworldness', required=False, action='store_true', help='compute smallwordness (very slow)', default=False)
    parser.add_argument('-snt', '--sigma_niter', type=int, help='Approximate number of rewiring per edge to compute the equivalent random graph.', default=100)
    parser.add_argument('-snr', '--sigma_nrand', type=int, help='Number of random graphs generated to compute the average clustering coefficient (Cr) and average shortest path length (Lr).', default=5)

    args = parser.parse_args()

    if args.output_dir == '':
        args.output_dir = args.input_data_dir

    output_df = analyze_networks(args.input_data_dir, args.matrix_name, args.minimum_connectivity_threshold, args.binned_thresholding,
                                 args.smallworldness, args.sigma_niter, args.sigma_nrand,
                                 args.connectivity_file_prefix)
    output_df.to_csv(os.path.join(args.output_dir, f'{args.matrix_name}_network_analysis.csv'), index=False)