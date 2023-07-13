import argparse
import os
import warnings
from pathlib import Path
from tqdm import tqdm
import pandas as pd
import numpy as np
import scipy.io as sio

from part2.network_analysis.network_analysis_tools import analyze_connectivity_graph
from part2.utils.utils import ensure_dir


def analyze_networks(data_dir:str, matrix_name: str = 'CM3D_z_norm', minimum_connectivity_threshold: float = 0.3,
                    binned_thresholding:bool = False,
                    compute_smallwordness:bool = False, sigma_niter:int = 100, sigma_nrand: int = 5,
                    save_graphs:bool = False,
                    connectivity_file_prefix:str = 'filtered_masked_transfer_',
                    control_folder_prefix:str = 'amc',
                    allow_multiple_connectivity_matrices_per_subject:bool=True,
                    restrict_to_subject_subdir:str=None) -> pd.DataFrame:
    """Analyze networks for all subjects in data_dir.

    Requires that the data_dir contains a subdirectory for each subject, and that each subject directory contains at least one .mat file with the connectivity matrix.

    Args:
        data_dir (str): Path to directory containing subject data.
        matrix_name (str, optional): Name of matrix to analyze (default: 'CM3D_z_norm').
        minimum_connectivity_threshold (float, optional): Minimum threshold to include for AUC computation (default: 0.3, i.e. [0.3-1]).
        connectivity_file_prefix (str, optional): Prefix of connectivity file (default: 'filtered_masked_transfer_').
        allow_multiple_connectivity_matrices_per_subject (bool, optional): Whether to allow multiple connectivity matrices per subject (default: True).
        restrict_to_subject_subdir (str, optional): If not None, only consider connectivity matrices in the subdirectory named restrict_to_subject_subdir with a supposed structure of subj/../restrict_to_subject_dir/X/CM.mat (default: None).
        save_graphs (bool, optional): Whether to save the undirected thresholded graphs computet from every connectivity file (default: False).
        Smallworldness (sigma) parameters:
                niter Approximate number of rewiring per edge to compute the equivalent random graph.
                nrand Number of random graphs generated to compute the average clustering coefficient (Cr) and average shortest path length (Lr).

    Returns:
        pd.DataFrame: Dataframe containing all metrics for all subjects.
    """
    # Constants
    index_to_region_mapping_name = 'allCodeBooks'

    # get list of subjects
    subjects = [subject for subject in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, subject))]

    # initialize output_df
    output_df = pd.DataFrame(columns=['subject', 'subject_type', 'subject_id', 'subject_timepoint', 'connectivity_file_name'])

    global_efficiencies_df = pd.DataFrame()

    # loop over subjects
    for subject in tqdm(subjects):
        # find the paths for all possible connectivity matrices
        connectivity_matrix_path_possibilities = [str(path) for path in Path(os.path.join(data_dir, subject))
                                                                .rglob(f'{connectivity_file_prefix}*.mat')]
        if restrict_to_subject_subdir:
            # ensure that the connectivity matrix is in the subject subdirectory named restrict_to_subject_subdir
            connectivity_matrix_path_possibilities = [path for path in connectivity_matrix_path_possibilities if os.path.dirname(os.path.dirname(path)).split('/')[-1] == restrict_to_subject_subdir]

        if not connectivity_matrix_path_possibilities:
            warnings.warn(f'Could not find connectivity matrix for subject {subject}. Skipping subject.')
            continue
        elif len(connectivity_matrix_path_possibilities) > 1 and not allow_multiple_connectivity_matrices_per_subject:
            warnings.warn(f'Multiple connectivity matrices found for subject {subject}, retaining only {os.path.basename(connectivity_matrix_path_possibilities[0])}')
            connectivity_matrix_path_possibilities = [connectivity_matrix_path_possibilities[0]]

        if subject.startswith(control_folder_prefix):
            subject_type = 'control'
            subject_id = subject.split('_')[1]
            subject_timepoint = subject.split('_')[3]
        else:
            subject_type = 'patient'
            subject_id = subject.split('_')[1]
            subject_timepoint = subject.split('_')[2]

        for connectivity_matrix_path in connectivity_matrix_path_possibilities:
            connectivity_matrix_object = sio.loadmat(connectivity_matrix_path)
            connectivity_matrix = connectivity_matrix_object[matrix_name]
            index_to_region_mapping = np.squeeze(connectivity_matrix_object[index_to_region_mapping_name].T)

            graphs, mean_degree_auc, median_degree_auc, \
                mean_betweenness_centrality_auc, median_betweenness_centrality_auc, \
                mean_clustering_coefficient_auc, median_clustering_coefficient_auc, \
                global_efficiency_auc, global_efficiencies, \
                order_parameter_auc, \
                modularity_auc, \
                overall_functional_connectivity, small_worldness_sigma_auc = analyze_connectivity_graph(connectivity_matrix, minimum_connectivity_threshold, index_to_region_mapping,
                                                                                                                                                                                                        binned_thresholding=binned_thresholding,
                                                                                                                                                                                      compute_smallwordness=compute_smallwordness,
                                                                                                                                                                                            sigma_niter=sigma_niter, sigma_nrand=sigma_nrand)

            if save_graphs:
                # transform keys of graphs dictionary to strings (matlab does not accept . in keys)
                graphs_str = {f't{int(threshold * 10)}': graph for threshold, graph in graphs.items()}
                graph_dir = os.path.join(data_dir, subject, 'undirected_thresholded_graphs')
                ensure_dir(graph_dir)
                sio.savemat(os.path.join(graph_dir, f'udt_graphs_{os.path.basename(connectivity_matrix_path)}'), graphs_str)

            output_df = pd.concat([output_df, pd.DataFrame({'subject': subject,
                                    'subject_type': subject_type,
                                    'subject_id': subject_id,
                                    'subject_timepoint': subject_timepoint,
                                    'connectivity_file_name': os.path.basename(connectivity_matrix_path),
                                    'mean_degree_auc': mean_degree_auc,
                                    'median_degree_auc': median_degree_auc,
                                    'mean_clustering_coefficient_auc': mean_clustering_coefficient_auc,
                                    'median_clustering_coefficient_auc': median_clustering_coefficient_auc,
                                    'mean_betweenness_centrality_auc': mean_betweenness_centrality_auc,
                                    'median_betweenness_centrality_auc': median_betweenness_centrality_auc,
                                    'global_efficiency_auc': global_efficiency_auc,
                                    'order_parameter_auc': order_parameter_auc,
                                    'modularity_auc': modularity_auc,
                                    'small_worldness_sigma_auc': small_worldness_sigma_auc,
                                    'overall_functional_connectivity': overall_functional_connectivity},
                                     index=[0])], ignore_index=True)

            global_efficiencies_subj_df = pd.DataFrame(global_efficiencies, index=[0]).melt(var_name='threshold', value_name='global_efficiency')
            global_efficiencies_subj_df['subject'] = subject
            global_efficiencies_subj_df['subject_type'] = subject_type
            global_efficiencies_subj_df['subject_id'] = subject_id
            global_efficiencies_subj_df['connectivity_file_name'] = os.path.basename(connectivity_matrix_path)
            global_efficiencies_subj_df['subject_timepoint'] = subject_timepoint
            global_efficiencies_df = pd.concat([global_efficiencies_df, global_efficiencies_subj_df], ignore_index=True)

    # add columns 'matrix_name', 'minimum_connectivity_threshold', 'connectivity_file_prefix' to output_df
    output_df['matrix_name'] = matrix_name
    output_df['sigma_parameters (niter, nrand)'] = f'({sigma_niter}, {sigma_nrand})'
    output_df['minimum_connectivity_threshold'] = minimum_connectivity_threshold
    output_df['binned_thesholding'] = binned_thresholding
    output_df['connectivity_file_prefix'] = connectivity_file_prefix

    return output_df, global_efficiencies_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze networks for all subjects in data_dir.')
    parser.add_argument('-i', '--input_data_dir', type=str, help='Path to directory containing subject data.', required=True)
    parser.add_argument('-r', '--restrict_to_subject_subdir', type=str, help='Restrict the search of CMs to a specific subdirectory in each subject directory, with a supposed structure of subj/../restrict_to_subject_dir/X/CM.mat (default: None).', default=None)
    parser.add_argument('-m', '--matrix_name', type=str, help='Name of matrix to analyze (default: CM3D_z_norm).', default='CM3D_z_norm')
    parser.add_argument('-t', '--minimum_connectivity_threshold', type=float, help='Minimum inclusive threshold to include for AUC computation (default: 0.3, i.e. [0.3-1]).', default=0.3)
    parser.add_argument('-b', '--binned_thresholding', required=False, action='store_true', help='use binned proportionnal threshold where each bin contains the links between pX and pX + 0.1')
    parser.add_argument('-p', '--connectivity_file_prefix', type=str, help='Prefix of connectivity file (default: filtered_masked_transfer_).', default='filtered_masked_transfer_')
    parser.add_argument('-c', '--control_folder_prefix', type=str, help='Prefix of control folder (default: amc).', default='amc')
    parser.add_argument('-o', '--output_dir', type=str, help='Path to output directory.', default='')
    parser.add_argument('-nm', '--do_not_allow_multiple_connectivity_matrices_per_subject', required=False, action='store_true', help='If set, will raise an error if multiple connectivity matrices are found for a subject.', default=False)
    parser.add_argument('-sg', '--save_graphs', required=False, action='store_true', help='If set, will save the graphs in matlab format.', default=False)
    # smallwordness arguments
    parser.add_argument('-ana_s', '--analytical_smallworldness', required=False, action='store_true', help='compute smallwordness analytically', default=False)
    parser.add_argument('-mc_s', '--montecarlo_smallworldness', required=False, action='store_true', help='compute smallwordness by monte carlo simulation (very slow)', default=False)
    parser.add_argument('-snt', '--sigma_niter', type=int, help='Approximate number of rewiring per edge to compute the equivalent random graph (Monte Carlo).', default=100)
    parser.add_argument('-snr', '--sigma_nrand', type=int, help='Number of random graphs generated to compute the average clustering coefficient (Cr) and average shortest path length (Lr) (Monte Carlo).', default=5)

    args = parser.parse_args()

    if args.output_dir == '':
        args.output_dir = args.input_data_dir

    smallworldness_flag = 0
    if args.analytical_smallworldness:
        smallworldness_flag = 1
    if args.montecarlo_smallworldness:
        smallworldness_flag = 2
    if args.analytical_smallworldness and args.montecarlo_smallworldness:
        raise ValueError('Cannot compute smallworldness analytically and by monte carlo simulation at the same time.')

    output_df, global_efficiencies_df = analyze_networks(args.input_data_dir, args.matrix_name, args.minimum_connectivity_threshold, args.binned_thresholding,
                                 smallworldness_flag, args.sigma_niter, args.sigma_nrand,
                                    args.save_graphs,
                                 args.connectivity_file_prefix,
                                    args.control_folder_prefix,
                                 allow_multiple_connectivity_matrices_per_subject=not args.do_not_allow_multiple_connectivity_matrices_per_subject,
                                 restrict_to_subject_subdir=args.restrict_to_subject_subdir)
    output_df = output_df.sort_values(by=['subject_type', 'subject_id', 'subject_timepoint', 'connectivity_file_name'])
    global_efficiencies_df = global_efficiencies_df.sort_values(by=['subject_type', 'subject_id', 'subject_timepoint',  'connectivity_file_name'])

    output_df.to_csv(os.path.join(args.output_dir, f'{args.matrix_name}_network_analysis.csv'), index=False)
    global_efficiencies_df.to_csv(os.path.join(args.output_dir, f'{args.matrix_name}_global_efficiencies.csv'), index=False)