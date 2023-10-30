import os
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
import scipy.io as sio

from part2.network_analysis.network_analysis_tools import thresholded_auc


def compute_small_worldness_auc(data_dir:str,
                                minimum_connectivity_threshold: float = 0.3,
                                restrict_to_prefix:str = 'udt_graphs_',
                                small_worldness_sigma_suffix:str='_sw_sigmas.mat',
                                control_folder_prefix:str = 'amc',
                                verbose:bool= True) -> pd.DataFrame:
    """Compute AUC for pre-computed smallworldness for all subjects in data_dir.

    Requires that the data_dir contains a subdirectory for each subject, and that each subject directory contains a "undirecred_thresholded_graphs" subdirectory with a file ending with small_worldness_sigma_suffix.

    Args:
        data_dir (str): Path to directory containing subject data.
        small_worldness_sigma_suffix (str, optional): Suffix of file containing pre-computed smallworldness. Defaults to '_sw_sigmas.mat'.
        control_folder_prefix (str, optional): Prefix of control subject folders. Defaults to 'amc'.
    Returns:
        pd.DataFrame: Dataframe containing all subjects and their smallworldness AUC (for everyone of their connectivity matrices).
    """
    # get list of subjects
    subjects = [subject for subject in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, subject))]

    # initialize output_df
    output_df = pd.DataFrame(columns=['subject', 'subject_type', 'subject_id', 'subject_timepoint', 'connectivity_file_name',
                                   'small_worldness_sigma_auc'])

    # loop over subjects
    for subject in tqdm(subjects):
        # check if subject has "undirected_thresholded_graphs" subdirectory
        if not os.path.isdir(os.path.join(data_dir, subject, 'undirected_thresholded_graphs')):
            if verbose:
                print(f'Skipping subject {subject} because no "undirected_thresholded_graphs" subdirectory was found.')
            continue

        if subject.startswith(control_folder_prefix):
            subject_type = 'control'
            subject_id = subject.split('_')[1]
            subject_timepoint = subject.split('_')[3]
        else:
            subject_type = 'patient'
            subject_id = subject.split('_')[1]
            subject_timepoint = subject.split('_')[2]

        # get list of small_worldness_sigma files
        small_worldness_sigma_files = [file for file in os.listdir(os.path.join(data_dir, subject, 'undirected_thresholded_graphs'))
                                       if (
                                           file.startswith(restrict_to_prefix) and
                                           file.endswith(small_worldness_sigma_suffix)
                                       )]

        for file in small_worldness_sigma_files:
            connectivity_matrix_name = file.split("udt_graphs_")[1].split(small_worldness_sigma_suffix)[0]
            small_worldness_sigmas_over_thresholds = np.rec.array(sio.loadmat(os.path.join(data_dir, subject, 'undirected_thresholded_graphs', file))['graph_sigmas'])
            sw_sigmas_dict = {}
            for string_threshold in small_worldness_sigmas_over_thresholds.dtype.names:
                threshold = float(string_threshold.split('t')[1]) / 10
                sw_sigmas_dict[threshold] = small_worldness_sigmas_over_thresholds[string_threshold][0][0][0][0]

            # transform into graph at multiple thresholds [0.1, 0.2, ..., 1]
            thresholds = np.arange(0.1, 1.1, 0.1)
            # only keep one decimal place
            thresholds = np.around(thresholds, decimals=1)
            small_worldness_sigma_auc = thresholded_auc(minimum_connectivity_threshold, thresholds, sw_sigmas_dict)

            output_df = pd.concat([output_df, pd.DataFrame({'subject': subject,
                                                            'subject_type': subject_type,
                                                            'subject_id': subject_id,
                                                            'subject_timepoint': subject_timepoint,
                                                            'connectivity_file_name': connectivity_matrix_name,
                                                            'small_worldness_sigma_auc': small_worldness_sigma_auc,
                                                            },
                                                           index=[0])], ignore_index=True)

    output_df['minimum_connectivity_threshold'] = minimum_connectivity_threshold

    return output_df



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze networks for all subjects in data_dir.')
    parser.add_argument('-i', '--input_data_dir', type=str, help='Path to directory containing subject data.', required=True)
    parser.add_argument('-t', '--minimum_connectivity_threshold', type=float, help='Minimum inclusive threshold to include for AUC computation (default: 0.3, i.e. [0.3-1]).', default=0.3)
    parser.add_argument('-p', '--restrict_to_prefix', type=str, help='Prefix of connectivity files to restrict to (default: udt_graphs_).', default='udt_graphs_')
    parser.add_argument('-s', '--small_worldness_sigma_suffix', type=str, help='Suffix of file containing pre-computed smallworldness (default: _sw_sigmas.mat).', default='_sw_sigmas.mat')
    parser.add_argument('-c', '--control_folder_prefix', type=str, help='Prefix of control folder (default: amc).', default='amc')
    parser.add_argument('-o', '--output_dir', type=str, help='Path to output directory.', default='')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print verbose output.', default=False)

    args = parser.parse_args()

    if args.output_dir == '':
        args.output_dir = args.input_data_dir

    output_df = compute_small_worldness_auc(args.input_data_dir,
                                            args.minimum_connectivity_threshold,
                                            args.restrict_to_prefix,
                                            args.small_worldness_sigma_suffix,
                                            control_folder_prefix=args.control_folder_prefix,
                                            verbose=args.verbose)

    output_df = output_df.sort_values(by=['subject_type', 'subject_id', 'subject_timepoint', 'connectivity_file_name'])

    output_df.to_csv(os.path.join(args.output_dir, 'small_worldness_auc.csv'), index=False)