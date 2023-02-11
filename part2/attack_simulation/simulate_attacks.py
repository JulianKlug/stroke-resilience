import os
from typing import List

import pandas as pd

from part2.preprocessing.mask_atlas_overlap_in_connectivity import detect_masked_regions
from part2.preprocessing.tools import remove_regions_from_connectivity_file
from part2.utils.utils import ensure_dir


def simulate_attacks(connectivity_path:str, atlas_path:str, mask_paths:List[str],
                        connectivity_name_start:list = ['CM'], timecourse_name_start:list = ['TCS'],
                        n_voxels_overlap_threshold:int=None, fractional_overlap_threshold:float=None, store_intermediates:bool=False,
                        output_dir:str='', logging:bool=True) -> None:
    """
    Simulate attacks on the network of a single subject.
    Gist: Given a set of masks of clinical lesions, remove the lesioned nodes covered by each mask from the network to simulate the attack.

    Arguments:
        connectivity_path {str} -- Path to the connectivity file
        atlas_path {str} -- Path to the atlas file
        mask_paths {list[str]} -- List of paths to the mask files
        connectivity_name_start {list} -- List of strings that the connectivity object names start with
        timecourse_name_start {list} -- List of strings that the timecourse object names start with
        n_voxels_overlap_threshold {int} -- Number of voxels that need to overlap between the mask and the atlas to be considered a lesion
        fractional_overlap_threshold {float} -- Fraction of voxels that need to overlap between the mask and the atlas to be considered a lesion
        store_intermediates {bool} -- Whether to store intermediary steps or not
        output_dir {str} -- Path to the output directory
        logging {bool} -- Whether to log the attack parameters or not

    Returns:
        None

    """

    if output_dir == '':
        output_dir = os.path.join(os.path.dirname(connectivity_path), 'attacks')
        ensure_dir(output_dir)

    if logging:
        log_df = pd.DataFrame()

    for mask_path in mask_paths:
        attack_id = os.path.basename(mask_path).split('_')[1]
        attack_dir = os.path.join(output_dir, attack_id)
        ensure_dir(attack_dir)

        masked_regions = detect_masked_regions(atlas_path, mask_path, n_voxels_overlap_threshold=n_voxels_overlap_threshold,
                                                 fractional_overlap_threshold=fractional_overlap_threshold, save_intermediary_steps=store_intermediates,
                                                 save_dir=attack_dir)

        if logging:
            n_masked_regions = len(masked_regions)
            log_df = log_df.append(
                {
                    'attack_id': attack_id,
                    'n_masked_regions': n_masked_regions,
                    'masked_regions': masked_regions,
                    'attack_mask_path': mask_path,
                },
                ignore_index=True,
            )

        attack_prefix = f'attack_{attack_id}_'

        remove_regions_from_connectivity_file(masked_regions, connectivity_path, connectivity_name_start, timecourse_name_start, save_dir=attack_dir,
                                                save_prefix=attack_prefix)

    if logging:
        # log all parameters
        log_df['atlas_path'] = atlas_path
        log_df['connectivity_path'] = connectivity_path
        log_df['connectivity_name_start'] = f'[{",".join(connectivity_name_start)}]'
        log_df['timecourse_name_start'] = f'[{",".join(timecourse_name_start)}]'
        log_df['n_voxels_overlap_threshold'] = n_voxels_overlap_threshold
        log_df['fractional_overlap_threshold'] = fractional_overlap_threshold
        # add time & date
        log_df['timestamp'] = pd.to_datetime('today').strftime('%Y-%m-%d %H:%M:%S')

        log_df.to_csv(os.path.join(output_dir, 'attack_log.csv'))



