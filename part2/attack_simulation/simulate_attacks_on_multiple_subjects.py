import glob
import os
import warnings

from tqdm import tqdm

from part2.attack_simulation.simulate_attacks import simulate_attacks


def simulate_attacks_on_multiple_subjects(subjects_data_dir:str, atlas_path:str, mask_dir:str,
                        connectivity_name_start:list = ['CM'], timecourse_name_start:list = ['TCS'],
                        lesion_mask_name_end:str = 'lesion_ifhfix.nii.gz', connectivity_prefix:str = 'noNaN_transfer_',
                        n_voxels_overlap_threshold:int=None, fractional_overlap_threshold:float=None, store_intermediates:bool=False,
                        output_dir:str='', logging:bool=True, verbose:bool=False) -> None:
    """
    Simulate multiple attacks on multiple subjects.
    Gist: Given a set of masks of clinical lesions, remove the lesioned nodes covered by each mask from the network to simulate the attack. This is repeated for all subjects.

    Args:
        subjects_data_dir {str} -- Path to the directory containing the subjects' data
        atlas_path {str} -- Path to the atlas file
        mask_dir {str} -- Path to the directory containing the lesion masks
        output_dir {str} -- Path to the output directory, if not specified, the output will be stored in the same directory as the subjects' data, under a new directory called 'attacks'
        lesion_mask_name_end {str} -- String that the lesion mask file names end with, default is 'lesion_ifhfix.nii.gz'
        connectivity_prefix {str} -- Prefix of the connectivity object names, default is 'noNaN_transfer_'
        connectivity_name_start {list} -- List of strings that the connectivity object names start with, default is ['CM']
        timecourse_name_start {list} -- List of strings that the timecourse object names start with, default is ['TCS']
        n_voxels_overlap_threshold {int} -- Number of voxels that need to overlap between the mask and the atlas to be considered a lesion
        fractional_overlap_threshold {float} -- Fraction of voxels that need to overlap between the mask and the atlas to be considered a lesion
        store_intermediates {bool} -- Whether to store intermediary steps or not
        logging {bool} -- Whether to log the attack parameters or not
        verbose {bool} -- Whether to print the progress or not

    Returns:
        None
    """

    # find all maks in mask_dir
    mask_paths = glob.glob(os.path.join(mask_dir, f'*{lesion_mask_name_end}'))

    if verbose:
        print(f'Found {len(mask_paths)} masks in {mask_dir}')

    # identify all subject directories in subjects_data_dir
    subjects_directories = [subject_dir for subject_dir in os.listdir(subjects_data_dir)
                            if not subject_dir.startswith('.') and os.path.isdir(os.path.join(subjects_data_dir, subject_dir))]
    if verbose:
        print(f'Found {len(subjects_directories)} subjects in {subjects_data_dir}')
        print(f'Simulating {len(mask_paths)} attacks on {len(subjects_directories)} subjects...')

    # iterate through all subjects directories
    for subject_dir in tqdm(subjects_directories):
        subject_path = os.path.join(subjects_data_dir, subject_dir)
        subject_output_dir = (
            os.path.join(output_dir, subject_dir) if output_dir else ''
        )
        # find the correct connectivity file in subject_path given the connectivity_prefix
        connectivity_path_possibilities = glob.glob(os.path.join(subject_path, f'{connectivity_prefix}*.mat'))

        if len(connectivity_path_possibilities) == 0:
            warnings.warn(f'No connectivity file found in {subject_path} with prefix {connectivity_prefix}. Skipping subject.')
            continue

        connectivity_path = connectivity_path_possibilities[0]

        # simulate attacks on subject
        simulate_attacks(connectivity_path, atlas_path, mask_paths,
                        connectivity_name_start, timecourse_name_start,
                        n_voxels_overlap_threshold, fractional_overlap_threshold, store_intermediates,
                        output_dir=subject_output_dir, logging=logging)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Simulate attacks on multiple subjects by removing the lesioned nodes from the network')
    parser.add_argument('--input_data_dir', '-i', required=True, type=str, help='path to the directory containing the subjects data')
    parser.add_argument('--atlas_path', '-a', required=True, type=str, help='path to the atlas file')
    parser.add_argument('--mask_dir', '-m', required=True, type=str, help='path to the directory containing the lesion masks')
    parser.add_argument('--n_voxels_overlap_threshold', '-nt', required=False, type=int,
                        help='number of voxels in overlap to be considered a masked region', default=None)
    parser.add_argument('--fractional_overlap_threshold', '-ft', required=False, type=float,
                        help='fraction of voxels in overlap to be considered a masked region', default=None)
    parser.add_argument('--store_intermediates', '-s', required=False, action='store_true', help='Store intermediate results (mask resampled to atlas space & only the masked regions in the atlas)')
    parser.add_argument('--output_dir', '-o', required=False, type=str, help='path to save the attacked connectivity files (as well as intermediate steps). Default is a subfolder of the subject folder.', default='')
    parser.add_argument('--connectivity_name_start', '-cns', nargs="*", required=False, help='list of strings that start the name of the connectivity matrix in the matlab file', default=['CM'])
    parser.add_argument('--timecourse_name_start', '-tcs', nargs="*", required=False, help='list of strings that start the name of the timecourse matrix in the matlab file', default=['TCS'])
    parser.add_argument('--lesion_mask_name_end', '-lme', required=False, type=str, help='string that the lesion mask file names end with', default='lesion_ifhfix.nii.gz')
    parser.add_argument('--connectivity_file_prefix', '-p', required=False, type=str, help='prefix of the connectivity file', default='noNaN_transfer_')
    parser.add_argument('--verbose', '-v', required=False, action='store_true', help='Verbose mode')
    parser.add_argument('--disable_logging', '-dl', required=False, action='store_true', help='Disable logging')

    args = parser.parse_args()

    simulate_attacks_on_multiple_subjects(
        subjects_data_dir=args.input_data_dir,
        atlas_path=args.atlas_path,
        mask_dir=args.mask_dir,
        connectivity_name_start=args.connectivity_name_start,
        lesion_mask_name_end=args.lesion_mask_name_end,
        connectivity_prefix=args.connectivity_file_prefix,
        n_voxels_overlap_threshold=args.n_voxels_overlap_threshold,
        fractional_overlap_threshold=args.fractional_overlap_threshold,
        store_intermediates=args.store_intermediates,
        output_dir=args.output_dir,
        logging=not args.disable_logging,
        verbose=args.verbose
    )
