import os
import os.path
import scipy.io as sio
import nibabel as nib
import numpy as np

from part2.preprocessing.mask_atlas_overlap_in_connectivity import remove_regions_from_original_connectivity_file
from part2.preprocessing.tools import remove_regions_from_connectivity


def filter_small_regions(data_directory:str, atlased_T1_name_end: str = '_AtlasBN_Atlas_246.nii',
                            masked_prefix:str = 'masked_',
                            connectivity_matrix_filename:str = 'transfer_exclude0_corrCMs_useAtlas1BN_Atlas_246_subband4_blockDriftCorrlinear_firstEV0_windsor1_95_regressOutMotion1_regressOutCSF1_regOutMat0_ROImeanScaling1_removeFirstNscans5_BootNone.mat',
                            connectivity_matrix_name_start: list = ['CM'], timecourse_matrix_name_start: list = ['TCS'],
                            save_prefix:str = 'filtered_', threshold: int = 2, verbose:bool=True):
    '''
    Remove regions with < [threshold] voxels for at least one subject in the dataset

    Prerequisite: connectivity matrices are already masked (presupposes structure of .mat after masking)

    Result: save connectivity matrix' with small regions removed

    Arguments:
        data_directory {str} -- path to directory containing all subjects' data
        atlased_T1_name_end {str} -- name of atlased T1 file (default: {'_AtlasBN_Atlas_246.nii'})
        masked_prefix {str} -- prefix of masked connectivity matrix file (default: {'masked_'})
        connectivity_matrix_name {str} -- name of masked connectivity matrix file (default: {'masked_transfer_exclude0_corrCMs_useAtlas1BN_Atlas_246_subband4_blockDriftCorrlinear_firstEV0_windsor1_95_regressOutMotion1_regressOutCSF1_regOutMat0_ROImeanScaling1_removeFirstNscans5_BootNone.mat'})
        connectivity_matrix_name_start {list} -- list of strings that start the name of the connectivity matrix in the matlab file (default: {['CM']})
        timecourse_matrix_name_start {list} -- list of strings that start the name of the timecourse matrix in the matlab file (default: {['TCS']})
        save_prefix {str} -- prefix to add to filtered connectivity matrix file (default: {'filtered_'})
        threshold {int} -- threshold (exclusive upper bound) for number of voxels in a region (default: {2})
        verbose {bool} -- print statements (default: {True})

    :return:
    '''

    masked_connectivity_matrix_name = masked_prefix + connectivity_matrix_filename
    connectivity_files_paths = []
    masked_connectivity_files_paths = []
    n_labels = 246

    # Find regions that are too small in at least one patient (i.e. have < [threshold] voxels)
    small_regions = []
    # - loop through main data directory
    for patient in os.listdir(data_directory):
        if not os.path.isdir(os.path.join(data_directory, patient)):
            continue
        # - loop through atlased T1s
        for file in os.listdir(os.path.join(data_directory, patient)):
            if file.endswith(atlased_T1_name_end):
                # - load atlased T1
                atlased_T1 = nib.load(os.path.join(data_directory, patient, file))
                # - loop through region labels
                for label in range(1, n_labels + 1):
                    # - find number of voxels in region
                    num_voxels = np.sum(atlased_T1.get_fdata() == label)
                    # - if number of voxels < threshold, add label to list of small regions
                    if num_voxels < threshold:
                        small_regions.append(label)

            if file == connectivity_matrix_filename:
                connectivity_files_paths.append(os.path.join(data_directory, patient, file))

            if file == masked_connectivity_matrix_name:
                masked_connectivity_files_paths.append(os.path.join(data_directory, patient, file))

    # retain only a collection of unique regions
    small_regions = list(set(small_regions))

    if verbose:
        print(f'List of small regions ({str(len(small_regions))}) using threshold of {threshold}:', small_regions)

    # Remove found small regions from connectivity matrices for all subjects
    for masked_connectivity_file in masked_connectivity_files_paths:
        remove_regions_from_masked_connectivity_file(small_regions, masked_connectivity_file, save_prefix=save_prefix,
                                                     connectivity_matrix_name_start=connectivity_matrix_name_start,
                                                     timecourse_matrix_name_start=timecourse_matrix_name_start
                                                     )

    for connectivity_file in connectivity_files_paths:
        remove_regions_from_original_connectivity_file(small_regions, connectivity_file, save_prefix=save_prefix,
                                                     connectivity_matrix_name_start=connectivity_matrix_name_start,
                                                     timecourse_matrix_name_start=timecourse_matrix_name_start
                                                     )

    return None


def remove_regions_from_masked_connectivity_file(regions_to_remove:list, connectivity_file_path:str,
                                          connectivity_matrix_name_start:list = ['CM'], timecourse_matrix_name_start:list=['TCS'],
                                          save_dir:str='', save_prefix:str='filtered_') -> None:
    """
    Method to remove regions from a connectivity matrix file (once it has been masked ie matfile structure modified) and save the filtered file
    :param regions_to_remove: list of regions to remove
    :param connectivity_file_path: path to connectivity matrix matlab (.mat) file
    :param connectivity_matrix_name_start: list of strings that start the name of the connectivity matrix in the matlab file
    :param timecourse_matrix_name_start: list of strings that start the name of the timecourse matrix in the matlab file
    :return: None
    """

    connectivity = sio.loadmat(connectivity_file_path)
    filtered_connectivity = connectivity.copy()
    initial_index_to_region_correspondence = np.squeeze(connectivity['allCodeBooks'])

    connectivity_matrix_keys = [key for key in connectivity.keys() if key.startswith(tuple(connectivity_matrix_name_start))]
    timecourse_matrix_keys = [key for key in connectivity.keys() if key.startswith(tuple(timecourse_matrix_name_start))]

    for connectivity_matrix_key in connectivity_matrix_keys:
        connectivity_matrix = connectivity[connectivity_matrix_key]
        filtered_connectivity_matrix, filtered_index_to_region_correspondence = remove_regions_from_connectivity(regions_to_remove,
                                                                                              connectivity_matrix,
                                                                                               list(initial_index_to_region_correspondence))
        filtered_connectivity[connectivity_matrix_key] = filtered_connectivity_matrix

    for timecourse_matrix_key in timecourse_matrix_keys:
        timecourse_matrix = connectivity[timecourse_matrix_key]
        filtered_timecourse_matrix, filtered_index_to_region_correspondence = remove_regions_from_connectivity(
            regions_to_remove,
            timecourse_matrix,
            list(initial_index_to_region_correspondence))
        filtered_connectivity[timecourse_matrix_key] = filtered_timecourse_matrix

    filtered_connectivity['allCodeBooks'] = filtered_index_to_region_correspondence

    if save_dir == '':
        save_dir = os.path.dirname(connectivity_file_path)

    sio.savemat(os.path.join(save_dir, save_prefix + os.path.basename(connectivity_file_path)), filtered_connectivity)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Filter out small regions from all connectivity matrix files')
    parser.add_argument('--data_directory', '-d', required=True, type=str, help='path to the mask file')
    parser.add_argument('--atlased_T1_name_end', '-a', required=False, type=str, help='end of atlas name to help find corresponding path to the atlas file',
                        default='_AtlasBN_Atlas_246.nii')
    parser.add_argument('--masked_file_prefix', '-m', required=False, type=str, help='prefix of masked connectivity matrix file', default='masked_')
    parser.add_argument('--connectivity_name', '-c', required=False, type=str, help='name of the masked connectivity file',
                        default='transfer_exclude0_corrCMs_useAtlas1BN_Atlas_246_subband4_blockDriftCorrlinear_firstEV0_windsor1_95_regressOutMotion1_regressOutCSF1_regOutMat0_ROImeanScaling1_removeFirstNscans5_BootNone.mat')
    parser.add_argument('--threshold', '-t', required=False, type=int,
                        help='exclusive upper bound of number of voxels to be considered a small region', default=None)
    parser.add_argument('--connectivity_name_start', '-cns', nargs="*", required=False, help='list of strings that start the name of the connectivity matrix in the matlab file', default=['CM'])
    parser.add_argument('--timecourse_name_start', '-tcs', nargs="*", required=False, help='list of strings that start the name of the timecourse matrix in the matlab file', default=['TCS'])
    parser.add_argument('--save_prefix', '-s', required=False, type=str, help='prefix for filtered connectivity file', default='filtered_')

    parser.add_argument('--verbose', '-v', required=False, action='store_true', help='Verbose mode')
    args = parser.parse_args()

    filter_small_regions(args.data_directory, atlased_T1_name_end=args.atlased_T1_name_end,
                         masked_prefix=args.masked_file_prefix, connectivity_matrix_filename=args.connectivity_name,
                         connectivity_matrix_name_start = args.connectivity_name_start,
                         timecourse_matrix_name_start = args.timecourse_name_start,
                         save_prefix=args.save_prefix, threshold=args.threshold, verbose=args.verbose
                         )