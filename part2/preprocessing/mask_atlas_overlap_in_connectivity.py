import os.path
import scipy.io as sio
import nibabel as nib
import numpy as np
from nilearn.image import resample_img

from part2.preprocessing.tools import remove_regions_from_connectivity


def detect_masked_regions(atlas_path:str, mask_path:str, n_voxels_overlap_threshold:int=None, fractional_overlap_threshold:float=None,
                          save_intermediary_steps:bool=False, save_dir:str='') -> list:
    """
    Method to detect masked regions in an atlas
    :param atlas_path: path to atlas
    :param mask_path: path to mask
    :param n_voxels_overlap_threshold: number of voxels in overlap to be considered a masked region
    :param fractional_overlap_threshold: fraction of voxels in overlap to be considered a masked region
    :param save_intermediary_steps: if true, intermediate resampled mask is saved
    :return: list of masked regions
    """

    # sanity check input arguments
    if n_voxels_overlap_threshold is None and fractional_overlap_threshold is None:
        raise ValueError('Must specify at least one overlap threshold')
    if n_voxels_overlap_threshold is not None and fractional_overlap_threshold is not None:
        raise ValueError('Must specify only one overlap threshold')

    atlas = nib.load(atlas_path)
    mask = nib.load(mask_path)
    atlas_data = np.squeeze(atlas.get_fdata())
    mask_data = np.squeeze(mask.get_fdata())

    if (save_intermediary_steps) & (save_dir == ''):
        save_dir = os.path.dirname(mask_path)

    # resample mask to atlas
    if mask_data.shape != atlas_data.shape:
        print('Resampling mask to atlas resolution')
        mask = resample_img(mask, target_affine=atlas.affine, target_shape=atlas.shape, interpolation='nearest')
        mask_data = np.squeeze(mask.get_fdata())

        if save_intermediary_steps:
            nib.save(mask, os.path.join(save_dir, 'resampled_' + os.path.basename(mask_path)))

    masked_regions = np.unique(atlas_data[mask_data > 0])

    masked_regions_over_threshold = []
    for region in masked_regions:
        if region == 0:
            continue

        # overlap threshold is specified in number of voxels overlapping with the mask
        if n_voxels_overlap_threshold is not None:
            if len(atlas_data[(mask_data > 0) & (atlas_data == region) & (atlas_data != 0)]) > n_voxels_overlap_threshold:
                masked_regions_over_threshold.append(region)

        # overlap threshold is specified as a fraction of voxels in the region overlapping with the mask
        if fractional_overlap_threshold is not None:
            fractional_region_mask_overlap = len(atlas_data[(mask_data > 0) & (atlas_data == region) & (atlas_data != 0)]) / len(atlas_data[(atlas_data == region) & (atlas_data != 0)])
            if fractional_region_mask_overlap > fractional_overlap_threshold:
                masked_regions_over_threshold.append(region)

    if save_intermediary_steps:
        # save a version of atlas only with masked regions
        masked_atlas_data = atlas_data.copy()
        masked_atlas_data[~np.isin(masked_atlas_data, masked_regions_over_threshold)] = 0
        masked_atlas = nib.Nifti1Image(masked_atlas_data, atlas.affine)
        nib.save(masked_atlas, os.path.join(save_dir, 'in_mask_' + os.path.basename(atlas_path)))


    return masked_regions_over_threshold


def remove_regions_from_original_connectivity_file(regions_to_remove:list, connectivity_file_path:str,
                                          connectivity_matrix_name_start:list = ['CM'], timecourse_matrix_name_start:list=['TCS'],
                                          save_dir:str='') -> None:
    """
    Method to remove regions from a connectivity matrix file and save the filtered file
    :param regions_to_remove: list of regions to remove
    :param connectivity_file_path: path to connectivity matrix matlab (.mat) file
    :param connectivity_matrix_name_start: list of strings that start the name of the connectivity matrix in the matlab file
    :param timecourse_matrix_name_start: list of strings that start the name of the timecourse matrix in the matlab file
    :return: None
    """

    connectivity = sio.loadmat(connectivity_file_path)
    filtered_connectivity = connectivity.copy()
    initial_index_to_region_correspondence = np.squeeze(connectivity['allCodeBooks'][0][0][0][0][0][0][0][0][0][0][0])

    connectivity_matrix_keys = [key for key in connectivity.keys() if key.startswith(tuple(connectivity_matrix_name_start))]
    timecourse_matrix_keys = [key for key in connectivity.keys() if key.startswith(tuple(timecourse_matrix_name_start))]

    for connectivity_matrix_key in connectivity_matrix_keys:
        connectivity_matrix = connectivity[connectivity_matrix_key][0][0][0]
        filtered_connectivity_matrix, filtered_index_to_region_correspondence = remove_regions_from_connectivity(regions_to_remove,
                                                                                              connectivity_matrix,
                                                                                               list(initial_index_to_region_correspondence))
        filtered_connectivity[connectivity_matrix_key] = filtered_connectivity_matrix

    for timecourse_matrix_key in timecourse_matrix_keys:
        timecourse_matrix = connectivity[timecourse_matrix_key][0][0][0][0][0]
        filtered_timecourse_matrix, filtered_index_to_region_correspondence = remove_regions_from_connectivity(
            regions_to_remove,
            timecourse_matrix,
            list(initial_index_to_region_correspondence))
        filtered_connectivity[timecourse_matrix_key] = filtered_timecourse_matrix

    filtered_connectivity['allCodeBooks'] = filtered_index_to_region_correspondence

    if save_dir == '':
        save_dir = os.path.dirname(connectivity_file_path)

    sio.savemat(os.path.join(save_dir, 'masked_' + os.path.basename(connectivity_file_path)), filtered_connectivity)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Mask regions from a connectivity matrix file')
    parser.add_argument('--mask_path', '-m', required=True, type=str, help='path to the mask file')
    parser.add_argument('--atlas_path', '-a', required=True, type=str, help='path to the atlas file')
    parser.add_argument('--connectivity_path', '-c', required=True, type=str, help='path to the atlas file')
    parser.add_argument('--n_voxels_overlap_threshold', '-nt', required=False, type=int,
                        help='number of voxels in overlap to be considered a masked region', default=None)
    parser.add_argument('--fractional_overlap_threshold', '-ft', required=False, type=float,
                        help='fraction of voxels in overlap to be considered a masked region', default=None)
    parser.add_argument('--store_intermediates', '-s', required=False, action='store_true', help='Store intermediate results (mask resampled to atlas space & only the masked regions in the atlas)')
    parser.add_argument('--output_dir', '-o', required=False, type=str, help='path to save the masked connectivity file (as well as intermediate steps)', default='')
    parser.add_argument('--connectivity_name_start', '-cns', nargs="*", required=False, help='list of strings that start the name of the connectivity matrix in the matlab file', default=['CM'])
    parser.add_argument('--timecourse_name_start', '-tcs', nargs="*", required=False, help='list of strings that start the name of the timecourse matrix in the matlab file', default=['TCS'])
    parser.add_argument('--verbose', '-v', required=False, action='store_true', help='Verbose mode')
    args = parser.parse_args()

    if args.output_dir == '':
        args.output_dir = os.path.dirname(args.connectivity_path)

    masked_regions = detect_masked_regions(args.atlas_path, args.mask_path, n_voxels_overlap_threshold=args.n_voxels_overlap_threshold,
                                             fractional_overlap_threshold=args.fractional_overlap_threshold, save_intermediary_steps=args.store_intermediates,
                                             save_dir=args.output_dir)

    if args.verbose:
        if args.n_voxels_overlap_threshold  is not None:
            print('Masked regions: {}'.format(masked_regions), 'at threshold: {}'.format(args.n_voxels_overlap_threshold))
        elif args.fractional_overlap_threshold is not None:
            print('Masked regions: {}'.format(masked_regions), 'at threshold: {}'.format(args.fractional_overlap_threshold))

    remove_regions_from_original_connectivity_file(masked_regions, args.connectivity_path, args.connectivity_name_start, args.timecourse_name_start, save_dir=args.output_dir)
