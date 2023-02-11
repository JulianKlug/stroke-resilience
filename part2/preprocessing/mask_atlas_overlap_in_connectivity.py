import os.path
import scipy.io as sio
import nibabel as nib
import numpy as np
from nilearn.image import resample_img

from part2.preprocessing.tools import remove_regions_from_connectivity, remove_regions_from_connectivity_file


def detect_masked_regions(atlas_path:str, mask_path:str, n_voxels_overlap_threshold:int=None, fractional_overlap_threshold:float=None,
                          save_intermediary_steps:bool=False, save_dir:str='', verbose:bool=False) -> list:
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
        if verbose:
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

    remove_regions_from_connectivity_file(masked_regions, args.connectivity_path, args.connectivity_name_start, args.timecourse_name_start, save_dir=args.output_dir)
