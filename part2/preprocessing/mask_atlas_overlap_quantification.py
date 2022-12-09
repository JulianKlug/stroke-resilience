import os.path
import pandas as pd
import nibabel as nib
import numpy as np
from nilearn.image import resample_img


def mask_atlas_overlap_quantification(atlas_path:str, mask_path:str,
                          save_intermediary_steps:bool=False, save_dir:str='', verbose:bool=False) -> pd.DataFrame:
    """
    Method to quantify masked regions in an atlas (number of voxels and fraction of total region volume)
    :param atlas_path: path to atlas
    :param mask_path: path to mask
    :param save_intermediary_steps: if true, intermediate resampled mask is saved
    :return: list of masked regions
    """

    atlas = nib.load(atlas_path)
    mask = nib.load(mask_path)
    atlas_data = np.squeeze(atlas.get_fdata())
    mask_data = np.squeeze(mask.get_fdata())

    if save_dir == '':
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

    masked_regions_quantification = []
    for region in masked_regions:
        n_voxels_overlap = len(atlas_data[(mask_data > 0) & (atlas_data == region)])
        fractional_region_mask_overlap = len(atlas_data[(mask_data > 0) & (atlas_data == region)]) / len(atlas_data[(atlas_data == region)])
        masked_regions_quantification.append([region, n_voxels_overlap, fractional_region_mask_overlap])

        if verbose:
            print(f'Region {region:.0f} has {n_voxels_overlap} voxels ({fractional_region_mask_overlap*100:.1f}%) overlap with the mask')

    masked_regions_quantification_df = pd.DataFrame(masked_regions_quantification, columns=['region', 'n_voxels_overlap', 'fractional_region_mask_overlap'])

    masked_regions_quantification_df.to_csv(os.path.join(save_dir, 'masked_regions_quantification.csv'))

    return masked_regions_quantification_df


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Quantify the regional overlap between an atlas and a mask')
    parser.add_argument('--mask_path', '-m', required=True, type=str, help='path to the mask file')
    parser.add_argument('--atlas_path', '-a', required=True, type=str, help='path to the atlas file')
    parser.add_argument('--store_intermediates', '-s', required=False, action='store_true',
                        help='Store intermediate results (mask resampled to atlas space)')
    parser.add_argument('--output_dir', '-o', required=False, type=str,
                        help='path to save the masked connectivity file (as well as intermediate steps)', default='')
    parser.add_argument('--verbose', '-v', required=False, action='store_true', help='Verbose mode')
    args = parser.parse_args()

    mask_atlas_overlap_quantification(args.atlas_path, args.mask_path, args.store_intermediates, args.output_dir, args.verbose)





