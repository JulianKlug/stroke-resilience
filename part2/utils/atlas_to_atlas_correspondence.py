import os.path
import pandas as pd
import nibabel as nib
import numpy as np
from scipy import stats
from nilearn.image import resample_img



def atlas_to_atlas_mapping(atlas1_path, atlas2_path, verbose=False, save_intermediary_steps=False):
    """
    Generate a correspondence (mapping) matrix between two atlases

    :param atlas1_path:
    :param atlas2_path:
    :return:
    """

    atlas1 = nib.load(atlas1_path)
    atlas2 = nib.load(atlas2_path)
    atlas1_data = np.squeeze(atlas1.get_fdata())
    atlas2_data = np.squeeze(atlas2.get_fdata())

    # resample atlas2 to atlas1
    if atlas2_data.shape != atlas1_data.shape:
        if verbose:
            print('Resampling atlas2 to atlas1 resolution')
        atlas2 = resample_img(atlas2, target_affine=atlas1.affine, target_shape=atlas1.shape, interpolation='nearest')
        atlas2_data = np.squeeze(atlas2.get_fdata())

        if save_intermediary_steps:
            nib.save(atlas2, os.path.join(os.path.dirname(atlas2_path), 'resampled_' + os.path.basename(atlas2_path)))

    # get unique regions in atlas2
    atlas1_regions = np.unique(atlas1_data[atlas1_data > 0])

    # for every region in atlas2, find the region in atlas1 with the highest overlap
    correspondence_df = pd.DataFrame(columns=['atlas1_region', 'atlas2_region', 'n_voxels_overlap', 'fractional_region_overlap'])
    for region in atlas1_regions:
        n_voxels_overlap = len(atlas2_data[(atlas2_data > 0) & (atlas1_data == region)])
        fractional_region_overlap = n_voxels_overlap / len(atlas1_data[(atlas1_data == region)])

        if n_voxels_overlap == 0:
            corresponding_atlas2_region = np.nan
        else:
            corresponding_atlas2_region = int(stats.mode(atlas2_data[(atlas2_data > 0) & (atlas1_data == region)])[0][0])

        correspondence_df = correspondence_df.append({'atlas1_region': region,
                                                      'atlas2_region': corresponding_atlas2_region,
                                                      'n_voxels_overlap': n_voxels_overlap,
                                                      'fractional_region_overlap': fractional_region_overlap},
                                                     ignore_index=True)


    atlas1_name = os.path.basename(atlas1_path).split('.')[0]
    atlas2_name = os.path.basename(atlas2_path).split('.')[0]
    correspondence_df['atlas1_name'] = atlas1_name
    correspondence_df['atlas2_name'] = atlas2_name
    return correspondence_df


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Generate mapping file between two atlases')
    parser.add_argument('--atlas1_path', '-a1', required=True, type=str, help='path to the 1st atlas file')
    parser.add_argument('--atlas2_path', '-a2', required=True, type=str, help='path to the 2nd atlas file')
    parser.add_argument('--store_intermediates', '-s', required=False, action='store_true',
                        help='Store intermediate results (mask resampled to atlas space)')
    parser.add_argument('--output_dir', '-o', required=False, type=str,
                        help='path to save the masked connectivity file', default='')
    parser.add_argument('--verbose', '-v', required=False, action='store_true', help='Verbose mode')
    args = parser.parse_args()

    if args.output_dir == '':
        args.output_dir = os.path.dirname(args.atlas1_path)

    correspondence_df = atlas_to_atlas_mapping(args.atlas1_path, args.atlas2_path, verbose=args.verbose,
                                                  save_intermediary_steps=args.store_intermediates)
    correspondence_df.to_csv(os.path.join(args.output_dir, 'atlas_correspondence.csv'), index=False)

