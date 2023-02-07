import os
import scipy.io as sio
import numpy as np

from part2.preprocessing.tools import remove_regions_from_connectivity_file


def filter_regions_with_nan(data_directory:str,
                            masked_prefix:str = 'masked_', filtered_prefix:str = 'filtered_',
                            connectivity_matrix_filename:str = 'transfer_exclude0_corrCMs_useAtlas1BN_Atlas_246_subband4_blockDriftCorrlinear_firstEV0_windsor1_95_regressOutMotion1_regressOutCSF1_regOutMat0_ROImeanScaling1_removeFirstNscans5_BootNone.mat',
                            connectivity_matrix_name_start: list = ['CM'], timecourse_matrix_name_start: list = ['TCS'],
                            save_prefix:str = 'noNaN_', verbose:bool = True):
    '''
    Remove regions with NaN values in connectivity matrix for at least one subject in the dataset

    Prerequisite: connectivity matrices are already filtered (presupposes structure of .mat after filtering)

    Result: save connectivity matrix' with NaN regions removed

    Arguments:
        data_directory {str} -- path to directory containing all subjects' data
        masked_prefix {str} -- prefix of masked connectivity matrix file (default: {'masked_'})
        filtered_prefix {str} -- prefix of filtered connectivity matrix file (default: {'filtered_'})
        connectivity_matrix_name {str} -- name of masked connectivity matrix file (default: {'masked_transfer_exclude0_corrCMs_useAtlas1BN_Atlas_246_subband4_blockDriftCorrlinear_firstEV0_windsor1_95_regressOutMotion1_regressOutCSF1_regOutMat0_ROImeanScaling1_removeFirstNscans5_BootNone.mat'})
        connectivity_matrix_name_start {list} -- list of strings that start the name of the connectivity matrix in the matlab file (default: {['CM']})
        timecourse_matrix_name_start {list} -- list of strings that start the name of the timecourse matrix in the matlab file (default: {['TCS']})
        save_prefix {str} -- prefix to add to filtered connectivity matrix file (default: {'filtered_'})
        verbose {bool} -- print statements (default: {True})

    :return:
    '''

    filtered_connectivity_matrix_name = filtered_prefix + masked_prefix + connectivity_matrix_filename
    filtered_connectivity_files_paths = []
    unmasked_connectivity_files_paths = []

    # Find regions that are have NaN values in at least one patient
    masked_regions_with_nan = []
    umasked_regions_with_nan = []
    # - loop through main data directory
    for patient in os.listdir(data_directory):
        if not os.path.isdir(os.path.join(data_directory, patient)):
            continue
        # - loop through filtered connectivity matrices
        for file in os.listdir(os.path.join(data_directory, patient)):
            if file == filtered_connectivity_matrix_name:
                # - load filtered connectivity matrix
                filtered_connectivity_file = sio.loadmat(os.path.join(data_directory, patient, file))
                filtered_CM = filtered_connectivity_file['CM3D']
                codeBook = filtered_connectivity_file['allCodeBooks']
                # - loop through region indexes
                for region_index in range(filtered_CM.shape[0]):
                    # - if region has only NaN values, find label and add it to list
                    if np.isnan(filtered_CM[region_index, :]).all():
                        masked_regions_with_nan.append(codeBook[:, region_index][0])
                        if verbose:
                            print(f'Region {codeBook[:, region_index][0]} has NaN values in {patient}')

                filtered_connectivity_files_paths.append(os.path.join(data_directory, patient, file))

            if file == connectivity_matrix_filename:
                # load connectivity matrix
                unmasked_connectivity_file = sio.loadmat(os.path.join(data_directory, patient, file))
                unmasked_CM = unmasked_connectivity_file['CM3D'][0][0][0]
                unmasked_codeBook = np.squeeze(unmasked_connectivity_file['allCodeBooks'][0][0][0][0][0][0][0][0][0][0][0])
                # - loop through region indexes
                for region_index in range(unmasked_CM.shape[0]):
                    # - if region has only NaN values, find label and add it to list
                    if np.isnan(unmasked_CM[region_index, :]).all():
                        umasked_regions_with_nan.append(unmasked_codeBook[region_index])
                        if verbose:
                            print(f'Unmasked data: Region {unmasked_codeBook[region_index]} has NaN values in {patient}')

                unmasked_connectivity_files_paths.append(os.path.join(data_directory, patient, file))

    # retain only a collection of unique regions
    regions_with_nan = list(set(masked_regions_with_nan))
    unmasked_regions_with_nan = list(set(umasked_regions_with_nan))

    if verbose:
        print(f'List of regions with NaN in masked data: ({str(len(regions_with_nan))}):', regions_with_nan)
        print(f'List of regions with NaN in unmasked data: ({str(len(unmasked_regions_with_nan))}):', unmasked_regions_with_nan)

    # Remove found small regions from connectivity matrices for all subjects
    for filtered_connectivity_file in filtered_connectivity_files_paths:
        remove_regions_from_connectivity_file(regions_with_nan, filtered_connectivity_file, save_prefix=save_prefix,
                                                     connectivity_matrix_name_start=connectivity_matrix_name_start,
                                                     timecourse_matrix_name_start=timecourse_matrix_name_start
                                                     )

    for unmasked_connectivity_file in unmasked_connectivity_files_paths:
        remove_regions_from_connectivity_file(unmasked_regions_with_nan, unmasked_connectivity_file, save_prefix=save_prefix,
                                                       connectivity_matrix_name_start=connectivity_matrix_name_start,
                                                       timecourse_matrix_name_start=timecourse_matrix_name_start
                                                       )

    return None


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Filter out regions with NaN values from all connectivity matrix files')
    parser.add_argument('--data_directory', '-d', required=True, type=str, help='path to the mask file')
    parser.add_argument('--masked_file_prefix', '-m', required=False, type=str, help='prefix of masked connectivity matrix file', default='masked_')
    parser.add_argument('--filtered_file_prefix', '-f', required=False, type=str, help='prefix of filtered connectivity matrix file', default='filtered_')
    parser.add_argument('--connectivity_name', '-c', required=False, type=str, help='name of the masked connectivity file',
                        default='transfer_exclude0_corrCMs_useAtlas1BN_Atlas_246_subband4_blockDriftCorrlinear_firstEV0_windsor1_95_regressOutMotion1_regressOutCSF1_regOutMat0_ROImeanScaling1_removeFirstNscans5_BootNone.mat')
    parser.add_argument('--connectivity_name_start', '-cns', nargs="*", required=False, help='list of strings that start the name of the connectivity matrix in the matlab file', default=['CM'])
    parser.add_argument('--timecourse_name_start', '-tcs', nargs="*", required=False, help='list of strings that start the name of the timecourse matrix in the matlab file', default=['TCS'])
    parser.add_argument('--save_prefix', '-s', required=False, type=str, help='prefix for filtered connectivity file', default='noNaN_')
    parser.add_argument('--verbose', '-v', required=False, action='store_true', help='Verbose mode')
    args = parser.parse_args()

    filter_regions_with_nan(args.data_directory,
                            masked_prefix=args.masked_file_prefix, filtered_prefix=args.filtered_file_prefix,
                            connectivity_matrix_filename=args.connectivity_name,
                            connectivity_matrix_name_start=args.connectivity_name_start,
                            timecourse_matrix_name_start=args.timecourse_name_start,
                            save_prefix=args.save_prefix,
                            verbose=args.verbose)
