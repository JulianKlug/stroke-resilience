import random
import numpy as np
import os
import os.path
import scipy.io as sio

def remove_regions_from_connectivity(regions_to_remove:list, connectivity_matrix:np.array, index_to_region_correspondence:list) -> np.array:
    """
    Method to remove regions from a connectivity matrix
    :param regions_to_remove: list of regions to remove
    :param connectivity_matrix: connectivity matrix
    :param index_to_region_correspondence: list of region names corresponding to connectivity matrix indices
    :return: connectivity matrix with regions removed
    """
    regions_to_remove_indices = [index_to_region_correspondence.index(region) for region in regions_to_remove
                                    if region in index_to_region_correspondence]
    filtered_connectivity_matrix = connectivity_matrix.copy()
    filtered_connectivity_matrix = np.delete(filtered_connectivity_matrix, regions_to_remove_indices, axis=0)
    filtered_connectivity_matrix = np.delete(filtered_connectivity_matrix, regions_to_remove_indices, axis=1)
    index_to_region_correspondence = [region for region in index_to_region_correspondence if region not in regions_to_remove]

    return filtered_connectivity_matrix, index_to_region_correspondence


def test_remove_regions_from_connectivity():
    """
    Test the function: remove_regions_from_original_connectivity_file
    """
    index_to_region_correspondence = [1,2,3,4,5]
    random.shuffle(index_to_region_correspondence)
    regions_to_remove_list = random.sample(index_to_region_correspondence, random.randrange(1, len(index_to_region_correspondence) -1))
    test_connectivity_matrix = np.random.random((5, 5))
    filtered_test_connectivity_matrix, filtered_index_to_region_correspondence = remove_regions_from_connectivity(regions_to_remove_list, test_connectivity_matrix, index_to_region_correspondence)

    assert filtered_test_connectivity_matrix.shape == (5 - len(regions_to_remove_list), 5 - len(regions_to_remove_list))
    #  assert that test_connectivity_matrix is equal to filtered_test_connectivity_matrix except for the rows and columns corresponding to the regions to remove
    index_to_remove_list = [index_to_region_correspondence.index(region) for region in regions_to_remove_list]
    remaining_indices = [index for index in range(len(index_to_region_correspondence)) if index not in index_to_remove_list]
    assert np.all(test_connectivity_matrix[remaining_indices, :][:, remaining_indices] == filtered_test_connectivity_matrix)


def remove_regions_from_connectivity_file(regions_to_remove:list, connectivity_file_path:str,
                                          connectivity_matrix_name_start:list = ['CM'], timecourse_matrix_name_start:list=['TCS'],
                                          save_dir:str='', save_prefix:str='masked_') -> None:
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

    initial_connectivity_matrix_state = True
    # check if the connectivity matrix is in its original state or if it has been reorganized
    try:
        initial_index_to_region_correspondence = np.squeeze(connectivity['allCodeBooks'][0][0][0][0][0][0][0][0][0][0][0])
    except IndexError:
        initial_connectivity_matrix_state = False
        initial_index_to_region_correspondence = np.squeeze(connectivity['allCodeBooks'])

    connectivity_matrix_keys = [key for key in connectivity.keys() if key.startswith(tuple(connectivity_matrix_name_start))]
    timecourse_matrix_keys = [key for key in connectivity.keys() if key.startswith(tuple(timecourse_matrix_name_start))]

    for connectivity_matrix_key in connectivity_matrix_keys:
        if initial_connectivity_matrix_state:
            connectivity_matrix = connectivity[connectivity_matrix_key][0][0][0]
        else:
            connectivity_matrix = connectivity[connectivity_matrix_key]

        filtered_connectivity_matrix, filtered_index_to_region_correspondence = remove_regions_from_connectivity(regions_to_remove,
                                                                                              connectivity_matrix,
                                                                                               list(initial_index_to_region_correspondence))
        filtered_connectivity[connectivity_matrix_key] = filtered_connectivity_matrix

    for timecourse_matrix_key in timecourse_matrix_keys:
        if initial_connectivity_matrix_state:
            timecourse_matrix = connectivity[timecourse_matrix_key][0][0][0][0][0]
        else:
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
