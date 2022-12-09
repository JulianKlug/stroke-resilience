import random
import numpy as np


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
