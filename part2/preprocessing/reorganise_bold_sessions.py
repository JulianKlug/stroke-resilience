import os
import shutil
import argparse

def reorganise_bold_sessions(subject_directory:str, n_first_images_to_remove:int=5, verbose:bool=False):
    """
    Reorganise all BOLD sessions from a subject to a single directory.
    The first n_first_images_to_remove images are removed from each session to account for stabilization of field fluctuations.
    :param subject_directory: path to the subject directory
    :param n_first_images_to_remove: number of first images to remove from each session (default: 5)
    :param verbose: if True, print information about the process
    :return: bool: True if successful, fails otherwise

    CLI Example:
    >> python reorganise_bold_sessions.py -s /path/to/subject_directory -v
    """

    new_directory_name = 'BOLD_all_sessions'
    bold_session_prefix = 'BOLD_'

    # get list of all subdirectories starting with Bold
    bold_subdirectories = [x[0] for x in os.walk(subject_directory)
                           if x[0].startswith(os.path.join(subject_directory, bold_session_prefix))
                           and not x[0].endswith(new_directory_name)]
    # ensure subdirectories are sorted by name (as name contains the session number)
    bold_subdirectories.sort()

    if not os.path.exists(os.path.join(subject_directory, new_directory_name)):
        os.mkdir(os.path.join(subject_directory, new_directory_name))

    for session_index, bold_subdirectory in enumerate(bold_subdirectories):
        if verbose:
            print(f'Original session index {bold_subdirectory.split("_")[-1]} -> New session index {session_index + 1}')
        # get list of all files in bold_subdirectory
        bold_files = [x for x in os.listdir(bold_subdirectory) if x.endswith('.nii')]
        # ensure files are sorted according to file name (which is named by time and image number)
        bold_files.sort()
        # remove first n images (to account for stabilization of fluctuations)
        bold_files = bold_files[n_first_images_to_remove:]
        # copy all files to BOLD_all_sessions directory
        for bold_file in bold_files:
            shutil.copy2(os.path.join(bold_subdirectory, bold_file),
                        os.path.join(subject_directory, new_directory_name))
            # rename new file to contain session number and respect new order of files
            new_file_index = str(len(os.listdir(os.path.join(subject_directory, new_directory_name)))).zfill(5)
            new_file_name = f'{bold_file.split("-")[0]}SESS{session_index + 1}-{bold_file.split("-")[1]}-{new_file_index}-{new_file_index}.nii'
            os.rename(os.path.join(subject_directory, new_directory_name, bold_file),
                        os.path.join(subject_directory, new_directory_name, new_file_name))


    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reorganise BOLD sessions from a subject to a single directory.')
    parser.add_argument('-s', '--subject_directory', type=str, required=True, help='path to the subject directory')
    parser.add_argument('-n', '--n_first_images_to_remove', type=int, default=5, help='number of first images to remove from each session')
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose output')
    args = parser.parse_args()

    reorganise_bold_sessions(args.subject_directory, args.n_first_images_to_remove, args.verbose)
    print('Reorganisation successful.')
    exit(0)