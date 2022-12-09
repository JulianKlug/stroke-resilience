import os
import numpy as np


def framewise_displacement(head_motion):
    """
    Method to calculate framewise Displacement (FD) as per Power et al., 2012

    Corresponding Matlab code (implementation in BRANT):
        motion_diff = diff(head_motion);
        FD = sum([abs(motion_diff(:, 1: 3)), 50 * abs(motion_diff(:, 4: 6))], 2);

    From Power et al., 2012:
    FDi = |Δdix| + |Δdiy| + |Δdiz| + |Δαi| + |Δβi| + |Δγi|, where Δdix = d(i − 1)x − dix, and similarly for the other rigid body parameters [dix diy diz αi βi γi]. Rotational displacements were converted from degrees to millimeters by calculating displace- ment on the surface of a sphere of radius 50 mm, which is approxi- mately the mean distance from the cerebral cortex to the center of the head.

    Parameters
    ----------
    head_motion : array
        movement parameters vector
    Returns
    -------
    FD : array
        Frame-wise displacement
    """

    motion_diff = np.diff(head_motion, axis=0, prepend=0)
    FD = np.sum(np.abs(motion_diff[:, 0:3]) + 50 * np.abs(motion_diff[:, 3:]), axis=1)

    return FD


def framewise_displacement_from_file(in_file:str, out_dir:str = ''):
    """
    Method to calculate framewise Displacement (FD) as per Power et al., 2012
    Parameters
    ----------
    in_file : string
        movement parameters vector file path (e.g. rp_***.txt from SPM)
    out_dir : string
        output directory

    Returns
    -------
    out_file : string
        Frame-wise displacement mat
        file path
    """

    head_motion = np.loadtxt(in_file)
    FD = framewise_displacement(head_motion)

    out_file_name = os.path.basename(in_file).replace('rp_', 'FD_')

    if out_dir == '':
        out_dir = os.path.dirname(in_file)

    np.savetxt(os.path.join(out_dir, out_file_name), FD)

    return os.path.join(out_dir, out_file_name)


def detect_excessive_framewise_displacement(in_file:str, displacement_threshold:float = 0.5,
                                            maximum_fraction_of_displaced_frames:float=0.6, verbose:bool=True) -> bool:
    """
    Method to detect excessive framewise displacement (FD) as per Power et al., 2012
    Parameters
    ----------
    in_file : string
        movement parameters vector file path (e.g. rp_***.txt from SPM)
    displacement_threshold : float
        displacement threshold in mm
    maximum_fraction_of_displaced_frames : float
        fraction of frames with displacement above threshold

    Returns
    -------
    bool:
        True if excessive displacement
    """

    head_motion = np.loadtxt(in_file)
    FD = framewise_displacement(head_motion)

    n_displaced_frames = np.sum(FD > displacement_threshold)
    fraction_displaced_frames = n_displaced_frames / len(FD)

    if fraction_displaced_frames > maximum_fraction_of_displaced_frames:
        if verbose:
            print(f'Excessive framewise displacement detected, with {n_displaced_frames} ({fraction_displaced_frames*100:.1f}%)'
                  f' of frames with displacement above {displacement_threshold}mm')
        return True
    else:
        if verbose:
            print(f'No excessive framewise displacement detected, with {n_displaced_frames} ({fraction_displaced_frames*100:.1f}%)'
                  f' of frames with displacement above {displacement_threshold}mm')
        return False


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Calculate framewise displacement')
    parser.add_argument('--compute_fd', '-fd', required=False, action='store_true', help='Compute framewise displacement')
    parser.add_argument('--in_file', '-i', required=True, type=str, help='input file')
    parser.add_argument('--out_dir', '-o', required=False, type=str, help='output directory', default='')
    parser.add_argument('--displacement_threshold', '-t', required=False, type=float, help='displacement threshold in mm', default=0.5)
    parser.add_argument('--maximum_fraction_of_displaced_frames', '-m', required=False, type=float, help='fraction of frames with displacement above threshold', default=0.6)
    parser.add_argument('--verbose', '-v', required=False, action='store_true', help='verbose')
    parser.set_defaults(compute_fd=False, verbose=True)
    args = parser.parse_args()

    detect_excessive_framewise_displacement(args.in_file, args.displacement_threshold, args.maximum_fraction_of_displaced_frames, args.verbose)

    if args.compute_fd:
        framewise_displacement_from_file(args.in_file, args.out_dir)