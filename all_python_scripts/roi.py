import sys
import os
import fnmatch
import shutil
import numpy as np
import subprocess
import nibabel as nib

def draw_seeds(seeds_file, out_file, mask_file):
    """
    read seeds_file (Nx3) and draw on a image
    """
    mask = nib.load(mask_file)
    seeds = np.loadtxt(seeds_file)
    
    T = mask.get_affine()
    S = np.linalg.inv(T[0:3, 0:3])
    offset = T[0:3, 3]
    radius = 3 * np.absolute(S[0,0])
    radius = radius.astype(int)

    (n_seeds, dummy) = seeds.shape
    D = np.zeros(mask.shape)

    for seed_id in range(n_seeds):
        vox_id = np.dot(S, seeds[seed_id,:] - offset)
        vox_id = tuple(np.around(vox_id).astype(int))
        cand = [(i, j, k) for i in range(vox_id[0] - radius, vox_id[0] + radius + 1) for j in range(vox_id[1] - radius, vox_id[1] + radius + 1) for k in range(vox_id[2] - radius, vox_id[2] + radius + 1) if (i-vox_id[0])**2 + (j-vox_id[1])**2 + (k-vox_id[2])**2 <= radius**2]
        for v in cand:
            D[v] = 1
            
    outimage = nib.Nifti1Image(D, header = mask.get_header(), affine = mask.get_affine())
    outimage.set_data_dtype('int16')
    outimage.to_filename(out_file)
    print '{0} saved. '.format(out_file)            
