import sys
import os
import fnmatch
import shutil
import numpy
import subprocess


def alignlabels(mydir, ref):
    """
    if the label map in a folder is not aligned to the Yeo map, align it.
    """

    # init_grp = "/usr/sci/scratch/weiliu/NYU_test_retest/Yeo2011_7Networks_MNI152_FreeSurferConformed3mm.nii.gz"
    all_files = [f for f in os.listdir(mydir) if fnmatch.fnmatch(f, '*.nii.gz') ]
    all_files.sort()

    alignlabels_set = set()
    for thisfile in all_files:
        alignlabels_set.add( subprocess.Popen(['/home/sci/weiliu/projects/group_mrf/build_uv/alignlabels', '-t', os.path.join(mydir, thisfile), '-r', ref, '-o', os.path.join(mydir, thisfile), '-k', '7']) )
        
    for p in alignlabels_set:
        if p.poll is None:
            p.wait()
                
                

if __name__ == '__main__':
    alignlabels(sys.argv[1], sys.argv[2])
