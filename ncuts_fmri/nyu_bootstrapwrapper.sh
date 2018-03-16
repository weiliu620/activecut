#!/usr/bin/env python 
import sys
from random import randint
def bootstrapsample(filedir, num_files):
    """
    Obtain one bootstrap sample from a number of fmri subject files.

    usage: bootstrapsample(filedir, num_files)
    """

        allsubfiles = [ thisfile for thisfile in os.listdir(filedir) if fnmatch.fnmatch(thisfile, '*.nii.gz') ]
        num_all_files = len(allsubfiles)
        allsubfiles.sort()
        for file_id in range(0, num_files):
            

if __name__ == '__main__':
    avgrain(sys.argv[1], sys.argv[2])
