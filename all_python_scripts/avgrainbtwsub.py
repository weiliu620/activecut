#!/usr/bin/env python 
import glob
import subprocess
import sys
import os
import fnmatch

def avgrain(refdir, targetdir):
    """
    Compute average rand index between a reference label map and each label map
    in a folder.
    
    usage: avggrain(refimage, targetdir)
    """

    avg = 0;
    numFiles = 0;
    allsubfiles = [ thisfile for thisfile in os.listdir(targetdir) if fnmatch.fnmatch(thisfile, 'sub*.nii.gz') ]
    allsubfiles.sort()
    for thisfile in allsubfiles:
        # rain = subprocess.check_output(['/home/sci/weiliu/projects/group_mrf_lemon/build_uv/compareclusterings', '-v', '0', '-r', os.path.join(refdir, thisfile), '-i', os.path.join(targetdir, thisfile)])
        # print('subject {0} Rand index: {1:.3f}'.format(thisfile, float(rain) ))


        rain = subprocess.Popen(['/home/sci/weiliu/projects/group_mrf_lemon/build_uv/compareclusterings', '-v', '0', '-r', os.path.join(refdir, thisfile), '-i', os.path.join(targetdir, thisfile)], stdout = subprocess.PIPE).communicate()[0]

        print('subject {0}: {1}'.format(thisfile, rain ))
        avg = avg + float(rain)
        numFiles = numFiles + 1

    avg = avg / numFiles;
    print ('{0} subjects. average Rand index: {1:.3f}'.format(numFiles, avg) )

if __name__ == '__main__':
    avgrain(sys.argv[1], sys.argv[2])
    
