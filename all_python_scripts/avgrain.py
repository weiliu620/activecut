#!/usr/bin/env python 
import subprocess
import sys
import fnmatch
import os

def avgrain(refimage, targetdir):
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
        # rain = subprocess.check_output(['/home/sci/weiliu/projects/group_mrf_lemon/build_uv/compareclusterings', '-v', '0', '-r', refimage, '-i', os.path.join(targetdir, thisfile) ])
        # print('subject {} Rand index: {:.3f}'.format(thisfile, float(rain) ))

        rain = subprocess.Popen(['/home/sci/weiliu/projects/group_mrf_lemon/build_uv/compareclusterings', '-v', '0', '-r', refimage, '-i', os.path.join(targetdir, thisfile) ], stdout = subprocess.PIPE).communicate()[0]
        print('subject {0}: {1}'.format(thisfile, rain ))

        avg = avg + float(rain)
        numFiles = numFiles + 1

    avg = avg / numFiles;
    print ('{0} subjects. average Rand index: {1}'.format(numFiles, avg) )

if __name__ == '__main__':
    avgrain(sys.argv[1], sys.argv[2])
    


