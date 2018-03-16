import sys
import os
import fnmatch
import subprocess
import shutil
import socket
import math
import time
import random

def collect_RI(ref_dir, image_dir, out_file):
    """
    compute Rand index btw each file in the 2 folder.
    """
    all_files = [f for f in os.listdir(ref_dir) if fnmatch.fnmatch(f, 'sub*.nii.gz')]
    all_files.sort()
    print 'total file number: {0}'.format(len(all_files))
    hostname = socket.gethostname()
    proc_set = set()
    for sub_file in all_files:
        proc_set.add(subprocess.Popen(['/home/sci/weiliu/projects/group_mrf_lemon/build_' + hostname + '/compareclusterings', '-v', '0', '-r', os.path.join(ref_dir, sub_file), '-i', os.path.join(image_dir, sub_file)], stdout = subprocess.PIPE, stderr = subprocess.PIPE) )

    myfile = open(out_file, 'w')
    for p in proc_set:
        if p.poll() is None:
            p.wait()

        (this_ri, mystderr) = p.communicate()
        myfile.write('{0}\n'.format(this_ri))

    myfile.close()

            
    
            
        
def collect_RI_ICA(ref_dir, image_dir, out_file):
    """
    Almost same with the above collect_RI, just slight change for the ICA between-session RI. 
    """
    all_files = [f for f in os.listdir(ref_dir) if fnmatch.fnmatch(f, '*.nii.gz')]
    all_files.sort()
    print 'total file number: {0}'.format(len(all_files))
    hostname = socket.gethostname()
    proc_set = set()
    for sub_file in all_files:
        proc_set.add(subprocess.Popen(['/home/sci/weiliu/projects/group_mrf_lemon/build_' + hostname + '/compareclusterings', '-v', '0', '-r', os.path.join(ref_dir, sub_file), '-i', os.path.join(image_dir, sub_file)], stdout = subprocess.PIPE, stderr = subprocess.PIPE) )

    myfile = open(out_file, 'w')
    for p in proc_set:
        if p.poll() is None:
            p.wait()

        (this_ri, mystderr) = p.communicate()
        myfile.write('{0}\n'.format(this_ri))

    myfile.close()
        
