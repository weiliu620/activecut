#!/usr/bin/env python

# add the path of other python scripts.

import os
import subprocess
import fnmatch
import sys
import random
import shutil
import socket

def bs_run(bs_start, bs_end, bs_dir):
    """
    bootstrap of groupmrf on read data.

    """
    bs_start = int(bs_start)
    bs_end = int(bs_end)

    fmri_smooth = os.path.join(bs_dir, 'fmri_smooth') 
    bs_ncuts = os.path.join(bs_dir, 'ncuts')

    session_smooth = './session1_smooth'
    
    host_name = socket.gethostname()

    if not os.path.exists(bs_dir):
        os.makedirs(bs_dir)

    if not os.path.exists( fmri_smooth ):
        os.makedirs( fmri_smooth )

    if not os.path.exists(bs_ncuts):
        os.makedirs(bs_ncuts)

    for bsid in range(bs_start, bs_end + 1):
        print('----->>>> Bootstrap {0} begin: <<<-----------'.format(bsid) )
        
        # resampling and normalizaion.
        resample(session_smooth, fmri_smooth, bsid)

        one_run(fmri_smooth, bsid, bs_dir)

    return

def resample(session_dir, normalized_dir, seed):
    """
    bootstrap resampling, and then project sampls to sphere.

    This is for both smooth and little smoothed data.
    """

    host_name = socket.gethostname()

    # subject folder
    all_sub_folders = [ thisfile for thisfile in os.listdir(session_dir)]
    all_sub_folders.sort()
    proc_set = set()    


    for sub_name in all_sub_folders:
        proc_set.add(subprocess.Popen(['/home/sci/weiliu/projects/group_mrf_lemon/build_' + host_name + '/fmriresample', '-i', os.path.join(session_dir, sub_name, 'func/lfo_res2standard.nii.gz'), '-o', os.path.join(normalized_dir, sub_name + '.nii.gz'), '-k', '197', '--seed', str(seed) ]) )

    # wait until all fmriresample finishes
    for p in proc_set:
        if p.poll() is None:
            p.wait()


    # Normalizatoin to sphere.
    all_subjects_files = os.listdir(normalized_dir)
    proc_set.clear()
    for sub_file in all_subjects_files:
        proc_set.add(subprocess.Popen(['/home/sci/weiliu/projects/group_mrf/build_' +host_name + '/projectfmri', '-i', os.path.join(normalized_dir, sub_file), '-o', os.path.join(normalized_dir, sub_file)]) )

    # wait until all fmriresample finishes
    for p in proc_set:
        if p.poll() is None:
            p.wait()


        
def one_run(fmri_smooth, bsid, bs_dir):

    mask_file = "/usr/sci/scratch/weiliu/NYU_test_retest/Yeo2011_7Networks_MNI152_FreeSurferConformed3mm.nii.gz"
    init_grp = "/usr/sci/scratch/weiliu/NYU_test_retest/Yeo2011_7Networks_MNI152_FreeSurferConformed3mm.nii.gz"

    bs_ncuts = os.path.join(bs_dir, 'ncuts')
    host_name = socket.gethostname()

    if not os.path.exists(bs_ncuts):
        os.makedirs(bs_ncuts)

    bsfix = '_' + str(bsid).zfill(3) # bootstrap id surfix.
    proc_set = set()


    # N-Cuts
    allsubfiles = os.listdir(fmri_smooth)
    ncuts_proc = set()
    for fmri_file in allsubfiles:
        sub_name = os.path.splitext(fmri_file)[0] 
        sub_name = os.path.splitext(sub_name)[0] 
        ncuts_sub_cmd =  'matlab -nodesktop  -nosplash -r "addpath ~/projects/ncuts_fmri; dbstop if error; nyu_segts(\'' + os.path.join(fmri_smooth, fmri_file) + '\', \'' + mask_file + '\', 7, 0.4, \'' + os.path.join(bs_ncuts, fmri_file) + '\',\'' + os.path.join(bs_ncuts, sub_name+'.mat') + '\'); quit"' 
        ncuts_proc.add(subprocess.Popen(ncuts_sub_cmd, shell = True) )

    for p in ncuts_proc:
        if p.poll() is None:
            p.wait()

    ncuts_grp_cmd =  'matlab -nodesktop  -nosplash -r "addpath ~/projects/ncuts_fmri; dbstop if error; nyu_groupcuts(\'' + bs_ncuts + '\', 0.4, \'' + mask_file + '\', 7, \'' + os.path.join(bs_ncuts, 'grp.nii.gz') + '\'); quit"'
    ncuts_grp_proc = subprocess.Popen(ncuts_grp_cmd, shell = True)
                     

    # wait until ncut group finishes. Align ncuts (sub and grp) label map to
    # Yeo map and save into gz file.
    if ncuts_grp_proc.poll() is None:
        ncuts_grp_proc.wait()

    all_ncuts_files = [ f for f in os.listdir(bs_ncuts) if fnmatch.fnmatch(f, '*.nii.gz') ]
    all_ncuts_files.sort()

    align_proc_set = set()
    for thisfile in all_ncuts_files:
        sub_name = os.path.splitext(thisfile)[0] 
        sub_name = os.path.splitext(sub_name)[0] 
        new_file = sub_name + bsfix + '.nii.gz'
        align_proc_set.add(subprocess.Popen(['/home/sci/weiliu/projects/group_mrf/build_' + host_name + '/alignlabels', '-t', os.path.join(bs_ncuts, thisfile), '-r', init_grp, '-o', os.path.join('ncuts', new_file), '-k', '7']) )

    for p in align_proc_set:
        if p.poll() is None:
            p.wait()

if __name__ == '__main__':
    bs_run(sys.argv[1], sys.argv[2], sys.argv[3])
