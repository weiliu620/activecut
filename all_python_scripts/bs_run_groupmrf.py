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

    fmri_little = os.path.join(bs_dir, 'fmri_little') 

    session_little = './session1_little'
    
    host_name = socket.gethostname()

    if not os.path.exists(bs_dir):
        os.makedirs(bs_dir)

    if not os.path.exists(fmri_little):
        os.makedirs(fmri_little)

    for bsid in range(bs_start, bs_end + 1):
        print('----->>>> Bootstrap {0} begin: <<<-----------'.format(bsid) )
        
        # resampling and normalizaion.
        resample(session_little, fmri_little, bsid)

        one_run(fmri_little, bsid, bs_dir)

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


        
def one_run(fmri_little, bsid, bs_dir):

    init_grp = "/usr/sci/scratch/weiliu/NYU_test_retest/Yeo2011_7Networks_MNI152_FreeSurferConformed3mm.nii.gz"
    cumusamples_dir = os.path.join(bs_dir, 'cumusamples')
    rsamples_dir = os.path.join(bs_dir, 'rsamples')
    mask_file = "/usr/sci/scratch/weiliu/NYU_test_retest/Yeo2011_7Networks_MNI152_FreeSurferConformed3mm.nii.gz"
    groupmrf_out_dir = "/usr/sci/scratch/weiliu/NYU_test_retest/bootstrap/groupmrf/out"

    bs_ncuts = os.path.join(bs_dir, 'ncuts')
    host_name = socket.gethostname()

    if not os.path.exists(cumusamples_dir):
        os.makedirs(cumusamples_dir)

    if not os.path.exists(rsamples_dir):
        os.makedirs(rsamples_dir)

    if not os.path.exists(bs_ncuts):
        os.makedirs(bs_ncuts)

    bsfix = '_' + str(bsid).zfill(3) # bootstrap id surfix.
    proc_set = set()


    groupmrf_proc = subprocess.Popen(['/home/sci/weiliu/projects/group_mrf_lemon/build_' + host_name + '/groupmrf',  '--nthreads', '100', '--sweepperthread', '1', '--burnin', '20', '--numSamples', '100', '--emiter', '50', '--alpha', '0.4', '--beta', '0.8', '--gamma', '1', '-k', '7', '--inittemp', '0.5', '--finaltemp', '0.5', '--initsame', '-i', init_grp, '-f', fmri_little, '--cumusampledir', cumusamples_dir, '--rsampledir', rsamples_dir, '-v', '1'])

    # wait until groupmrf finishes, then change groupmrf output filenames
    if groupmrf_proc.poll() is None:
        groupmrf_proc.wait()

    for file in os.listdir(cumusamples_dir):
        sub_name = os.path.splitext(file)[0]
        sub_name = os.path.splitext(sub_name)[0]
        new_file = sub_name + bsfix + '.nii.gz'

        # os.rename(os.path.join(cumusamples_dir, file), os.path.join(groupmrf_out_dir, new_file) )
        sample2label_set = set()
        sample2label_set.add( subprocess.Popen(['/home/sci/weiliu/projects/group_mrf_lemon/build_' + host_name + '/sample2label', '-i', os.path.join(cumusamples_dir, file), '-o', os.path.join(groupmrf_out_dir, new_file)]) )

    for p in sample2label_set:
        if p.poll() is None:
            p.wait()


if __name__ == '__main__':
    bs_run(sys.argv[1], sys.argv[2], sys.argv[3])
