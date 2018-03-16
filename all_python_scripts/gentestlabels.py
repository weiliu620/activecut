import sys
import os
import fnmatch
import subprocess
import shutil
import socket
import math
import time
import random
        
def genallsamples3(alpha, outdir, num_samples):
    """
    draw samples of test subject label map, given group label map. Save sampels in separate files.
    """
    groupmap = "/usr/sci/scratch/weiliu/NYU_test_retest/Yeo2011_7Networks_MNI152_FreeSurferConformed3mm.nii.gz"
    beta = "1"
    burnin = "10000"
    initmap = "/usr/sci/scratch/weiliu/NYU_test_retest/Yeo2011_7Networks_MNI152_FreeSurferConformed3mm.nii.gz"

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    proc_set = set()
    hostname = socket.gethostname()

    # sampling in blocks. with-in block is parallel, between-block is done in
    # serial.
    block_size = 150
    num_blocks = int(math.ceil(num_samples / block_size));
    for block_id in range(num_blocks):
        start_sid = block_id * block_size;
        if block_id < num_blocks - 1:
            end_sid = (block_id + 1) * block_size - 1
        else:
            end_sid = num_samples - 1

        print 'Samling {0} to {1}'.format(start_sid, end_sid)
        for sid in range(start_sid, end_sid + 1):
            proc_set.add(subprocess.Popen(['/home/sci/weiliu/projects/group_mrf_lemon/build_' + hostname + '/gentestlabel2', '-k', '7', '-i', groupmap, '-b', burnin, '--beta', beta, '--alpha', alpha, '-t', initmap, '-o', os.path.join(outdir, 'yt' + str(sid).zfill(4) + '.nii.gz'), '--seed', str(sid), '-v', '0'] ) )

        for p in proc_set:
            if p.poll() is None:
                p.wait()

def genallsamples4(alpha, outdir, num_samples, initpercent, bestfit_map, burnin):
    """
    draw samples of test subject label map, given group label map. Save sampels in separate files. Save with genallsamples3 but use a process pool for implementation.
    """
    groupmap = "/usr/sci/scratch/weiliu/NYU_test_retest/Yeo2011_7Networks_MNI152_FreeSurferConformed3mm.nii.gz"
    beta = "1"

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    proc_set = set()
    hostname = socket.gethostname()
    burnin = str(burnin)

    RAND_START = 0
    RAND_END = initpercent[0]
    BESTFIT_START = initpercent[0]
    BESTFIT_END = initpercent[0] + initpercent[1]
    GRP_START = initpercent[0] + initpercent[1]
    GRP_END = 1

    # sampling in blocks. with-in block is parallel, between-block is done in
    # serial.

    sid_set = range(num_samples)
    proc_set = set()
    pool_size = 160
    while True:
        new_proc_set = proc_set.copy()
        for p in proc_set:
            if p.poll() is not None:
                new_proc_set.remove(p)

        proc_set = new_proc_set.copy()

        while len(proc_set) < pool_size and len(sid_set) > 0:
            sid = sid_set.pop(0)
            myinitchoice = random.random()
            if myinitchoice > RAND_START and myinitchoice < RAND_END:
                # init with random map.
                proc_set.add(subprocess.Popen(['/home/sci/weiliu/projects/group_mrf_lemon/build_' + hostname + '/gentestlabel2', '-k', '7', '-i', groupmap, '-b', burnin, '--beta', beta, '--alpha', alpha, '-o', os.path.join(outdir, 'yt' + str(sid).zfill(4) + '.nii.gz'), '--seed', str(sid), '-v', '0'], stdout=subprocess.PIPE ) )

            elif myinitchoice > BESTFIT_START and myinitchoice < BESTFIT_END:
                # init with bestfit map.
                proc_set.add(subprocess.Popen(['/home/sci/weiliu/projects/group_mrf_lemon/build_' + hostname + '/gentestlabel2', '-k', '7', '-i', groupmap, '-b', burnin, '--beta', beta, '--alpha', alpha, '-t', bestfit_map, '-o', os.path.join(outdir, 'yt' + str(sid).zfill(4) + '.nii.gz'), '--seed', str(sid), '-v', '0'], stdout=subprocess.PIPE ) )
            else:
                # init with group map.
                proc_set.add(subprocess.Popen(['/home/sci/weiliu/projects/group_mrf_lemon/build_' + hostname + '/gentestlabel2', '-k', '7', '-i', groupmap, '-b', burnin, '--beta', beta, '--alpha', alpha, '-t', groupmap, '-o', os.path.join(outdir, 'yt' + str(sid).zfill(4) + '.nii.gz'), '--seed', str(sid), '-v', '0'], stdout=subprocess.PIPE ) )

        if len(proc_set) == 0:
            break
        time.sleep(2)

def compll(indir, fmri_file):
    """
    compute log likelihood given all test sub samples in a dir.
    """
    ll = 0;
    allsamples = [f for f in os.listdir(indir) if fnmatch.fnmatch(f, 'yt*.nii.gz')]
    allsamples.sort()
    num_samples = len(allsamples)
    hostname = socket.gethostname()

    proc_set = set()
    pool_size = 300
    while True:
        new_proc_set = proc_set.copy()
        for p in proc_set:
            if p.poll() is not None:
                (this_ll, mystderr) = p.communicate()
                # print 'this_ll = {0}, stderr = {1}'.format(this_ll, mystderr)
                this_ll = float(this_ll)
                ll = ll + this_ll
                new_proc_set.remove(p)
                
        proc_set = new_proc_set.copy()

        while len(proc_set) < pool_size and len(allsamples) > 0:
            fname = allsamples.pop(0)
            # print 'adding process for {0}'.format(fname)
            proc_set.add(subprocess.Popen(['/home/sci/weiliu/projects/group_mrf_lemon/build_' + hostname + '/margll', '-i', os.path.join(indir, fname), '-f', fmri_file], stdout=subprocess.PIPE ) )

        if len(proc_set) == 0:
            break
      
        time.sleep(0.5)

    ll = ll / num_samples
    return ll

def testallapha(alpha, initpercent, burnin, sample_dir):
    """
    generate samples for different alpha, and compute log-likelihood. Use test subject label map to initialize the sampling.
    """
    
    # sample_dir = "/usr/sci/scratch/weiliu/NYU_test_retest/testalpha_bayes/sample_dir"
    bestfit_dir = "/usr/sci/scratch/weiliu/NYU_test_retest/testalpha_bayes/bestfit_labels"
    test_fmri_dir = "/usr/sci/scratch/weiliu/NYU_test_retest/littlesmooth/s1/allfunc_normalized/"
    
    num_samples = 160
    all_test_fmri = [f for f in os.listdir(test_fmri_dir) if fnmatch.fnmatch(f, '*.nii.gz')]
    all_test_fmri.sort()
        
    ll_set = list()
    ll = 0
    for test_fmri in all_test_fmri:
        genallsamples4(str(alpha), sample_dir, num_samples, initpercent, os.path.join(bestfit_dir, test_fmri), burnin)
        this_ll = compll(sample_dir, os.path.join(test_fmri_dir, test_fmri))
            
        print '    test sub: {0}. ll = {1:E}'.format(test_fmri, this_ll)
        ll = ll + this_ll

    ll = ll / len(all_test_fmri)
    print 'test sub log likelihood: {0:E}\n'.format(ll)
                
def bestfit_labelmap(test_fmri_dir, outdir):
    """
    Compute a bestfit (in the sense of likelihood) label map for each subject fmri files.
    """
    all_test_fmri = [f for f in os.listdir(test_fmri_dir) if fnmatch.fnmatch(f, '*.nii.gz')]
    all_test_fmri.sort()
    for test_fmri in all_test_fmri:
        # get a kmeans segmentation of test subject.
        subprocess.call(['/home/sci/weiliu/projects/group_mrf/build_kourosh/subkmeans', '-k', '7', '-f', os.path.join(test_fmri_dir, test_fmri), '-o', os.path.join(outdir, 'test_sub_kmeans_seg.nii.gz'), '-m', '/usr/sci/scratch/weiliu/NYU_test_retest/Yeo2011_7Networks_MNI152_FreeSurferConformed3mm.nii.gz'])
        # compute test subject labbel map as initial map.
        subprocess.call(['/home/sci/weiliu/projects/mcem_lemon/build_kourosh/mcem_lemon', '--inittemp', '0.1', '--finaltemp', '0.1', '-n', '1', '--beta', '0', '-k', '7', '-i', os.path.join(outdir, 'test_sub_kmeans_seg.nii.gz'), '-f', os.path.join(test_fmri_dir, test_fmri), '--cumusample', os.path.join(outdir, 'test_sub_cumusample.nii.gz')])
        subprocess.call(['/home/sci/weiliu/projects/group_mrf_lemon/build_kourosh/sample2label', '-i', os.path.join(outdir, 'test_sub_cumusample.nii.gz'), '-o', os.path.join(outdir, test_fmri)])

        os.remove(os.path.join(outdir, 'test_sub_cumusample.nii.gz'))
        os.remove(os.path.join(outdir, 'test_sub_kmeans_seg.nii.gz'))
