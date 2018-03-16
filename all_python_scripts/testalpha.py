#!/usr/bin/env python 
import subprocess
import sys
import fnmatch
import os
import shlex

def testalpha(fmri_dir, out_dir):

    alpha_all = [0.0, 0.3, 0.5, 0.8, 1.0]
    initgrpimage = "/usr/sci/scratch/weiliu/NYU_test_retest/Yeo2011_7Networks_MNI152_FreeSurferConformed3mm.nii.gz"

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    groupmrf_set = set()
    all_files = set()
    for alpha in alpha_all:
        cumusamples_dir = os.path.join(out_dir, 'alpha_' + str(alpha), 'cumusamples')
        rsamples_dir = os.path.join(out_dir, 'alpha_' + str(alpha), 'rsamples')
        if not os.path.exists(cumusamples_dir):
            os.makedirs(cumusamples_dir)
        if not os.path.exists(rsamples_dir):
            os.makedirs(rsamples_dir)

        mystdoutfile = open(os.path.join(out_dir, 'stdout_' + str(alpha) + '.txt'), 'w')
        mystderrfile = open(os.path.join(out_dir, 'stderr_' + str(alpha) + '.txt'), 'w')

        all_files.add(mystdoutfile)
        all_files.add(mystderrfile)
            
        cmd = ['/home/sci/weiliu/projects/group_mrf_lemon/build_uv/groupmrf',  '--nthreads',  '100', '--sweepperthread', '1', '--burnin', '20', '--numSamples', '100', '--emiter', '30', '--alpha', str(alpha), '--beta', '0.8', '--gamma', '1', '-k', '7', '--inittemp', '1.0', '--finaltemp', '1.0',  '--initsame', '-i', initgrpimage, '-f', fmri_dir, '--cumusampledir', cumusamples_dir, '--rsampledir', rsamples_dir, '--verbose', '2']

        groupmrf_set.add(subprocess.Popen(cmd, stdout = mystdoutfile, stderr = mystderrfile) )

    for p in groupmrf_set:
        if p.poll() is None:
            p.wait()

                
    for file in all_files:
        file.flush()
        file.close()
        
if __name__ == '__main__':
    testalpha(sys.argv[1], sys.argv[2])        

    
