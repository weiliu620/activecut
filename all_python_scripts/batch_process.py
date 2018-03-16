
import os
import subprocess
import fnmatch
import sys
import random
import shutil

def batch_process(session_dir, fwhm):
    """
    fcon batch_process script, and save to session_dir
    """

    session_dir_orig = '/usr/sci/scratch/weiliu/nyu_trt_zip/unzipped/session1'
    shutil.rmtree(session_dir)
    shutil.copytree(session_dir_orig, session_dir)

    scripts_dir = './'
    anat_name = 'mprage_skullstripped'
    func_name = 'lfo'
    standard_brain = './templates/MNI152_T1_3mm_brain.nii.gz'

    all_sub_folders = os.listdir(session_dir)
    all_sub_folders.sort()

    proc_set = set()
    for sub_name in all_sub_folders:
        proc_set.add(subprocess.Popen("./1_anatpreproc.sh" +" " + sub_name + " " + os.path.join(session_dir, sub_name) + " " +  anat_name, shell = True ) )

    for p in proc_set:
        if p.poll() is None:
            p.wait()

    proc_set.clear()
    for sub_name in all_sub_folders:
        proc_set.add( subprocess.Popen("./2_funcpreproc.sh" + " " +  sub_name + " " + os.path.join(session_dir, sub_name) + " " + func_name + " " + "0" + " " + "196" + " " + "2" + " " + fwhm, shell = True) )

    for p in proc_set:
        if p.poll() is None:
            p.wait()


    proc_set.clear()
    for sub_name in all_sub_folders:
        proc_set.add( subprocess.Popen("./3_registration.sh" + " " + sub_name + " " + os.path.join(session_dir, sub_name) + " " + anat_name + " " + standard_brain, shell = True) )

    for p in proc_set:
        if p.poll() is None:
            p.wait()




    proc_set.clear()
    for sub_name in all_sub_folders:
        proc_set.add( subprocess.Popen("./4_segment.sh " + sub_name + " " +  os.path.join(session_dir, sub_name) + " " + anat_name + " " + func_name +  " ./tissuepriors/3mm/", shell = True) )


    for p in proc_set:
        if p.poll() is None:
            p.wait()


    proc_set.clear()
    for sub_name in all_sub_folders:
        proc_set.add( subprocess.Popen("./5_nuisance.sh " + sub_name + " " + os.path.join(session_dir, sub_name) + " " + func_name + " 2 " + "197 " + "templates/nuisance.fsf", shell = True) )

    for p in proc_set:
        if p.poll() is None:
            p.wait()



if __name__ == '__main__':
    batch_process(sys.argv[1], sys.argv[2])
