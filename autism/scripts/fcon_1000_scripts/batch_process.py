import os
import subprocess
import fnmatch
import sys
import random
import shutil

def batch_process(site_dir):
    """
    fcon batch_process script, and save to session_dir
    """
    anat_name = 'mprage'
    func_name = 'rest'
    standard_brain = './templates/MNI152_T1_3mm_brain.nii.gz'

    all_sub_folders = os.listdir(site_dir)
    all_sub_folders.sort()

    proc_set = set()
    for sub_name in all_sub_folders:
        print "./1_anatpreproc.sh" +" " + sub_name + " " + os.path.join(site_dir, sub_name, 'session_1') + " " +  anat_name
        proc_set.add(subprocess.Popen("./1_anatpreproc.sh" +" " + sub_name + " " + os.path.join(site_dir, sub_name, 'session_1') + " " +  anat_name, shell = True ) )

    for p in proc_set:
        if p.poll() is None:
            p.wait()

    proc_set.clear()
    for sub_name in all_sub_folders:
        proc_set.add( subprocess.Popen("./2_funcpreproc.sh" + " " +  sub_name + " " + os.path.join(site_dir, sub_name, 'session_1') + " " + func_name + " " + "0" + " " + "239" + " " + "2", shell = True) )

    for p in proc_set:
        if p.poll() is None:
            p.wait()

    proc_set.clear()
    for sub_name in all_sub_folders:
        proc_set.add( subprocess.Popen("./3_registration.sh" + " " + sub_name + " " + os.path.join(site_dir, sub_name, 'session_1') + " " + anat_name + " " + standard_brain, shell = True) )

    for p in proc_set:
        if p.poll() is None:
            p.wait()

    proc_set.clear()
    for sub_name in all_sub_folders:
        print "./4_segment.sh " + sub_name + " " +  os.path.join(site_dir, sub_name, 'session_1') + " " + anat_name + " " + func_name +  " ./tissuepriors/3mm"
        proc_set.add( subprocess.Popen("./4_segment.sh " + sub_name + " " +  os.path.join(site_dir, sub_name, 'session_1') + " " + anat_name + " " + func_name +  " /home/sci/weiliu/projects/autism/scripts/fcon_1000_scripts/tissuepriors/3mm", shell = True) )

    for p in proc_set:
        if p.poll() is None:
            p.wait()


    proc_set.clear()
    for sub_name in all_sub_folders:
        proc_set.add( subprocess.Popen("./5_nuisance.sh " + sub_name + " " + os.path.join(site_dir, sub_name, 'session_1') + " " + func_name + " 2 " + "240 " + "/home/sci/weiliu/projects/autism/scripts/fcon_1000_scripts/templates/nuisance.fsf", shell = True) )

    for p in proc_set:
        if p.poll() is None:
            p.wait()

if __name__ == '__main__':
    batch_process(sys.argv[1])
