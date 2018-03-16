import os
import subprocess
import fnmatch
import sys
import random
import shutil

def batch_process(sub_dir):
    """
    fcon batch_process script. This is different from the batch_process.py only that it takes single subject's dir as input. 
    """
    anat_name = 'mprage'
    func_name = 'rest'
    standard_brain = './templates/MNI152_T1_3mm_brain.nii.gz'

    # get the subject id from the subject dir.
    sub_name = os.path.split(os.path.abspath(sub_dir))[1]
    site_dir = os.path.split(os.path.abspath(sub_dir))[0]
    print sub_name
    print sub_dir
    
    print "./1_anatpreproc.sh" +" " + sub_name + " " + os.path.join(site_dir, sub_name, 'session_1') + " " +  anat_name
    subprocess.call("./1_anatpreproc.sh" +" " + sub_name + " " + os.path.join(site_dir, sub_name, 'session_1') + " " +  anat_name, shell = True )
    
    subprocess.call("./2_funcpreproc.sh" + " " +  sub_name + " " + os.path.join(site_dir, sub_name, 'session_1') + " " + func_name + " " + "0" + " " + "239" + " " + "2", shell = True)

    subprocess.call("./3_registration.sh" + " " + sub_name + " " + os.path.join(site_dir, sub_name, 'session_1') + " " + anat_name + " " + standard_brain, shell = True) 
    
    # print "./4_segment.sh " + sub_name + " " +  os.path.join(site_dir, sub_name, 'session_1') + " " + anat_name + " " + func_name +  " ./tissuepriors/3mm"
    subprocess.call("./4_segment.sh " + sub_name + " " +  os.path.join(site_dir, sub_name, 'session_1') + " " + anat_name + " " + func_name +  " /home/sci/weiliu/projects/autism/scripts/fcon_1000_scripts/tissuepriors/3mm", shell = True) 
    subprocess.call("./5_nuisance.sh " + sub_name + " " + os.path.join(site_dir, sub_name, 'session_1') + " " + func_name + " 2 " + "240 " + "/home/sci/weiliu/projects/autism/scripts/fcon_1000_scripts/templates/nuisance.fsf", shell = True) 

if __name__ == '__main__':
    batch_process(sys.argv[1])
