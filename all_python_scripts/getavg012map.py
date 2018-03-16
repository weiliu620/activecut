import sys
import os
import fnmatch
import shutil
import numpy
import subprocess
import nibabel as nib
from nipy.labs.viz import plot_map, mni_sform, coord_transform
import  matplotlib

def getavg012map(testname, out_dir):
    """
    Get the 0-1-2 variance map for a test.
    """
    
    nyu_root_dir = '/usr/sci/scratch/weiliu/NYU_test_retest'

    # just choose a folder to know all subject filenames.
    all_sub_files = [ f for f in os.listdir(os.path.join(nyu_root_dir, 'preprocessed_nosmooth/session1_stage2b/test1b/cumusamples') ) if fnmatch.fnmatch(f, 'sub*') ]
    all_sub_files.sort()   

    # groupmrf session folder. 
    S1_dir = os.path.join(nyu_root_dir, 'preprocessed_nosmooth/session1_stage2b', testname, 'labelmap')
    S2_dir = os.path.join(nyu_root_dir, 'preprocessed_nosmooth/session2_stage2b', testname, 'labelmap')
    S3_dir = os.path.join(nyu_root_dir, 'preprocessed_nosmooth/session3_stage2b', testname, 'labelmap')

    getmethodmap(all_sub_files, S1_dir, S2_dir, S3_dir, os.path.join(out_dir, 'groupmrf' + testname + '.nii.gz') )

    # kmeans session folder
    S1_dir = os.path.join(nyu_root_dir, 'preprocessed_smoothed/session1_stage2b/subinit')
    S2_dir = os.path.join(nyu_root_dir, 'preprocessed_smoothed/session2_stage2b/subinit')
    S3_dir = os.path.join(nyu_root_dir, 'preprocessed_smoothed/session3_stage2b/subinit')
    getmethodmap(all_sub_files, S1_dir, S2_dir, S3_dir, os.path.join(out_dir, 'kmeans.nii.gz') )

    # ncuts session folder
    S1_dir = os.path.join(nyu_root_dir, 'preprocessed_smoothed/session1_stage2b/ncuts')
    S2_dir = os.path.join(nyu_root_dir, 'preprocessed_smoothed/session2_stage2b/ncuts')
    S3_dir = os.path.join(nyu_root_dir, 'preprocessed_smoothed/session3_stage2b/ncuts')
    getmethodmap(all_sub_files, S1_dir, S2_dir, S3_dir, os.path.join(out_dir, 'ncuts.nii.gz') )




def getmethodmap(all_sub_files, S1, S2, S3, out_file):
    """
    Get 0-1-2 map for one method. The final results is added 1, so inbrrain and outbrain is different.
    """
    maskfile="/usr/sci/scratch/weiliu/NYU_test_retest/Yeo2011_7Networks_MNI152_FreeSurferConformed3mm.nii.gz"
    subprocess.call(['/home/sci/weiliu/packages/fsl5/bin/fslmaths', maskfile, '-thr', '100', out_file])

    # each map has max value 2, and there are 25 such maps adding together.
    max_var = 2 * len(all_sub_files)
    for sub_file in all_sub_files:
        subprocess.call(['/home/sci/weiliu/projects/group_mrf_lemon/build_wukong/get012', '--s1', os.path.join(S1, sub_file), '--s2', os.path.join(S2, sub_file), '--s3', os.path.join(S3, sub_file), '-m', maskfile, '-o', 'this012map.nii.gz'])
        subprocess.call(['/home/sci/weiliu/projects/group_mrf_lemon/build_wukong/copyheader', '-s', maskfile, '-i', 'this012map.nii.gz', '-o', 'this012map.nii.gz'])
        subprocess.call(['/home/sci/weiliu/packages/fsl5/bin/fslmaths', out_file, '-add', 'this012map.nii.gz', out_file])
        
    subprocess.call(['/home/sci/weiliu/packages/fsl5/bin/fslmaths', out_file, '-div', str(max_var),  out_file])
    subprocess.call(['/home/sci/weiliu/packages/fsl5/bin/fslmaths', maskfile, '-bin', 'binmaskfile.nii.gz'])
    subprocess.call(['/home/sci/weiliu/packages/fsl5/bin/fslmaths', out_file, '-add', 'binmaskfile.nii.gz', out_file])
        
    os.remove('this012map.nii.gz')
    os.remove('binmaskfile.nii.gz')
        
def draw012map(in_file, myslicer, my_cut_coords, out_file):
    """
    Draw the 012 variance map and save into a png/pdf file.
    """

    # mycmap = 'RdBu'
    mycmap = 'RdYlBu'
    # mycmap = 'GnBu'
    matplotlib.pyplot.ioff()
    t3 = "/usr/sci/scratch/weiliu/NYU_test_retest/Yeo2011_7Networks_MNI152_FreeSurferConformed3mm.nii.gz"
    t3_img = nib.load(t3)

    img = nib.load(in_file)
    matplotlib.use('Agg')
    h = matplotlib.pyplot.figure()

    h.clear()
    plot_map(img.get_data(), t3_img.get_affine(), slicer = myslicer, cut_coords = my_cut_coords, figure = h,  threshold = 0.99, vmin = 1, vmax = 1.8, cmap = mycmap, interpolation = 'nearest', draw_cross = True, annotate = False)
    h.savefig(out_file)

    subprocess.call(['/usr/bin/convert', out_file,  '-crop', '800x315+0+140', out_file ]) 
    
def draw012mapwrapper(var_dir):
    """
    Draw 012 map for all 3 methods. Each method has two slices at x, y, z direction respectively. Run this fun after getavg012map.
    """
    # all_cuts = [(0, 0, 0), (10, 10, 10), (-5, 30, 30)]
    # all_cuts = [(-3, -78, 6), (-3, 15, 48), (-6, -66, 27)] // before July 22
    all_cuts = [(-3, -78, 6), (-3, -50, 48), (-6, -66, 27)] 

    all_files = [f for f in os.listdir(var_dir) if fnmatch.fnmatch(f, '*.nii.gz') ]
    for thisfile in all_files:
        methodname = os.path.splitext(thisfile)[0]
        methodname = os.path.splitext(methodname)[0]
        for i, thiscut in enumerate(all_cuts):
            draw012map(os.path.join(var_dir, thisfile), 'ortho', thiscut, os.path.join(var_dir, methodname + str(i).zfill(2) + '.png'))
                   

