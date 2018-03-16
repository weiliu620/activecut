import sys
import os
import fnmatch
import shutil
import numpy
import subprocess
import nibabel as nib
from nipy import load_image
import  matplotlib
from nipy.labs.viz import plot_map, mni_sform, coord_transform


    
def get_comp_mean(comp, sub_name, label_map_dir):
    """
    return mean of the map over all bootstrap samples.
    """

    all_sub_files = [ f for f in os.listdir(label_map_dir) if fnmatch.fnmatch(f, sub_name + '*')]
    all_sub_files.sort()

    firstimage = nib.load(os.path.join(label_map_dir, all_sub_files[0]))
    mean = numpy.zeros(firstimage.shape)
    mean = numpy.float32(mean)

    print "get_comp_mean(), comp {0}, total file: {1}. From {2} to {3}".format(comp, len(all_sub_files), all_sub_files[0], all_sub_files[-1])

    comp = float(comp)
    for sf in all_sub_files:
        ci = nib.load(os.path.join(label_map_dir, sf))
        cd = ci.get_data()
        cd[cd > comp] = 0
        cd[cd < comp] = 0
        cd[cd == comp] = 1
        mean = mean + cd
        
    mean = mean / len(all_sub_files)
    return mean

def save_mean_map(sub_name, in_dir, out_dir, method_name):
    """
    Save the mean map of the group map, over all bootstrap samples. This is not for the subject mean map. 
    """
    t3 = "/usr/sci/scratch/weiliu/NYU_test_retest/Yeo2011_7Networks_MNI152_FreeSurferConformed3mm.nii.gz"
    t3_img = load_image(t3)
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    mycolmap = "autumn"
    # mycolmap = matplotlib.pyplot.cm.get_cmap('autumn')
    # alphas = numpy.abs(numpy.linspace(0.5, 1.0, mycolmap.N) )
    # mycolmap._lut[:-3, -1] = alphas
    
    myinterpo = 'nearest'
    matplotlib.pyplot.ioff()

    cut_pos = [(-3, -78, 6), # visual
               (-3, -18, 51), # motor
               (-3, 15, 48),  # dorsal attention.
               (6, 18, -4), # salience
               (0, 45, -24), # brain stem?
               (6, 18, 48), # executive control
               (-6, -66, 27)] # default
               
    for comp, thispos in enumerate(cut_pos):
        mean_vol = get_comp_mean(comp + 1, sub_name, in_dir)
        h = matplotlib.pyplot.figure()
        h.clear()
        print comp, thispos
        plot_map(mean_vol + 1, t3_img.affine, slicer = 'ortho', cut_coords=thispos, figure = h,  threshold = 1.0, vmin=1.0, vmax = 2, cmap = mycolmap, interpolation = myinterpo, draw_cross = True, annotate = False)
        
        out_file = os.path.join(out_dir, sub_name + '_' + method_name + '_mean_' + str(comp + 1).zfill(2) + '.png')
        print "save_mean_map(): for comp {0}, file saved to {1}".format(comp, out_file)
        h.savefig( out_file )
        subprocess.call(['/usr/bin/convert', out_file,  '-crop', '800x315+0+140', out_file ])

def save_mean_map2(sub_name, in_dir, out_dir, method_name):
    """
    This seems for the subject mean map.
    """

    t3 = "/usr/sci/scratch/weiliu/NYU_test_retest/Yeo2011_7Networks_MNI152_FreeSurferConformed3mm.nii.gz"
    t3_img = load_image(t3)
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    mycolmap = "autumn"
    # mycolmap = matplotlib.pyplot.cm.get_cmap('autumn')
    # alphas = numpy.abs(numpy.linspace(0.5, 1.0, mycolmap.N) )
    # mycolmap._lut[:-3, -1] = alphas
    
    myinterpo = 'nearest'
    matplotlib.pyplot.ioff()

    direc = ['z', 'z', 'z', 'x', 'y', 'y', 'z']
    # the second coordinates -10 is a dummy cut, just for get the correct size
    # of the png image.
    cut_pos = [(6, -10), (51, -10), (48, -10), (6, -10), (45, -10), (18, -10), (27, -10)]
               
    for comp, thispos in enumerate(cut_pos):
        mean_vol = get_comp_mean(comp + 1, sub_name, in_dir)
        h = matplotlib.pyplot.figure()
        h.clear()
        plot_map(mean_vol + 1, t3_img.affine, slicer = direc[comp], cut_coords=thispos, figure = h,  threshold = 1, vmin=1, vmax = 2, cmap = mycolmap, interpolation = myinterpo, draw_cross = True, annotate = False)
        
        out_file = os.path.join(out_dir, method_name + '_' +  sub_name + '_' + str(comp + 1).zfill(2) + '.png')
        print "save_mean_map(): for comp {0}, file saved to {1}".format(comp + 1, out_file)        
        h.savefig( out_file )
        subprocess.call(['/usr/bin/convert', out_file,  '-crop', '800x510+400+45', out_file ])

def get_all_var(comp, label_map_dir):
    """
    average variance of the binary component map, across all bootstrap samples and all subjects.
    """
    
    comp = float(comp)
    all_sub_names = [ f for f in os.listdir('/usr/sci/scratch/weiliu/NYU_test_retest/bootstrap/session1_smooth') ]

    all_sub_names.sort()
    all_files = [ f for f in os.listdir(label_map_dir) if fnmatch.fnmatch(f, all_sub_names[0] + '*')]
    firstimage = nib.load(os.path.join(label_map_dir, all_files[0]))
    mean = numpy.zeros(firstimage.shape)
    sub_var  = numpy.zeros(firstimage.shape)
    all_var = numpy.zeros(firstimage.shape)

    for subject in all_sub_names:
        all_files = [ f for f in os.listdir(label_map_dir) if fnmatch.fnmatch(f, subject + '*')]
        all_files.sort()

        mean = numpy.zeros(firstimage.shape)
        for sf in all_files:
            ci = nib.load(os.path.join(label_map_dir, sf))
            cd = ci.get_data()
            cd[cd > comp] = 0
            cd[cd < comp] = 0
            cd[cd == comp] = 1
            mean = mean + cd
            
        mean = mean / len(all_files)

        sub_var  = numpy.zeros(firstimage.shape)
        for sf in all_files:
            ci = nib.load(os.path.join(label_map_dir, sf))
            cd = ci.get_data()
            cd[cd > comp] = 0
            cd[cd < comp] = 0
            cd[cd == comp] = 1
            sub_var = sub_var + (cd - mean)**2

        sub_var = sub_var / len(all_files)

        all_var = all_var + sub_var

    all_var = all_var / len(all_sub_names)

    return all_var

def save_var_map(in_dir, out_dir, method_name):
    """
    Save average variance map.
    """

    cut_pos = [(-3, -78, 6), # visual
               (-3, -18, 51), # motor
               (-3, 15, 48),  # dorsal attention.
               (6, 18, -4), # salience
               (0, 45, -24), # brain stem?
               (6, 18, 48), # executive control
               (-6, -66, 27)] # default


    t3 = "/usr/sci/scratch/weiliu/NYU_test_retest/Yeo2011_7Networks_MNI152_FreeSurferConformed3mm.nii.gz"
    t3_img = load_image(t3)

    matplotlib.pyplot.ioff()

    for comp in range(7):
        var_vol = get_all_var(comp+1, in_dir)
        h = matplotlib.pyplot.figure()
        h.clear()
        # before July 26, the vmin is 1.05 and vmax is 1.25 so the display range
        # is [0.05, 0.25]. after that, we change it to vmin = 1.0 and vmax =
        # 1.15, the display range will be [0, 0.15].
        plot_map(var_vol + 1, t3_img.affine, slicer = 'ortho', cut_coords=cut_pos[comp], figure = h,  threshold = 1.05, vmin=1.0, vmax = 1.15, cmap = 'autumn', interpolation = 'nearest', draw_cross = True, annotate = False)


        out_file = os.path.join(out_dir, method_name + str(comp+1).zfill(2) + '.png') 
        h.savefig(out_file)

        subprocess.call(['/usr/bin/convert', out_file,  '-crop', '800x315+0+140',out_file ])
    

# Jan 5, 2014. To add a t test between K-Means/N-Cuts and HMRF for a specific
# component, I need to compute the variance volume first, and do a t test.
def var_ttest(hmrf_in_dir, method_in_dir, method_name, comp_id):
    """
    Compute the variance volume for the HMRF and the other method, and compare them by a T test.
    
    hmrf_in_dir: the hmrf label map dir.
    method_in_dir: the other method's label map dir.
    method name: a string. Can be 'kmeans' or 'ncuts'.
    comp_id: 1-based component id. An integer. 
    """
    var_hmrf = get_all_var(comp_id, hmrf_in_dir)
    var_control = get_all_var(comp_id, method_in_dir)

    # for now I don't bother choosing a subset of voxels for T test. Just assume
    # all voxels are samples. This may not be the optimal since some voxels are
    # not for this particular component. 
    scipy.stats.ttest_rel(var_hmrf, var_control)
    

def get_comp_meanvar(comp, set_name, label_dir):
    """
    Compute the mean and var map of for one component. Used for both group and subjects.

    Jan 2014: since I didn't find the old func for computing the group var, just write one.
    grp_name: the string for the file name of the group file. Should be just "grp".
    
    comp: comp id.
    set_name: for group, will be 'grp', for subject, will be subject name.
    label_dir: group bootstrap sample label map.

    return: (mean, var) volume
    """
    all_files = [ f for f in os.listdir(label_dir) if fnmatch.fnmatch(f, set_name + '*.nii.gz')]
    all_files.sort()

    firstimage = nib.load(os.path.join(label_dir, all_files[0]))
    mean = numpy.zeros(firstimage.shape)
    mean = numpy.float32(mean)

    comp = float(comp)
    # compute the mean map first.
    for f in all_files:
        ci = nib.load(os.path.join(label_dir, f))
        cd = ci.get_data()
        cd[cd > comp] = 0
        cd[cd < comp] = 0
        cd[cd == comp] = 1
        mean = mean + cd
        
    mean = mean / len(all_files)

    # compute the var. 
    var  = numpy.zeros(firstimage.shape)
    for f in all_files:
        ci = nib.load(os.path.join(label_dir, f))
        cd = ci.get_data()
        cd[cd > comp] = 0
        cd[cd < comp] = 0
        cd[cd == comp] = 1
        var = var + (cd - mean)**2

    var = var / len(all_files)

    return (mean, var)
    

def save_grp_meanvar_17(grp_name, in_dir, out_dir, method_name):
    """
    Jan 2014, for 17 networks, 
    Save the mean and var map of the group map, over all bootstrap samples. This is not for the subject mean map. 
    """
    t3 = "/usr/sci/scratch/weiliu/NYU_test_retest/Yeo2011_7Networks_MNI152_FreeSurferConformed3mm.nii.gz"
    t3_img = load_image(t3)
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    mycolmap = "autumn"
    # mycolmap = matplotlib.pyplot.cm.get_cmap('autumn')
    # alphas = numpy.abs(numpy.linspace(0.5, 1.0, mycolmap.N) )
    # mycolmap._lut[:-3, -1] = alphas
    
    myinterpo = 'nearest'
    matplotlib.pyplot.ioff()

    cut_pos = [(-30, -93, 6), # 1, visual
                (-3, -78, 6), # 2, visual
                (-3, -24, 60), # 3
                (54, -12, 6), # 4,
                (48, -66, 6), # 5,
                (12, -57, 66), # 6
                (-42, 6, 3), # 7
                (27, 48, 30), #8
                (33, -15, -30), # 9
                (-3, 33, -21), #10
                (3, -63, 51), # 11
                (45, 15, 27), # 12
                (48, -51, 51), # 13
                (57, -39, 6), # 14
                (12, -51, 2), # 15
                (6, 51, 3), #16
                (-3, 51, 30)] # 17
               
               
    for comp, thispos in enumerate(cut_pos):
        (mean_vol, var_vol) = get_comp_meanvar(comp + 1, grp_name, in_dir)

        h = matplotlib.pyplot.figure()
        h.clear()
        plot_map(mean_vol + 1, t3_img.affine, slicer = 'ortho', cut_coords=thispos, figure = h,  threshold = 1.0, vmin=1.0, vmax = 2, cmap = mycolmap, interpolation = myinterpo, draw_cross = True, annotate = False)
        
        out_file = os.path.join(out_dir, grp_name + '_' + method_name + '_mean_' + str(comp + 1).zfill(2) + '.png')
        print "save_mean_map(): for comp {0}, file saved to {1}".format(comp + 1, out_file)
        h.savefig( out_file )
        subprocess.call(['/usr/bin/convert', out_file,  '-crop', '800x315+0+140', out_file ])

        h.clear()
        plot_map(var_vol, t3_img.affine, slicer = 'ortho', cut_coords=thispos, figure = h,  threshold = 0.05, vmin = 0.05, vmax = 0.15, cmap = mycolmap, interpolation = myinterpo, draw_cross = True, annotate = False)
        
        out_file = os.path.join(out_dir, grp_name + '_' + method_name + '_var_' + str(comp + 1).zfill(2) + '.png')
        print "save_mean_map(): for comp {0}, file saved to {1}".format(comp + 1, out_file)
        h.savefig( out_file )
        subprocess.call(['/usr/bin/convert', out_file,  '-crop', '800x315+0+140', out_file ])
        
def save_sub_mean_17(sub_name, in_dir, out_dir, method_name):
    """
    Save the mean subject map over all bootstrap samples and over all subjects. For 17 networks. 
    """

    t3 = "/usr/sci/scratch/weiliu/NYU_test_retest/Yeo2011_7Networks_MNI152_FreeSurferConformed3mm.nii.gz"
    t3_img = load_image(t3)
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    mycolmap = "autumn"
    myinterpo = 'nearest'
    matplotlib.pyplot.ioff()

    direc = ['z', 'z', 'z', 'y', 'x',
             'z', 'y', 'z', 'x', 'x',
             'x', 'y', 'z', 'z', 'x',
             'x', 'x']
    # the second coordinates -10 is a dummy cut, just for get the correct size
    # of the png image.
    cut_pos = [(6, -10), #1z
                (6, -10), # 2z
                (60, -10), # 3z
                (-12, -10), #4y
                (48, -10), #5x
                (66, -10), #6z
                (-6, -10), #7y
                (30, -10), #8z
                (33, -10), #9x
                (-3, -10), #10x
                (3, -10), #11
                (15, -10), #12y
                (51, -10), #13z
                (6, -10), #14z
                (12, -10), #15 
                (6, -10), #16 
                (-3, -10)] # 17 
               
    for comp, thispos in enumerate(cut_pos):
        mean_vol = get_comp_mean(comp + 1, sub_name, in_dir)
        h = matplotlib.pyplot.figure()
        h.clear()
        plot_map(mean_vol + 1, t3_img.affine, slicer = direc[comp], cut_coords=thispos, figure = h,  threshold = 1, vmin=1, vmax = 2, cmap = mycolmap, interpolation = myinterpo, draw_cross = True, annotate = False)
        
        out_file = os.path.join(out_dir, method_name + '_' +  sub_name + '_' + str(comp + 1).zfill(2) + '.png')
        print "save_mean_map(): for comp {}, seed {} = {}, file saved to {}".format(comp + 1, direc[comp], thispos, out_file)        
        h.savefig( out_file )
        subprocess.call(['/usr/bin/convert', out_file,  '-crop', '800x510+400+45', out_file ])



def save_var_17(in_dir, out_dir, method_name):
    """
    wrapper for saving subject variance map, over all subjects and bootstrap samples. (17 networks)
    """

    cut_pos = [(-30, -93, 6), # 1, visual
                (-3, -78, 6), # 2, visual
                (-3, -24, 60), # 3
                (54, -12, 6), # 4,
                (48, -66, 6), # 5,
                (12, -57, 66), # 6
                (-42, 6, 3), # 7
                (27, 48, 30), #8
                (33, -15, -30), # 9
                (-3, 33, -21), #10
                (3, -63, 51), # 11
                (45, 15, 27), # 12
                (48, -51, 51), # 13
                (57, -39, 6), # 14
                (12, -51, 2), # 15
                (6, 51, 3), #16
                (-3, 51, 30)] # 17


    t3 = "/usr/sci/scratch/weiliu/NYU_test_retest/Yeo2011_7Networks_MNI152_FreeSurferConformed3mm.nii.gz"
    t3_img = load_image(t3)

    matplotlib.pyplot.ioff()

    for comp, thispos in enumerate(cut_pos):
        var_vol = get_all_var(comp+1, in_dir)
        h = matplotlib.pyplot.figure()
        h.clear()
        # before July 26, the vmin is 1.05 and vmax is 1.25 so the display range
        # is [0.05, 0.25]. after that, we change it to vmin = 1.0 and vmax =
        # 1.15, the display range will be [0, 0.15].
        plot_map(var_vol, t3_img.affine, slicer = 'ortho', cut_coords = thispos, figure = h,  threshold = 0.05, vmin=0.05, vmax = 0.15, cmap = 'autumn', interpolation = 'nearest', draw_cross = True, annotate = False)


        out_file = os.path.join(out_dir, method_name + str(comp+1).zfill(2) + '.png') 
        h.savefig(out_file)

        subprocess.call(['/usr/bin/convert', out_file,  '-crop', '800x315+0+140',out_file ])
        
