#!/usr/bin/env python

# compute mean and variance map 

import os
import fnmatch
import sys
import shutil
import nibabel as nib
import numpy

def kmeansvar(comp, sub_name, label_map_dir, out_filename):
    """
    The consistency across all subjects and bootstrap samples.
    """

    comp = float(comp)
    all_sub_files = [ f for f in os.listdir(label_map_dir) if fnmatch.fnmatch(f, sub_name + '*')]
    all_sub_files.sort()

    firstimage = nib.load(os.path.join(label_map_dir, all_sub_files[0]))
    mean = numpy.zeros(firstimage.shape)
    mean = numpy.float32(mean)


    for sf in all_sub_files:
        print "computing mean, at {0}".format(sf)
        ci = nib.load(os.path.join(label_map_dir, sf))
        cd = ci.get_data()
        cd[cd > comp] = 0
        cd[cd < comp] = 0
        mean = mean + cd
        
    mean = mean / len(all_sub_files)

    # consistency
    co  = numpy.zeros(firstimage.shape)
    co = numpy.float32(co)
    
    for sf in all_sub_files:
        ci = nib.load(os.path.join(label_map_dir, sf))
        cd = ci.get_data()
        cd[cd > comp] = 0
        cd[cd < comp] = 0
        co = co + (cd - mean)**2

    co = co / len(all_sub_files)


    outimg = nib.Nifti1Image(co, numpy.eye(4))
    outimg.to_filename(out_filename)



def kmeans_mean(comp, sub_name, label_map_dir, out_filename):
    """
    """

    comp = float(comp)
    all_sub_files = [ f for f in os.listdir(label_map_dir) if fnmatch.fnmatch(f, sub_name + '*')]
    all_sub_files.sort()

    firstimage = nib.load(os.path.join(label_map_dir, all_sub_files[0]))
    mean = numpy.zeros(firstimage.shape)
    mean = numpy.float32(mean)


    for sf in all_sub_files:
        print "computing mean, at {0}".format(sf)
        ci = nib.load(os.path.join(label_map_dir, sf))
        cd = ci.get_data()
        cd[cd > comp] = 0
        cd[cd < comp] = 0
        mean = mean + cd
        
    mean = mean / len(all_sub_files)

    mean = mean + 1

    outimg = nib.Nifti1Image(mean, numpy.eye(4))
    outimg.to_filename(out_filename)


def groupmrf_mean(comp, sub_name, label_map_dir, out_filename, num_mcsamples):
    """
    """

    comp = int(comp)
    all_sub_files = [ f for f in os.listdir(label_map_dir) if fnmatch.fnmatch(f, sub_name + '*')]
    all_sub_files.sort()

    firstimage = nib.load(os.path.join(label_map_dir, all_sub_files[0]))
    mean = numpy.zeros(firstimage.shape[0:3])
    mean = numpy.float32(mean)

    for sf in all_sub_files:
        print "computing mean, at {0}".format(sf)
        sample_4d = nib.load(os.path.join(label_map_dir, sf))
        list_3d = nib.funcs.four_to_three(sample_4d)
        mean = mean + list_3d[comp-1].get_data()
        
    mean = mean / len(all_sub_files)

    mean = mean / num_mcsamples

    mean = mean + 1

    outimg = nib.Nifti1Image(mean, numpy.eye(4))
    outimg.to_filename(out_filename)


def groupmrf_var(comp, sub_name, label_map_dir, out_filename):
    """
    """

    comp = int(comp)
    all_sub_files = [ f for f in os.listdir(label_map_dir) if fnmatch.fnmatch(f, sub_name+'*')]
    all_sub_files.sort()

    firstimage = nib.load(os.path.join(label_map_dir, all_sub_files[0]))
    mean = numpy.zeros(firstimage.shape[0:3])
    mean = numpy.float32(mean)


    for sf in all_sub_files:
        print "computing mean, at {0}".format(sf)
        sample_4d = nib.load(os.path.join(label_map_dir, sf))
        list_3d = nib.funcs.four_to_three(sample_4d)
        mean = mean + list_3d[comp-1].get_data()
        
    mean = mean / len(all_sub_files)

    var_data = numpy.zeros(firstimage.shape[0:3])
    for sf in all_sub_files:
        print "computing var, at {0}".format(sf)
        sample_4d = nib.load(os.path.join(label_map_dir, sf))
        list_3d = nib.funcs.four_to_three(sample_4d)
        var_data = var_data + (list_3d[comp-1].get_data() - mean)**2

    var_data = var_data / len(all_sub_files)
    
    outimg = nib.Nifti1Image(var_data, numpy.eye(4))
    outimg.to_filename(out_filename)        
    



if __name__ == '__main__':
    compvar(sys.argv[1], sys.argv[2], sys.argv[3])
