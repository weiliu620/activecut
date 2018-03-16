#!/usr/bin/bash

# This script use NCuts to segment the NYU dataset. First each subject is
# segmented, then the group segmentation is obtained again by NCuts.

# normalized fmri file. Un-normalized should also be fine since the segts
# routine will normalize all time series to zero mean and 1 variance. However,
# un-mormalized fmri data are scattered in subject folder, which make this
# script hard to process. So we just use normalized fmri, which is in a single
# folder.
fmri_dir=$1

# mask file should be from the Yeo's paper.
mask_file=$2

# output directory. Need to make the directory before running the script.
out_dir=$3

NCLUSTERS=$4

for fmri_file in `ls  $fmri_dir`
do
    out_file=$out_dir/${fmri_file}
    subname=`basename ${fmri_file} .nii.gz`
    aff_name=${out_dir}/${subname}.mat # affinity matrix file.
    matlab -nodesktop  -nosplash -r "addpath ~/projects/ncuts_fmri; dbstop if error; nyu_segts('${fmri_dir}/${fmri_file}', '${mask_file}', ${NCLUSTERS}, 0.4, '${out_file}', '${aff_name}'); quit" 
done

# Give a file name for group segmentation.
grpsegfile=${out_dir}/grp.nii.gz
matlab -nodesktop  -nosplash -r "addpath ~/projects/ncuts_fmri; dbstop if error; nyu_groupcuts('$out_dir', 0.4, '${mask_file}', ${NCLUSTERS}, '${grpsegfile}'); quit" 



