#!/usr/bin/bash

fmri_dir=/usr/sci/scratch/weiliu/ADHD200_mypreprocess/Pittsburgh_stage2/smoothed/
mask_file=/usr/sci/scratch/weiliu/ADHD200_mypreprocess/Pittsburgh_stage2/mask/avg_4mm_thr0.3_bin.nii.gz
out_dir=/usr/sci/scratch/weiliu/ADHD200_mypreprocess/Pittsburgh_stage2/ncuts/
for fmri_file in `ls  $fmri_dir`; do
    fmri_file_base=`basename $fmri_file .gz`
    fmri_file_full=$fmri_dir$fmri_file
    out_file=$out_dir$fmri_file_base
    matlab -nodesktop  -nosplash -r "addpath ~/projects/ncuts_fmri; segts('$fmri_file_full', '$mask_file', 20, 0.4, '$out_file'); quit" 
done
