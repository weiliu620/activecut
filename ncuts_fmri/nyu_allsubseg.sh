#!/usr/bin/bash

# This script use NCuts to segment the NYU dataset's all subject fMRI data.

# normalized fmri file. Un-normalized should also be fine since the segts
# routine will normalize all time series to zero mean and 1 variance. However,
# un-mormalized fmri data are scattered in subject folder, which make this
# script hard to process. So we just use normalized fmri, which is in a single
# folder.
fmri_dir=$1

# mask file should be from the Yeo's paper. as Jan 19, 2013. Both mask file and
# fmri file must be nii format rather than nii.gz file. The script will unzip
# the gz file and get the temporary nii file.
mask_file=$2

# output directory. Need to make the directory before running the script.
out_dir=$3


# define macros.
mask_file_nii='tmp_mask.nii'

if [ ${mask_file: -3} == ".gz" ]
then
    gunzip -c ${mask_file} > ${mask_file_nii}
elif [ ${mask_file: -4} == ".nii" ]
then
    mask_file_nii=${mask_file}
else
    echo "input file must be nii or gz format."
fi


for fmri_file in `ls  $fmri_dir`
do
    if [ ${fmri_file: -3} == ".gz" ]
    then
    	fmri_file_nii=`basename ${fmri_file} .gz`
    	gunzip -c ${fmri_dir}/${fmri_file} > ${fmri_file_nii}
    elif [ ${fmri_file: -4} == ".nii" ]
    then
    	fmri_file_nii=${fmri_file}
    else
    echo "input file must be nii or gz format."
    fi
    out_file=$out_dir/${fmri_file_nii}

    matlab -nodesktop  -nosplash -r "nyu_segts('${fmri_file_nii}', '${mask_file_nii}', 7, 0.4, '${out_file}'); quit" 
    rm ${fmri_file_nii}
done

rm ${mask_file_nii}

