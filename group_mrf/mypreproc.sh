#!/bin/bash

# Usage: myscript subjhome subjout verboseLevel

# filename convention: 
# t -- truncated
# r -- motion corrected.
# s -- slice timing corrected.
# w -- warpped to MNI space.

subjhome=$1
subjout=$2
verbose=$3
cd ${subjout}

# pickout series 12 and merge into a single nii file.
fslmerge -t fmri ${subjhome}/Resting/f*-001[2,6]-00???-000???-01.hdr

# retrieve T1 and T2 image...\n"
fslchfiletype NIFTI_GZ ${subjhome}/MPRAGE/s*.hdr t1
fslchfiletype NIFTI_GZ ${subjhome}/T2/s*.hdr t2


# remove first 5 time points.
nvols=`fslnvols fmri` # number of time points
echo  -e "Now remove first five time points...\n"
fslroi fmri tfmri 5 $nvols-5

nvols=`fslnvols tfmri` # number of time points
ref_vol_id=`expr $nvols "/" 2` # mid volume number. assume nvols is even.
fslroi tfmri mid_vol $ref_vol_id 1  # extract mid volume.

# Motion correction. Save affine registration parameters Tx6 matrix for
# regression use in the future. The parameter file has same name with output
# file.
mcflirt -in tfmri -o rtfmri -refvol ${ref_vol_id} -plots

# Slice timing correction.
slicetimer -i rtfmri -o srtfmri  -r 2 --ocustom=/home/sci/weiliu/bin/slice_timing_order.txt

# Remove scalp
bet  t1 t1_brain -R -f 0.1
bet t2 t2_brain -R -f 0.1 
bet mid_vol mid_vol_brain -R -f 0.1

# Register scalp removed T2 to T1. Also apply same transformation to normal T2 so it aligns to T1.
flirt -in t2_brain -ref t1_brain -out rt2_brain -omat t21.mat -dof 6
flirt -in t2 -ref t1 -out rt2 -applyxfm -init t21.mat

# register fmri's middle volume to T2 (also aligned to T1 space). Only use it's
# xfm matrix in future.
flirt -in mid_vol_brain -ref rt2_brain  -omat mid_rt2.mat -dof 6

# Linear Register T1 to MNI152 standard image (as initial value for warping).
flirt -in t1_brain -ref /usr/local/fsl/data/standard/MNI152_T1_2mm_brain  -omat t1_MNI.mat -dof 12

# Compute warping from T1 to MNI space.
fnirt --in=t1 --ref=/usr/local/fsl/data/standard/MNI152_T1_2mm --iout=wt1 --aff=t1_MNI.mat --cout=t1_MNI_warp.mat --config=T1_2_MNI152_2mm

# also transform T2 to MNI space.
applywarp --in=rt2 --ref=/usr/local/fsl/data/standard/MNI152_T1_2mm_brain --out=wrt2 --warp=t1_MNI_warp.mat


# apply (fmri->T1 -> MNI). Refer to:
# http://www.fmrib.ox.ac.uk/fsl/fnirt/index.html#fnirt_examples. premat is from
# fmri to T1 registration.
applywarp --ref=/usr/local/fsl/data/standard/MNI152_T1_2mm_brain --in=srtfmri --warp=t1_MNI_warp.mat --premat=mid_rt2.mat --out=wsrtfmri
fslmaths wsrtfmri -subsamp2 kwsrtfmri

# regression on white matter
/home/sci/weiliu/projects/group_mrf/regfilter -m ~/dataset/templates/MNI152_T1_4mm_brain_mask.nii.gz -r ~/dataset/templates/avg152T1_white_4mm_mask100.nii.gz -i kwsrtfmri.nii.gz -o kwsrtfmri_reg1.nii.gz -b beta_w.nii.gz -d yes 

# regression on CSF
/home/sci/weiliu/projects/group_mrf/regfilter -m ~/dataset/templates/MNI152_T1_4mm_brain_mask.nii.gz -r ~/dataset/templates/avg152T1_csf_4mm_mask100.nii.gz  -i kwsrtfmri_reg1.nii.gz -o kwsrtfmri_reg2.nii.gz -b beta_c.nii.gz -d yes 

# regression on motion parameters
/home/sci/weiliu/projects/group_mrf/regmotion -m ~/dataset/templates/MNI152_T1_4mm_brain_mask.nii.gz -i kwsrtfmri_reg2.nii.gz -p rtfmri.par -o kwsrtfmri_reg.nii.gz -b beta_m.nii.gz 

# Spatial smoothing.
fslmaths kwsrtfmri_reg  -s 3 skwsrtfmri_reg

# Band pass filtering on smoothed.
3dFourier -retrend -lowpass 0.1 -highpass 0.01 -prefix fskwsrtfmri_reg.nii.gz skwsrtfmri_reg.nii.gz

# Band pass filtering on non-smoothed 
3dFourier -retrend -lowpass 0.1 -highpass 0.01 -prefix fkwsrtfmri_reg.nii.gz kwsrtfmri_reg.nii.gz



if [ "$verbose" == 1 ]
then
    echo "verbose = 1. keep temp file."
else
    rm fmri.nii.gz
    rm mid_vol.nii.gz
    rm tfmri.nii.gz
    rm rtfmri.nii.gz
    rm srtfmri.nii.gz
    rm t1_brain.nii.gz
    rm t2_brain.nii.gz
    rm mid_vol_brain.nii.gz
    rm t1_MNI_warp.mat.nii.gz
    rm rt2.nii.gz
    rm rt2_brain.nii.gz
    rm wsrtfmri.nii.gz
    rm kwsrtfmri.nii.gz
    rm kwsrtfmri_reg1.nii.gz
    rm kwsrtfmri_reg2.nii.gz
    rm kwsrtfmri_reg.nii.gz
    rm skwsrtfmri_reg.nii.gz
    rm beta_w.nii.gz
    rm beta_c.nii.gz
    rm beta_m.nii.gz
fi
