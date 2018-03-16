#!/usr/bin/bash
NUMSUBS=16

# create mask.
/home/sci/weiliu/projects/group_mrf/createmask -X 64 -Y 64 -Z 16 -R 64 -r 0 --mask /home/sci/weiliu/projects/group_mrf/syn_results/mask.nii.gz -v 1

# generate group label map.
/home/sci/weiliu/projects/group_mrf/creategrpmap --beta 2 --numclusters 4 --scan 500 --mask /home/sci/weiliu/projects/group_mrf/syn_results/mask.nii.gz -v 1 --true ~/projects/group_mrf/syn_results/truegrpmap.nii.gz

# Given true group map, generate individual sub map.
for i in $(seq -w 1 $NUMSUBS); do
    /home/sci/weiliu/projects/group_mrf/gensubmap -a 0.7 -b 2 -s 300 -g ~/projects/group_mrf/syn_results/truegrpmap.nii.gz -o /home/sci/weiliu/projects/group_mrf/syn_results/sub${i}.nii.gz --seed ${i} &
done

wait 

# extract mean time series, given label map and fmri.
for i in $(seq -w 1 $NUMSUBS); do
    /home/sci/weiliu/projects/group_mrf/xtrctmeants -k 8 -f /home/sci/weiliu/dataset/allfmri/niftiDATA_Subject0${i}_Condition001_sphere.nii --labelmap /home/sci/weiliu/projects/group_mrf/init/best_8_clusters.nii -o  ~/dataset/mean_time_series/sub${i}.txt &
done

# wait

# generate fmri image based on label map and mean time series.
for i in $(seq -w 1 $NUMSUBS); do
    /home/sci/weiliu/projects/group_mrf/createfmri -k 25 --obs ~/dataset/synthetic_fmri/sub${i}.nii.gz -m ~/dataset/mean_time_series/sub${i}.txt.xtrd -p 238 -v 1 -t /home/sci/weiliu/projects/group_mrf/syn_results/sub${i}.nii.gz &
done

##############################################################
######## End of generating fmri data #########################
##############################################################

# project fmri data on to sphere.
for i in $(seq -w 1 $NUMSUBS); do
/home/sci/weiliu/projects/group_mrf/projectfmri --in ~/dataset/synthetic_fmri/sub${i}.nii.gz --out ~/dataset/synthetic_fmri_onsphere/sub${i}.nii.gz &
done

# Compute initial labelmap by K-means.
for i in $(seq -w 1 10); do
/home/sci/weiliu/projects/group_mrf/initmrf --fmripath  ~/dataset/synthetic_fmri_onsphere/ -m /home/sci/weiliu/projects/group_mrf/syn_results/mask.nii.gz --kmeansiter 20 --numClusters 4 --seed ${i} --outgrouplabel /home/sci/weiliu/projects/group_mrf/init/grouplabel${i}.nii.gz -v 1 > /home/sci/weiliu/projects/group_mrf/init/log${i} &
done


# run groupmrf on synthetic data (given initial group label map as input)

# Test 1. No beta term. small mc sample size.
/home/sci/weiliu/projects/group_mrf/groupmrf --seed 0 -v 1 -b 20 --inittemp 1 --finaltemp 0.2 -n 5 --emiter 30 --betag 0.1 --betaz 0.1 --alpha 0.7 -k 4  -f ~/dataset/synthetic_fmri_onsphere -g ~/dataset/groupmrf_output/groupmap1.nii.gz -i /home/sci/weiliu/projects/group_mrf/init/best_4_clusters_syn.nii.gz --sampleprefix ~/dataset/groupmrf_output/test1_sample_sub --subbasename ~/dataset/groupmrf_output/test1_submap > ~/dataset/groupmrf_output/log01.txt &

# Test 2. with beta and alpha term. mc sample 30.
/home/sci/weiliu/projects/group_mrf/groupmrf --seed 0 -v 1 -b 30 --inittemp 1 --finaltemp 0.2 -n 20 --emiter 30 --betag 2 --betaz 2 --alpha 0.7 -k 4  -f ~/dataset/synthetic_fmri_onsphere -g ~/dataset/groupmrf_output/groupmap2.nii.gz -i /home/sci/weiliu/projects/group_mrf/init/best_4_clusters_syn.nii.gz --sampleprefix ~/dataset/groupmrf_output/test2_sample_sub --subbasename ~/dataset/groupmrf_output/test2_submap > ~/dataset/groupmrf_output/log02.txt &

# Test 3. Same with test 2 but small MC sample size. To see if small sample size have any effect.
/home/sci/weiliu/projects/group_mrf/groupmrf --seed 0 -v 1 -b 30 --inittemp 1 --finaltemp 0.2 -n 1 --emiter 30 --betag 2 --betaz 2 --alpha 0.7 -k 4  -f ~/dataset/synthetic_fmri_onsphere -g ~/dataset/groupmrf_output/groupmap3.nii.gz -i /home/sci/weiliu/projects/group_mrf/init/best_4_clusters_syn.nii.gz --sampleprefix ~/dataset/groupmrf_output/test3_sample_sub --subbasename ~/dataset/groupmrf_output/test3_submap > ~/dataset/groupmrf_output/log03.txt &


# Test 4. Same with Test 2.
/home/sci/weiliu/projects/group_mrf/groupmrf --seed 0 -v 1 -b 30 --inittemp 1 --finaltemp 0.2 -n 20 --emiter 30 --betag 2 --betaz 2 --alpha 0.7 -k 4  -f ~/dataset/synthetic_fmri_onsphere -g ~/dataset/groupmrf_output/groupmap4.nii.gz -i /home/sci/weiliu/projects/group_mrf/init/best_4_clusters_syn.nii.gz --sampleprefix ~/dataset/groupmrf_output/test4_sample_sub --subbasename ~/dataset/groupmrf_output/test4_sumap > ~/dataset/groupmrf_output/log04.txt &

# Test 5. no group beta, but has sub level beta and alpha.
/home/sci/weiliu/projects/group_mrf/groupmrf --seed 0 -v 1 -b 20 --inittemp 1 --finaltemp 0.2 -n 5 --emiter 30 --betag 0.1 --betaz 2 --alpha 0.7 -k 4  -f ~/dataset/synthetic_fmri_onsphere -g ~/dataset/groupmrf_output/groupmap5.nii.gz -i /home/sci/weiliu/projects/group_mrf/init/best_4_clusters_syn.nii.gz --sampleprefix ~/dataset/groupmrf_output/test5_sample_sub --subbasename ~/dataset/groupmrf_output/test5_sumap > ~/dataset/groupmrf_output/log05.txt &

# Test 6. Same with Test 2, but init'd with true group label map.
/home/sci/weiliu/projects/group_mrf/groupmrf --seed 0 -v 1 -b 30 --inittemp 1 --finaltemp 0.2 -n 20 --emiter 30 --betag 2 --betaz 2 --alpha 0.7 -k 4  -f ~/dataset/synthetic_fmri_onsphere -g ~/dataset/groupmrf_output/groupmap6.nii.gz -i /home/sci/weiliu/projects/group_mrf/syn_results/truegrpmap.nii.gz --sampleprefix ~/dataset/groupmrf_output/test6_sample_sub --subbasename ~/dataset/groupmrf_output/test6_submap > ~/dataset/groupmrf_output/log06.txt &

# Test 7. Same with Test 1 but init with true grp label map.
/home/sci/weiliu/projects/group_mrf/groupmrf --seed 0 -v 1 -b 20 --inittemp 1 --finaltemp 0.2 -n 5 --emiter 30 --betag 0.1 --betaz 0.1 --alpha 0.7 -k 4  -f ~/dataset/synthetic_fmri_onsphere -g ~/dataset/groupmrf_output/groupmap7.nii.gz -i /home/sci/weiliu/projects/group_mrf/syn_results/truegrpmap.nii.gz --sampleprefix ~/dataset/groupmrf_output/test7_sample_sub --subbasename ~/dataset/groupmrf_output/test7_submap > ~/dataset/groupmrf_output/log07.txt &

# Test 8. Same with Test 3 but init with grp label map obtained from test 1. This is like bootstrap method?
/home/sci/weiliu/projects/group_mrf/groupmrf --seed 0 -v 1 -b 20 --inittemp 1 --finaltemp 0.2 -n 1 --emiter 30 --betag 2 --betaz 2 --alpha 0.7 -k 4  -f ~/dataset/synthetic_fmri_onsphere -g ~/dataset/groupmrf_output/groupmap8.nii.gz -i /home/sci/weiliu/dataset/groupmrf_output/groupmap1.nii.gz --sampleprefix ~/dataset/groupmrf_output/test8_sample_sub --subbasename ~/dataset/groupmrf_output/test8_submap > ~/dataset/groupmrf_output/log08.txt &

# Test 1a. same with test1 but with smaller sample size. This is to see if smaller same size get better or worse resutls.
/home/sci/weiliu/projects/group_mrf/groupmrf --seed 0 -v 1 -b 20 --inittemp 1 --finaltemp 0.2 -n 1 --emiter 30 --betag 0.1 --betaz 0.1 --alpha 0.7 -k 4  -f ~/dataset/synthetic_fmri_onsphere -g ~/dataset/groupmrf_output/groupmap1a.nii.gz -i /home/sci/weiliu/projects/group_mrf/init/best_4_clusters_syn.nii.gz --sampleprefix ~/dataset/groupmrf_output/test1a_sample_sub --subbasename ~/dataset/groupmrf_output/test1a_submap > ~/dataset/groupmrf_output/log01a.txt &

# Test 1b. same with test1 but with bigger sample size. To see why the samples are abnormal.
/home/sci/weiliu/projects/group_mrf/groupmrf --seed 0 -v 1 -b 20 --inittemp 1 --finaltemp 0.2 -n 10 --emiter 30 --betag 0.1 --betaz 0.1 --alpha 0.7 -k 4  -f ~/dataset/synthetic_fmri_onsphere -g ~/dataset/groupmrf_output/groupmap1b.nii.gz -i /home/sci/weiliu/projects/group_mrf/init/best_4_clusters_syn.nii.gz --sampleprefix ~/dataset/groupmrf_output/test1b_sample_sub --subbasename ~/dataset/groupmrf_output/test1b_submap > ~/dataset/groupmrf_output/log01b.txt &

# Test 1c. same with test1 but with different temperature. for debugging.
/home/sci/weiliu/projects/group_mrf/groupmrf --seed 0 -v 1 -b 20 --inittemp 1 --finaltemp 0.1 -n 5 --emiter 30 --betag 0.1 --betaz 0.1 --alpha 0.7 -k 4  -f ~/dataset/synthetic_fmri_onsphere -g ~/dataset/groupmrf_output/groupmap1c.nii.gz -i /home/sci/weiliu/projects/group_mrf/init/best_4_clusters_syn.nii.gz --sampleprefix ~/dataset/groupmrf_output/test1c_sample_sub --subbasename ~/dataset/groupmrf_output/test1c_submap > ~/dataset/groupmrf_output/log01c.txt &





############### vmm test #######################

# After upgrading vmmseg with MRF.
/home/sci/weiliu/projects/group_mrf/vmmseg -b 20 -n 2 --emiter 30 --seed 0 --numClusters 6 --betag 0.1 -f ~/dataset/allfmri -i /home/sci/weiliu/dataset/allfmri_mask_and_init/6_clusters --grpsamples /home/sci/weiliu/dataset/allfmri_vmmseg_output/grpSample -g /home/sci/weiliu/dataset/allfmri_vmmseg_output/outGrpLabel.nii.gz -v 1

# test2: 8 clusters
/home/sci/weiliu/projects/group_mrf/vmmseg -b 20 -n 2 --emiter 30 --seed 0 --numClusters 8 --betag 0.1 -f ~/dataset/allfmri -i /home/sci/weiliu/dataset/allfmri_mask_and_init/8_clusters --grpsamples /home/sci/weiliu/dataset/allfmri_vmmseg_output/grpSample_test2.nii.gz -g /home/sci/weiliu/dataset/allfmri_vmmseg_output/outGrpLabel_test2.nii.gz -v 1   > /home/sci/weiliu/dataset/allfmri_vmmseg_output/test2.log &

# Test 3. Still 8 clusters but 50 iteration.
/home/sci/weiliu/projects/group_mrf/vmmseg -b 20 -n 2 --emiter 50 --seed 0 --numClusters 8 --betag 0.1 -f ~/dataset/allfmri -i /home/sci/weiliu/dataset/allfmri_mask_and_init/8_clusters --grpsamples /home/sci/weiliu/dataset/allfmri_vmmseg_output/grpSample_test3.nii.gz -g /home/sci/weiliu/dataset/allfmri_vmmseg_output/outGrpLabel_test3.nii.gz -v 1  > /home/sci/weiliu/dataset/allfmri_vmmseg_output/test3.log &

# Test 4: 10 clusters.
/home/sci/weiliu/projects/group_mrf/vmmseg -b 20 -n 2 --emiter 30 --seed 0 --numClusters 10 --betag 0.1 -f ~/dataset/allfmri -i /home/sci/weiliu/dataset/allfmri_mask_and_init/10_clusters --grpsamples /home/sci/weiliu/dataset/allfmri_vmmseg_output/grpSample_test4.nii.gz -g /home/sci/weiliu/dataset/allfmri_vmmseg_output/outGrpLabel_test4.nii.gz -v 1 > /home/sci/weiliu/dataset/allfmri_vmmseg_output/test4.log &

# Test 5: 10 clusters.
/home/sci/weiliu/projects/group_mrf/vmmseg -b 20 -n 2 --emiter 30 --seed 0 --numClusters 12 --betag 0.1 -f ~/dataset/allfmri -i /home/sci/weiliu/dataset/allfmri_mask_and_init/12_clusters --grpsamples /home/sci/weiliu/dataset/allfmri_vmmseg_output/grpSample_test5.nii.gz -g /home/sci/weiliu/dataset/allfmri_vmmseg_output/outGrpLabel_test5.nii.gz -v 1 > /home/sci/weiliu/dataset/allfmri_vmmseg_output/test5.log &

# Test 6: 25 clusters.
/home/sci/weiliu/projects/group_mrf/vmmseg -b 20 -n 2 --emiter 30 --seed 0 --numClusters 25 --betag 0.1 -f ~/dataset/allfmri -i /home/sci/weiliu/dataset/allfmri_mask_and_init/25_clusters --grpsamples /home/sci/weiliu/dataset/allfmri_vmmseg_output/grpSample_test6.nii.gz -g /home/sci/weiliu/dataset/allfmri_vmmseg_output/outGrpLabel_test6.nii.gz -v 1 > /home/sci/weiliu/dataset/allfmri_vmmseg_output/test6.log &

###################################################
# run vmmseg on unsmoothed data.

# test 7. 
/home/sci/weiliu/projects/group_mrf/vmmseg -b 20 -n 2 --emiter 30 --seed 0 --numClusters 8 --betag 0.1 -f ~/dataset/allfmri_unsmoothed -i /home/sci/weiliu/dataset/allfmri_mask_and_init/8_clusters --grpsamples /home/sci/weiliu/dataset/allfmri_unsmoothed_vmmseg_output/grpSample_test7.nii.gz -g /home/sci/weiliu/dataset/allfmri_unsmoothed_vmmseg_output/outGrpLabel_test7.nii.gz -v 1 > /home/sci/weiliu/dataset/allfmri_unsmoothed_vmmseg_output/test7.log &

# test 7a. Different with test7 only on iterations of EM.
/home/sci/weiliu/projects/group_mrf/vmmseg -b 20 -n 2 --emiter 50 --seed 0 --numClusters 8 --betag 0.1 -f ~/dataset/allfmri_unsmoothed -i /home/sci/weiliu/dataset/allfmri_mask_and_init/8_clusters --grpsamples /home/sci/weiliu/dataset/allfmri_unsmoothed_vmmseg_output/grpSample_test7a.nii.gz -g /home/sci/weiliu/dataset/allfmri_unsmoothed_vmmseg_output/outGrpLabel_test7a.nii.gz -v 1 > /home/sci/weiliu/dataset/allfmri_unsmoothed_vmmseg_output/test7a.log &


# test 8. 
/home/sci/weiliu/projects/group_mrf/vmmseg -b 20 -n 2 --emiter 30 --seed 0 --numClusters 8 --betag 1 --gamma 20 -f ~/dataset/allfmri_unsmoothed -i /home/sci/weiliu/dataset/allfmri_mask_and_init/8_clusters --grpsamples /home/sci/weiliu/dataset/allfmri_unsmoothed_vmmseg_output/grpSample_test8.nii.gz -g /home/sci/weiliu/dataset/allfmri_unsmoothed_vmmseg_output/outGrpLabel_test8.nii.gz -v 1 > /home/sci/weiliu/dataset/allfmri_unsmoothed_vmmseg_output/test8.log &
