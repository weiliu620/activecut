segts('/usr/sci/scratch/weiliu/conn_all_smooth/results/preprocessing/niftiDATA_Subject001_Condition001.nii',...
    '/home/sci/weiliu/dataset/R_data/mask_and_init/grey_3mm_thr0.4_clean.nii.gz',...
    20, 0.4, '/home/sci/weiliu/dataset/R_test/ncuts/sub1.nii');


M = segts('~/dataset/R_data/sub01_trunc.nii.gz',...
    '/home/sci/weiliu/dataset/R_data/mask_trunc.nii',...
    20, 0.6, '/home/sci/weiliu/dataset/R_test/ncuts/sub1_trunc.nii');

addpath ~/packages/nifti_to_matlab
mask_nii = load_untouch_nii('/home/sci/weiliu/dataset/R_data/grey_3mm_thr0.4_clean.nii');
mask_nii.img(:,:,1:20) = 0;
mask_nii.img(:,:,26:46) = 0;
save_untouch_nii(mask_nii, '/home/sci/weiliu/dataset/R_data/mask_trunc.nii');





