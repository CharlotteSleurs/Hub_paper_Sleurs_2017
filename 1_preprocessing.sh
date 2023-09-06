 #!/bin/bash

# Use:  ./preprocessing.sh C01

# convert niis (from dcm2nii) to .mif
mrconvert OriginalFiles/$1_DKI25d40d7001000.nii -fslgrad OriginalFiles/$1_DKI25d40d7001000.bvec OriginalFiles/$1_DKI25d40d7001000.bval MergedFiles2/$1_DKI1.mif
mrconvert OriginalFiles/$1_DKI75d2800.nii -fslgrad OriginalFiles/$1_DKI75d2800.bvec OriginalFiles/$1_DKI75d2800.bval MergedFiles2/$1_DKI2.mif
mrcat -axis 3 MergedFiles2/$1_DKI1.mif MergedFiles2/$1_DKI2.mif MergedFiles2/$1_DKItmp.mif
rm MergedFiles2/$1_DKI1.mif MergedFiles2/$1_DKI2.mif

# remove ugly B0 in acquired image
mrconvert -coord 3 0:73,75:end MergedFiles2/$1_DKItmp.mif MergedFiles2/$1_DKI.mif
rm MergedFiles2/$1_DKItmp.mif

# create APPA (combi of anterior-posterior and reversed phase B0 image)
mrconvert -coord 3 0 OriginalFiles/$1_DKI25d40d7001000.nii - | mrcat -axis 3 - OriginalFiles/$1_DKIrev6d.nii MergedFiles2/$1_APPAtmp.mif
mrpad MergedFiles2/$1_APPAtmp.mif MergedFiles2/$1_APPA.mif -axis 1 6 6 # pad your APPAimage in case that subjects' images reach the borders of the FOV
rm MergedFiles2/$1_APPAtmp.mif

# denoising
dwidenoise MergedFiles2/$1_DKI.mif MergedFiles2/$1_DKI_denoised.mif

# gibbsringing
mrdegibbs MergedFiles2/$1_DKI_denoised.mif MergedFiles2/$1_DKI_denoised_gibbstmp.mif
mrpad MergedFiles2/$1_DKI_denoised_gibbstmp.mif MergedFiles2/$1_DKI_denoised_gibbs.mif -axis 1 6 6 # also pad your dwi-image for topup in next step, to have similar sizes of these images
rm MergedFiles2/$1_DKI_denoised_gibbstmp.mif

# dwipreproc
dwipreproc MergedFiles2/$1_DKI_denoised_gibbs.mif MergedFiles2/$1_DKI_denoised_gibbs_preproc.mif -pe_dir AP -rpe_pair -se_epi MergedFiles2/$1_APPA.mif -eddy_options "--data_is_shelled --repol "
#mrcrop -mask MergedFiles2/$1_DKI.mif MergedFiles2/$1_DKI_denoised_gibbs_preproc.mif MergedFiles2/$1_DKI_denoised_gibbs_preproc_origdims.mif # crop to original size again (after the padding initially)

# dwibiascorrect
mask_dir=MergedFiles2/Masks
mkdir $mask_dir
tonorm_dir=MergedFiles2/images2normalise
mkdir $tonorm_dir
norm_dir=MergedFiles2/normalised

mask=$mask_dir'/'$1_dwi_bin.mif
dwi2mask MergedFiles2/$1_DKI_denoised_gibbs_preproc.mif $mask
dwibiascorrect -ants -mask $mask MergedFiles2/$1_DKI_denoised_gibbs_preproc.mif $tonorm_dir'/'$1_DKI_denoised_gibbs_preproc_bias.mif

# TEST dwipreprocparallel
# GNU_par_nodes=my_cluster_uzl.txt
# input=parallel_input_pp.txt
# f_in  = MergedFiles2/$1_DKI_denoised_gibbs.mif; % Full path input *.mif file name.
# f_out = MergedFiles2/$1_DKI_denoised_gibbs_preproc.mif; % Full path output file name.
# f_b0s = MergedFiles2/$1_APPA.mif
# f_option = -rpe_pair -pe_dir AP -eddy_options "--repol " -se_epi ' f_b0s ' '];
# fprintf(input,'%s\t%s\t%s\t%s\n',f_in,f_b0s,f_out,f_option);
# parallel --will-cite --colsep ''\t'' --joblog ' Dir_Main_ex '/' Dir_log '/parallel_pp.log --nice 19 --cleanup --transferfile {1} --transferfile {2} --return {3} eval
#    MRtrix3 '/' 'dwipreproc {1} {3} {4} :::: input_for_parallel_pp];


