##!/bin/bash

# Use:  ./ACT_tracto.sh, only if all DWI preprocessing steps and T1 segmentations are performed! 
# Define directories

#./ACT_tracto.sh C01; ./ACT_tracto.sh C02; ./ACT_tracto.sh C03; ./ACT_tracto.sh C04; ./ACT_tracto.sh C05; ./ACT_tracto.sh C06; ./ACT_tracto.sh C07;
#./ACT_tracto.sh C08; ./ACT_tracto.sh C09; ./ACT_tracto.sh C10; ./ACT_tracto.sh C11; ./ACT_tracto.sh C12; ./ACT_tracto.sh C13; ./ACT_tracto.sh C14;
#./ACT_tracto.sh C15; ./ACT_tracto.sh C16; ./ACT_tracto.sh C17; ./ACT_tracto.sh C18; ./ACT_tracto.sh C19; ./ACT_tracto.sh C20; ./ACT_tracto.sh C21;
#./ACT_tracto.sh D01; ./ACT_tracto.sh D02; ./ACT_tracto.sh D03; ./ACT_tracto.sh D04; ./ACT_tracto.sh D05; ./ACT_tracto.sh D06; ./ACT_tracto.sh D07;
#./ACT_tracto.sh D08; ./ACT_tracto.sh D09; ./ACT_tracto.sh D10; ./ACT_tracto.sh D11; ./ACT_tracto.sh D12; ./ACT_tracto.sh D13; ./ACT_tracto.sh D14;
#./ACT_tracto.sh D15; ./ACT_tracto.sh D16; ./ACT_tracto.sh D17; ./ACT_tracto.sh D18; ./ACT_tracto.sh D19; ./ACT_tracto.sh D20; ./ACT_tracto.sh D21;

# T1
T1_betted=/Volumes/DriveChar/Fossa_Posteriors/T1/WIP3DTFE/standard_oriented/CAT_segmentations/segmentation_native/AAL_native_registrations'/'$1_T1_brain.nii
T1_wm_segm=/Volumes/DriveChar/Fossa_Posteriors/T1/WIP3DTFE/standard_oriented/CAT_segmentations/segmentation_native/mri/p2$1_WIP3DTFE_reor2std.nii
T1_betted_mask=/Volumes/DriveChar/Fossa_Posteriors/T1/WIP3DTFE/standard_oriented/CAT_segmentations/segmentation_native/AAL_native_registrations'/'$1_T1_brain_bin04pixdim.nii
fivetissueimage=/Volumes/DriveChar/Fossa_Posteriors/T1/WIP3DTFE/standard_oriented/CAT_5ttgen_segmentations/$1_CAT_5ttgen.nii
AALsubject=/Volumes/DriveChar/Fossa_Posteriors/T1/WIP3DTFE/standard_oriented/CAT_segmentations/segmentation_native/AAL_native_registrations/$1_AAL2native_labels_fnirt.nii

#gunzip /Users/csleur0/Google_drive/Fossa_Posteriors/T1/WIP3DTFE/standard_oriented/CAT_segmentations/segmentation_native/AAL_native_registrations/*.gz

# DWImasks
mask_dir=MergedFiles2/Masks
maskdil_dir=MergedFiles2/Masks_dil
#mkdir $maskdil_dir
# DWI_preprocessed - denoised, gibbs, dwipreproc, dwibiascorrected - images
tonorm_dir=MergedFiles2/images2normalise
# DWI_normalised
norm_dir=MergedFiles2/normalised
mask=$mask_dir'/'$1_dwi_bin.mif
# FODdir & tckdir
FODdir=MergedFiles2/normalised/FODs
FODaverageResponse=MergedFiles2/normalised/FODs_averageResp; mkdir $FODaverageResponse
#mkdir $FODdir
tckdir=MergedFiles2/normalised/tcks
tckwtsAverageResponse=MergedFiles2/normalised/tckwtsAverageResponse; mkdir $tckwtsAverageResponse
#mkdir $tckdir

# T1_DWI_space
T1_DWI_space=MergedFiles2/anat2DWI
#mkdir $T1_DWI_space

###############################################################################

# dwiintensitynorm (normalise dwi intensities across group)
# dwiintensitynorm $tonorm_dir $mask_dir $norm_dir fa_template.mif wm_mask.mif

    #(( Register T1-betted image rigidly to DWI B0 (to average DWI performs better, see line 44)
    #mrconvert -coord 3 0 $norm_dir'/'$2_DKI_denoised_gibbs_preproc_bias.mif $norm_dir'/'$2_B0.nii
    #flirt -in $norm_dir'/'$2_B0.nii -ref $T1_betted -dof 6 -omat tmp.mat
    #flirt -in $norm_dir'/'$2_B0.nii -ref $T1_betted -dof 6 -cost bbr -wmseg T1_wm_segm -init tmp.nii.gz -omat diff_to_structural-bbr.mat -schedule $FSLDIR/etc/flirtsch/bbr.sch
    #transformconvert diff_to_structural-bbr.mat $norm_dir'/'$2_B0.nii $T1_betted flirt_import diff_to_structural-bbr-mrtrixformat.txt
    #mrtransform $T1_betted -linear diff_to_structural-bbr-mrtrixformat.txt $T1_DWI_space'/'$1_anat2dwi.mif -inverse))

# Registration to average DWI
#dwiextract -no_bzero $norm_dir'/'$1_DKI_denoised_gibbs_preproc_bias.mif $norm_dir'/'$1_noB0.nii
#mrmath -axis 3 $norm_dir'/'$1_noB0.nii mean $norm_dir'/'$1_noB0_average.nii
#flirt -in $norm_dir'/'$1_noB0_average.nii -ref $T1_betted -dof 6 -omat $1_tmp.mat
#flirt -in $norm_dir'/'$1_noB0_average.nii -ref $T1_betted -dof 6 -cost bbr -wmseg $T1_wm_segm -init $1_tmp.mat -omat $1_diff_to_structural-bbr.mat -schedule $FSLDIR/etc/flirtsch/bbr.sch
#transformconvert $1_diff_to_structural-bbr.mat $norm_dir'/'$1_noB0_average.nii $T1_betted flirt_import $1_diff_to_structural-bbr-mrtrixformat.txt
##mrtransform $T1_betted -linear $1_diff_to_structural-bbr-mrtrixformat.txt $T1_DWI_space'/'$1_anat2dwi.mif -inverse -force
###mrtransform $fivetissueimage -linear $1_diff_to_structural-bbr-mrtrixformat.txt $T1_DWI_space'/'$1_5tt2dwi.mif -inverse -force
#mrtransform $AALsubject -linear $1_diff_to_structural-bbr-mrtrixformat.txt -interp nearest $T1_DWI_space'/'$1_AAL2dwi.mif -inverse -force

###
#rm $1_diff_to_structural-bbr-mrtrixformat.mat; rm $1_tmp.mat; rm $norm_dir'/'$1_noB0.nii;
#rm $norm_dir'/'$1_noB0_average.nii

# dilate mask-image for ACT tractography
#maskfilter $mask dilate $maskdil_dir'/'$1_dwi_bin_dilate.mif -force

# Generate tracks using DWI & T1
#dwi2response msmt_5tt $norm_dir'/'$1_DKI_denoised_gibbs_preproc_bias.mif $T1_DWI_space'/'$1_5tt2dwi.mif $FODdir/$1_wm $FODdir/$1_gm $FODdir/$1_csf -force
#dwi2fod -mask $maskdil_dir/$1_dwi_bin_dilate.mif msmt_csd $norm_dir'/'$1_DKI_denoised_gibbs_preproc_bias.mif $FODdir/$1_wm $FODdir/$1_wm_odf.mif $FODdir/$1_gm $FODdir/$1_gm_odf.mif $FODdir/$1_csf $FODdir/$1_csf_odf.mif -force
#
#rm -rf *-tmp*

#tckgen
#5tt2gmwmi $T1_DWI_space'/'$1_5tt2dwi.mif $T1_DWI_space'/'$1_gmwmi.mif -force
#tckgen -seed_gmwmi $T1_DWI_space'/'$1_gmwmi.mif -select 10000000 -act $T1_DWI_space'/'$1_5tt2dwi.mif $FODdir/$1_wm_odf.mif $tckdir/$1_tracks.tck -cutoff 0.05 -backtrack
##sift1 if you want to check visually
#tcksift -output_at_counts 1000 $tckdir/$1_tracks.tck $FODdir/$1_wm_odf.mif $tckdir/$1_tracks1000.tck

## Connectome constructions
##sift2 for trackweights for connectome construction
#tcksift2 -act $T1_DWI_space'/'$1_5tt2dwi.mif $tckdir/$1_tracks.tck $FODdir/$1_wm_odf.mif $tckdir/$1_tckweights
#tck2connectome $tckdir/$1_tracks.tck $T1_DWI_space'/'$1_AAL2dwi.mif $tckdir/$1_connectome.csv -assignment_reverse_search 10 -zero_diagonal -symmetric -stat_edge sum -tck_weights_in $tckdir/$1_tckweights
## Use FA for connectomes
#echo "Tensor calculation Subject" $1
#dwi2tensor $norm_dir'/'$1_DKI_denoised_gibbs_preproc_bias.mif $FODdir/$1_tensor.mif -mask $maskdil_dir'/'$1_dwi_bin_dilate.mif -force
#tensor2metric $FODdir/$1_tensor.mif -fa $FODdir/$1_fa.mif -adc $FODdir/$1_md.mif -mask $maskdil_dir'/'$1_dwi_bin_dilate.mif -force
#tcksample -stat_tck median $tckdir/$1_tracks.tck $FODdir/$1_fa.mif $FODdir/$1_medianFA_streamlines.txt -force
#tck2connectome $tckdir/$1_tracks.tck $T1_DWI_space'/'$1_AAL2dwi.mif $tckdir/$1_famedianconnectome.csv -assignment_reverse_search 10 -zero_diagonal -symmetric -stat_edge mean -scale_file $FODdir/$1_medianFA_streamlines.txt -force

#((IF USING AVERAGE response functions (necessary for fixelbased, discutable for connectome construction):
#average_response $FODdir/*_wm.txt $FODdir/average_wm.txt
#average_response $FODdir/*_gm.txt $FODdir/average_gm.txt
#average_response $FODdir/*_csf.txt $FODdir/average_csf.txt
dwi2fod -mask $maskdil_dir/$1_dwi_bin_dilate.mif msmt_csd $norm_dir'/'$1_DKI_denoised_gibbs_preproc_bias.mif $FODdir/average_wm $FODaverageResponse/$1_wm_odf.mif $FODdir/average_gm $FODaverageResponse/$1_gm_odf.mif $FODdir/average_csf $FODaverageResponse/$1_csf_odf.mif -force
mtnormalise $FODaverageResponse/$1_wm_odf.mif $FODaverageResponse/$1_wmfod_norm.mif $FODaverageResponse/$1_gm_odf.mif $FODaverageResponse/gmfod_norm.mif $FODaverageResponse/$1_csf_odf.mif $FODaverageResponse/csffod_norm.mif -mask $maskdil_dir/$1_dwi_bin_dilate.mif -force
tcksift2 -act $T1_DWI_space'/'$1_5tt2dwi.mif $tckdir/$1_tracks.tck $FODaverageResponse/$1_wmfod_norm.mif $tckwtsAverageResponse/$1_tckweights_AR -out_mu $tckwtsAverageResponse/$1_tckweights_AR_mu -force
tck2connectome $tckdir/$1_tracks.tck $T1_DWI_space'/'$1_AAL2dwi.mif $tckwtsAverageResponse/$1_sift_connectome_AR.csv -assignment_reverse_search 10 -zero_diagonal -symmetric -stat_edge sum -tck_weights_in $tckwtsAverageResponse/$1_tckweights_AR -force

## Use NBS for statistics
# NBS algorithm
connectomestats nbs_files.txt nbse designMatrix.txt contrastPatLess.txt nbs_
# NBSE algorithm
connectomestats -threshold 1.5 nbs_files.txt nbs designMatrix.txt contrastPatLess.txt nbse_
