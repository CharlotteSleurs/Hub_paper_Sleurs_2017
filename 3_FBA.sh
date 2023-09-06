#!/bin/sh

#  3_FBA.sh
#  
#
#  Created by Charlotte Sleurs on 4/12/17.
#

#./3_FBA.sh C01; ./3_FBA.sh C02; ./3_FBA.sh C03; ./3_FBA.sh C04; ./3_FBA.sh C05; ./3_FBA.sh C06;
#./3_FBA.sh C07;
#./3_FBA.sh C08; ./3_FBA.sh C09; ./3_FBA.sh C10; ./3_FBA.sh C11; ./3_FBA.sh C12; ./3_FBA.sh C13;
#./3_FBA.sh C14;
#./3_FBA.sh C15; ./3_FBA.sh C16; ./3_FBA.sh C17; ./3_FBA.sh C18; ./3_FBA.sh C19; ./3_FBA.sh C20;
#./3_FBA.sh C21; ./3_FBA.sh D01; ./3_FBA.sh D02; ./3_FBA.sh D03; ./3_FBA.sh D04; ./3_FBA.sh D05;
#./3_FBA.sh D06; ./3_FBA.sh D07;
#./3_FBA.sh D08; ./3_FBA.sh D09; ./3_FBA.sh D10; ./3_FBA.sh D11; ./3_FBA.sh D12; ./3_FBA.sh D13;
#./3_FBA.sh D14;
#./3_FBA.sh D15; ./3_FBA.sh D16; ./3_FBA.sh D17; ./3_FBA.sh D18; ./3_FBA.sh D19; ./3_FBA.sh D20;
#./3_FBA.sh D21;


T1_betted=/Volumes/DriveChar/Fossa_Posteriors/T1/WIP3DTFE/standard_oriented/CAT_segmentations/segmentation_native/AAL_native_registrations'/'$1_T1_brain.nii
T1_wm_segm=/Volumes/DriveChar/Fossa_Posteriors/T1/WIP3DTFE/standard_oriented/CAT_segmentations/segmentation_native/mri/p2$1_WIP3DTFE_reor2std.nii
T1_betted_mask=/Volumes/DriveChar/Fossa_Posteriors/T1/WIP3DTFE/standard_oriented/CAT_segmentations/segmentation_native/AAL_native_registrations'/'$1_T1_brain_bin04pixdim.nii

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
FODaverageResponse=MergedFiles2/normalised/FODs_averageResp
#mkdir $FODdir
tckdir=MergedFiles2/normalised/tcks
tckwtsAverageResponse=MergedFiles2/normalised/tckwtsAverageResponse;
#mkdir $tckwtsAverageResponse
#mkdir $tckdir
# T1_DWI_space
T1_DWI_space=MergedFiles2/anat2DWI
#mkdir $T1_DWI_space
# FODtemplate
FODtemp=MergedFiles2/normalised/FODs_averageResp/FOD_template
# ADC-maps & FA
ADC=/Volumes/DriveChar/Fossa_Posteriors/DWI/MergedFiles2/normalised/tensor
ADCtemp=/Volumes/DriveChar/Fossa_Posteriors/DWI/MergedFiles2/normalised/tensor/subj2template
#mkdir $ADCtemp
# FBA folder
FBA=/Volumes/DriveChar/Fossa_Posteriors/DWI/MergedFiles2/normalised/FBA
fbastats=/Volumes/DriveChar/Fossa_Posteriors/DWI/MergedFiles2/normalised/FBA/stats
#mkdir $FBA;
mkdir $fbastats

################################################################
# Assuming that you followed preprocessing dwidenoise, mrdegibbbs, dwipreproc, dwibiascorrect,
# dwiintensitynorm, dwi2fod (averageresp), mtnormalise fods:
# Create FOD template
#cp $FODaverageResponse/*wmfod_norm* $FODtemp/
#population_template $FODtemp/ -mask_dir $maskdil_dir $FODtemp/wm_fod_template.mif

# Register FODs to template
#mrregister $FODaverageResponse/$1_wmfod_norm.mif -mask1 $mask_dir/$1_dwi_bin.mif $FODtemp/wm_fod_template.mif -nl_warp $FODtemp/$1_subject2template_warp.mif $FODtemp/$1_template2subject_warp.mif
#
## ( If you also want to compare voxelparams eg FA/ADC, just apply
#mrtransform $ADC/$1_adc.mif -warp $FODtemp/$1_subject2template_warp.mif $ADCtemp/$1_adc_templspace.mif
#mrtransform $ADC/$1_fa.mif -warp $FODtemp/$1_subject2template_warp.mif $ADCtemp/$1_fa_templspace.mif  )

# Apply warps to masks
#mrtransform $mask_dir/$1_dwi_bin.mif -warp $FODtemp/$1_subject2template_warp.mif -interp nearest $ADCtemp/$1_dwi_bin_templspace.mif

# Compute the intersection of all warped masks
#mrmath $ADCtemp/*bin_templspace.mif min $ADCtemp/mask_intersection.mif

# Compute a white matter analysis voxel & fixel mask
#mrconvert $FODtemp/wm_fod_template.mif -coord 3 0 - | mrthreshold - $ADCtemp/voxel_mask.mif
# Create FixelMask :
#fod2fixel -mask $ADCtemp/voxel_mask.mif -fmls_peak_value 0.1 $FODtemp/wm_fod_template.mif $FBA/fixel_mask -force

# Warp FODs to template
# mrtransform $FODaverageResponse/$1_wmfod_norm.mif -warp $FODtemp/$1_subject2template_warp.mif -noreorientation $FBA/$1_fod_in_template_space.mif

# Create AFD fixels
#fod2fixel $FBA/$1_fod_in_template_space.mif -mask $ADCtemp/voxel_mask.mif $FBA/$1_fixeltemp -afd $1_fd.mif -force
#fixelreorient $FBA/$1_fixeltemp $FODtemp/$1_subject2template_warp.mif $FBA/$1_fixeltemp -force
#fixelcorrespondence $FBA/$1_fixeltemp/$1_fd.mif $FBA/fixel_mask $FBA/fd $1_fd.mif -force
###
##### Create FC fixel
#warp2metric $FODtemp/$1_subject2template_warp.mif -fc $FBA/fixel_mask $FBA/fc $1_fc.mif -force
#mkdir $FBA/fc_log
#cp $FBA/fc/index.mif $FBA/fc_log;
#cp $FBA/fc/directions.mif $FBA/fc_log;
#mrcalc $FBA/fc/$1_fc.mif -log $FBA/fc_log/$1_fc_log.mif
###
#mkdir $FBA/fdc
#cp $FBA/fc/index.mif $FBA/fdc/;
#cp $FBA/fc/directions.mif $FBA/fdc/;
#mrcalc $FBA/fd/$1_fd.mif $FBA/fc/$1_fc.mif -mult $FBA/fdc/$1_fdc.mif
##
# Wholebrain tractography
#tckgen -angle 22.5 -maxlen 250 -minlen 10 -power 1.0 $FODtemp/wm_fod_template.mif -seed_image $ADCtemp/voxel_mask.mif -mask $ADCtemp/voxel_mask.mif -select 20000000 $FBA/tracks_20_million.tck
#
#tcksift $FBA/tracks_20_million.tck $FODtemp/wm_fod_template.mif $FBA/tracks_2_million_sift.tck -term_number 2000000

# FBA statistics
cd $FBA
#cat fd_files.txt | tr "\r" "\n" > fd_files.txt
#cat logfc_files.txt | tr "\r" "\n" > fc_filesUN.txt; rm fc_files.txt; mv fc_filesUN.txt logfc_files.txt;
#cat fdc_files.txt | tr "\r" "\n" > fdc_filesUN.txt; rm fdc_files.txt; mv fdc_filesUN.txt fdc_files.txt;
#fixelcfestats fd fd_files.txt designMatrix.txt contrastPatLess.txt tracks_2_million_sift.tck stats_fd_PatLess -mask fixel_mask/directions.mif
fixelcfestats fd fd_files.txt designMatrix.txt contrastConLess.txt tracks_2_million_sift.tck stats_fd_ConLess -mask fixel_mask/directions.mif
fixelcfestats fc_log logfc_files.txt designMatrix.txt contrastPatLess.txt $FBA/tracks_2_million_sift.tck stats_log_fc_PatLess -mask fixel_mask/directions.mif
fixelcfestats fc_log logfc_files.txt designMatrix.txt contrastConLess.txt $FBA/tracks_2_million_sift.tck stats_log_fc_ConLess -mask fixel_mask/directions.mif
fixelcfestats fdc fdc_files.txt designMatrix.txt contrastPatLess.txt $FBA/tracks_2_million_sift.tck stats_fdc_PatLess -mask fixel_mask/directions.mif
fixelcfestats fdc fdc_files.txt designMatrix.txt contrastConLess.txt $FBA/tracks_2_million_sift.tck stats_fdc_ConLess -mask fixel_mask/directions.mif

fixelcfestats fd fd_patientsOnly.txt design_patientsOnly.txt contrastRTLess.txt tracks_2_million_sift.tck stats_fd_RTLess -mask fixel_mask/directions.mif
fixelcfestats fd fd_patientsOnly.txt design_patientsOnly.txt contrastPatLess.txt tracks_2_million_sift.tck stats_fd_noRTLess -mask fixel_mask/directions.mif
