#!/bin/bash

################################################################
# MUST FOLLOW BIDS ARCHITECTURE:
# sub
# 	-anat
# 	 	-*T1w.nii.gz
# 	-dwi
# 	 	-*.bvec
# 	 	-*.bval
# 	 	-*.json
# 	 	-*dwi.nii.gz
################################################################

################################################################
# For the command "mri_aparc2aseg", there can be an error which is due to the way multithreading is handled. Just rerun the command manually for the subjects which did not have an output hcpmmp1.mgz
# You can find these subjects by typing "find . -name *.mgz" in a terminal in the directory that contains all your subjects folders
################################################################

################################################################
# Find number of available cpus
NPROC=$(nproc)
################################################################

############################### STEP 1 ###############################
#             Convert data to .mif format and denoise                #
######################################################################

# Also consider doing Gibbs denoising (using mrdegibbs). Check your diffusion data for ringing artifacts before deciding whether to use it
for_each -nthreads $NPROC  -info */dwi : mrconvert IN/*dwi.nii.gz IN/dwi.mif
for_each -nthreads $NPROC -info */dwi : mrconvert IN/dwi.mif -fslgrad IN/*.bvec IN/*.bval IN/dwi_header.mif

for_each -nthreads $NPROC  -info */dwi : dwidenoise IN/dwi_header.mif IN/dwi_den.mif -noise IN/noise.mif
for_each -nthreads $NPROC  -info */dwi : mrdegibbs IN/dwi_den.mif IN/dwi_den_unr.mif

# Extract the b0 images from the diffusion data acquired in the AP direction
for_each -nthreads $NPROC  -info */dwi : dwiextract IN/dwi_den.mif - -bzero \| mrmath - mean IN/mean_b0_AP.mif -axis 3

# Runs the dwipreproc command, which is a wrapper for eddy and topup.
for_each -nthreads 8 -info */dwi : dwifslpreproc IN/dwi_den.mif IN/dwi_den_preproc.mif -pe_dir AP -rpe_none -readout_time 0.0342002 -eddy_options " --slm=linear --data_is_shelled" -nthreads 8

# Performs bias field correction. Needs ANTs to be installed in order to use the "ants" option (use "fsl" otherwise)
for_each -nthreads $NPROC  -info */dwi : dwibiascorrect ants IN/dwi_den_preproc.mif IN/dwi_den_preproc_unbiased.mif -bias IN/bias.mif

# Create a mask for future processing steps
for_each -nthreads $NPROC  -info */dwi : dwi2mask IN/dwi_den_preproc_unbiased.mif IN/mask.mif

########################### STEP 2 ###################################
#             Basis function for each tissue type                    #
######################################################################

# Create a basis function from the subject's DWI data. The "dhollander" function is best used for multi-shell acquisitions; it will estimate different basis functions for each tissue type. For single-shell acquisition, use the "tournier" function instead
for_each -nthreads $NPROC  -info */dwi : dwi2response dhollander IN/dwi_den_preproc_unbiased.mif IN/wm.txt IN/gm.txt IN/csf.txt -voxels IN/voxels.mif

# Performs multishell-multitissue constrained spherical deconvolution, using the basis functions estimated above
for_each -nthreads $NPROC  -info */dwi : dwi2fod msmt_csd IN/dwi_den_preproc_unbiased.mif -mask IN/mask.mif IN/wm.txt IN/wmfod.mif IN/gm.txt IN/gmfod.mif IN/csf.txt IN/csffod.mif

# Creates an image of the fiber orientation densities overlaid onto the estimated tissues (Blue=WM; Green=GM; Red=CSF)
# You should see FOD's mostly within the white matter. These can be viewed later with the command "mrview vf.mif -odf.load_sh wmfod.mif"
for_each -nthreads $NPROC  -info */dwi : mrconvert -coord 3 0 IN/wmfod.mif - \| mrcat IN/csffod.mif IN/gmfod.mif - IN/vf.mif

# Now normalize the FODs to enable comparison between subjects
for_each -nthreads $NPROC  -info */dwi : mtnormalise IN/wmfod.mif IN/wmfod_norm.mif IN/gmfod.mif IN/gmfod_norm.mif IN/csffod.mif IN/csffod_norm.mif -mask IN/mask.mif

########################### STEP 3 ###################################
#            Create a GM/WM boundary for seed analysis               #
######################################################################

# Convert the anatomical image to .mif format, and then extract all five tissue catagories (1=GM; 2=Subcortical GM; 3=WM; 4=CSF; 5=Pathological tissue)
for_each -nthreads $NPROC  -info */dwi : mrconvert IN/../anat/*T1w.nii.gz IN/T1.mif
for_each -nthreads 8 -info */dwi : 5ttgen fsl IN/T1.mif IN/5tt_nocoreg.mif -nthreads 8

# The following series of commands will take the average of the b0 images (which have the best contrast), convert them and the 5tt image to NIFTI format, and use it for coregistration.
for_each -nthreads $NPROC  -info */dwi : dwiextract IN/dwi_den_preproc_unbiased.mif - -bzero \| mrmath - mean IN/mean_b0_processed.mif -axis 3
for_each -nthreads $NPROC  -info */dwi : mrconvert IN/mean_b0_processed.mif IN/mean_b0_processed.nii.gz
for_each -nthreads $NPROC  -info */dwi : mrconvert IN/5tt_nocoreg.mif IN/5tt_nocoreg.nii.gz

# Uses FSL commands fslroi and flirt to create a transformation matrix for regisitration between the tissue map and the b0 images
for_each -nthreads $NPROC  -info */dwi : fslroi IN/5tt_nocoreg.nii.gz IN/5tt_vol0.nii.gz 0 1 #Extract the first volume of the 5tt dataset (since flirt can only use 3D images, not 4D images)
for_each -nthreads $NPROC  -info */dwi : flirt -in IN/mean_b0_processed.nii.gz -ref IN/5tt_vol0.nii.gz -interp nearestneighbour -dof 6 -omat IN/diff2struct_fsl.mat
for_each -nthreads $NPROC  -info */dwi : transformconvert IN/diff2struct_fsl.mat IN/mean_b0_processed.nii.gz IN/5tt_nocoreg.nii.gz flirt_import IN/diff2struct_mrtrix.txt
for_each -nthreads $NPROC  -info */dwi : mrtransform IN/5tt_nocoreg.mif -linear IN/diff2struct_mrtrix.txt -inverse IN/5tt_coreg.mif

#Create a seed region along the GM/WM boundary
for_each -nthreads $NPROC  -info */dwi : 5tt2gmwmi IN/5tt_coreg.mif IN/gmwmSeed_coreg.mif

########################## STEP 4 ###################################
#                 Run the streamline analysis                        #
######################################################################

# Create streamlines
# Note that the "right" number of streamlines is still up for debate. Last I read from the MRtrix documentation,
# They recommend about 100 million tracks. Here I use 10 million, if only to save time. Read their papers and then make a decision
for_each -nthreads 8 -info */dwi : tckgen -act IN/5tt_coreg.mif -backtrack -seed_gmwmi IN/gmwmSeed_coreg.mif -nthreads 8 -maxlength 250 -cutoff 0.06 -select 10000k IN/wmfod_norm.mif IN/tracks_10M.tck

# Extract a subset of tracks (here, 200 thousand) for ease of visualization
# tckedit tracks_10M.tck -number 200k smallerTracks_200k.tck

# Reduce the number of streamlines with tcksift
for_each -nthreads 8 -info */dwi : tcksift2 -act IN/5tt_coreg.mif -out_mu IN/sift_mu.txt -out_coeffs IN/sift_coeffs.txt -nthreads 8 IN/tracks_10M.tck IN/wmfod_norm.mif IN/sift_1M.txt


########################### RECON-ALL ################################
######################################################################

mkdir FS

ls .| grep ^sub- > subjList.txt

for sub in `cat subjList.txt`; do
	cp ./${sub}/anat/${sub}_T1w.nii.gz ./FS
done

cd FS

gunzip *.gz

SUBJECTS_DIR=`pwd`
ls *.nii | parallel --jobs $NPROC recon-all -s {.}_recon -i {} -all
rm *.nii


# Move into freesurfer's subjects directory
echo <password> | sudo -S mv *T1w_recon $SUBJECTS_DIR
# Grant permission to write the file
echo <password> | sudo -S chmod -R a+w $FREESURFER_HOME

######################################################################
# Glasser annotation
######################################################################

ls . | grep ^sub- > subjList.txt

for sub in `cat subjList.txt`; do
	mri_surf2surf --srcsubject fsaverage --trgsubject ${sub}_T1w_recon --hemi lh --sval-annot $SUBJECTS_DIR/fsaverage/label/lh.hcpmmp1.annot --tval $SUBJECTS_DIR/${sub}_T1w_recon/label/lh.hcpmmp1.annot
	mri_surf2surf --srcsubject fsaverage --trgsubject ${sub}_T1w_recon --hemi rh --sval-annot $SUBJECTS_DIR/fsaverage/label/rh.hcpmmp1.annot --tval $SUBJECTS_DIR/${sub}_T1w_recon/label/rh.hcpmmp1.annot

	cd ./${sub}/dwi
	mri_aparc2aseg --old-ribbon --s ${sub}_T1w_recon --annot hcpmmp1 --o hcpmmp1.mgz --nthreads $NPROC
	cd ../..
done

######################################################################
# SC Generation (Streamline density)
######################################################################

for_each -nthreads $NPROC -info */dwi : mrconvert –datatype uint32  IN/hcpmmp1.mgz  IN/hcpmmp1.mif

# Replace the random integers of the hcpmmp1.mif file with integers
# that start at 1 and increase by 1.
for_each -nthreads $NPROC -info */dwi : labelconvert  IN/hcpmmp1.mif $MRtrix3/labelconvert/hcpmmp1_original.txt $MRtrix3/labelconvert/hcpmmp1_ordered.txt  IN/hcpmmp1_parcels_nocoreg.mif

# Register the ordered atlas-based volumetric parcellation to diffusion space.
for_each -nthreads $NPROC -info */dwi : mrtransform  IN/hcpmmp1_parcels_nocoreg.mif –linear  IN/diff2struct_mrtrix.txt –inverse –datatype uint32  IN/hcpmmp1_parcels_coreg.mif

# Create a whole-brain connectome, representing the streamlines between each parcellation pair in the atlas (in this case, 84x84). The "symmetric" option will make the lower diagonal the same as the upper diagonal, 
# and the "scale_invnodevol" option will scale the connectome by the inverse of the size of the node

for_each -nthreads $NPROC -info * : tck2connectome -symmetric -zero_diagonal -scale_invnodevol -tck_weights_in IN/dwi/sift_1M.txt IN/dwi/tracks_10M.tck IN/dwi/hcpmmp1_parcels_coreg.mif IN/dwi/IN_hcpmmp1_parcels_coreg.csv -out_assignment IN/dwi/assignments_IN_hcpmmp1_parcels_coreg.csv
