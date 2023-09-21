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

ls . | grep ^sub- > subjList.txt

for sub in `cat subjList.txt`; do
	echo "${sub}";

	cd ./${sub}/dwi

	############################### STEP 1 ###############################
	#             Convert data to .mif format and denoise                #
	######################################################################

	# Also consider doing Gibbs denoising (using mrdegibbs). Check your diffusion data for ringing artifacts before deciding whether to use it
	mrconvert ${sub}_dwi.nii.gz dwi.mif
	mrconvert dwi.mif -fslgrad *.bvec *.bval dwi_header.mif

	dwidenoise dwi_header.mif dwi_den.mif -noise noise.mif
	mrdegibbs dwi_den.mif dwi_den_unr.mif

	# Extract the b0 images from the diffusion data acquired in the AP direction
	dwiextract dwi_den.mif - -bzero | mrmath - mean mean_b0_AP.mif -axis 3

	# Runs the dwipreproc command, which is a wrapper for eddy and topup. 
	# CAMCAN has LAS (check mrinfo dwi.nii: data strides -1, 2, 3, 4) and PhaseEncodingDirection: "j-" --> AP
	dwifslpreproc dwi_den.mif dwi_den_preproc.mif -pe_dir AP -rpe_none -readout_time 0.0342002 -eddy_options " --slm=linear --data_is_shelled" -nthreads 8

	# Performs bias field correction. Needs ANTs to be installed in order to use the "ants" option (use "fsl" otherwise)
	dwibiascorrect ants dwi_den_preproc.mif dwi_den_preproc_unbiased.mif -bias bias.mif

	# Create a mask for future processing steps
	dwi2mask dwi_den_preproc_unbiased.mif mask.mif

	########################### STEP 2 ###################################
	#             Basis function for each tissue type                    #
	######################################################################

	# Create a basis function from the subject's DWI data. The "dhollander" function is best used for multi-shell acquisitions; it will estimate different basis functions for each tissue type. For single-shell acquisition, use the "tournier" function instead
	dwi2response dhollander dwi_den_preproc_unbiased.mif wm.txt gm.txt csf.txt -voxels voxels.mif

	# Performs multishell-multitissue constrained spherical deconvolution, using the basis functions estimated above
	dwi2fod msmt_csd dwi_den_preproc_unbiased.mif -mask mask.mif wm.txt wmfod.mif gm.txt gmfod.mif csf.txt csffod.mif

	# Creates an image of the fiber orientation densities overlaid onto the estimated tissues (Blue=WM; Green=GM; Red=CSF)
	# You should see FOD's mostly within the white matter. These can be viewed later with the command "mrview vf.mif -odf.load_sh wmfod.mif"
	mrconvert -coord 3 0 wmfod.mif - | mrcat csffod.mif gmfod.mif - vf.mif

	# Now normalize the FODs to enable comparison between subjects
	mtnormalise wmfod.mif wmfod_norm.mif gmfod.mif gmfod_norm.mif csffod.mif csffod_norm.mif -mask mask.mif

	########################### STEP 3 ###################################
	#            Create a GM/WM boundary for seed analysis               #
	######################################################################

	# Convert the anatomical image to .mif format, and then extract all five tissue catagories (1=GM; 2=Subcortical GM; 3=WM; 4=CSF; 5=Pathological tissue)
	mrconvert ../anat/*T1w.nii.gz T1.mif
	5ttgen fsl T1.mif 5tt_nocoreg.mif

	# The following series of commands will take the average of the b0 images (which have the best contrast), convert them and the 5tt image to NIFTI format, and use it for coregistration.
	dwiextract dwi_den_preproc_unbiased.mif - -bzero | mrmath - mean mean_b0_processed.mif -axis 3
	mrconvert mean_b0_processed.mif mean_b0_processed.nii.gz
	mrconvert 5tt_nocoreg.mif 5tt_nocoreg.nii.gz

	# Uses FSL commands fslroi and flirt to create a transformation matrix for regisitration between the tissue map and the b0 images
	fslroi 5tt_nocoreg.nii.gz 5tt_vol0.nii.gz 0 1 #Extract the first volume of the 5tt dataset (since flirt can only use 3D images, not 4D images)
	flirt -in mean_b0_processed.nii.gz -ref 5tt_vol0.nii.gz -interp nearestneighbour -dof 6 -omat diff2struct_fsl.mat
	transformconvert diff2struct_fsl.mat mean_b0_processed.nii.gz 5tt_nocoreg.nii.gz flirt_import diff2struct_mrtrix.txt
	mrtransform 5tt_nocoreg.mif -linear diff2struct_mrtrix.txt -inverse 5tt_coreg.mif

	#Create a seed region along the GM/WM boundary
	5tt2gmwmi 5tt_coreg.mif gmwmSeed_coreg.mif

	########################### STEP 4 ###################################
	#                 Run the streamline analysis                        #
	######################################################################

	# Create streamlines
	# Note that the "right" number of streamlines is still up for debate. Last I read from the MRtrix documentation,
	# They recommend about 100 million tracks. Here I use 10 million, if only to save time. Read their papers and then make a decision
	tckgen -act 5tt_coreg.mif -backtrack -seed_gmwmi gmwmSeed_coreg.mif -nthreads 8 -maxlength 250 -cutoff 0.06 -select 10000k wmfod_norm.mif tracks_10M.tck

	# Extract a subset of tracks (here, 200 thousand) for ease of visualization
	# tckedit tracks_10M.tck -number 200k smallerTracks_200k.tck

	# Reduce the number of streamlines with tcksift
	tcksift2 -act 5tt_coreg.mif -out_mu sift_mu.txt -out_coeffs sift_coeffs.txt -nthreads 8 tracks_10M.tck wmfod_norm.mif sift_1M.txt

	cd ../..;
done


########################### RECON-ALL ################################
######################################################################
NPROC=$(nproc --all)
SUBJECTS_DIR=`pwd`

recon-all -s ${sub}_T1w_recon -i ./${sub}/anat/*T1w.nii.gz -all -nthreads $NPROC

# Move into freesurfer's subjects directory
mv *T1w_recon $SUBJECTS_DIR

######################################################################
# Glasser annotation & SC Generation (Streamline density)
######################################################################


ls . | grep ^sub- > subjList.txt

for sub in `cat subjList.txt`; do
	mri_surf2surf --srcsubject fsaverage --trgsubject ${sub}_T1w_recon --hemi lh --sval-annot $SUBJECTS_DIR/fsaverage/label/lh.hcpmmp1.annot --tval $SUBJECTS_DIR/${sub}_T1w_recon/label/lh.hcpmmp1.annot
	mri_surf2surf --srcsubject fsaverage --trgsubject ${sub}_T1w_recon --hemi rh --sval-annot $SUBJECTS_DIR/fsaverage/label/rh.hcpmmp1.annot --tval $SUBJECTS_DIR/${sub}_T1w_recon/label/rh.hcpmmp1.annot

	cd ./${sub}/dwi
  ################################################################
  # For the command "mri_aparc2aseg", there can be an error which is due to the way multithreading is handled. 
  # Just rerun the command manually until it works or use a single thread
  ################################################################
	mri_aparc2aseg --old-ribbon --s ${sub}_T1w_recon --annot hcpmmp1 --o hcpmmp1.mgz --nthreads $NPROC/2

	mrconvert –datatype uint32 hcpmmp1.mgz hcpmmp1.mif

	# Replace the random integers of the hcpmmp1.mif file with integers
	# that start at 1 and increase by 1.
	labelconvert hcpmmp1.mif $MRtrix3/labelconvert/hcpmmp1_original.txt $MRtrix3/labelconvert/hcpmmp1_ordered.txt hcpmmp1_parcels_nocoreg.mif

	# Register the ordered atlas-based volumetric parcellation to diffusion space.
	mrtransform hcpmmp1_parcels_nocoreg.mif –linear diff2struct_mrtrix.txt –inverse –datatype uint32 hcpmmp1_parcels_coreg.mif

	# Create a whole-brain connectome, representing the streamlines between each parcellation pair in the atlas (in this case, 379x379). The "symmetric" option will make the lower diagonal the same as the upper diagonal, and the "scale_invnodevol" option will scale the connectome by the inverse of the size of the node
	tck2connectome -symmetric -zero_diagonal -scale_invnodevol -tck_weights_in sift_1M.txt tracks_10M.tck hcpmmp1_parcels_coreg.mif ${sub}_hcpmmp1_parcels_coreg.csv -out_assignment assignments_${sub}_hcpmmp1_parcels_coreg.csv
	
	cd ../..
done

