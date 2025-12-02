import os, sys, argparse

import numpy as np
import pandas as pd

import nibabel as nib

import random


## --- Add arguments
parser = argparse.ArgumentParser(description="Perform ICA component cleaning.")
parser.add_argument("--working_directory",type=str,default="none")
parser.add_argument("--subject",type=str,default="none")
parser.add_argument("--session",type=str,default="none")
parser.add_argument("--direction",type=str,default="none")
parser.add_argument("--tr",type=float,default="none")
args=parser.parse_args()


## --- Parse Arguments
work_dir=args.working_directory
subject=args.subject
session=args.session
direction=args.direction
tr=args.tr


###################
## INITIAL SETUP ##
###################

## --- Set path to data
ica_dir=os.path.join(work_dir,"data","ICA",subject,"melodic_out_REST%s_%s" % (session, direction))
fmri_dir=os.path.join(work_dir,"data","preprocessed","rest",subject+"_3T_rfMRI_REST"+session+"_preproc",subject,"MNINonLinear","Results","rfMRI_REST"+session+"_"+direction)
parc_dir=os.path.join(work_dir,"parcellations")

## --- Load images
ica_img=nib.load(os.path.join(ica_dir,"melodic_IC.nii.gz"))
ica_thr_img=nib.load(os.path.join(ica_dir,"melodic_IC_thresh50.nii.gz"))
mask_img=nib.load(os.path.join(ica_dir,"mask.nii.gz"))
gm_mask_img=nib.load(os.path.join(parc_dir,subject,"MMP1_MNI_REST"+session+"_"+direction+"_binarized.nii.gz"))

## --- Extract data from images
ica_data=ica_img.get_fdata()
ica_thr_data=ica_thr_img.get_fdata()
mask_data=mask_img.get_fdata()
gm_mask_data=gm_mask_img.get_fdata()

## --- Read in classification variables
mp=np.loadtxt(os.path.join(fmri_dir,"Movement_Regressors.txt"))
FT=np.loadtxt(os.path.join(ica_dir,"melodic_FTmix"))
melmix=np.loadtxt(os.path.join(ica_dir,"melodic_mix"))

## --- Create list of components
n_components=ica_data.shape[3]
components=list(range(n_components))

## --- Old python division function for classification 
def old_div(a, b):
    import numbers
    if isinstance(a, numbers.Integral) and isinstance(b, numbers.Integral):
        return a // b
    else:
        return a / b

## --- Cross correlations between columns of two matrices
def cross_correlation(a, b):
    """Cross Correlations between columns of two matrices"""
    assert a.ndim == b.ndim == 2
    _, ncols_a = a.shape
    # nb variables in columns rather than rows hence transpose
    # extract just the cross terms between cols in a and cols in b
    return np.corrcoef(a.T, b.T)[:ncols_a, ncols_a:]

## --- Set high frequency content classification threshold (this is fairly conservative, original had 0.35)
thr_HFC = 0.30

## --- Set motion parameter classification threshold
thr_MP = 0.30

## --- Set skewness distribution classification threshold
thr_SKEW = 0.1

## --- Set grey matter fraction classification threshold
thr_FRACT = 0.48


##################################
## RUN FREQUENCY CLASSIFICATION ##
##################################

## --- Determine sample frequency
Fs = old_div(1, tr)

## --- Determine Nyquist-frequency
Ny = old_div(Fs, 2)

## --- Determine which frequencies are associated with every row in the mleodic_FTmix file
f = Ny * (np.array(list(range(1, FT.shape[0] + 1)))) / (FT.shape[0])

## --- Only include frequencies higher than 0.01Hz
fincl = np.squeeze(np.array(np.where(f > 0.01)))
FT = FT[fincl, :]
f = f[fincl]

## --- Set frequency range to [0-1]
f_norm = old_div((f - 0.01), (Ny - 0.01))

## --- For every IC; get the cumulative sum as a fraction of the total sum
fcumsum_fract = old_div(np.cumsum(FT, axis=0), np.sum(FT, axis=0))

## --- Determine the index of the frequency with the fractional cumulative sum closest to 0.5
idx_cutoff = np.argmin(np.abs(fcumsum_fract - 0.5), axis=0)

## --- Now get the fractions associated with those indices index, these are the final feature scores
HFC = f_norm[idx_cutoff]

## --- Extract which components do not pass the high frequency content threshold, i.e., ICs with too much high frequency content
HFC_components = np.where(HFC > thr_HFC)[0]


###############################
## RUN MOTION CLASSIFICATION ##
###############################

## --- Get shape from motion parameters file
_, nparams = mp.shape

## --- Get forward shifted motion parameter version
mp_1fw = np.vstack((
    np.zeros(nparams),
    mp[:-1]
))

## --- Get backward shifted motion parameter version
mp_1bw = np.vstack((
    mp[1:],
    np.zeros(nparams)
))

## --- Concatenate original, forward shifted, and backward shifted motion parameters
mp_model = np.hstack((mp, mp_1fw, mp_1bw))

## --- Determine the maximum correlation between motion parameters and IC time-series
nsplits=10000
nmixrows, nmixcols=melmix.shape
nrows_to_choose=int(round(0.9 * nmixrows))

## --- Max correlations for multiple splits of the dataset (for a robust estimate)
max_correls=np.empty((nsplits, nmixcols))
for i in range(nsplits):
    # Select a random subset of 90% of the dataset rows (*without* replacement)
    chosen_rows=random.sample(population=range(nmixrows), k=nrows_to_choose)
    correl_nonsquared = cross_correlation(melmix[chosen_rows],mp_model[chosen_rows])
    correl_squared = cross_correlation(melmix[chosen_rows]**2,mp_model[chosen_rows]**2)
    correl_both = np.hstack((correl_squared, correl_nonsquared))
    max_correls[i] = np.abs(correl_both).max(axis=1)

## --- Feature score is the mean of the maximum correlation over all the random splits (each split will have 3 correlations, which we got the maximum from in the loop above)
mean_correls=np.nanmean(max_correls,axis=0)
MP_components=np.where(mean_correls > thr_MP)[0]

## --- Add noise components
HFC_MP_components=np.concatenate((HFC_components,MP_components),axis=0)
HFC_MP_components_unique=np.unique(HFC_MP_components)

## --- Set temporary good components, so we can loop through a smaller number (i.e., the ones that passed previous classifications)
good_components_tmp = list(set(components) - set(HFC_MP_components_unique))
good_components_tmp.sort(reverse=False)


#################################
## RUN SKEWNESS CLASSIFICATION ##
#################################

## --- Create empty list to store Pearson's skewness coefficients
skew_list=[]

## --- Flatten mask data for indexing
mask_data_flat=mask_data.flatten()

## Loop through components that haven't been classified as bad
for comp in good_components_tmp:
    comp_data=ica_data[:,:,:,comp]
    comp_data_flat=comp_data.flatten()
    brain_idx=np.where(mask_data_flat > 0)[0].tolist()
    comp_data_flat_brain=comp_data_flat[brain_idx]
    mean=np.mean(comp_data_flat_brain)
    median=np.median(comp_data_flat_brain)
    s=np.std(comp_data_flat_brain)
    skew=(3 * (mean-median)) / s
    skew_list.append(skew)

## --- Get index of bad components
good_components_tmp_array=np.array(good_components_tmp)
skew_array=np.array(skew_list)
SKEW_components=good_components_tmp_array[np.where(skew_array < thr_SKEW)[0].tolist()]
HFC_MP_SKEW_components=np.concatenate((HFC_MP_components,SKEW_components),axis=0)
HFC_MP_SKEW_components_unique=np.unique(HFC_MP_SKEW_components)

## --- Create second temporary good components list
good_components_tmp2=list(set(components) - set(HFC_MP_SKEW_components))
good_components_tmp2.sort(reverse=False)


#############################################
## RUN GREY MATTER FRACTION CLASSIFICATION ##
#############################################

## --- Create empty list to store fractions
fraction_list=[]

## --- Flatten grey matter mask
gm_mask_flat=gm_mask_data.flatten()

## --- Loop through remaining good components
for comp in good_components_tmp2:
    comp_data=ica_thr_data[:,:,:,comp]
    comp_data_flat=comp_data.flatten()
    nz_idx=np.where(comp_data_flat != 0)[0].tolist()
    total_voxels=comp_data_flat[nz_idx].shape[0]
    total_mean=np.mean(comp_data_flat[nz_idx])
    total_summary=total_voxels*total_mean
    gm_idx=np.where(gm_mask_flat != 0)[0].tolist()
    gm_voxels=np.array(list(set(gm_idx) & set(nz_idx))).shape[0]
    gm_mean=np.mean(comp_data_flat[list(set(gm_idx) & set(nz_idx))])
    gm_summary=gm_voxels*gm_mean
    gm_fraction=gm_summary/total_summary
    #gm_fraction=gm_voxels/total_voxels
    fraction_list.append(gm_fraction)

## --- Get index of bad components
fraction_array=np.array(fraction_list)
good_components_tmp2_array=np.array(good_components_tmp2)
FRACT_components=good_components_tmp2_array[np.where(fraction_array <= thr_FRACT)[0].tolist()]
HFC_MP_SKEW_FRACT_components=np.concatenate((HFC_MP_SKEW_components,FRACT_components),axis=0)
HFC_MP_SKEW_FRACT_components_unique=np.unique(HFC_MP_SKEW_FRACT_components)

## --- Create final good components list
good_components_final=list(set(components) - set(HFC_MP_SKEW_FRACT_components_unique))
good_components_final.sort(reverse=False)


#####################
## REFINE AND SAVE ##
#####################

## --- Save good components
np.savetxt(os.path.join(ica_dir,"good_components.csv"),good_components_final,delimiter=",")

## --- Extract only good components from the thresholded ICA image
ica_thr_data_good=ica_thr_data[:,:,:,good_components_final]

## --- Remove values below zero, joint structure-function model won't take negative values
ica_thr_data_good=ica_thr_data_good.clip(min=0)
ica_thr_img_good=nib.Nifti1Image(ica_thr_data_good,ica_thr_img.affine)
nib.save(ica_thr_img_good,os.path.join(ica_dir,"melodic_IC_classified.nii.gz"))
