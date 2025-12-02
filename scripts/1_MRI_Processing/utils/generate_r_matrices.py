import os, sys, argparse

import numpy as np
import pandas as pd

import nibabel as nib
from nilearn.maskers import NiftiLabelsMasker


###################
## ADD ARGUMENTS ##
###################
parser = argparse.ArgumentParser(description="Generate R matrices for multi-modal network modeling.")
parser.add_argument("--working_directory",type=str,default="none")
parser.add_argument("--subject",type=str,default="none")
parser.add_argument("--session",type=str,default="none")
parser.add_argument("--direction",type=str,default="none")
args=parser.parse_args()


#######################
## EXTRACT ARGUMENTS ##
#######################
work_dir=args.working_directory
subject=args.subject
session=args.session
direction=args.direction


#########################
## GET DIRECTORY PATHS ##
#########################
ica_dir=os.path.join(work_dir,"data","ICA",subject,"melodic_out_REST%s_%s" % (session, direction))
parc_dir=os.path.join(work_dir,"data","parcellations",subject)

ica_img=nib.load(os.path.join(ica_dir,"melodic_IC_classified.nii.gz"))
parc_img=nib.load(os.path.join(parc_dir,"MMP1_MNI_REST%s_%s.nii.gz" % (session, direction)))

ica_data=ica_img.get_fdata()
parc_data=parc_img.get_fdata().astype(int)

region_labels=np.unique(parc_data)
region_labels=region_labels[region_labels!=0]
region_labels=region_labels[region_labels!=361]

n_components=ica_data.shape[3]
voxel_number_region=[]

for comp_idx in range(n_components):
    comp_data=ica_data[..., comp_idx]
    for region in region_labels:
        region_mask = parc_data == region
        nonzero_voxels = np.sum((comp_data != 0) & region_mask)
        voxel_number_region.append({
            'Component': comp_idx,
            'Region': region,
            'Nonzero_Voxels': nonzero_voxels
        })

voxel_number_region=pd.DataFrame(voxel_number_region)

mask_data=np.copy(ica_data)
mask_data[mask_data>0]=1
mask_img=nib.Nifti1Image(mask_data,ica_img.affine)

lut=pd.read_csv(os.path.join(work_dir,"parcellations","hcpmmp1_ordered.txt"),escapechar='#',skipinitialspace=True,delim_whitespace=True)

masker=NiftiLabelsMasker(
    labels_img=parc_img,
    lut=lut,
    standardize=False,
    strategy="sum")

r_matrix_sum=masker.fit_transform(ica_img)
r_matrix_sum=np.transpose(r_matrix_sum)
np.savetxt(os.path.join(ica_dir,"R_matrix_sum.csv"),r_matrix_sum,delimiter=",")


r_matrix_subset_avg=np.copy(r_matrix_sum)
for comp_idx in range(n_components):
    tmp=r_matrix_subset_avg[:,comp_idx]
    voxel_sums=voxel_number_region.loc[voxel_number_region['Component']==comp_idx, 'Nonzero_Voxels'].to_numpy()
    tmp=tmp / voxel_sums
    tmp=np.nan_to_num(tmp,nan=0,posinf=0,neginf=0)
    r_matrix_subset_avg[:,comp_idx]=tmp

np.savetxt(os.path.join(ica_dir,"R_matrix_subset_average.csv"),r_matrix_subset_avg,delimiter=",")


masker=NiftiLabelsMasker(
    labels_img=parc_img,
    lut=lut,
    standardize=False,
    strategy="mean")

r_matrix_avg=masker.fit_transform(ica_img)
r_matrix_avg=np.transpose(r_matrix_avg)
np.savetxt(os.path.join(ica_dir,"R_matrix_average.csv"),r_matrix_avg,delimiter=",")

