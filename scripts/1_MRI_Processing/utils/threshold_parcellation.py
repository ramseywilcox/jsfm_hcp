import os, sys, argparse

import numpy as np
import nibabel as nib

parser = argparse.ArgumentParser(description="Threshold parcellation.")
parser.add_argument("--working_directory",type=str,default="none")
parser.add_argument("--subject",type=str,default="none")
args=parser.parse_args()

work_dir=args.working_directory
subject=args.subject

parc_dir=os.path.join(work_dir,"data","parcellations",subject)

parc_img=nib.load(os.path.join(parc_dir,"MMP1_MNI.nii.gz"))
parc_data=parc_img.get_fdata()

parc_data_thresh=np.where((parc_data >= 1) & (parc_data <= 360), parc_data, 0)
parc_data_thresh=parc_data_thresh.astype(np.int16)

parc_img_thresh=nib.Nifti1Image(parc_data_thresh, affine=parc_img.affine)
nib.save(parc_img_thresh, os.path.join(parc_dir,"MMP1_MNI.nii.gz"))
