#!/bin/bash
#
#SBATCH --time=10:59:00
#SBATCH --mem=32GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -J run_parc_reg
#SBATCH --output /work/out_batch/out_parc_reg-%A_%a.out
#SBATCH --error /work/out_batch/parc_reg-%A_%a.err
#

module purge
module load fsl
module load connectome-workbench

subjs=($@)
subjs=$1
S=${subjs[${SLURM_ARRAY_TASK_ID}]}

HDIR=$WORK/projects/HCP
ADIR=$HDIR/data/preprocessed/anat
ODIR=$HDIR/data/tractography
PDIR=$HDIR/data/parcellations

mkdir ${PDIR}/${S}

unzip ${ADIR}/${S}_3T_Structural_preproc.zip -d ${ADIR}/

LH_MIDTHICKNESS=${ADIR}/${S}/MNINonLinear/fsaverage_LR32k/${S}.L.midthickness.32k_fs_LR.surf.gii
RH_MIDTHICKNESS=${ADIR}/${S}/MNINonLinear/fsaverage_LR32k/${S}.R.midthickness.32k_fs_LR.surf.gii

LH_WHITE=${ADIR}/${S}/MNINonLinear/fsaverage_LR32k/${S}.L.white.32k_fs_LR.surf.gii
RH_WHITE=${ADIR}/${S}/MNINonLinear/fsaverage_LR32k/${S}.R.white.32k_fs_LR.surf.gii

LH_PIAL=${ADIR}/${S}/MNINonLinear/fsaverage_LR32k/${S}.L.pial.32k_fs_LR.surf.gii
RH_PIAL=${ADIR}/${S}/MNINonLinear/fsaverage_LR32k/${S}.R.pial.32k_fs_LR.surf.gii

MNI_REF=${ADIR}/${S}/MNINonLinear/T1w_restore_brain.nii.gz
NATIVE_REF=${ADIR}/${S}/T1w/T1w_acpc_dc_restore_brain.nii.gz

WARP_IMG=${ADIR}/${S}/MNINonLinear/xfms/standard2acpc_dc.nii.gz

wb_command -label-to-volume-mapping ${PDIR}/MMP1_MNI_lh.label.gii ${LH_MIDTHICKNESS} ${MNI_REF} ${PDIR}/${S}/MMP1_MNI_lh.nii.gz -ribbon-constrained ${LH_WHITE} ${LH_PIAL}
wb_command -label-to-volume-mapping ${PDIR}/MMP1_MNI_rh.label.gii ${RH_MIDTHICKNESS} ${MNI_REF} ${PDIR}/${S}/MMP1_MNI_rh.nii.gz -ribbon-constrained ${RH_WHITE} ${RH_PIAL}
fslmaths ${PDIR}/${S}/MMP1_MNI_lh.nii.gz -add ${PDIR}/${S}/MMP1_MNI_rh.nii.gz ${PDIR}/${S}/MMP1_MNI.nii.gz

module load anaconda
conda activate $WORK/mri_conda_env

python ${HDIR}/scripts/1_MRI_Processing/utils/threshold_parcellation.py --working_directory=${HDIR} \
	--subject=${S}

conda deactivate

applywarp \
    -i ${PDIR}/${S}/MMP1_MNI.nii.gz \
    -r ${NATIVE_REF} \
    -o ${PDIR}/${S}/MMP1_T1w.nii.gz \
    -w ${WARP_IMG} \
    --interp=nn

fslmaths ${PDIR}/${S}/MMP1_T1w.nii.gz -thr 0 -mul 1 -odt int ${PDIR}/${S}/MMP1_T1w.nii.gz

module load mrtrix3
mrconvert ${PDIR}/${S}/MMP1_T1w.nii.gz ${PDIR}/${S}/MMP1_T1w.mif -force

