#!/bin/bash
#
#SBATCH --time=1-10:59:00
#SBATCH --mem=20GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -J run_ica_classify
#SBATCH --output /work/out_batch/out_ica_classify-%A_%a.out
#SBATCH --error /work/out_batch/ica_classify-%A_%a.err
#


module purge
module load fsl
module load freesurfer
module load anaconda

conda activate $WORK/mri_conda_env

SESSIONS=$(seq 1 2)
DIRECTIONS=("LR RL")

HDIR=$WORK/projects/HCP
IDIR=$HDIR/data/preprocessed/rest
ICADIR=$HDIR/data/ICA
PDIR=$HDIR/data/parcellations

subjs=($@)
subjs=$1
S=${subjs[${SLURM_ARRAY_TASK_ID}]}

if [ -f ${ICADIR}/${S}/melodic_out_REST1_LR/melodic_IC_thresh50.nii.gz ] &&
   [ -f ${ICADIR}/${S}/melodic_out_REST2_LR/melodic_IC_thresh50.nii.gz ] &&
   [ -f ${ICADIR}/${S}/melodic_out_REST1_RL/melodic_IC_thresh50.nii.gz ] &&
   [ -f ${ICADIR}/${S}/melodic_out_REST2_RL/melodic_IC_thresh50.nii.gz ]; then

    for SES in $SESSIONS; do

        for DIR in $DIRECTIONS; do

            flirt -in ${PDIR}/${S}/MMP1_MNI.nii.gz \
                -ref ${IDIR}/${S}_3T_rfMRI_REST${SES}_preproc/${S}/MNINonLinear/Results/rfMRI_REST${SES}_${DIR}/rfMRI_REST${SES}_${DIR}_SBRef.nii.gz \
                -out ${PDIR}/${S}/MMP1_MNI_REST${SES}_${DIR}.nii.gz -dof 6 -interp nearestneighbour

            mri_binarize --i ${PDIR}/${S}/MMP1_MNI_REST${SES}_${DIR}.nii.gz --min 0.00001 --o ${PDIR}/${S}/MMP1_MNI_REST${SES}_${DIR}_binarized.nii.gz

            python ${HDIR}/scripts/1_MRI_Processing/utils/ICA_Component_Classification.py --working_directory=${HDIR} \
                --subject=${S} \
                --session=${SES} \
                --direction=${DIR} \
                --tr=0.72

        done

    done

fi
