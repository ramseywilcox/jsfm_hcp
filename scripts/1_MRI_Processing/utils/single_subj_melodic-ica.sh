#!/bin/bash
#
#SBATCH --time=1-10:59:00
#SBATCH --mem=64GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -J run_ica
#SBATCH --output /work/out_batch/out_HCPica-%A_%a.out
#SBATCH --error /work/out_batch/HCPica-%A_%a.err
#

module purge
module load fsl

RUNS=$(seq 1 2)
DIRECTIONS=("LR RL")

HDIR=$WORK/projects/HCP
IDIR=$HDIR/data/preprocessed/rest
ODIR=$HDIR/data/ICA

SUBJS=($@)
SUBJS=$1
S=${SUBJS[${SLURM_ARRAY_TASK_ID}]}

if [ -f ${IDIR}/${S}_3T_rfMRI_REST1_preproc/${S}/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_smoothed_HP_detrended.nii.gz ] && \
   [ -f ${IDIR}/${S}_3T_rfMRI_REST2_preproc/${S}/MNINonLinear/Results/rfMRI_REST2_RL/rfMRI_REST2_LR_smoothed_HP_detrended.nii.gz ] && \
   [ -f ${IDIR}/${S}_3T_rfMRI_REST1_preproc/${S}/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_RL_smoothed_HP_detrended.nii.gz ] && \
   [ -f ${IDIR}/${S}_3T_rfMRI_REST2_preproc/${S}/MNINonLinear/Results/rfMRI_REST2_RL/rfMRI_REST2_RL_smoothed_HP_detrended.nii.gz ] && \
   [ ! -f ${ODIR}/${S}/melodic_out_REST1_LR/melodic_IC.nii.gz ] && \
   [ ! -f ${ODIR}/${S}/melodic_out_REST1_RL/melodic_IC.nii.gz ] && \
   [ ! -f ${ODIR}/${S}/melodic_out_REST2_LR/melodic_IC.nii.gz ] && \
   [ ! -f ${ODIR}/${S}/melodic_out_REST2_RL/melodic_IC.nii.gz ]; then

    for RUN in $RUNS; do

        for D in $DIRECTIONS; do

            mkdir ${ODIR}/${S}/melodic_out_REST${RUN}_${D}
	        melodic -i ${IDIR}/${S}_3T_rfMRI_REST${RUN}_preproc/${S}/MNINonLinear/Results/rfMRI_REST${RUN}_${D}/rfMRI_REST${RUN}_${D}_smoothed_HP_detrended.nii.gz \
	                -o ${ODIR}/${S}/melodic_out_REST${RUN}_${D} --Ostats --Oorig

        done

    done

fi
