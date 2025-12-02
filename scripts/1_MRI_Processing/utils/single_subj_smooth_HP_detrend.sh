#!/bin/bash
#
#SBATCH --time=1-10:59:00
#SBATCH --mem=32GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -J run_ica
#SBATCH --output /work/out_batch/out_HCPsmooth-%A_%a.out
#SBATCH --error /work/out_batch/HCPsmooth-%A_%a.err
#

module purge
module load fsl
module load afni

RUNS=$(seq 1 2)
DIRECTIONS=("LR RL")

HDIR=$WORK/projects/HCP
IDIR=$HDIR/data/preprocessed/rest

SUBJS=($@)
SUBJS=$1
S=${SUBJS[${SLURM_ARRAY_TASK_ID}]}

unzip ${IDIR}/${S}_3T_rfMRI_REST1_preproc.zip -d ${IDIR}/${S}_3T_rfMRI_REST1_preproc
unzip ${IDIR}/${S}_3T_rfMRI_REST2_preproc.zip -d ${IDIR}/${S}_3T_rfMRI_REST2_preproc

if [ -f ${IDIR}/${S}_3T_rfMRI_REST1_preproc/${S}/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR.nii.gz ] && \
   [ -f ${IDIR}/${S}_3T_rfMRI_REST2_preproc/${S}/MNINonLinear/Results/rfMRI_REST2_LR/rfMRI_REST2_LR.nii.gz ] && \
   [ -f ${IDIR}/${S}_3T_rfMRI_REST1_preproc/${S}/MNINonLinear/Results/rfMRI_REST1_RL/rfMRI_REST1_RL.nii.gz ] && \
   [ -f ${IDIR}/${S}_3T_rfMRI_REST2_preproc/${S}/MNINonLinear/Results/rfMRI_REST2_RL/rfMRI_REST2_RL.nii.gz ] && \
   [ ! -f ${IDIR}/${S}_3T_rfMRI_REST1_preproc/${S}/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_smoothed_HP_detrended.nii.gz ] && \
   [ ! -f ${IDIR}/${S}_3T_rfMRI_REST2_preproc/${S}/MNINonLinear/Results/rfMRI_REST2_RL/rfMRI_REST2_LR_smoothed_HP_detrended.nii.gz ] && \
   [ ! -f ${IDIR}/${S}_3T_rfMRI_REST1_preproc/${S}/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_RL_smoothed_HP_detrended.nii.gz ] && \
   [ ! -f ${IDIR}/${S}_3T_rfMRI_REST2_preproc/${S}/MNINonLinear/Results/rfMRI_REST2_RL/rfMRI_REST2_RL_smoothed_HP_detrended.nii.gz ]; the

    for RUN in $RUNS; do

        for D in $DIRECTIONS; do

            x=$(fslstats ${IDIR}/${S}_3T_rfMRI_REST${RUN}_preproc/${S}/MNINonLinear/Results/rfMRI_REST${RUN}_${D}/rfMRI_REST${RUN}_${D}.nii.gz -p 2 -p 98)

            x=$(echo $x | awk '{print $2}')

            stat1=$(bc <<< "${x} * 0.10")

            fslmaths ${IDIR}/${S}_3T_rfMRI_REST${RUN}_preproc/${S}/MNINonLinear/Results/rfMRI_REST${RUN}_${D}/rfMRI_REST${RUN}_${D}.nii.gz -thr ${stat1} -Tmin -bin \
                     ${IDIR}/${S}_3T_rfMRI_REST${RUN}_preproc/${S}/MNINonLinear/Results/rfMRI_REST${RUN}_${D}/pre_thr_REST${RUN}_${D}_mask.nii.gz -odt char

            stat2=$(fslstats ${IDIR}/${S}_3T_rfMRI_REST${RUN}_preproc/${S}/MNINonLinear/Results/rfMRI_REST${RUN}_${D}/rfMRI_REST${RUN}_${D}.nii.gz \
                     -k ${IDIR}/${S}_3T_rfMRI_REST${RUN}_preproc/${S}/MNINonLinear/Results/rfMRI_REST${RUN}_${D}/pre_thr_REST${RUN}_${D}_mask.nii.gz -p 50)

            stat3=$(bc <<< "${stat2} * 0.75")

            fslmaths ${IDIR}/${S}_3T_rfMRI_REST${RUN}_preproc/${S}/MNINonLinear/Results/rfMRI_REST${RUN}_${D}/rfMRI_REST${RUN}_${D}.nii.gz \
                    -mas ${IDIR}/${S}_3T_rfMRI_REST${RUN}_preproc/${S}/MNINonLinear/Results/rfMRI_REST${RUN}_${D}/pre_thr_REST${RUN}_${D}_mask.nii.gz \
                    ${IDIR}/${S}_3T_rfMRI_REST${RUN}_preproc/${S}/MNINonLinear/Results/rfMRI_REST${RUN}_${D}/pre_thr_REST${RUN}_${D}_func.nii.gz

            fslmaths ${IDIR}/${S}_3T_rfMRI_REST${RUN}_preproc/${S}/MNINonLinear/Results/rfMRI_REST${RUN}_${D}/pre_thr_REST${RUN}_${D}_func.nii.gz \
                    -Tmean ${IDIR}/${S}_3T_rfMRI_REST${RUN}_preproc/${S}/MNINonLinear/Results/rfMRI_REST${RUN}_${D}/pre_thr_REST${RUN}_${D}_func_mean.nii.gz

            susan ${IDIR}/${S}_3T_rfMRI_REST${RUN}_preproc/${S}/MNINonLinear/Results/rfMRI_REST${RUN}_${D}/pre_thr_REST${RUN}_${D}_func.nii.gz ${stat3} 1.698514 3 1 1 \
                  ${IDIR}/${S}_3T_rfMRI_REST${RUN}_preproc/${S}/MNINonLinear/Results/rfMRI_REST${RUN}_${D}/pre_thr_REST${RUN}_${D}_func_mean.nii.gz \
                  ${stat2} ${HDIR}/preprocessed/${S}_3T_rfMRI_REST${RUN}_preproc/${S}/MNINonLinear/Results/rfMRI_REST${RUN}_${D}/rfMRI_REST${RUN}_${D}_smoothed.nii.gz

            3dTproject -input ${IDIR}/${S}_3T_rfMRI_REST${RUN}_preproc/${S}/MNINonLinear/Results/rfMRI_REST${RUN}_${D}/rfMRI_REST${RUN}_${D}_smoothed.nii.gz \
                       -prefix ${IDIR}/${S}_3T_rfMRI_REST${RUN}_preproc/${S}/MNINonLinear/Results/rfMRI_REST${RUN}_${D}/rfMRI_REST${RUN}_${D}_smoothed_HP_detrended.nii.gz \
                       -polort 1 \
                       -stopband 0 0.01 \
                       -automask \
                       -overwrite

        done

    done

else

	rm -r ${IDIR}/${S}_3T_rfMRI_REST1_preproc
	rm -r ${IDIR}/${S}_3T_rfMRI_REST2_preproc

fi

