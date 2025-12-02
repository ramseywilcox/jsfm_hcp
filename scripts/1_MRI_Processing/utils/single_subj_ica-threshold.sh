#!/bin/bash
#
#SBATCH --time=1-10:59:00
#SBATCH --mem=64GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -J run_ica
#SBATCH --output /work/out_batch/out_HCPica_thresh-%A_%a.out
#SBATCH --error /work/out_batch/HCPica_thresh-%A_%a.err
#

module purge
module load fsl

RUNS=$(seq 1 2)
DIRECTIONS=("LR RL")

HDIR=$WORK/projects/HCP
ODIR=$HDIR/data/ICA

SUBJS=($@)
SUBJS=$1
S=${SUBJS[${SLURM_ARRAY_TASK_ID}]}

if [ -f ${ODIR}/${S}/melodic_out_REST1_LR/melodic_IC.nii.gz ] && \
   [ -f ${ODIR}/${S}/melodic_out_REST1_RL/melodic_IC.nii.gz ] && \
   [ -f ${ODIR}/${S}/melodic_out_REST2_LR/melodic_IC.nii.gz ] && \
   [ -f ${ODIR}/${S}/melodic_out_REST2_RL/melodic_IC.nii.gz ]; then

    for RUN in $RUNS; do

        for D in $DIRECTIONS; do

			if [ ! -f ${ODIR}/${S}/melodic_out_REST${RUN}_${D}/melodic_IC_thresh50.nii.gz ]; then

						MAPLIST=$( \ls ${ODIR}/${S}/melodic_out_REST${RUN}_${D}/stats | grep 'probmap*' )

						for MAP in $MAPLIST; do
							echo " running ${MAP} "
							fslmaths ${ODIR}/${S}/melodic_out_REST${RUN}_${D}/stats/${MAP} -thr 0.5 -bin ${ODIR}/${S}/melodic_out_REST${RUN}_${D}/stats/Thresh50_${MAP}
						done

						mkdir ${ODIR}/${S}/melodic_out_REST${RUN}_${D}/temp

						NCOMP=$(echo ${MAPLIST} | wc -w)

						for COMP in $(seq 1 ${NCOMP}); do
							echo " running ${COMP} "
							ROIIDX=$(bc <<< "${COMP} - 1")
							fslroi ${ODIR}/${S}/melodic_out_REST${RUN}_${D}/melodic_IC.nii.gz ${ODIR}/${S}/melodic_out_REST${RUN}_${D}/temp/comp${COMP}.nii.gz ${ROIIDX} 1
							fslmaths ${ODIR}/${S}/melodic_out_REST${RUN}_${D}/temp/comp${COMP}.nii.gz -mas ${ODIR}/${S}/melodic_out_REST${RUN}_${D}/stats/Thresh50_probmap_${COMP}.nii.gz ${ODIR}/${S}/melodic_out_REST${RUN}_${D}/temp/comp${COMP}_thresh50.nii.gz
						done

						for COMP in $(seq 1 ${NCOMP}); do
							echo ${ODIR}/${S}/melodic_out_REST${RUN}_${D}/temp/comp${COMP}_thresh50.nii.gz >> ${ODIR}/${S}/melodic_out_REST${RUN}_${D}/temp/file_array.csv
						done

						fslmerge -t ${ODIR}/${S}/melodic_out_REST${RUN}_${D}/melodic_IC_thresh50.nii.gz `cat ${ODIR}/${S}/melodic_out_REST${RUN}_${D}/temp/file_array.csv`

						rm -r ${ODIR}/${S}/melodic_out_REST${RUN}_${D}/temp

			fi

        done

    done

fi
