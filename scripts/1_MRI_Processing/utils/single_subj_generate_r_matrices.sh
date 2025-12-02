#!/bin/bash
#
#SBATCH --time=10:59:00
#SBATCH --mem=32GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -J run_r_matrix
#SBATCH --output /work/out_batch/out_r_matrix-%A_%a.out
#SBATCH --error /work/out_batch/r_matrix-%A_%a.err
#

module purge
module load anaconda

conda activate $WORK/mri_conda_env

SESSIONS=$(seq 1 2)
DIRECTIONS=("LR RL")

HDIR=$WORK/projects/HCP
IDIR=${HDIR}/data/ICA

subjs=($@)
subjs=$1
S=${subjs[${SLURM_ARRAY_TASK_ID}]}

for SES in $SESSIONS; do

    for DIR in $DIRECTIONS; do

        if [ ! -f ${IDIR}/${S}/melodic_out_REST${SES}_${DIR}/R_matrix_sum.csv ] &&
           [ -f ${IDIR}/${S}/melodic_out_REST${SES}_${DIR}/melodic_IC_classified.nii.gz ]; then

            python ${HDIR}/scripts/1_MRI_Processing/utils/generate_r_matrices.py --working_directory=${HDIR} \
                --subject=${S} \
                --session=${SES} \
                --direction=${DIR}

        fi

    done

done
