#!/bin/bash
#
#SBATCH --time=10:59:00
#SBATCH --mem=24GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -J run_dist
#SBATCH --output /work/out_batch/out_dist-%A_%a.out
#SBATCH --error /work/out_batch/dist-%A_%a.err
#

module purge
module load mrtrix3

subjs=($@)
subjs=$1
S=${subjs[${SLURM_ARRAY_TASK_ID}]}

HDIR=$WORK/projects/HCP
ODIR=$HDIR/data/tractography
PDIR=$HDIR/data/parcellations

if [ ! -f ${ODIR}/${S}/MMP1_structural_distome.csv ]; then

    tck2connectome -symmetric -zero_diagonal \
        ${ODIR}/${S}/tracks_50M.tck \
        ${PDIR}/${S}/MMP1_T1w_dwi.mif \
        ${ODIR}/${S}/MMP1_structural_distome.csv \
        -scale_length \
        -stat_edge mean \
        -force

fi
