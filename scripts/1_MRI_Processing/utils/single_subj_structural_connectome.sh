#!/bin/bash
#
#SBATCH --time=10:59:00
#SBATCH --mem=32GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -J run_conn
#SBATCH --output /work/out_batch/out_conn-%A_%a.out
#SBATCH --error /work/out_batch/conn-%A_%a.err
#

module purge
module load fsl
module load mrtrix3

subjs=($@)
subjs=$1
S=${subjs[${SLURM_ARRAY_TASK_ID}]}

HDIR=$WORK/projects/HCP
ODIR=$HDIR/data/tractography
DWDIR=$HDIR/data/preprocessed/dwi
PDIR=$HDIR/data/parcellations

if [ -f ${ODIR}/${S}/tracks_50M.tck ]; then

    flirt -in ${PDIR}/${S}/MMP1_T1w.nii.gz -ref ${DWDIR}/${S}/T1w/T1w_acpc_dc_restore_1.25.nii.gz \
        -out ${PDIR}/${S}/MMP1_T1w_dwi.nii.gz -dof 6 -interp nearestneighbour

    mrconvert ${PDIR}/${S}/MMP1_T1w_dwi.nii.gz ${PDIR}/${S}/MMP1_T1w_dwi.mif -force

    tck2connectome -symmetric -zero_diagonal \
        -tck_weights_in ${ODIR}/${S}/sift_weights.txt \
        ${ODIR}/${S}/tracks_50M.tck \
        ${PDIR}/${S}/MMP1_T1w_dwi.mif \
        ${ODIR}/${S}/MMP1_structural_connectome.csv \
        -out_assignment ${ODIR}/${S}/MMP1_assignments_structural_connectome.csv \
        -force

fi
