#!/bin/bash
#
#SBATCH --time=1-10:59:00
#SBATCH --mem=32GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -J run_trax
#SBATCH --output /work/out_batch/out_trax-%A_%a.out
#SBATCH --error /work/out_batch/trax-%A_%a.err
#

module purge
module load mrtrix3

subjs=($@)
subjs=$1
S=${subjs[${SLURM_ARRAY_TASK_ID}]}

HDIR=$WORK/projects/HCP
ADIR=$HDIR/data/preprocessed/anat
ODIR=$HDIR/data/tractography

if [ -f ${ODIR}/${S}/wmfod_norm.mif ]; then

    tckgen \
        -force \
        -algorithm iFOD2 \
        -select 50M \
        -cutoff 0.06 \
        -minlength 5.0 \
        -maxlength 275.0 \
        -act ${ODIR}/${S}/5TT.mif -crop_at_gmwmi -backtrack \
        -max_attempts_per_seed 1000 \
        -seed_dynamic ${ODIR}/${S}/wmfod_norm.mif \
        -output_seeds ${ODIR}/${S}/seeds.txt \
        ${ODIR}/${S}/wmfod_norm.mif ${ODIR}/${S}/tracks_50M.tck -nthreads 32 -force

    tckedit ${ODIR}/${S}/tracks_50M.tck -number 200k ${ODIR}/${S}/tracks_200k.tck

    #mrview ${ODIR}/${S}/DWI_biascorrected.mif -tractography.load ${ODIR}/${S}/tracks_200k.tck

    tcksift2 -act ${ODIR}/${S}/5TT.mif \
            -out_mu ${ODIR}/${S}/sift_mu.txt \
            -out_coeffs ${ODIR}/${S}/sift_coeffs.txt \
            -nthreads 32 \
            ${ODIR}/${S}/tracks_50M.tck \
            ${ODIR}/${S}/wmfod_norm.mif \
            ${ODIR}/${S}/sift_weights.txt -force

fi
