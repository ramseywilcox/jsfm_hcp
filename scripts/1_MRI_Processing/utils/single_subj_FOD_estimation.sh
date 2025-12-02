#!/bin/bash
#
#SBATCH --time=1-10:59:00
#SBATCH --nodes=1
#SBATCH --mem=32GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -J run_FOD
#SBATCH --output /work/out_batch/out_FOD-%A_%a.out
#SBATCH --error /work/out_batch/FOD-%A_%a.err
#

module purge
module load mrtrix3

subjs=($@)
subjs=$1
S=${subjs[${SLURM_ARRAY_TASK_ID}]}

HDIR=$WORK/projects/HCP
ODIR=$HDIR/data/tractography

if [ -f ${ODIR}/${S}/DWI_biascorrected.mif ]; then

    cd ${ODIR}/${S}

    dwi2fod msmt_csd ${ODIR}/${S}/DWI_biascorrected.mif \
        ${ODIR}/wm_group_average_response.txt ${ODIR}/${S}/wmfod.mif \
        ${ODIR}/gm_group_average_response.txt ${ODIR}/${S}/gmfod.mif \
        ${ODIR}/csf_group_average_response.txt ${ODIR}/${S}/csffod.mif \
        -mask ${ODIR}/${S}/DWI_biascorrected_mask.mif \
        -nthreads 32 \
        -force

    mtnormalise \
        ${ODIR}/${S}/wmfod.mif ${ODIR}/${S}/wmfod_norm.mif \
        ${ODIR}/${S}/gmfod.mif ${ODIR}/${S}/gmfod_norm.mif \
        ${ODIR}/${S}/csffod.mif ${ODIR}/${S}/csffod_norm.mif \
        -mask ${ODIR}/${S}/DWI_biascorrected_mask.mif \
        -check_norm ${ODIR}/${S}/mtnormalise_norm.mif \
        -check_mask ${ODIR}/${S}/mtnormalise_mask.mif \
        -nthreads 32 \
        -force

    mrconvert -coord 3 0 ${ODIR}/${S}/wmfod_norm.mif - | mrcat ${ODIR}/${S}/csffod_norm.mif ${ODIR}/${S}/gmfod_norm.mif - ${ODIR}/${S}/vf.mif \
        -nthreads 32 \
        -force

fi
