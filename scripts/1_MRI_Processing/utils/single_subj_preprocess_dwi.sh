#!/bin/bash
#
#SBATCH --time=1-10:59:00
#SBATCH --nodes=1
#SBATCH --mem=32GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -J run_DWIprocess
#SBATCH --output /work/out_batch/out_DWIprocess-%A_%a.out
#SBATCH --error /work/out_batch/DWIprocess-%A_%a.err
#

module purge
module load ants
module load mrtrix3
module load freesurfer

subjs=($@)
subjs=$1
S=${subjs[${SLURM_ARRAY_TASK_ID}]}

HDIR=$WORK/projects/HCP
DWDIR=$HDIR/data/preprocessed/dwi
ADIR=$HDIR/data/preprocessed/anat
ODIR=$HDIR/data/tractography

unzip ${DWDIR}/${S}_3T_Diffusion_preproc.zip -d ${DWDIR}/

if [ -f ${DWDIR}/${S}/T1w/Diffusion/data.nii.gz ]; then

    mkdir $ODIR/$S
    cd ${ODIR}/${S}

    mrconvert ${DWDIR}/${S}/T1w/Diffusion/data.nii.gz ${ODIR}/${S}/DWI.mif \
        -fslgrad ${DWDIR}/${S}/T1w/Diffusion/bvecs ${DWDIR}/${S}/T1w/Diffusion/bvals \
        -datatype float32 \
        -stride 0,0,0,1 \
        -info \
        -nthreads 32 \
        -force

    dwibiascorrect ants ${ODIR}/${S}/DWI.mif ${ODIR}/${S}/DWI_biascorrected.mif \
        -bias ${ODIR}/${S}/bias_ants_field.mif \
        -info \
        -nthreads 32 \
        -force

    dwi2response dhollander ${ODIR}/${S}/DWI_biascorrected.mif ${ODIR}/${S}/wm_response.txt ${ODIR}/${S}/gm_response.txt ${ODIR}/${S}/csf_response.txt \
        -voxels ${ODIR}/${S}/RF_voxels.mif \
        -info \
        -nthreads 32 \
        -force

    dwiextract ${ODIR}/${S}/DWI_biascorrected.mif - -bzero | mrmath - mean ${ODIR}/${S}/meanb0.mif \
        -axis 3 \
        -info \
        -nthreads 32 \
        -force

    #mrview ${HDIR}/out/meanb0.mif -overlay.load ${HDIR}/out/RF_voxels.mif -overlay.opacity 0.5 -force -info

    dwi2mask ${ODIR}/${S}/DWI_biascorrected.mif ${ODIR}/${S}/DWI_biascorrected_mask.mif -nthreads 32 -force

    5ttgen freesurfer ${ADIR}/${S}/T1w/aparc+aseg.nii.gz ${ODIR}/${S}/5TT.mif -nthreads 32 -force

    rm -r ${ODIR}/${S}/5ttgen-tmp-* ## sometimes mrtrix3 doesn't delete

else

    rm -r ${DWDIR}/${S}

fi
