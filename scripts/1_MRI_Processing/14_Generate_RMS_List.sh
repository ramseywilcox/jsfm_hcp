HDIR=$WORK/projects/HCP
IDIR=$HDIR/data/preprocessed/rest

SESSIONS=$(seq 1 2)
DIRECTIONS=("LR RL")

subjs=( $( cat ${HDIR}/subject_list.txt ) )


for S in ${subjs[@]}; do

    for SES in ${SESSIONS}; do

        for DIR in ${DIRECTIONS}; do

            rms_mean=$(cat ${IDIR}/${S}_3T_rfMRI_REST${SES}_preproc_unzipped/${S}/MNINonLinear/Results/rfMRI_REST${SES}_${DIR}/Movement_RelativeRMS_mean.txt)
            echo ${S},${SES},${DIR},${rms_mean} >> ${HDIR}/data/mri_qaqc/rms_list.csv

        done

    done

done
