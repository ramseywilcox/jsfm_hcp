#!/bin/bash

HCP=$WORK/projects/HCP

subjs=( $( cat ${HCP}/subject_list2.txt ) )

len=$(expr ${#subjs[@]} )
echo " Spawning ${#subjs[@]} "

proc=$HCP/scripts/1_MRI_Processing/utils/single_subj_parcellation_registration.sh

sbatch --array=0-$len $proc ${subjs[@]}
