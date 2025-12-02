#!/bin/bash

HCP=$WORK/projects/HCP

subjs=( $( cat ${HCP}/subject_list.txt ) )

len=$(expr ${#subjs[@]} )
echo " Spawning ${#subjs[@]} "

proc=$WORK/projects/HCP/scripts/2_Joint_Structure_Function_Modeling/BNIF/utils/single_subj_BNIF.sh

sbatch --array=0-$len $proc ${subjs[@]}
