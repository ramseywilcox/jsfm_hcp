#!/bin/bash
#
#SBATCH --time=4-10:59:00
#SBATCH --mem=100gb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -J run_bnif
#SBATCH --output /work/out_batch/out_bnif-%A_%a.out
#SBATCH --error /work/out_batch/bnif-%A_%a.err
#

module purge
module load ibm-ilog-cplex/12.10
module load matlab/r2019b

subjs=($@)
subjs=$1
S=${subjs[${SLURM_ARRAY_TASK_ID}]}

HDIR=$WORK/projects/HCP
SDIR=$HDIR/scripts
ICADIR=$HDIR/data/ICA
TDIR=$HDIR/data/tractography
ODIR=$HDIR/data/jsfm

cd $SDIR/BNIF

SESSIONS=$(seq 1 2)
DIRECTIONS=("LR RL")

if [ -f ${ICADIR}/${S}/melodic_out_REST1_LR/R_matrix_sum.csv ] && \
 [ -f ${ICADIR}/${S}/melodic_out_REST1_RL/R_matrix_sum.csv ] && \
 [ -f ${ICADIR}/${S}/melodic_out_REST2_LR/R_matrix_sum.csv ] && \
 [ -f ${ICADIR}/${S}/melodic_out_REST2_RL/R_matrix_sum.csv ]; then

mkdir $ODIR/$S

for SES in $SESSIONS; do
	for DIR in $DIRECTIONS; do
			if [ ! -f ${ODIR}/${S}/REST${SES}_${DIR}/MMP1_maxed_structure_function_network.csv ]; then
    				mkdir $ODIR/$S/REST${SES}_${DIR}
				SODIR=$ODIR/$S/REST${SES}_${DIR}
				SICADIR=$ICADIR/$S/melodic_out_REST${SES}_${DIR}
    				v="run_BNIF('$S','$SICADIR','$TDIR','$SODIR'); exit"
    				matlab -nodisplay -nosplash -nodesktop -r $v
			else
				echo "Subject jsfm network already generated"
			fi
	done
done

fi

