HDIR=$WORK/projects/HCP
ODIR=$HDIR/data/tractography

module load mrtrix3tissue/5.2

responsemean ${ODIR}/*/wm_response.txt ${ODIR}/wm_group_average_response.txt
responsemean ${ODIR}/*/gm_response.txt ${ODIR}/gm_group_average_response.txt
responsemean ${ODIR}/*/csf_response.txt ${ODIR}/csf_group_average_response.txt
