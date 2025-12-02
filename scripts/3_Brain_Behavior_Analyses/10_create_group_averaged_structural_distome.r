library(data.table)
rdir=getwd()

# --- Get subject list
subjects<-read.csv(paste0(rdir,"/data/mri_qaqc/good_subjects.csv"))
subjects<-subjects$Subject

# --- Create group averaged file and save
mat_list<-list()
i=1
for(sub in subjects) {
  mat<-read.csv(paste0(rdir,"/data/tractography/",sub,"/MMP1_structural_distome.csv"),header=F,sep=",")
  mat_list[[i]]<-mat
  i<-i+1
}

matrix<-Reduce("+",mat_list) / length(mat_list)
fwrite(matrix,paste0(rdir,"/results/MMP1_structural_distome_GROUP_AVERAGED.csv"),row.names=F,col.names=F,sep=",")