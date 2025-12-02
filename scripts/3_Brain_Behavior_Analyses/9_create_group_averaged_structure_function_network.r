library(data.table)

# --- Set working directory
rdir=getwd()

# --- Get subject list
subjects<-read.csv(paste0(rdir,"/data/mri_qaqc/good_subjects.csv"))
subjects<-subjects$Subject

# --- Create group averaged file and save
mat_list<-list()
i=1
for(sub in subjects) {
  print(i)
  mat<-read.csv(paste0(rdir,"/data/jsfm/",sub,"/MMP1_maxed_averaged_structure_function_network.csv"),header=F,sep=" ")
  mat_list[[i]]<-mat
  i<-i+1
}

matrix<-Reduce("+",mat_list) / length(mat_list)
fwrite(matrix,paste0(rdir,"/results/MMP1_maxed_averaged_structure_function_network_GROUP_AVERAGED.csv"),row.names=F,col.names=F,sep=",")