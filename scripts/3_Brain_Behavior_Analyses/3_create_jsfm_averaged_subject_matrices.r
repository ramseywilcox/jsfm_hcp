library(tidyr)
library(data.table)

# --- Set working directory
rdir<-getwd()

# --- Get subject list
subjects<-list.files(paste0(rdir,"/data/jsfm"))

# --- Set sessions and runs
sessions=c(1,2)
runs=c("LR","RL")

# --- Create averaged subject files and save
for(sub in subjects) {
  matrix_list=list()
  if(file.exists(paste0(rdir,"/data/jsfm/",sub,"/REST1_LR/MMP1_maxed_structure_function_network.csv")) &&
     file.exists(paste0(rdir,"/data/jsfm/",sub,"/REST1_RL/MMP1_maxed_structure_function_network.csv")) &&
     file.exists(paste0(rdir,"/data/jsfm/",sub,"/REST2_LR/MMP1_maxed_structure_function_network.csv")) &&
     file.exists(paste0(rdir,"/data/jsfm/",sub,"/REST2_RL/MMP1_maxed_structure_function_network.csv"))) {
    print(sub)
    i=1
    sub_num_zeros<-data.frame()
    for(ses in sessions) {
      for(run in runs) {
        tmp<-read.table(paste0(rdir,"/data/jsfm/",sub,"/REST",ses,"_",run,"/MMP1_maxed_structure_function_network.csv"),header=F,sep=",")
        matrix_list[[count]]<-tmp
        i<-i+1
      }
    }
    matrix<-Reduce("+",matrix_list) / length(matrix_list)
    write.table(matrix,paste0(rdir,"/data/jsfm/",sub,"/MMP1_maxed_averaged_structure_function_network.csv"),row.names=F,col.names=F)
  }
}
