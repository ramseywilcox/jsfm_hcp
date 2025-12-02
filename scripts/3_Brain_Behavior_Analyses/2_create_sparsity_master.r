library(tidyr)
library(data.table)

# --- Set working directoy
rdir<-getwd()

# --- Get subject list
subjects<-list.files(paste0(rdir,"/data/jsfm"))

# --- Set sessions and runs
sessions=c(1,2)
runs=c("LR","RL")

# --- Create sparsity master file and save
group_data<-list()
i=1
for(sub in subjects) {
  matrix_list=list()
  if(file.exists(paste0(rdir,"/data/jsfm/",sub,"/REST1_LR/MMP1_maxed_structure_function_network.csv")) &&
     file.exists(paste0(rdir,"/data/jsfm/",sub,"/REST1_RL/MMP1_maxed_structure_function_network.csv")) &&
     file.exists(paste0(rdir,"/data/jsfm/",sub,"/REST2_LR/MMP1_maxed_structure_function_network.csv")) &&
     file.exists(paste0(rdir,"/data/jsfm/",sub,"/REST2_RL/MMP1_maxed_structure_function_network.csv"))) {
    print(sub)
    sub_data<-data.frame()
    for(ses in sessions) {
      for(run in runs) {
        tmp<-read.table(paste0(rdir,"/data/jsfm/",sub,"/REST",ses,"_",run,"/MMP1_maxed_structure_function_network.csv"),header=F,sep=",")
        n_zeros<-table(tmp==0)[[2]]-360
        sub_tmp<-data.frame(Subject=sub,Ses=ses,Run=run,Sparsity=n_zeros)
        sub_data<-rbind(sub_data,sub_tmp)
      }
    }
    group_data[[i]]<-sub_data
    i=i+1
  }
}

data<-do.call(rbind,group_data)
data$Sparsity_scaled<-scale(data$Sparsity)
fwrite(data,paste0(rdir,"/data/mri_qaqc/master_sparsity_data.csv"),row.names=F,sep=",")
