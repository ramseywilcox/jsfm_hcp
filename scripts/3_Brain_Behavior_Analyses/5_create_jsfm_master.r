library(tidyr)
library(data.table)

# --- Set working directory
rdir<-getwd()

# --- Get subject list
subjects<-read.csv(paste0(rdir,"/data/mri_qaqc/good_subjects.csv"))
subjects<-subjects$Subject

# --- Create JSFM master file and save
data_list=list()
i=1
for(sub in subjects) {
  if(file.exists(paste0(rdir,"/data/jsfm/",sub,"/MMP1_maxed_averaged_structure_function_network.csv"))) {
    print(i)
    tmp<-read.table(paste0(rdir,"/data/jsfm/",sub,"/MMP1_maxed_averaged_structure_function_network.csv"),header=F,sep=" ")
    num_feats<-(nrow(tmp) * (nrow(tmp)-1)) / 2
    upper_tri<-data.frame(Subject=sub,Headers=paste0("V",seq(1,num_feats,by=1)),Values=tmp[upper.tri(tmp)])
    data_list[[i]]<-upper_tri
    i=i+1
  }
}

data<-do.call(rbind,data_list)
data<-spread(key=Headers,value=Values,data=data)
fwrite(data,paste0(rdir,"/data/jsfm/master_files/master_maxed_averaged_structure_function_network.csv"),row.names=F, sep=",")