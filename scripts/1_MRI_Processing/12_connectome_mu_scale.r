rdir="/work/barbey/rwilcox5/projects/HCP"

subjects<-list.files(paste0(rdir,"/data/tractography"))
subjects<-subjects[grep("^[0-9]",subjects)]

for(sub in subjects) {
  if(file.exists(paste0(rdir,"/data/tractography/",sub,"/MMP1_structural_connectome.csv"))) {
    print(paste("Running", sub))
    mat<-read.csv(paste0(rdir,"/data/tractography/",sub,"/MMP1_structural_connectome.csv"),header=F)
    mu<-read.table(text=readLines(paste0(rdir,"/data/tractography/",sub,"/sift_mu.txt"),warn=F),header=F)
    mat_scaled<-mat*mu$V1
    write.table(mat_scaled,paste0(rdir,"/data/tractography/",sub,"/MMP1_structural_connectome_mu_scaled.csv"),row.names=F,col.names=F,sep=",")
  }
}
  
