library(tidyr)

# --- Set working directory
rdir=getwd()

# --- Set threshold parameter
n_folds=5

# --- Map connections and save
feats<-read.csv(paste0(rdir,"/results/stable_significant_features_folds-",n_folds,".csv"))
map<-read.csv(paste0(rdir,"/data/parcellations/network_connection_header_map.csv"))

sig_conns<-merge(feats,map,by.x="Var1",by.y="header",all=F)
sig_conns<-sig_conns[,'connection',drop=F]
sig_conns<-separate(sig_conns,col=connection,into = c("Node1","Node2"))
sig_conns$Node1<-as.numeric(sig_conns$Node1)
sig_conns$Node2<-as.numeric(sig_conns$Node2)

write.table(sig_conns,paste0(rdir,"/results/stable_significant_connections_folds-",n_folds,".csv"),row.names=F,col.names=F)
