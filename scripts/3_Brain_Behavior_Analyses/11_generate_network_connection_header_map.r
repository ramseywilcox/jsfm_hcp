# --- Set working directory
rdir<-getwd()

# -- Read in network assignment map
net_map<-read.csv(paste0(rdir,"/data/parcellations/cortex_parcel_network_assignments.txt"),header=F)
names(net_map)<-"Network"
net_map$Network[which(net_map$Network==1)]<-"V1"
net_map$Network[which(net_map$Network==2)]<-"V2"
net_map$Network[which(net_map$Network==3)]<-"SMN"
net_map$Network[which(net_map$Network==4)]<-"CON"
net_map$Network[which(net_map$Network==5)]<-"DAN"
net_map$Network[which(net_map$Network==6)]<-"LAN"
net_map$Network[which(net_map$Network==7)]<-"FPN"
net_map$Network[which(net_map$Network==8)]<-"AUD"
net_map$Network[which(net_map$Network==9)]<-"DMN"
net_map$Network[which(net_map$Network==10)]<-"PMM"
net_map$Network[which(net_map$Network==11)]<-"VMM"
net_map$Network[which(net_map$Network==12)]<-"OA"
net_map$Node<-seq(1,360,by=1)

# -- Create connection to header map and save
region_sequence<-rep(1:360,each=360)
row_mat<-matrix(data=region_sequence,nrow=360,ncol=360,byrow=T)
col_mat<-matrix(data=region_sequence,nrow=360,ncol=360,byrow=F)
row_idx<-row_mat[upper.tri(row_mat)]
col_idx<-col_mat[upper.tri(col_mat)]
conn_map<-data.frame(node_1=col_idx,node_2=row_idx,connection=paste(col_idx,row_idx,sep="_"),header=paste0("V",seq(1,length(row_idx),by=1)))

net_conn_map<-merge(net_map,conn_map,by.x="Node",by.y="node_1")
names(net_conn_map)<-c("node1","network1","node2","connection","header")
net_conn_map<-merge(net_map,net_conn_map,by.x="Node",by.y="node2")
names(net_conn_map)<-c("node2","network2","node1","network1","connection","header")
net_conn_map<-net_conn_map[,c(3,4,1,2,5,6)]
net_conn_map$header<-sub("V","",net_conn_map$header)
net_conn_map<-net_conn_map[order(as.numeric(net_conn_map$header)),]
net_conn_map$header<-paste0("V",net_conn_map$header)

write.csv(net_conn_map,paste0(rdir,"/data/parcellations/network_connection_header_map.csv"),row.names = F)
