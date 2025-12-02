# --- Set working directory
rdir<-getwd()

# --- Read in data and merge
data<-read.csv(paste0(rdir,"/data/graph_metrics/master_small_worldness_telesford.csv"))
g_scores<-read.csv(paste0(rdir,"/data/behavioral/g_scores.csv"))
data<-merge(g_scores,data,by="Subject",all=F)

# --- Fix skewness
data$C_real <- log(data$C_real)

# --- Run correlations
cor.test(data$g,data$Sigma)
cor.test(data$g,data$C_real)
cor.test(data$g,data$L_real)