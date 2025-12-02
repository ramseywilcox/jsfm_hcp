library(data.table)
library(dplyr)
library(tidyr)

# --- Set working directory
rdir='/home/rwilcox5/Projects/HCP'

n_folds=5 # # of folds predictive connections need to have survived

# --- Load and setup data
dist_data<-as.data.frame(fread(paste0(rdir,"/data/tractography/master_files/master_structural_distome.csv")))
haufe_data<-as.data.frame(fread(paste0(rdir,"/results/haufe_coefficients.csv")))
sig_feats_data<-as.data.frame(fread(paste0(rdir,"/results/stable_significant_features_folds-",n_folds,".csv")))
sig_feats<-sig_feats_data$Var1

# --- Subset features
haufe_data<-haufe_data[which(haufe_data$Feats %in% sig_feats),]
haufe_data$iter<-NULL
dist_data<-dist_data[,which(names(dist_data) %in% sig_feats)]

# --- Summarize features
haufe_summary <- haufe_data %>% group_by(Feats) %>% summarise(haufe_mean = mean(Haufe, na.rm=T))
haufe_summary <- as.data.frame(haufe_summary)
dist_data_long <- gather(key = Feats, value = Values, dist_data)
dist_summary <- dist_data_long %>% group_by(Feats) %>% summarise(dist_mean = mean(Values, na.rm=T))
dist_summary <- as.data.frame(dist_summary)

# --- Setup for model
data <- merge(dist_summary, haufe_summary, by = "Feats")
data$sign[data$haufe_mean > 0] <- 0
data$sign[data$haufe_mean < 0] <- 1
data$sign<-factor(data$sign, levels = c(0,1))
data_pos <- data[which(data$haufe_mean > 0),]
data_neg <- data[which(data$haufe_mean < 0),]

# --- Model 1
mod1<-glm(sign ~ dist_mean, data=data, family = "binomial")
summary(mod1)
wilcox.test(data$dist_mean[which(data$sign==1)],data$dist_mean[which(data$sign==0)])
mean(data$dist_mean[which(data$sign==1)])
mean(data$dist_mean[which(data$sign==0)])
sd(data$dist_mean[which(data$sign==1)]) / sqrt(length(data$dist_mean[which(data$sign==1)]))
sd(data$dist_mean[which(data$sign==0)]) / sqrt(length(data$dist_mean[which(data$sign==0)]))

# --- Model 2
mod2<-glm(haufe_mean ~ dist_mean*sign, data=data, family="gaussian")
summary(mod2)
cor.test(data_pos$dist_mean,data_pos$haufe_mean)
cor.test(data_neg$dist_mean,data_neg$haufe_mean)
