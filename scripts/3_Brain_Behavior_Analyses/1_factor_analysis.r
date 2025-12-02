library(psych)
library(lavaan)

# --- Set function
cfi_bf <- function(x){
  return((1-((x$stats$STATISTIC-x$stats$dof))/(x$stats$null.chisq-x$stats$null.dof)))
}

# --- Set working directory
rdir<-getwd()

# --- Load and setup data
data_beh<-read.csv(paste0(rdir,"/data/behavioral/unrestricted_data.csv"))

cols_keep<-c("Subject","PicVocab_Unadj","ReadEng_Unadj","CardSort_Unadj",
             "Flanker_Unadj","ProcSpeed_Unadj","VSPLOT_TC","PMAT24_A_CR",
             "PicSeq_Unadj","ListSort_Unadj","IWRD_TOT","PMAT_Compl",
             "NEO.FFI_Compl","Non.TB_Compl","VisProc_Compl","SCPT_Compl",
             "IWRD_Compl","VSPLOT_Compl","MMSE_Compl","MMSE_Score",
             "X3T_RS.fMRI_PctCompl","X3T_dMRI_PctCompl","NEORAW_01")

cog_tests<-c("PicVocab_Unadj","ReadEng_Unadj","PicSeq_Unadj","Flanker_Unadj",
             "CardSort_Unadj","ProcSpeed_Unadj","PMAT24_A_CR","VSPLOT_TC","IWRD_TOT","ListSort_Unadj")

data_beh<-data_beh[,which(names(data_beh) %in% cols_keep),drop=F]

# --- Remove subjects with missing neuropsych data
data_beh<-data_beh[which(data_beh$PMAT_Compl=="true" &
                           data_beh$NEO.FFI_Compl=="true" &
                           data_beh$Non.TB_Compl=="true" &
                           data_beh$VisProc_Compl=="true" &
                           data_beh$SCPT_Compl=="true" &
                           data_beh$IWRD_Compl=="true" &
                           data_beh$VSPLOT_Compl=="true" &
                           data_beh$MMSE_Compl=="true"),]

# --- Remove subjects with low MMST
data_beh<-data_beh[which(data_beh$MMSE_Score > 26),,drop=F]

# --- Remove subjects with incomplete behavioral data
data_beh<-data_beh[complete.cases(data_beh),,drop=F]

# --- Get only variables to be analyzed
cogdf<-data_beh[,which(names(data_beh) %in%  cog_tests),drop=F]
cogdf<-as.data.frame(lapply(cogdf, scale)) # scale variables

# --- Run Horn's parallel analysis
fa.parallel(cogdf,plot=F,fa="fa")

# --- Run factor analysis
bfm <- omega(cogdf,nfactors=4,fm="minres",sl=TRUE,n.obs=NA,rotate="oblimin",Phi = NULL,option="equal",covar=FALSE)

# --- Get fit indices
CFI<-cfi_bf(bfm)
RMSEA<-bfm$schmid$RMSEA[1]
fit_inds<-data.frame(CFI=CFI,RMSEA=RMSEA)
write.csv(fit_inds,paste0(rdir,"/results/bifactor_model_fit_indices.csv"),row.names=F)

# --- Get scores and save
bfm_scores<-as.data.frame(factor.scores(cogdf,bfm$schmid$sl[,1:5])$scores)
bfm_scores<-cbind(data_beh$Subject,bfm_scores)
names(bfm_scores)<-c("Subject","g","cry","spd","vis","mem")
write.csv(bfm_scores,paste0(rdir,"/data/behavioral/factor_scores.csv"),row.names = F)

g_scores<-bfm_scores[,c(1,2),drop=F]
write.csv(g_scores,paste0(rdir,"/data/behavioral/g_scores.csv"),row.names=F)
