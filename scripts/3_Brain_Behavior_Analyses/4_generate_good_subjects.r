# --- Set working directory
rdir<-getwd()

# --- Set movement and sparsity thresholds
rms_thresh<-0.20
sparsity_thresh<-3.5
mmse_thresh<-26

# --- Load and setup data
data_sparsity<-read.csv(paste0(rdir,"/data/mri_qaqc/master_sparsity_data.csv"))
data_rms<-read.csv(paste0(rdir,"/data/mri_qaqc/resting_state_rms_movement_data.csv"))
data_beh<-read.csv("/home/rwilcox5/Projects/HCP/data/behavioral/unrestricted_data.csv")

cols_keep<-c("Subject","PicVocab_Unadj","ReadEng_Unadj","CardSort_Unadj","Flanker_Unadj","ProcSpeed_Unadj","VSPLOT_TC","PMAT24_A_CR","PicSeq_Unadj","ListSort_Unadj","IWRD_TOT",
             "PMAT_Compl","NEO.FFI_Compl","Non.TB_Compl","VisProc_Compl","SCPT_Compl","IWRD_Compl","VSPLOT_Compl","MMSE_Compl","MMSE_Score","X3T_RS.fMRI_PctCompl","X3T_dMRI_PctCompl","NEORAW_01")
data_beh<-data_beh[,which(names(data_beh) %in% cols_keep),drop=F]

# --- Subset on data completeness
data_beh<-data_beh[which(data_beh$PMAT_Compl=="true" &
                           data_beh$NEO.FFI_Compl=="true" &
                           data_beh$Non.TB_Compl=="true" &
                           data_beh$VisProc_Compl=="true" &
                           data_beh$SCPT_Compl=="true" &
                           data_beh$IWRD_Compl=="true" &
                           data_beh$VSPLOT_Compl=="true" &
                           data_beh$MMSE_Compl=="true" &
                           data_beh$X3T_RS.fMRI_PctCompl==100 &
                           data_beh$X3T_dMRI_PctCompl==100),,drop=F]

# --- Subset on mini-mental state exam
data_beh<-data_beh[which(data_beh$MMSE_Score > mmse_thresh),,drop=F]
data_beh<-data_beh[complete.cases(data_beh),,drop=F]

# --- Get good and bad subjects
bad_subs_rms<-data_rms$Subject[which(data_rms$rms > rms_thresh)]
bad_subs_sparsity<-data_sparsity$Subject[which(data_sparsity$Sparsity_scaled > sparsity_thresh)]
bad_subs<-unique(c(bad_subs_rms,bad_subs_sparsity))
good_subs_rms<-data_rms$Subject[which(!(data_rms$Subject %in% bad_subs_rms))]
good_subs_sparsity<-data_sparsity$Subject[which(!(data_sparsity$Subject %in% bad_subs_sparsity))]
good_subs<-unique(intersect(good_subs_rms,good_subs_sparsity))

## --- Subset and save
good_subjects<-data_beh$Subject[which(data_beh$Subject %in% good_subs)]
good_subjects<-data.frame(Subject=good_subjects)
write.csv(good_subjects,paste0(rdir,"/data/mri_qaqc/good_subjects.csv"),row.names=F)
