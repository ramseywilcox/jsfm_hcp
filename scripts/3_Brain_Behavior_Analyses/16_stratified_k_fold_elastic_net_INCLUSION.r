library(caret)
library(e1071)
library(data.table)

# --- Set functions
remove_skewed_columns <- function(df, threshold = 2, type=3) {
  skew_vals <- sapply(df, e1071::skewness, na.rm = TRUE, type=type)
  keep_cols <- names(skew_vals[abs(skew_vals) <= threshold])
  return(keep_cols)
}

cv_r_squared_from_folds <- function(obs, pred, fold) {
  
  fold_ids <- unique(fold)
  ss_res_total <- 0
  ss_tot_total <- 0
  
  for (f in fold_ids) {
    test_idx <- which(fold == f)
    train_idx <- which(fold != f)
    
    obs_fold <- obs[test_idx]
    pred_fold <- pred[test_idx]
    train_mean <- mean(obs[train_idx])
    
    ss_res <- sum((obs_fold - pred_fold)^2)
    ss_tot <- sum((obs_fold - train_mean)^2)
    
    ss_res_total <- ss_res_total + ss_res
    ss_tot_total <- ss_tot_total + ss_tot
  }
  
  r2_cv <- 1 - (ss_res_total / ss_tot_total)
  return(r2_cv)
}

# --- Set working directory
rdir=getwd()

# --- Load and setup data
g_scores<-read.csv(paste0(rdir,"/data/behavioral/g_scores.csv"))
data<-fread(paste0(rdir,"/data/jsfm_sum/master_files/master_maxed_averaged_structure_function_network.csv"))
data<-merge(g_scores,data,by="Subject",all=F)
data$Subject<-NULL
net_map<-read.csv(paste0(rdir,"/data/parcellations/network_connection_header_map.csv"))

# --- Set networks
networks <- c("V1","V2","SMN","CON","DAN","LAN","FPN","AUD","DMN","PMM","VMM","OA")

# --- Remove skewed features
keep_cols<-remove_skewed_columns(data[,-1],threshold=3,type=3)
data<-data[,which(names(data) %in% c("g",keep_cols)),drop=F]

# --- Set model parameters
tune.grid<-expand.grid(
  alpha=0.001,
  lambda=10^seq(-4, 1, length.out=100)
)

alpha<-0.1 # p-value threshold for feature seleection

n_bins<-4 # number of bins for stratified sampling

k=5 # number of outer folds

inner_k=3 # number of inner folds


# --- Run model
set.seed(1)
performance_metrics<-data.frame()
strata<-cut(data$g, breaks = quantile(data$g, probs = seq(0,1,length.out = n_bins+1)),
            include.lowest=T, labels=F)
folds<-createFolds(strata, k=k, list=T, returnTrain=F)
for(net in networks) {
  preds<-data.frame()
  print(paste0("Running ",net))
  cols_keep<-net_map$header[which(net_map$network1==net & net_map$network2==net)]
  data_sub<-data[,which(names(data) %in% c("g",cols_keep))]
  
  for(i in seq_along(folds)) {
    skip_outer<-FALSE
    
    print(paste0("Running fold ",i))
    fold=folds[[i]]
    train<-data_sub[-fold,]
    test<-data_sub[fold,]
    test_obv<-test$g
    
    cor_df<-data.frame()
    for(z in 2:length(train)) {
      r<-cor.test(train$g,train[,z],alternative="two.sided",method="pearson")
      r_p<-r$p.value
      if(r_p<=alpha && !is.na(r_p)) {
        cor_df_tmp<-data.frame(node=names(train)[z],R=r$estimate,r_p=r_p)
        cor_df<-rbind(cor_df,cor_df_tmp)
      }
    }
    if(nrow(cor_df) >= 2) {  
      train<-train[,which(names(train) %in% c("g",cor_df$node))]
      test<-test[,which(names(test) %in% cor_df$node)]
      
      inner_strata<-cut(train$g, breaks = quantile(train$g, probs = seq(0,1,length.out = n_bins+1)),
                        include.lowest=T, labels=F)
      inner_folds_test_tmp<-createFolds(inner_strata,k=inner_k,list=T,returnTrain=F)
      inner_folds_test<-list()
      for(z in seq_along(inner_folds_test_tmp)) {
        inner_folds_test[[z]]<-inner_folds_test_tmp[[z]]
      }
      
      inner_folds_train<-list()
      idx_seq<-seq(1,nrow(train))
      for(z in seq_along(inner_folds_test)){
        inner_fold_test<-inner_folds_test[[z]]
        train_idx<-idx_seq[!(idx_seq %in% inner_fold_test)]
        inner_folds_train[[z]]<-train_idx
      }
      means<-as.numeric(as.data.frame(lapply(train[,-1],mean)))
      stds<-as.numeric(as.data.frame(lapply(train[,-1],sd)))
      train[,-1]<-sweep(train[,-1],2,means,"-")
      train[,-1]<-sweep(train[,-1],2,stds,"/")
      test<-sweep(test,2,means,"-")
      test<-sweep(test,2,stds,"/")
      
      model<-train(
        g ~ ., data=train, method="glmnet",
        trControl=trainControl(index=inner_folds_train,indexOut=inner_folds_test),
        tuneGrid=tune.grid,
        metric="RMSE"
      )
      prediction<-model %>% predict(test)
      preds_tmp<-data.frame(predictions=prediction,observations=test_obv)
      preds<-rbind(preds,preds_tmp)
      
    } else {
      skip_outer<-TRUE
      break
    }
  }
  if(skip_outer) {
    next
  }
  write.csv(preds,paste0(rdir,"/results/lesion/inclusion/predictions_network",net,".csv"),row.names=F)
  
  r = cor(preds$observations,preds$predictions,method="pearson")
  r2 = cv_r_squared_from_folds(preds$observations,preds$predictions,fold_info$count)
  nRMSD=sqrt(1-r2)
  tmp_metrics<-data.frame(network=net,r=r,r2=r2,nRMSD=nRMSD)
  performance_metrics<-rbind(performance_metrics,tmp_metrics)
}

# --- Save performance metrics
write.csv(performance_metrics,paste0(rdir,"/results/lesion/inclusion/performance_metrics.csv"),row.names=F)

# --- Save fold info
fold_info<-data.frame()
for(z in seq_along(folds)) {
  fold<-folds[[z]]
  tmp<-data.frame(idx=fold,count=z)
  fold_info<-rbind(fold_info,tmp)
}
fold_info<-fold_info[order(fold_info$idx,decreasing = F),]
write.csv(fold_info,paste0(rdir,"/results/lesion/inclusion/fold_info.csv"),row.names=F)
