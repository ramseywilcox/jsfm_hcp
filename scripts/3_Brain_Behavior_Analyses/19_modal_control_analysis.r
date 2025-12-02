library(ggplot2)
library(ggrepel)

# --- Set function
mdmr <- function(D, X, n_perm = 1000, seed = 123) {
  set.seed(seed)
  
  n <- nrow(D)
  
  # double-center the distance matrix
  D2 <- as.matrix(D)^2
  H_center <- diag(n) - matrix(1, n, n) / n
  G <- -0.5 * H_center %*% D2 %*% H_center  # Gower's centered matrix
  
  # construct model matrix
  X_full <- model.matrix(~ ., data = as.data.frame(X))
  Q <- ncol(X_full)
  
  # calculate hat matrix for full model
  H <- X_full %*% solve(t(X_full) %*% X_full) %*% t(X_full)
  
  # calculate residual projection matrix
  I_H <- diag(n) - H
  
  # calculate pseudo F-statistic for the full model
  tr_HG <- sum(diag(H %*% G))
  tr_E <- sum(diag(I_H %*% G))
  F_obs <- (tr_HG / (Q - 1)) / (tr_E / (n - Q))
  
  # permutation test
  F_perm <- numeric(n_perm)
  for (i in 1:n_perm) {
    print(i)
    perm_idx <- sample(n)
    G_perm <- G[perm_idx, perm_idx]
    tr_HG_perm <- sum(diag(H %*% G_perm))
    tr_E_perm <- sum(diag(I_H %*% G_perm))
    F_perm[i] <- (tr_HG_perm / (Q - 1)) / (tr_E_perm / (n - Q))
  }
  
  # compute p-value
  p_val <- (sum(F_perm >= F_obs)+1) / n_perm

  return(list(
    F_stat = F_obs,
    p_value = p_val,
    F_distribution = F_perm
  ))
}

# --- Set working directory
rdir<-getwd()

# --- Load and setup data
data<-read.csv(paste0(rdir,"/data/graph_metrics/modal_control_continuous_structural_connectome_mu_scaled.csv"))
g_scores<-read.csv(paste0(rdir,"/data/behavioral/g_scores.csv"))
comm<-read.csv(paste0(rdir,"/data/parcellations/cortex_parcel_network_assignments.txt"),header=F)
comm$Node<-paste0("Node_",seq(1,360,by=1))
names(comm)<-c("Network","Node")

# --- Get top 75% nodes
data_means<-data
data_means$Subject<-NULL
data_means<-as.data.frame(colMeans(data_means))
data_means$Node<-row.names(data_means)
names(data_means)<-c("mean","Node")
q<-quantile(data_means$mean, probs = seq(0,1,0.25))
top_nodes<-data_means$Node[which(data_means$mean > q[[2]])]
data<-data[,which(names(data) %in% c("Subject",top_nodes))]
data<-merge(g_scores,data,by="Subject",all=F)
data$Subject<-NULL

# --- Run MDMR
data_x<-data[,"g",drop=F]
data_y<-data[,-which(names(data) %in% "g")]
dist_mat<-as.matrix( as.dist(1 - cor(t(data_y), method="pearson")) )
mod<-mdmr(dist_mat, data_x,n_perm=1000)
