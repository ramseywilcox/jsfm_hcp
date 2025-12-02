# --- Set working directory
rdir='/home/rwilcox5/Projects/HCP'

# --- Set # of permutations
n_perms<-1000

# --- Load and setup data
g_scores<-read.csv(paste0(rdir,"/data/behavioral/g_scores.csv"))
good_subs<-read.csv(paste0(rdir,"/data/mri_qaqc/good_subjects.csv"))
g_scores<-g_scores[which(g_scores$Subject %in% good_subs$Subject),]
g_scores<-g_scores$g

# --- Generate permutations
set.seed(123)
permutations <- as.data.frame(replicate(n_perms, sample(g_scores), simplify=T))
names(permutations) <- paste0("perm_",1:n_perms)
write.csv(permutations,paste0(rdir,"/results/permutation/g_scores_permuted.csv"),row.names=F)
