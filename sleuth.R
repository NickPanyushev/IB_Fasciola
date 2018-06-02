library(sleuth)
library(data.table)
library(corrplot)
library(psych)
library(ggplot2)

# path to kallisto data
base_dir <- "~/ncbi/kallisto_b"
sample_id <- dir(file.path(base_dir))
s2c <- read.table(file.path("~/ncbi/info.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample, condition)                 
kal_dirs <- file.path(base_dir, sample_id)
s2c <- dplyr::mutate(s2c, path = kal_dirs)

# create sleuth model 
so <- sleuth_prep(s2c, extra_bootstrap_summary = T)

so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')

# model - likehood ratio test
so <- sleuth_lrt(so, 'reduced', 'full')

# show our model
models(so)

# matrix of differentially expressed TEs
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval < 0.05)
head(sleuth_significant, 20)

# significant vs non significant TEs
table(sleuth_table[,"qval"] < 0.05)

# how different is the expression of the most significant TEs
mean(sleuth_significant$test_stat[1:20] / sd(sleuth_significant$test_stat))

# some statistical tests 
shapiro.test(sleuth_significant$test_stat)
wilcox.test(sleuth_significant$test_stat, (sleuth_significant$test_stat[1:20]))
# write results
write.csv(sleuth_significant, file = "~/ncbi/sleuth_res/sleuth_significant.csv")
# table of results
kallisto_table(so, use_filtered=T, normalized=T, include_covariates=T)

# matrix of results
matrix = sleuth_to_matrix(so, "obs_norm", "tpm")
write.csv(matrix, file = "~/ncbi/sleuth_res/res_matrix.csv")

# correlation matrix
filtered.cor <- function(x){    
  num_var <- sapply(x, function(x) is.numeric(x))    
  cor_mat <- cor(x[, num_var], method = 'spearman')   
  diag(cor_mat) <- 1    
  return(cor_mat)}
data2 = data1[c(1,5,2,8,7,4,3,6)]
str(data1)
data1 = as.data.frame(matrix)
data2 = data1[c(3,6,7,4,2,8,1,5)]
setnames(data2, old=c( "data.ad1","data.ad2",  "data.juv1","data.juv2" ,"data.met1", "data.met2", "data.egg1" ,"data.egg2"), new=c( "ad1","ad2",  "juv1","juv2" ,"met1", "met2", "egg1" ,"egg2"))


cor_matrix = filtered.cor(data2)

# p-value of correlation
r = corr.test(data1, method = "spearman" )
corr.p(r$r, 8)

# vizualization of correlation matrix
corrplot.mixed(cor_matrix)

write.csv(cor_matrix, file = "~/ncbi/sleuth_res/cor_matrix.csv")

# boxplot of expression of particular TE
plot_bootstrap(so, "rnd-4_family-47#LINE/Penelope", units = "est_counts", color_by = "condition", divide_groups=F ) + theme_bw() 


# bootstrap values used in boxplots 
penelope_bootstrap <- get_bootstrap_summary(so, "rnd-4_family-47#LINE/Penelope", units = "est_counts")
write.csv(penelope_bootstrap, file = "~/ncbi/sleuth_res/penelope_bootstrap.csv")

# heatmap of 20 most significant TEs
transcripts = sleuth_significant$target_id[1:20]
plot_transcript_heatmap(so, transcripts, units = "tpm", trans = "log")+theme(axis.text.x = element_text(size=20))

# pca
p <- plot_pca(so, pc_x=1L, pc_y=2L, units="est_counts",text_labels = T, color_by = "condition", use_filtered = T) + theme_bw()
  
# % of variance in PC
plot_pc_variance(so,units='est_counts') + theme_bw()

# loadings of PC
plot_loadings(so, use_filtered = TRUE, sample = NULL, pc_input = "PC1",
              units = "est_counts", pc_count=5, scale = F,
              pca_loading_abs = TRUE) + theme_bw() + theme(axis.text.x=element_text(angle=70, hjust=1))


#  PCA values
mat = sleuth:::spread_abundance_by(
  abund = so$obs_norm_filt,
  var = "est_counts",
  which_order = so$sample_to_covariates$sample)
pca_res <- prcomp(mat,scale. = F,center = F)


#  tpm and est_counts for a particular TE
tmp <- so$obs_raw %>% dplyr::filter(target_id == "rnd-4_family-47#LINE/Penelope")
tmp <- dplyr::full_join(so$sample_to_covariates, tmp, by = 'sample')
tmp
