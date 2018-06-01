library(sleuth)
library(data.table)
library(corrplot)
library(psych)
library(ggplot2)

base_dir <- "~/ncbi/kallisto_b"
sample_id <- dir(file.path(base_dir))
s2c <- read.table(file.path("~/ncbi/info.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample, condition)                 
kal_dirs <- file.path(base_dir, sample_id)
s2c <- dplyr::mutate(s2c, path = kal_dirs)

so <- sleuth_prep(s2c, extra_bootstrap_summary = T)

so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')

# модель - likehood ratio test
so <- sleuth_lrt(so, 'reduced', 'full')

models(so)

# таблица наиболее диф экспресс транспазонов (20 шт)
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval < 0.05)
head(sleuth_significant, 20)

# сколько значимых и незначимых транспазонов
table(sleuth_table[,"qval"] < 0.05)

# на сколько сигм отличается диф экспрессия значимых транспазонов
mean(sleuth_significant$test_stat[1:20] / sd(sleuth_significant$test_stat))

# стат тесты, достоверно ли отличается диф экспрессия 20 топовых транспазонов по сравнению со всеми
shapiro.test(sleuth_significant$test_stat)
wilcox.test(sleuth_significant$test_stat, (sleuth_significant$test_stat[1:20]))
# запись в файл
write.csv(sleuth_significant, file = "~/ncbi/sleuth_res/sleuth_significant.csv")
# представление результатов в виде таблицы
kallisto_table(so, use_filtered=T, normalized=T, include_covariates=T)

# представление результатов в виде удобной матрицы
matrix = sleuth_to_matrix(so, "obs_norm", "tpm")
write.csv(matrix, file = "~/ncbi/sleuth_res/res_matrix.csv")


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

# корреляционная матрица
cor_matrix = filtered.cor(data2)

# другой корреляционный тест - узнаем p-value
r = corr.test(data1, method = "spearman" )
corr.p(r$r, 8)

# визуализация матрицы
corrplot.mixed(cor_matrix)

# запись матрицы в файл
write.csv(cor_matrix, file = "~/ncbi/sleuth_res/cor_matrix.csv")

# боксплот по данной транспазоне по всем условиям

plot_bootstrap(so, "rnd-4_family-47#LINE/Penelope", units = "est_counts", color_by = "condition", divide_groups=F ) + theme_bw() 


# файл с учетом бутстрепа. Максимальное значение, минимальное, среднее и тп. Поэтому оценка общей экспрессии производится с поправкой на значения каждого бутстрепа
penelope_bootstrap <- get_bootstrap_summary(so, "rnd-4_family-47#LINE/Penelope", units = "est_counts")
write.csv(penelope_bootstrap, file = "~/ncbi/sleuth_res/penelope_bootstrap.csv")

# heatmap для 20 топовых транспазонов
transcripts = sleuth_significant$target_id[1:20]
plot_transcript_heatmap(so, transcripts, units = "tpm", trans = "log")

# pca

p <- plot_pca(so, pc_x=1L, pc_y=2L, units="est_counts",text_labels = T, color_by = "condition", use_filtered = T) + theme_bw()
  
  
  
#text_high <- textGrob("PC2=35%", gp=gpar(fontsize=13, fontface="bold"))
#text_low <- textGrob("PC1=55%", gp=gpar(fontsize=13, fontface="bold"))
  
#+theme(axis.text.y=element_text(size=10),
#        axis.text.x = element_text(size=10),
#        legend.position = "bottom",
#        legend.text = element_text(size=10),
#        legend.title = element_blank(), 
#        panel.background = element_rect(color='black', fill='white'),
#        panel.grid.major = element_line(size = 0.5, linetype="dotted", colour = "gray"), 
#        panel.grid.minor = element_line(size = 0.25, linetype = 'dotted', colour = "gray"))+annotation_custom(text_low,xmin=100000,xmax=5,ymin=-1,ymax=-50000)+annotation_custom(text_high,xmin=100000,xmax=5,ymin=-1,ymax=-60000)

# c tpm и est_counts PCA немного разный

gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)

# сколько процентов занимает каждая компонента
plot_pc_variance(so,units='est_counts') + theme_bw()

# какие транспазоны участвуют в формировании каждой компоненты
plot_loadings(so, use_filtered = TRUE, sample = NULL, pc_input = "PC1",
              units = "est_counts", pc_count=5, scale = F,
              pca_loading_abs = TRUE) + theme_bw() + theme(axis.text.x=element_text(angle=70, hjust=1))


# извлекаем значения PCA
mat = sleuth:::spread_abundance_by(
  abund = so$obs_norm_filt,
  var = "est_counts",
  which_order = so$sample_to_covariates$sample)
pca_res <- prcomp(mat,scale. = F,center = F)

# то же самое 
pcs <- sleuth:::as_df(pca_res$rotation[,c(1L,2L,3L)])
pcs$sample <- rownames(pcs)


# значения tpm и est_counts для каждого транспазона
tmp <- so$obs_raw %>% dplyr::filter(target_id == "rnd-4_family-47#LINE/Penelope")
tmp <- dplyr::full_join(so$sample_to_covariates, tmp, by = 'sample')
tmp
