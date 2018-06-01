# This file is written based on TEdiff.R, part of TEtools (https://github.com/l-modolo/TEtools). 
# You can see the source code here: https://github.com/l-modolo/TEtools/blob/master/TEdiff.R
#
# TEtools is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# TEtools is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with TEtools.  If not, see <http://www.gnu.org/licenses/>.


# parse input arguments
args = commandArgs(trailingOnly=FALSE)
scriptPath <- dirname(sub("--file=","",args[grep("--file",args)]))
scriptPath = paste0(scriptPath,"/TEdiff")
print(scriptPath)
.libPaths(c(.libPaths(), scriptPath))

args = commandArgs(trailingOnly=TRUE)
for(i in 2:length(args))
{
    eval(parse(text=sub("--", "", args[i])))
}

help = FALSE
if(!exists("FDR_level")){help = TRUE; print("FDR_level")}
if(!exists("count_column")){help = TRUE; print("count_column")}
if(!exists("count_file")){help = TRUE; print("count_file")}
if(!exists("experiment_formula")){help = TRUE; print("experiment_formula")}
if(!exists("sample_names")){help = TRUE; print("sample_names")}

if (!is.null(outdir)) {
  # Create the directory where we put all the html support files
  dir.create(outdir,FALSE)
}

output_figures = list()

if(help==TRUE)
{
    print("TEdiff.R --args --FDR_level=0.05 --count_column=2 --count_file=\\\"count.txt\\\" --experiment_formula=\\\"population:type\\\" --sample_names=\\\"population1:type1,population1:type2,population2:type1,population2:type2\\\"" )
    quit("no")
}
if(exists("version"))
{
    print("1.0.0")
}
count_column = count_column-1
print(FDR_level)
print(count_column)
print(count_file)
print(experiment_formula)
print(sample_names)

#count_file = "all_same_3.csv"
#count_column = 2
#experiment_formula = "sample:replicant"
#sample_names = "adult:1,adult:2,egg:1,egg:2,juv:1,juv:2,met:1,met:2"
#outdir = "~/ncbi/DeSeq2_res"


# we get the counts
counts=read.csv(count_file, header=F, sep=' ')
rownames(counts) = counts[,1]
counts_information = counts[, c(1:count_column) ]
counts_total = counts[, dim(counts)[2]]
counts = counts[, -c(1:count_column, dim(counts)[2]) ]

# we get the variables
variable_names = strsplit(experiment_formula, "[:]")[[1]]
sample_names = strsplit(sample_names,",")[[1]]
variables = strsplit(sample_names[1], "[:]")[[1]]
variable_number = length(variables)

for( i in c(2:length(sample_names)))
{
    variable = strsplit(sample_names[i], "[:]")[[1]]
    variables = cbind(variables, variable)
}
variables = t(variables)
rownames(variables) = sample_names
colnames(variables) =  variable_names
variables = as.data.frame(variables)

names(counts) = sample_names
counts = as.data.frame(counts)
counts = counts[rowSums(counts)!=0,]

print(head(counts))
print(variables)

# we run DeSeq
suppressMessages(require(genefilter, quietly = TRUE))
suppressMessages(require(DESeq2, quietly = TRUE))
suppressMessages(require(gplots, quietly = TRUE))
suppressMessages(require(ggplot2, quietly = TRUE))
suppressMessages(require(RColorBrewer, quietly = TRUE))

# LRT test
TE = DESeqDataSetFromMatrix(countData = counts,
                            colData = variables,
                            design = formula(paste0("~", variable_names[1]))
                            )

TE = DESeq(TE, test='LRT', full=~sample, reduced=~1)

dds <- estimateSizeFactors(TE)
dds <- estimateDispersions(dds)
dds <- nbinomLRT(dds, full=~sample, reduced = ~ 1)
res <- results(dds)
res_sorted <- res[order(-res$stat),] 
# list of TEs sorted by significance
write.csv(res_sorted, file = paste0(outdir, "/dds_sorted", ".csv"))

mcols(res, use.names=TRUE)

# PCA 
ntop = 500
rld = rlogTransformation(TE, blind=T)
rv = rowVars(assay(rld))
select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
pca = prcomp(t(assay(rld)[select, ]))


fac = factor(apply(as.data.frame(colData(rld)[, variable_names, drop = FALSE]),
    1, paste, collapse = " : "))
x = ggplot(data=as.data.frame(pca$x), aes(PC1, PC2, color=fac)) +
        geom_point(size = 6) +
        xlab(paste0("PC1: ",100*summary(pca)[6]$importance[2,][1],"% of variance")) +
        ylab(paste0("PC1: ",100*summary(pca)[6]$importance[2,][2],"% of variance")) +
        theme_bw() +
        guides(color=guide_legend(title="stages"))+
  theme(axis.text.x = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),
        axis.text.y = element_text(colour="black",size=15,angle=0,hjust=1,vjust=0,face="plain"),
        axis.title.x = element_text(colour="black",size=20),
        axis.title.y = element_text(colour="black",size=20),
        legend.text=element_text(size=20))
# }
ggsave(file=paste0(outdir, "/PCA.pdf"),x, width=20, height=20, units="in", dpi=1200)
ggsave(file=paste0(outdir, "/PCA.png"),x , width = 8, height = 8, units="in", dpi=100)

# PCA loadings
loadings <- as.data.frame(pca$rotation)
aload <- abs(pca$rotation)

# proportions of loadings
proportions <- as.data.frame(sweep(aload, 2, colSums(aload), "/"))

# matrix of TEs sorted by their effect on PC1
prop_pc1 <- proportions[order(-proportions$PC1),]
pc1 <- prop_pc1['PC1']
write.csv(pc1, file = paste0(outdir, "/pc1", ".csv"))

# plot of PCA1 loadings
ten_pc1 <- as.data.frame(pc1[1:10,1])
ten_pc1$TEs <- rownames(pc1)[1:10]
legend_ord <- levels(with(ten_pc1, reorder(TEs, -ten_pc1$`pc1[1:10, 1]`)))

x <-  ggplot(ten_pc1, aes(y = ten_pc1$`pc1[1:10, 1]`, x = reorder(TEs, -ten_pc1$`pc1[1:10, 1]`), fill=TEs))+
  geom_bar(stat = 'identity', width = .4,  position=position_dodge())+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=1,vjust=0,face="plain"))+
  scale_fill_discrete(breaks=legend_ord)+
  labs(title = "PC1 loadings", y = "proportional contribution")
ggsave(file=paste0(outdir, "/plot_pca1_load", "~", variable_names[1], ".pdf"), x, width = 20, height = 20)
ggsave(file=paste0(outdir, "/plot_pca1_load", "~", variable_names[1], ".png"), x, width = 8, height = 8, units="in", dpi=100)

x <- barplot(prop_pc1[, 1], col = 'pink')
# ggsave(file=paste0(outdir, "/pca1_100", "~", variable_names[1], ".pdf"), x)
# ggsave(file=paste0(outdir, "/pca1_100", "~", variable_names[1], ".png"), x)

# matrix of TEs sorted by their effect on PC2
prop_pc2 <- proportions[order(-proportions$PC2),]
pc2 <- prop_pc2['PC2']
write.csv(pc2, file = paste0(outdir, "/pc2", ".csv"))
ten_pc2 <- as.data.frame(pc2[1:10, 1])
ten_pc2$TEs <- rownames(pc2)[1:10]
legend_ord <- levels(with(ten_pc2, reorder(TEs, -ten_pc2$`pc2[1:10, 1]`)))

# plot of PCA2 loadings
x <-  ggplot(ten_pc2, aes(y = ten_pc2$`pc2[1:10, 1]`, x = reorder(TEs, -ten_pc2$`pc2[1:10, 1]`), fill=TEs))+
  geom_bar(stat = 'identity', width = .4,  position=position_dodge())+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=1,vjust=0,face="plain"))+
  scale_fill_discrete(breaks=legend_ord)+
  labs(title = "PC2 loadings", y = "proportional contribution")
ggsave(file=paste0(outdir, "/plot_pca2_load", "~", variable_names[1], ".pdf"), x, width = 20, height = 20)
ggsave(file=paste0(outdir, "/plot_pca2_load", "~", variable_names[1], ".png"), x, width = 8, height = 8, units="in", dpi=100)

x <- barplot(prop_pc2[, 2], col = 'pink')



# euclidean distance between samples
TE_vsd = varianceStabilizingTransformation(TE)
TE_row = order(rowMeans(counts(TE,normalized=TRUE)),decreasing=TRUE)
old_i = 1
output_figures_heatmap = c()
TE_step = min(length(TE_row), 30)
for(i in seq(from=TE_step, to=length(TE_row), by = TE_step))
{
    select = order(rownames(TE),decreasing=FALSE)[old_i:i]
    hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
    old_i = i
}

distsRL = dist(t(assay(rld)))
mat <- as.matrix(distsRL)
hc <- hclust(distsRL)

pdf(paste0(outdir, "/Sample-to-sample distances~", variable_names[1], ".pdf") , height=20,width=20)
heatmap.2(mat, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col = rev(hmcol), margin=c(13, 13))
x = dev.off()
png(paste0(outdir, "/Sample-to-sample distances~", variable_names[1], ".png") , height=800,width=800)
heatmap.2(mat, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col = rev(hmcol), margin=c(13, 13), cexRow=2.5, cexCol=2.5)
x = dev.off()












