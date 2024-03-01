suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(EnhancedVolcano))
suppressPackageStartupMessages(library(apeglm))
suppressPackageStartupMessages(library("IHW"))
suppressPackageStartupMessages(library(hexbin))
suppressPackageStartupMessages(library(xlsx))
suppressPackageStartupMessages(library(ashr))
suppressPackageStartupMessages(library(apeglm))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(optparse))

defaultW <- getOption("warn")
options(warn = -1)

option_list = list(
  make_option(c('-t', '--sampletable'), action='store', help='Sample Table with metadata'),
  make_option(c('-c', '--count'), action='store', type='character', help='Read count for filtering the matrix'),
  make_option(c('-a', '--alpha'), action='store', type='character', help='Alpha threshold used to extract the result table from a DESeq analysis'),
  make_option(c('-H', '--Height'), action='store', type='integer', help='Plots height'),
  make_option(c('-w', '--width'), action='store', type='integer', help='Plots width'),
  make_option(c('-o', '--output'), action='store', type='character', help='Output folder')
)

opt = parse_args(OptionParser(option_list=option_list))

## set working directory

#setwd("count/")

if (dir.exists(opt$output)){

	files <- list.files(opt$output)

	if (length(files) > 0){

		stop("[Error] Results Directory is not empty specify another one")

	}

} else {

	dir.create(opt$output)

}


indir <- file.path(getwd(),"count")

#preprocessing
sampleTable <- read.table(opt$sampletable, sep="\t", header=TRUE)
sampleTable$sampleName <- factor(sampleTable$sampleName)
sampleTable$condition <- factor(sampleTable$condition)
sampleTable$group <- factor(sampleTable$group) 

#load 
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = indir, design = ~ condition)
ddsHTSeq <- estimateSizeFactors(ddsHTSeq)

## Filter low counts gene, filter out genes (row) with no count
smallestGroupSize <- min(c(length(which(sampleTable$condition == "TUMOR")),length(which(sampleTable$condition == "CONTROL"))))
keep <- rowSums(counts(ddsHTSeq) >= opt$count) >= smallestGroupSize
dds <- ddsHTSeq[keep,]
dds$condition <- relevel(dds$condition, ref = "CONTROL")

# Visualize Boxplot of normalized counts
pdf(file.path(opt$output,"boxplot.counts.pdf"), height=opt$Height, width=opt$width)
rs = rowSums(counts(dds))
par(mfrow=c(1,2)) # plots two plots 
boxplot(log2(counts(dds)[rs > 0,] + 1), cex.lab = 2, cex.axis = 2,xaxt = "n")# not normalized
tick <- seq_along(as.character(dds$group))
axis(1, at = tick, labels = F)
text(tick, par("usr")[3] - 2.5, as.character(dds$group), srt = 90, xpd = T, offset =5)
boxplot(log2(counts(dds, normalized=TRUE)[rs > 0,] + 1), cex.lab = 2, cex.axis = 2,xaxt = "n")# not normalized
tick <- seq_along(as.character(dds$group))
axis(1, at = tick, labels = F)
text(tick, par("usr")[3] - 2.5, as.character(dds$group), srt = 90, xpd = T, offset =5)
dev.off()

# Transform the raw count data

vst <- vst(dds, blind = TRUE)
rld <- rlog(dds, blind = TRUE)

write.table(assay(rld), file.path(opt$output,"normalizedCountRlog.txt"), quote = FALSE, sep = "\t", row.names = TRUE)
write.csv(assay(rld), file.path(opt$output,"normalizedCountRlog.csv"), quote = FALSE, row.names = TRUE)
write.table(assay(vst), file.path(opt$output,"normalizedCountVst.txt"), quote = FALSE, sep = "\t")
write.csv(assay(vst), file.path(opt$output,"normalizedCountVst.csv"), quote = FALSE, row.names = TRUE)

## Scatter plot sample 1 2 
df <- bind_rows(
	as_data_frame(log2(counts(dds, normalized=TRUE)[, c(1,2)]+1)) %>%
		mutate(transformation = "log2(x + 1)"),
	as_data_frame(assay(vst)[, c(1,2)]) %>% mutate(transformation = "vst"),
	as_data_frame(assay(rld)[, c(1,2)]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

p <- ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) + coord_fixed() + facet_grid( . ~ transformation)  

ggsave(file.path(opt$output,"transformedcounts.2samples.pdf"), height=opt$Height, width=opt$width)

# Effects of transformations on the variance
pdf(file.path(opt$output, "transformationvsvariance.vst.pdf"), height=opt$Height, width=opt$width)
meanSdPlot(assay(vst))
dev.off()

pdf(file.path(opt$output, "transformationvsvariance.rlog.pdf"), height=opt$Height, width=opt$width)
meanSdPlot(assay(rld))
dev.off()

## sample distance rlog
sampleDists <- dist(t(assay(rld)))

## Plot Heatmap with sample-to-sample distance
sampleDistMatrix <- as.matrix( sampleDists )
annotation <- data.frame(Phase=as.character(dds$condition))
rownames(annotation) <- as.character(dds$group)

cols <- colorRampPalette(brewer.pal(9, "Set1"))
mycolors <- cols(length(unique(annotation$pts)))
names(mycolors) <- unique(annotation$pts)
annotation_colors = list(
	Phase = c(TUMOR="red", CONTROL="green"))

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)


pdf(file.path(opt$output,"heatmap_dist_sample.rlog.pdf"), height=opt$Height, width=opt$width)
pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors,annotation = annotation,annotation_colors = annotation_colors)
dev.off()

## sample distance vst
sampleDists <- dist(t(assay(vst)))

## Plot Heatmap with sample-to-sample distance
sampleDistMatrix <- as.matrix( sampleDists )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf(file.path(opt$output,"heatmap_dist_sample.vst.pdf"), height=opt$Height, width=opt$width)
pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors,annotation=annotation,annotation_colors = annotation_colors)
dev.off()

## Plot PCA in order to see sample-to-sample distance

## PCA
rld_pca <- plotPCA(rld, intgroup = c("condition"), returnData = TRUE)
percentVar <- round(100 * attr(rld_pca, "percentVar"))
pca_rld <- ggplot(rld_pca, aes(x = PC1, y = PC2, color = condition)) +
	geom_point(size =2) +
	xlab(paste0("PC1: ", percentVar[1], "% variance")) +
	ylab(paste0("PC2: ", percentVar[2], "% variance")) +
	coord_fixed() +
	ggtitle("PCA") + 	
	geom_label_repel(aes(label = name),
									 box.padding   = 0.25, 
									 point.padding = 0.35,
									 segment.color = 'grey50',
									 max.overlaps = 60) 

ggsave(pca_rld, filename = file.path(opt$output,"2DPCA.rlog.pdf"),height=opt$Height, width=opt$width)


## PCA
vst_pca <- plotPCA(vst, intgroup = c("condition"), returnData = TRUE)
percentVar <- round(100 * attr(vst_pca, "percentVar"))
vst_pca$MF <- annotation$MF
pca_vst <- ggplot(vst_pca, aes(x = PC1, y = PC2,color = condition)) +
	geom_point(size =2) +
	xlab(paste0("PC1: ", percentVar[1], "% variance")) +
	ylab(paste0("PC2: ", percentVar[2], "% variance")) +
	coord_fixed() +
	ggtitle("PCA") +
	geom_label_repel(aes(label = name),
									 box.padding   = 0.25, 
									 point.padding = 0.35,
									 segment.color = 'grey50',
									 max.overlaps = 60) 
ggsave(pca_vst, filename = file.path(opt$output,"2DPCA.vst.pdf"),height=opt$Height, width=opt$width)

## Differential Expression Analysis 
design(dds) <- ~ condition
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "TUMOR","CONTROL"), alpha = opt$alpha)

## export MA
pdf(file.path(opt$output,"MA_plot.pdf"), height = opt$Height, width = opt$width)
plotMA(res, ylim = c(-10,10), main = "MA plot")
dev.off()

## export dispersion plot plot that is important to evaluate when performing QC on RNA-seq data is the plot of dispersion versus the mean of normalized counts. 
#For a good dataset, we expect the dispersion to decrease as the mean of normalized counts increases for each gene. 
pdf(file.path(opt$output,"dispersion_plot.pdf"), height = opt$Height, width = opt$width)
plotDispEsts(dds, ylim = c(1e-6, 1e1))
dev.off()

pdf(paste0(file.path(opt$output,"hist_pvalue.padj.pdf")), height = opt$Height, width = opt$width)
par(mfrow = c(1,2))
hist(res$padj, breaks=20, col="grey")
hist(res$pvalue, breaks=20, col="grey")
dev.off()


ResAnnotation <- function (restable){

	#AnnotationDBI
	ens.str <- substr(rownames(restable), 1, 15)
	restable$symbol <- mapIds(org.Hs.eg.db,
													keys=rownames(restable),
													column="SYMBOL",
													keytype="ENSEMBL",
													multiVals = "first")
	restable$entrez <- mapIds(org.Hs.eg.db,
													keys=rownames(restable),
													column="ENTREZID",
													keytype="ENSEMBL")
	restable$alias <- mapIds(org.Hs.eg.db,
												 keys=rownames(restable),
												 column="ALIAS",
												 keytype="ENSEMBL")
	#Annotation with biomatr
	#ensembl <- useEnsembl(biomart = "genes",dataset = "hsapiens_gene_ensembl")
	#restable$ensembl <- sapply( strsplit( rownames(restable), split="\\+" ), "[", 1 )
	#ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl")
	#genemap <- getBM( attributes = c("ensembl_gene_id", "external_gene_name"),
	#									filters = "ensembl_gene_id",
	#									values = restable$ensembl,
	#									mart = ensembl,
	#									useCache = FALSE)
	#idx <- match(restable$ensembl, genemap$ensembl_gene_id, genemap$external_gene_name)
	#restable$ext <- genemap$external_gene_name[ idx ]


	return(restable)

}

res <- ResAnnotation(res)

pdf(file.path(opt$output,"VolcanoPlot.results.pdf"), height = opt$Height, width = opt$width)
EnhancedVolcano(res,
								lab = res$symbol,
								x = 'log2FoldChange',
								y = 'padj',
								FCcutoff = 1.5,
								pCutoff = 0.05,
								title = "Volcano Plot")
dev.off()



## lfcShrink function to shrink the log2 fold changes for the comparison of dex Tumor vs Control samples
resLFC <- lfcShrink(dds, coef="condition_TUMOR_vs_CONTROL", type="apeglm", res = res)
resNorm <- lfcShrink(dds, coef="condition_TUMOR_vs_CONTROL", type="normal", res = res)
resAsh <- lfcShrink(dds, coef="condition_TUMOR_vs_CONTROL", type="ashr", res = res)

#Independent hypothesis weighting
resIHW <- results(dds, filterFun=ihw, alpha = 0.05)
metadata(resIHW)$ihwResult

## Export MA plot
pdf(file.path(opt$output,"multipleMAplot.pdf"),height = opt$Height, width = opt$width)
par(mfrow=c(2,2), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-10,10)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
plotMA(resIHW, xlim=xlim, ylim=ylim, main="Independent Hypotesis Weighting")
dev.off()

#resApeglm
resApeglm <- ResAnnotation(resLFC)

pdf(file.path(opt$output,"VolcanoPlot.apeglm.pdf"), height = opt$Height, width = opt$width)
EnhancedVolcano(resApeglm,
								lab = resApeglm$symbol,
								x = 'log2FoldChange',
								y = 'padj',
								FCcutoff = 1.5,
								pCutoff = 0.05,
								title = "Volcano Plot")
dev.off()


#resNorm
resNorm <- ResAnnotation(resNorm)

pdf(file.path(opt$output,"VolcanoPlot.Normal.pdf"), height = opt$Height, width = opt$width)
EnhancedVolcano(resNorm,
								lab = resNorm$symbol,
								x = 'log2FoldChange',
								y = 'padj',
								FCcutoff = 1.5,
								pCutoff = 0.05,
								title = "Volcano Plot")
dev.off()


#resAsh
resAsh <- ResAnnotation(resAsh)

pdf(file.path(opt$output,"VolcanoPlot.Normal.pdf"), height = opt$Height, width = opt$width)
EnhancedVolcano(resAsh,
								lab = resAsh$symbol,
								x = 'log2FoldChange',
								y = 'padj',
								FCcutoff = 1.5,
								pCutoff = 0.05,
								title = "Volcano Plot")
dev.off()


## two option for colorPalette
heat_color <-  brewer.pal(6, "YlOrRd")
redgreen <- c("green", "black", "red") 
pal <- colorRampPalette(redgreen)(500)

## Gene clustering top 20 genes vst
topVarGenes20 <- head(order(rowVars(assay(vst)), decreasing = TRUE), 20)
mat20  <- assay(vst)[ topVarGenes20, ]

#export
pdf(file.path(opt$output,"heatmap.top20.vst.pdf"),height = opt$Height, width = opt$width)
pheatmap(mat20, annotation = annotation, annotation_colors = annotation_colors, color = pal, scale = "row")
dev.off()

df_20 <- as.data.frame(res[which(rownames(res) %in% rownames(mat20)),])
write.table(df_20, file.path(opt$output,"top20deg.vst.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
write.csv(df_20, file.path(opt$output,"top20deg.vst.csv"), row.names = TRUE, quote = FALSE)
write.xlsx(df_20, file.path(opt$output,"top20deg.vst.xlsx"), row.names = TRUE)

## Gene clustering top 50 genes vst
topVarGenes50 <- head(order(rowVars(assay(vst)), decreasing = TRUE), 50)
mat50  <- assay(vst)[ topVarGenes50, ]

df_50 <- as.data.frame(res[which(rownames(res) %in% rownames(mat50)),])
write.table(df_50, file.path(opt$output,"top50deg.vst.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
write.csv(df_50, file.path(opt$output,"top50deg.vst.csv"), row.names = TRUE, quote = FALSE)
write.xlsx(df_50, file.path(opt$output,"top50deg.vst.xlsx"), row.names = TRUE)

#export
pdf(file.path(opt$output,"heatmap_top50.vst.pdf"),height = opt$Height, width = opt$width)
pheatmap(mat50, annotation = annotation, annotation_colors = annotation_colors, color = pal, scale = "row")
dev.off()

## Gene clustering top 100 genes vst
topVarGenes100 <- head(order(rowVars(assay(vst)), decreasing = TRUE), 100)
mat100  <- assay(vst)[ topVarGenes100, ]

df_100 <- as.data.frame(res[which(rownames(res) %in% rownames(mat100)),])
write.table(df_100, file.path(opt$output,"top100deg.vst.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
write.csv(df_100, file.path(opt$output,"top100deg.vst.csv"), row.names = TRUE, quote = FALSE)
write.xlsx(df_100, file.path(opt$output,"top100deg.vst.xlsx"), row.names = TRUE)

#export
pdf(file.path(opt$output,"heatmap_top100.vst.pdf"),,height = opt$Height, width = opt$width)
pheatmap(mat100, annotation = annotation, annotation_colors = annotation_colors, color = pal, scale = "row")
dev.off()

## Gene clustering top 1000 genes vst
topVarGenes1000 <- head(order(rowVars(assay(vst)), decreasing = TRUE), 1000)
mat1000  <- assay(vst)[ topVarGenes1000, ]

df_1000 <- as.data.frame(res[which(rownames(res) %in% rownames(mat1000)),])
write.table(df_1000, file.path(opt$output,"top1000deg.vst.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
write.csv(df_1000, file.path(opt$output,"top1000deg.vst.csv"), row.names = TRUE, quote = FALSE)
write.xlsx(df_1000, file.path(opt$output,"top1000deg.vst.xlsx"), row.names = TRUE)

#export
pdf(file.path(opt$output,"heatmap_top1000.vst.pdf"), ,height = opt$Height, width = opt$width)
pheatmap(mat1000, annotation = annotation, annotation_colors = annotation_colors, fontsize_row = 8, color = pal,show_rownames = FALSE, scale = "row")
dev.off()

## Gene clustering top 20 genes rld
topVarGenes20 <- head(order(rowVars(assay(rld)), decreasing = TRUE), 20)
mat20  <- assay(rld)[ topVarGenes20, ]

#export
pdf(file.path(opt$output,"heatmap.top20.rlog.pdf"),height = opt$Height, width = opt$width)
pheatmap(mat20, annotation = annotation, annotation_colors = annotation_colors, color = pal, scale = "row")
dev.off()

df_20 <- as.data.frame(res[which(rownames(res) %in% rownames(mat20)),])
write.table(df_20, file.path(opt$output,"top20deg.rlog.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
write.csv(df_20, file.path(opt$output,"top20deg.rlog.csv"), row.names = TRUE, quote = FALSE)
write.xlsx(df_20, file.path(opt$output,"top20deg.rlog.xlsx"), row.names = TRUE)

## Gene clustering top 50 genes rld
topVarGenes50 <- head(order(rowVars(assay(rld)), decreasing = TRUE), 50)
mat50  <- assay(rld)[ topVarGenes50, ]

df_50 <- as.data.frame(res[which(rownames(res) %in% rownames(mat50)),])
write.table(df_50, file.path(opt$output,"top50deg.rlog.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
write.csv(df_50, file.path(opt$output,"top50deg.rlog.csv"), row.names = TRUE, quote = FALSE)
write.xlsx(df_50, file.path(opt$output,"top50deg.rlog.xlsx"), row.names = TRUE)

#export
pdf(file.path(opt$output,"heatmap_top50.rlog.pdf"),height = opt$Height, width = opt$width)
pheatmap(mat50, annotation = annotation, annotation_colors = annotation_colors, color = pal, scale = "row")
dev.off()

## Gene clustering top 100 genes rld
topVarGenes100 <- head(order(rowVars(assay(rld)), decreasing = TRUE), 100)
mat100  <- assay(rld)[ topVarGenes100, ]

df_100 <- as.data.frame(res[which(rownames(res) %in% rownames(mat100)),])
write.table(df_100, file.path(opt$output,"top100deg.rlog.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
write.csv(df_100, file.path(opt$output,"top100deg.rlog.csv"), row.names = TRUE, quote = FALSE)
write.xlsx(df_100, file.path(opt$output,"top100deg.rlog.xlsx"), row.names = TRUE)

#export
pdf(file.path(opt$output,"heatmap_top100.rlog.pdf"),height = opt$Height, width = opt$width)
pheatmap(mat100, annotation = annotation, annotation_colors = annotation_colors, color = pal, scale = "row")
dev.off()

## Gene clustering top 1000 genes rld
topVarGenes1000 <- head(order(rowVars(assay(rld)), decreasing = TRUE), 1000)
mat1000  <- assay(rld)[ topVarGenes1000, ]

df_1000 <- as.data.frame(res[which(rownames(rld) %in% rownames(mat1000)),])
write.table(df_1000, file.path(opt$output,"top1000deg.rlog.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
write.csv(df_1000, file.path(opt$output,"top1000deg.rlog.csv"), row.names = TRUE, quote = FALSE)
write.xlsx(df_1000, file.path(opt$output,"top1000deg.rlog.xlsx"), row.names = TRUE)

#export
pdf(file.path(opt$output,"heatmap_top1000.rlog.pdf"),height = opt$Height, width = opt$width)
pheatmap(mat1000, annotation = annotation, annotation_colors = annotation_colors, fontsize_row = 8, color = pal,show_rownames = FALSE, scale = "row")
dev.off()

#Export the results
res_ <- as.data.frame(res[order(res$padj, decreasing = FALSE),])
write.table(res_, file.path(opt$output, "results.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
write.xlsx(res_, file.path(opt$output,"results.xlsx"),row.names=FALSE)
write.csv(res_, file.path(opt$output,"results.csv"), row.names = FALSE, quote = FALSE)

#Export the results apeglm
resApeglm_ <- as.data.frame(resApeglm[order(resApeglm$padj, decreasing = FALSE),])
write.table(resApeglm_, file.path(opt$output, "results.apeglm.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
write.xlsx(resApeglm_, file.path(opt$output,"results.apeglm.xlsx"),row.names=FALSE)
write.csv(resApeglm_, file.path(opt$output,"results.apeglm.csv"), row.names = FALSE, quote = FALSE)

#Export the results Norm
resNorm_ <- as.data.frame(resNorm[order(resNorm$padj, decreasing = FALSE),])
write.table(resNorm_, file.path(opt$output, "results.norm.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
write.xlsx(resNorm_, file.path(opt$output,"results.norm.xlsx"),row.names=FALSE)
write.csv(resNorm_, file.path(opt$output,"results.norm.csv"), row.names = FALSE, quote = FALSE)

#Export the results Ash
resAsh_ <- as.data.frame(resAsh[order(resAsh$padj, decreasing = FALSE),])
write.table(resAsh_, file.path(opt$output, "results.ash.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
write.xlsx(resAsh_, file.path(opt$output,"results.ash.xlsx"),row.names=FALSE)
write.csv(resAsh_, file.path(opt$output,"results.ash.csv"), row.names = FALSE, quote = FALSE)

























