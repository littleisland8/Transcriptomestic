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
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(PCAtools))

defaultW <- getOption("warn")
options(warn = -1)

options(java.parameters = "-Xmx8000m")

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

#Define path and list of featureCounts read counts files.
indir <- file.path(getwd(),"count/featureCounts/STAR")

files_counts <- list.files(indir,"*.tsv$", full.names = T)

#Load first file geneid and read count
countData <- data.frame(fread(files_counts[1]))[c(1,7)]

# Loop and read the 7th column remaining files.
for(i in 2:length(files_counts)) {
  countData = cbind(countData, data.frame(fread(files_counts[i]))[7])
}

#preprocessing
sampleTable <- read.table(opt$sampletable, sep="\t", header=TRUE)

#Rename columns
names <- sampleTable$sampleName
colnames(countData) <- c("GeneID", names)
rownames(countData) <- countData$GeneID
countData <- countData[,c(2:ncol(countData))]

# export table countdata
write.table(countData, file.path(opt$output,"countdata_featurecounts.txt"), quote = FALSE, sep = "\t", row.names = TRUE)

#Countdata generation
#Convert to matrix
countdata <- countData
countdata <- as.matrix(countdata)

#Assign condition
condition <- factor(sampleTable$condition)

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))

# Build a DESeqDataSet from countData with DESeqDataSetFromMatrix
dds_featurecounts <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds_featurecounts <- estimateSizeFactors(dds_featurecounts)

## Filter low counts gene, filter out genes (row) with no count
smallestGroupSize <- min(c(length(which(sampleTable$condition == "TUMOR")),length(which(sampleTable$condition == "CONTROL"))))
keep <- rowSums(counts(dds_featurecounts) >= as.numeric(opt$count)) >= smallestGroupSize
dds <- dds_featurecounts[keep,]
dds$condition <- relevel(dds$condition, ref = "CONTROL")
dds$group <- sampleTable$group

# Visualize Boxplot of normalized counts
statusCol <- as.numeric(factor(sampleTable$condition)) + 1
pdf(file.path(opt$output,"boxplot.counts.pdf"), height=as.numeric(opt$Height), width=as.numeric(opt$width))
rs = rowSums(counts(dds))
par(mfrow=c(1,2)) # plots two plots 
boxplot(log2(counts(dds)[rs > 0,] + 1), cex.lab = 2, cex.axis = 2,xaxt = "n", col=statusCol)# not normalized
tick <- seq_along(as.character(dds$group))
axis(1, at = tick, labels = F)
text(tick, par("usr")[3] - 2.5, as.character(dds$group), srt = 90, xpd = T, offset =5)
boxplot(log2(counts(dds, normalized=TRUE)[rs > 0,] + 1), cex.lab = 2, cex.axis = 2,xaxt = "n", col=statusCol)# not normalized
tick <- seq_along(as.character(dds$group))
axis(1, at = tick, labels = F)
text(tick, par("usr")[3] - 2.5, as.character(dds$group), srt = 90, xpd = T, offset =5)
dev.off()

# Transform the raw count data

vst <- vst(dds, blind = TRUE)
rlog <- rlog(dds, blind = TRUE)

write.table(assay(rlog), file.path(opt$output,"normalizedCountRlog.txt"), quote = FALSE, sep = "\t", row.names = TRUE)
write.csv(assay(rlog), file.path(opt$output,"normalizedCountRlog.csv"), quote = FALSE, row.names = TRUE)
write.table(assay(vst), file.path(opt$output,"normalizedCountVst.txt"), quote = FALSE, sep = "\t")
write.csv(assay(vst), file.path(opt$output,"normalizedCountVst.csv"), quote = FALSE, row.names = TRUE)

## Scatter plot sample 1 2 
df <- bind_rows(
	as_data_frame(log2(counts(dds, normalized=TRUE)[, c(1,2)]+1)) %>%
		mutate(transformation = "log2(x + 1)"),
	as_data_frame(assay(vst)[, c(1,2)]) %>% mutate(transformation = "vst"),
	as_data_frame(assay(rlog)[, c(1,2)]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

p <- ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) + coord_fixed() + facet_grid( . ~ transformation)  

ggsave(file.path(opt$output,"transformedcounts.2samples.pdf"), height=as.numeric(opt$Height), width=as.numeric(opt$width))

# Effects of transformations on the variance
pdf(file.path(opt$output, "transformationvsvariance.vst.pdf"), height=as.numeric(opt$Height), width=as.numeric(opt$width))
meanSdPlot(assay(vst))
dev.off()

pdf(file.path(opt$output, "transformationvsvariance.rlog.pdf"), height=as.numeric(opt$Height), width=as.numeric(opt$width))
meanSdPlot(assay(rlog))
dev.off()

## sample distance rlog
sampleDists <- dist(t(assay(rlog)))

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


pdf(file.path(opt$output,"heatmap_dist_sample.rlog.pdf"), height=as.numeric(opt$Height), width=as.numeric(opt$width))
pheatmap(sampleDists,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors,annotation = annotation,annotation_colors = annotation_colors)
dev.off()

## Plot clustering dendogram rlog
pdf(file.path(opt$output,"cluster_dendo.rlog.pdf"), height=as.numeric(opt$Height), width=as.numeric(opt$width))
plot(hclust(sampleDists, method = "ward.D2"))
dev.off()

## sample distance vst
sampleDists <- dist(t(assay(vst)))

## Plot Heatmap with sample-to-sample distance
sampleDistMatrix <- as.matrix( sampleDists )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf(file.path(opt$output,"heatmap_dist_sample.vst.pdf"), height=as.numeric(opt$Height), width=as.numeric(opt$width))
pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors,annotation=annotation,annotation_colors = annotation_colors)
dev.off()

## Plot clustering dendogram vst
pdf(file.path(opt$output,"cluster_dendo.vst.pdf"), height=opt$Height, width=opt$width)
plot(hclust(sampleDists, method = "ward.D2"))
dev.off()


## Plot PCA in order to see sample-to-sample distance

## PCA
rlog_pca <- plotPCA(rlog, intgroup = c("condition"), returnData = TRUE)
percentVar <- round(100 * attr(rlog_pca, "percentVar"))
pca_rlog <- ggplot(rlog_pca, aes(x = PC1, y = PC2, color = condition)) +
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

ggsave(pca_rlog, filename = file.path(opt$output,"2DPCA.rlog.pdf"),height=as.numeric(opt$Height), width=as.numeric(opt$width))

## PCA PCAplot
## Generate p
rlog.output <- assay(rlog)
p <- pca(rlog.output, metadata = coldata, removeVar = 0.1)
#Determine optimum number of PCs to retain
#Let's perform Horn's parallel analysis 
horn <- parallelPCA(rlog.output)
#Let's perform elbow's parallel analysis
elbow <- findElbowPoint(p$variance)

#Taking these values, we can produce a new scree plot and mark these
# change number of PC for each experiment
pdf(file.path(opt$output,"screeplot.rlog.pdf"), height=as.numeric(opt$Height), width=as.numeric(opt$width))
screeplot(p,
          components = getComponents(p, 1:11),
          vline = c(horn$n, elbow)) +
  geom_label(aes(x = horn$n + 1, y = 50,
                 label = 'Horn\'s', vjust = -1, size = 8)) +
  geom_label(aes(x = elbow + 1, y = 50,
                 label = 'Elbow method', vjust = 1, size = 8))
dev.off()

# plotting
pdf(file.path(opt$output,"PCA.biplot.rlog.pdf"), height=as.numeric(opt$Height), width=as.numeric(opt$width))
biplot(p, title = "Principal Component Analysis (PCA)")
dev.off()


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
ggsave(pca_vst, filename = file.path(opt$output,"2DPCA.vst.pdf"),height=as.numeric(opt$Height), width=as.numeric(opt$width))

## PCA PCAplot
## Generate p
vst.output <- assay(vst)
p <- pca(vst.output, metadata = coldata, removeVar = 0.1)
#Determine optimum number of PCs to retain
#Let's perform Horn's parallel analysis 
horn <- parallelPCA(vst.output)
#Let's perform elbow's parallel analysis
elbow <- findElbowPoint(p$variance)

#Taking these values, we can produce a new scree plot and mark these
# change number of PC for each experiment
pdf(file.path(opt$output,"screeplot.vst.pdf"), height=as.numeric(opt$Height), width=as.numeric(opt$width))
screeplot(p,
          components = getComponents(p, 1:11),
          vline = c(horn$n, elbow)) +
  geom_label(aes(x = horn$n + 1, y = 50,
                 label = 'Horn\'s', vjust = -1, size = 8)) +
  geom_label(aes(x = elbow + 1, y = 50,
                 label = 'Elbow method', vjust = 1, size = 8))
dev.off()

# plotting
pdf(file.path(opt$output,"PCA.biplot.vst.pdf"), height=as.numeric(opt$Height), width=as.numeric(opt$width))
biplot(p, title = "Principal Component Analysis (PCA)")
dev.off()


## Differential Expression Analysis 
design(dds) <- ~ condition
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "TUMOR","CONTROL"), alpha = as.numeric(opt$alpha))

## export MA
pdf(file.path(opt$output,"MA_plot.pdf"), height = as.numeric(opt$Height), width = as.numeric(opt$width))
plotMA(res, ylim = c(-10,10), main = "MA plot")
dev.off()

## export dispersion plot plot that is important to evaluate when performing QC on RNA-seq data is the plot of dispersion versus the mean of normalized counts. 
#For a good dataset, we expect the dispersion to decrease as the mean of normalized counts increases for each gene. 
pdf(file.path(opt$output,"dispersion_plot.pdf"), height = as.numeric(opt$Height), width = as.numeric(opt$width))
plotDispEsts(dds, ylim = c(1e-6, 1e1))
dev.off()

pdf(paste0(file.path(opt$output,"hist_pvalue.padj.pdf")), height = as.numeric(opt$Height), width = as.numeric(opt$width))
par(mfrow = c(1,2))
hist(res$padj, breaks=20, col="grey")
hist(res$pvalue, breaks=20, col="grey")
dev.off()

# outlayer plot
pdf(file.path(opt$output,"Outlayer.pdf"), height = as.numeric(opt$Height), width = as.numeric(opt$width))
par(mar=c(8,5,3,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2, main = "Outlayer", col=statusCol)
dev.off()


ResAnnotation <- function (restable){

	#AnnotationDBI
  restable <- res
	ens.str <- substr(rownames(restable), 1, 15)
	restable$symbol <- mapIds(org.Hs.eg.db,
													keys=ens.str,
													column="SYMBOL",
													keytype="ENSEMBL",
													multiVals = "first")
	restable$entrez <- mapIds(org.Hs.eg.db,
													keys=ens.str,
													column="ENTREZID",
													keytype="ENSEMBL")
	restable$alias <- mapIds(org.Hs.eg.db,
												 keys=ens.str,
												 column="ALIAS",
												 keytype="ENSEMBL")
	restable$EnsemblID <- ens.str
	
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

pdf(file.path(opt$output,"VolcanoPlot.results.pdf"), height = as.numeric(opt$Height), width = as.numeric(opt$width))
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
pdf(file.path(opt$output,"multipleMAplot.pdf"),height = as.numeric(opt$Height), width = as.numeric(opt$width))
par(mfrow=c(2,2), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-10,10)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
plotMA(resIHW, xlim=xlim, ylim=ylim, main="Independent Hypotesis Weighting")
dev.off()

#resApeglm
resApeglm <- ResAnnotation(resLFC)

pdf(file.path(opt$output,"VolcanoPlot.apeglm.pdf"), height = as.numeric(opt$Height), width = as.numeric(opt$width))
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

pdf(file.path(opt$output,"VolcanoPlot.Normal.pdf"), height = as.numeric(opt$Height), width = as.numeric(opt$width))
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

pdf(file.path(opt$output,"VolcanoPlot.Normal.pdf"), height = as.numeric(opt$Height), width = as.numeric(opt$width))
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
pdf(file.path(opt$output,"heatmap.top20.vst.pdf"),height = as.numeric(opt$Height), width = as.numeric(opt$width))
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
pdf(file.path(opt$output,"heatmap_top50.vst.pdf"),height = as.numeric(opt$Height), width = as.numeric(opt$width))
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
pdf(file.path(opt$output,"heatmap_top100.vst.pdf"),,height = as.numeric(opt$Height), width = as.numeric(opt$width))
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
pdf(file.path(opt$output,"heatmap_top1000.vst.pdf"), ,height = as.numeric(opt$Height), width = as.numeric(opt$width))
pheatmap(mat1000, annotation = annotation, annotation_colors = annotation_colors, fontsize_row = 8, color = pal,show_rownames = FALSE, scale = "row")
dev.off()

## Gene clustering top 20 genes rlog
topVarGenes20 <- head(order(rowVars(assay(rlog)), decreasing = TRUE), 20)
mat20  <- assay(rlog)[ topVarGenes20, ]

#export
pdf(file.path(opt$output,"heatmap.top20.rlog.pdf"),height = as.numeric(opt$Height), width = as.numeric(opt$width))
pheatmap(mat20, annotation = annotation, annotation_colors = annotation_colors, color = pal, scale = "row")
dev.off()

df_20 <- as.data.frame(res[which(rownames(res) %in% rownames(mat20)),])
write.table(df_20, file.path(opt$output,"top20deg.rlog.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
write.csv(df_20, file.path(opt$output,"top20deg.rlog.csv"), row.names = TRUE, quote = FALSE)
write.xlsx(df_20, file.path(opt$output,"top20deg.rlog.xlsx"), row.names = TRUE)

## Gene clustering top 50 genes rlog
topVarGenes50 <- head(order(rowVars(assay(rlog)), decreasing = TRUE), 50)
mat50  <- assay(rlog)[ topVarGenes50, ]

df_50 <- as.data.frame(res[which(rownames(res) %in% rownames(mat50)),])
write.table(df_50, file.path(opt$output,"top50deg.rlog.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
write.csv(df_50, file.path(opt$output,"top50deg.rlog.csv"), row.names = TRUE, quote = FALSE)
write.xlsx(df_50, file.path(opt$output,"top50deg.rlog.xlsx"), row.names = TRUE)

#export
pdf(file.path(opt$output,"heatmap_top50.rlog.pdf"),height = as.numeric(opt$Height), width = as.numeric(opt$width))
pheatmap(mat50, annotation = annotation, annotation_colors = annotation_colors, color = pal, scale = "row")
dev.off()

## Gene clustering top 100 genes rlog
topVarGenes100 <- head(order(rowVars(assay(rlog)), decreasing = TRUE), 100)
mat100  <- assay(rlog)[ topVarGenes100, ]

df_100 <- as.data.frame(res[which(rownames(res) %in% rownames(mat100)),])
write.table(df_100, file.path(opt$output,"top100deg.rlog.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
write.csv(df_100, file.path(opt$output,"top100deg.rlog.csv"), row.names = TRUE, quote = FALSE)
write.xlsx(df_100, file.path(opt$output,"top100deg.rlog.xlsx"), row.names = TRUE)

#export
pdf(file.path(opt$output,"heatmap_top100.rlog.pdf"),height = as.numeric(opt$Height), width = as.numeric(opt$width))
pheatmap(mat100, annotation = annotation, annotation_colors = annotation_colors, color = pal, scale = "row")
dev.off()

## Gene clustering top 1000 genes rlog
topVarGenes1000 <- head(order(rowVars(assay(rlog)), decreasing = TRUE), 1000)
mat1000  <- assay(rlog)[ topVarGenes1000, ]

df_1000 <- as.data.frame(res[which(rownames(rlog) %in% rownames(mat1000)),])
write.table(df_1000, file.path(opt$output,"top1000deg.rlog.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
write.csv(df_1000, file.path(opt$output,"top1000deg.rlog.csv"), row.names = TRUE, quote = FALSE)
write.xlsx(df_1000, file.path(opt$output,"top1000deg.rlog.xlsx"), row.names = TRUE)

#export
pdf(file.path(opt$output,"heatmap_top1000.rlog.pdf"),height = as.numeric(opt$Height), width = as.numeric(opt$width))
pheatmap(mat1000, annotation = annotation, annotation_colors = annotation_colors, fontsize_row = 8, color = pal,show_rownames = FALSE, scale = "row")
dev.off()

#Export the results
res_ <- as.data.frame(res[order(res$padj, decreasing = FALSE),])
#Add column of FC
res_$FC <-  2**res_$log2FoldChange
#Add column of absolute FC
res_$absFC <- 2**abs(res_$log2FoldChange)
#Relocate FC and absolute FC columns
res_ <- res_ %>% relocate(FC, .before = log2FoldChange)
res_ <- res_ %>% relocate(absFC, .before = log2FoldChange)

write.table(res_, file.path(opt$output, "results.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
write.xlsx(res_, file.path(opt$output,"results.xlsx"),row.names=FALSE)
write.csv(res_, file.path(opt$output,"results.csv"), row.names = FALSE, quote = FALSE)

#Export the filtered results
res_filter <- res_[(which(res_$pvalue <= 0.05 & abs(res_$log2FoldChange)>=0.58)),]
write.table(res_filter, file.path(opt$output, "results.filtered.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
write.xlsx(res_filter, file.path(opt$output,"results.filtered.xlsx"),row.names=FALSE)
write.csv(res_filter, file.path(opt$output,"results.filtered.csv"), row.names = FALSE, quote = FALSE)

#Heatmap DEG res_
res_filter <- merge(as.data.frame(res_filter), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
norm_counts_DEG <- res_filter[c(1,14:ncol(res_filter))]
rownames(norm_counts_DEG) <- norm_counts_DEG$Row.names
norm_counts_DEG$Row.names <- NULL

norm_counts_DEG_mat <- as.matrix(norm_counts_DEG)
scaledata <- t(scale(t(norm_counts_DEG_mat)))
hc <- hclust(as.dist(1-cor(scaledata, method="spearman")), method="complete")
sampleTree = as.dendrogram(hc, method="average")
#plot(sampleTree,
#    main = "Sample Clustering",
#     ylab = "Height")

hr <- hclust(as.dist(1-cor(t(scaledata), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
geneTree = as.dendrogram(hr, method="average")
#plot(geneTree,
#     leaflab = "none",             
#     main = "Gene Clustering",
#     ylab = "Height")

pdf(file.path(opt$output,"heatmap_res_DEG.pdf"),as.numeric(opt$Height), width = as.numeric(opt$width))
nrow_heat <- nrow(norm_counts_DEG_mat)
heatmap.2(norm_counts_DEG_mat,
          Rowv=as.dendrogram(hr), 
          Colv=as.dendrogram(hc),
          col=redgreen(100),
          scale="row",
          margins = c(7, 7),
          cexCol = 0.7,
          labRow = F,
          main = paste0("DEG results heatmap n=", nrow_heat),
          trace = "none")
dev.off()

# graph for first 100 top regulated genes for res_ vst
mypalette <- brewer.pal(11, "RdYlBu")
morecols <- colorRampPalette(mypalette)
ann_colors <- list(CellType= c(CONTROL="orange", TUMOR="purple"))
topVarGenes_res <- res_filter$EnsemblID[1:100]
rownames(vst) <- gsub("\\..*","",rownames(vst))
mat  <- assay(vst)[ topVarGenes_res, ]
mat  <- mat - rowMeans(mat)
annocol <- as.data.frame(colData(vst)[,1,drop=FALSE])
pdf(file.path(opt$output,"Top100.vst.res.DEGs.pdf"), as.numeric(opt$Height), width = as.numeric(opt$width))
pheatmap::pheatmap(mat, annotation_col = annocol, cellwidth = 20, cellheight = 10,
                   cutree_cols = 2, cutree_rows = 2, col=rev(morecols(50)), 
                   treeheight_row = 150, annotation_colors = ann_colors,
                   main = "Top 100 most significant regulated genes", border_color = "black" )
dev.off()


# graph for first 100 top regulated genes for res_ rlog
#mypalette <- brewer.pal(11, "RdYlBu")
#morecols <- colorRampPalette(mypalette)
#ann_colors <- list(CellType= c(CONTROL="orange", TUMOR="purple"))
#topVarGenes_res <- res_filter$EnsemblID[1:100]
rownames(rlog) <- gsub("\\..*","",rownames(rlog))
mat  <- assay(rlog)[ topVarGenes_res, ]
mat  <- mat - rowMeans(mat)
annocol <- as.data.frame(colData(rlog)[,1,drop=FALSE])
pdf(file.path(opt$output,"Top100.rlog.res.DEGs.pdf"), as.numeric(opt$Height), width = as.numeric(opt$width))
pheatmap::pheatmap(mat, annotation_col = annocol, cellwidth = 20, cellheight = 10,
                   cutree_cols = 2, cutree_rows = 2, col=rev(morecols(50)), 
                   treeheight_row = 150, annotation_colors = ann_colors,
                   main = "Top 100 most significant regulated genes", border_color = "black" )
dev.off()

#Export the results apeglm
resApeglm_ <- as.data.frame(resApeglm[order(resApeglm$padj, decreasing = FALSE),])
#Add column of FC
resApeglm_$FC <-  2**resApeglm_$log2FoldChange
#Add column of absolute FC
resApeglm_$absFC <- 2**abs(resApeglm_$log2FoldChange)
#Relocate FC and absolute FC columns
resApeglm_ <- resApeglm_ %>% relocate(FC, .before = log2FoldChange)
resApeglm_ <- resApeglm_ %>% relocate(absFC, .before = log2FoldChange)
write.table(resApeglm_, file.path(opt$output, "results.apeglm.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
write.xlsx(resApeglm_, file.path(opt$output,"results.apeglm.xlsx"),row.names=FALSE)
write.csv(resApeglm_, file.path(opt$output,"results.apeglm.csv"), row.names = FALSE, quote = FALSE)

#Export the filtered results apeglm
resApeglm_filter <- resApeglm_[(which(resApeglm_$pvalue <= 0.05 & abs(resApeglm_$log2FoldChange)>=0.58)),]
write.table(resApeglm_filter, file.path(opt$output, "results.filtered.apeglm.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
write.xlsx(resApeglm_filter, file.path(opt$output,"results.filtered.apeglm.xlsx"),row.names=FALSE)
write.csv(resApeglm_filter, file.path(opt$output,"results.filtered.apeglm.csv"), row.names = FALSE, quote = FALSE)

#Heatmap DEG apeglm
resApeglm_filter <- merge(as.data.frame(resApeglm_filter), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
norm_counts_DEG <- resApeglm_filter[c(1,14:ncol(resApeglm_filter))]
rownames(norm_counts_DEG) <- norm_counts_DEG$Row.names
norm_counts_DEG$Row.names <- NULL

norm_counts_DEG_mat <- as.matrix(norm_counts_DEG)
scaledata <- t(scale(t(norm_counts_DEG_mat)))
hc <- hclust(as.dist(1-cor(scaledata, method="spearman")), method="complete")
sampleTree = as.dendrogram(hc, method="average")
#plot(sampleTree,
#     main = "Sample Clustering",
#    ylab = "Height")

hr <- hclust(as.dist(1-cor(t(scaledata), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
geneTree = as.dendrogram(hr, method="average")
#plot(geneTree,
#     leaflab = "none",             
#     main = "Gene Clustering",
#     ylab = "Height")

pdf(file.path(opt$output,"heatmap_resApeglm_DEG.pdf"),as.numeric(opt$Height), width = as.numeric(opt$width))
nrow_heat <- nrow(norm_counts_DEG_mat)
heatmap.2(norm_counts_DEG_mat,
          Rowv=as.dendrogram(hr), 
          Colv=as.dendrogram(hc),
          col=redgreen(100),
          scale="row",
          margins = c(7, 7),
          cexCol = 0.7,
          labRow = F,
          main = paste0("DEG Apeglm heatmap n=", nrow_heat),
          trace = "none")
dev.off()

# graph for first 100 top regulated genes for res_apeglm vst
mypalette <- brewer.pal(11, "RdYlBu")
morecols <- colorRampPalette(mypalette)
ann_colors <- list(CellType= c(CONTROL="orange", TUMOR="purple"))
topVarGenes_res <- resApeglm_filter$EnsemblID[1:100]
rownames(vst) <- gsub("\\..*","",rownames(vst))
mat  <- assay(vst)[ topVarGenes_res, ]
mat  <- mat - rowMeans(mat)
annocol <- as.data.frame(colData(vst)[,1,drop=FALSE])
pdf(file.path(opt$output,"Top100.vst.res.apeglm.DEGs.pdf"), as.numeric(opt$Height), width = as.numeric(opt$width))
pheatmap::pheatmap(mat, annotation_col = annocol, cellwidth = 20, cellheight = 10,
                   cutree_cols = 2, cutree_rows = 2, col=rev(morecols(50)), 
                   treeheight_row = 150, annotation_colors = ann_colors, 
                   main = "Top 100 most significant regulated genes", border_color = "black" )
dev.off()

# graph for first 100 top regulated genes for res_apeglm rlog
#mypalette <- brewer.pal(11, "RdYlBu")
#morecols <- colorRampPalette(mypalette)
#ann_colors <- list(CellType= c(CONTROL="orange", TUMOR="purple"))
#topVarGenes_res <- resApeglm_filter$EnsemblID[1:100]
rownames(rlog) <- gsub("\\..*","",rownames(rlog))
mat  <- assay(rlog)[ topVarGenes_res, ]
mat  <- mat - rowMeans(mat)
annocol <- as.data.frame(colData(rlog)[,1,drop=FALSE])
pdf(file.path(opt$output,"Top100.rlog.res.apeglm.DEGs.pdf"), as.numeric(opt$Height), width = as.numeric(opt$width))
pheatmap::pheatmap(mat, annotation_col = annocol, cellwidth = 20, cellheight = 10,
                   cutree_cols = 2, cutree_rows = 2, col=rev(morecols(50)), 
                   treeheight_row = 150, annotation_colors = ann_colors, 
                   main = "Top 100 most significant regulated genes", border_color = "black" )
dev.off()

#Export the results Norm
resNorm_ <- as.data.frame(resNorm[order(resNorm$padj, decreasing = FALSE),])
#Add column of FC
resNorm_$FC <-  2**resNorm_$log2FoldChange
#Add column of absolute FC
resNorm_$absFC <- 2**abs(resNorm_$log2FoldChange)
#Relocate FC and absolute FC columns
resNorm_ <- resNorm_ %>% relocate(FC, .before = log2FoldChange)
resNorm_ <- resNorm_ %>% relocate(absFC, .before = log2FoldChange)
write.table(resNorm_, file.path(opt$output, "results.norm.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
write.xlsx(resNorm_, file.path(opt$output,"results.norm.xlsx"),row.names=FALSE)
write.csv(resNorm_, file.path(opt$output,"results.norm.csv"), row.names = FALSE, quote = FALSE)

#Export the filtered results Norm
resNorm_filter <- resNorm_[(which(resNorm_$pvalue <= 0.05 & abs(resNorm_$log2FoldChange)>=0.58)),]
write.table(resNorm_filter, file.path(opt$output, "results.filtered.norm.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
write.xlsx(resNorm_filter, file.path(opt$output,"results.filtered.norm.xlsx"),row.names=FALSE)
write.csv(resNorm_filter, file.path(opt$output,"results.filtered.norm.csv"), row.names = FALSE, quote = FALSE)

#Heatmap DEG norm
resNorm_filter <- merge(as.data.frame(resNorm_filter), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
norm_counts_DEG <- resNorm_filter[c(1,14:ncol(resNorm_filter))]
rownames(norm_counts_DEG) <- norm_counts_DEG$Row.names
norm_counts_DEG$Row.names <- NULL

norm_counts_DEG_mat <- as.matrix(norm_counts_DEG)
scaledata <- t(scale(t(norm_counts_DEG_mat)))
hc <- hclust(as.dist(1-cor(scaledata, method="spearman")), method="complete")
sampleTree = as.dendrogram(hc, method="average")
#plot(sampleTree,
#     main = "Sample Clustering",
#     ylab = "Height")

hr <- hclust(as.dist(1-cor(t(scaledata), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
geneTree = as.dendrogram(hr, method="average")
#plot(geneTree,
#     leaflab = "none",             
#     main = "Gene Clustering",
#     ylab = "Height")

pdf(file.path(opt$output,"heatmap_resnormDEG.pdf"),as.numeric(opt$Height), width = as.numeric(opt$width))
nrow_heat <- nrow(norm_counts_DEG_mat)
heatmap.2(norm_counts_DEG_mat,
          Rowv=as.dendrogram(hr), 
          Colv=as.dendrogram(hc),
          col=redgreen(100),
          scale="row",
          margins = c(7, 7),
          cexCol = 0.7,
          labRow = F,
          main = paste0("DEG norm heatmap n=", nrow_heat),
          trace = "none")
dev.off()

# graph for first 100 top regulated genes for res_norm vst
mypalette <- brewer.pal(11, "RdYlBu")
morecols <- colorRampPalette(mypalette)
ann_colors <- list(CellType= c(CONTROL="orange", TUMOR="purple"))
topVarGenes_res <- resNorm_filter$EnsemblID[1:100]
rownames(vst) <- gsub("\\..*","",rownames(vst))
mat  <- assay(vst)[ topVarGenes_res, ]
mat  <- mat - rowMeans(mat)
annocol <- as.data.frame(colData(vst)[,1,drop=FALSE])
pdf(file.path(opt$output,"Top100.vst.res.norm.DEGs.pdf"), as.numeric(opt$Height), width = as.numeric(opt$width))
pheatmap::pheatmap(mat, annotation_col = annocol, cellwidth = 20, cellheight = 10,
                   cutree_cols = 2, cutree_rows = 2, col=rev(morecols(50)), 
                   treeheight_row = 150, annotation_colors = ann_colors, 
                   main = "Top 100 most significant regulated genes", border_color = "black" )
dev.off()

# graph for first 100 top regulated genes for res_norm rlog
#mypalette <- brewer.pal(11, "RdYlBu")
#morecols <- colorRampPalette(mypalette)
#ann_colors <- list(CellType= c(CONTROL="orange", TUMOR="purple"))
#topVarGenes_res <- resNorm_filter$EnsemblID[1:100]
rownames(rlog) <- gsub("\\..*","",rownames(rlog))
mat  <- assay(rlog)[ topVarGenes_res, ]
mat  <- mat - rowMeans(mat)
annocol <- as.data.frame(colData(rlog)[,1,drop=FALSE])
pdf(file.path(opt$output,"Top100.rlog.res.norm.DEGs.pdf"), as.numeric(opt$Height), width = as.numeric(opt$width))
pheatmap::pheatmap(mat, annotation_col = annocol, cellwidth = 20, cellheight = 10,
                   cutree_cols = 2, cutree_rows = 2, col=rev(morecols(50)), 
                   treeheight_row = 150, annotation_colors = ann_colors, 
                   main = "Top 100 most significant regulated genes", border_color = "black" )
dev.off()

#Export the results Ash
resAsh_ <- as.data.frame(resAsh[order(resAsh$padj, decreasing = FALSE),])
#Add column of FC
resAsh_$FC <-  2**resAsh_$log2FoldChange
#Add column of absolute FC
resAsh_$absFC <- 2**abs(resAsh_$log2FoldChange)
#Relocate FC and absolute FC columns
resAsh_ <- resAsh_ %>% relocate(FC, .before = log2FoldChange)
resAsh_ <- resAsh_ %>% relocate(absFC, .before = log2FoldChange)
write.table(resAsh_, file.path(opt$output, "results.ash.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
write.xlsx(resAsh_, file.path(opt$output,"results.ash.xlsx"),row.names=FALSE)
write.csv(resAsh_, file.path(opt$output,"results.ash.csv"), row.names = FALSE, quote = FALSE)

#Export the filtered results Ash
resAsh_filter <- resAsh_[(which(resAsh_$pvalue <= 0.05 & abs(resAsh_$log2FoldChange)>=0.58)),]
write.table(resAsh_filter, file.path(opt$output, "results.filtered.ash.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
write.xlsx(resAsh_filter, file.path(opt$output,"results.filtered.ash.xlsx"),row.names=FALSE)
write.csv(resAsh_filter, file.path(opt$output,"results.filtered.ash.csv"), row.names = FALSE, quote = FALSE)

#Heatmap DEG Ash
resAsh_filter <- merge(as.data.frame(resAsh_filter), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
norm_counts_DEG <- resAsh_filter[c(1,14:ncol(resAsh_filter))]
rownames(norm_counts_DEG) <- norm_counts_DEG$Row.names
norm_counts_DEG$Row.names <- NULL

norm_counts_DEG_mat <- as.matrix(norm_counts_DEG)
scaledata <- t(scale(t(norm_counts_DEG_mat)))
hc <- hclust(as.dist(1-cor(scaledata, method="spearman")), method="complete")
sampleTree = as.dendrogram(hc, method="average")
#plot(sampleTree,
#     main = "Sample Clustering",
#     ylab = "Height")

hr <- hclust(as.dist(1-cor(t(scaledata), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
geneTree = as.dendrogram(hr, method="average")
#plot(geneTree,
#     leaflab = "none",             
#     main = "Gene Clustering",
#     ylab = "Height")

pdf(file.path(opt$output,"heatmap_Ash_DEG.pdf"),as.numeric(opt$Height), width = as.numeric(opt$width))
nrow_heat <- nrow(norm_counts_DEG_mat)
heatmap.2(norm_counts_DEG_mat,
          Rowv=as.dendrogram(hr), 
          Colv=as.dendrogram(hc),
          col=redgreen(100),
          scale="row",
          margins = c(7, 7),
          cexCol = 0.7,
          labRow = F,
          main = paste0("DEG Ash heatmap n=", nrow_heat),
          trace = "none")
dev.off()

# graph for first 100 top regulated genes for res_ash vst
mypalette <- brewer.pal(11, "RdYlBu")
morecols <- colorRampPalette(mypalette)
ann_colors <- list(CellType= c(CONTROL="orange", TUMOR="purple"))
topVarGenes_res <- resAsh_filter$EnsemblID[1:100]
rownames(vst) <- gsub("\\..*","",rownames(vst))
mat  <- assay(vst)[ topVarGenes_res, ]
mat  <- mat - rowMeans(mat)
annocol <- as.data.frame(colData(vst)[,1,drop=FALSE])
pdf(file.path(opt$output,"Top100.vst.res.ash.DEGs.pdf"), as.numeric(opt$Height), width = as.numeric(opt$width))
pheatmap::pheatmap(mat, annotation_col = annocol, cellwidth = 20, cellheight = 10,
                   cutree_cols = 2, cutree_rows = 2, col=rev(morecols(50)), 
                   treeheight_row = 150, annotation_colors = ann_colors, 
                   main = "Top 100 most significant regulated genes", border_color = "black" )
dev.off()

# graph for first 100 top regulated genes for res_ash rlog
#mypalette <- brewer.pal(11, "RdYlBu")
#morecols <- colorRampPalette(mypalette)
#ann_colors <- list(CellType= c(CONTROL="orange", TUMOR="purple"))
#topVarGenes_res <- resAsh_filter$EnsemblID[1:100]
rownames(rlog) <- gsub("\\..*","",rownames(rlog))
mat  <- assay(rlog)[ topVarGenes_res, ]
mat  <- mat - rowMeans(mat)
annocol <- as.data.frame(colData(rlog)[,1,drop=FALSE])
pdf(file.path(opt$output,"Top100.rlog.res.ash.DEGs.pdf"), as.numeric(opt$Height), width = as.numeric(opt$width))
pheatmap::pheatmap(mat, annotation_col = annocol, cellwidth = 20, cellheight = 10,
                   cutree_cols = 2, cutree_rows = 2, col=rev(morecols(50)), 
                   treeheight_row = 150, annotation_colors = ann_colors, 
                   main = "Top 100 most significant regulated genes", border_color = "black" )
dev.off()