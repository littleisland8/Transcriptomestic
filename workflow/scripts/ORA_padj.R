.libPaths(new="/home/simone.romagnoli/R/x86_64-pc-linux-gnu-library/4.3/")
library(ReactomePA)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(ggplot2)
library(enrichplot)
library(meshes)
library(MeSHDbi)
library(msigdbr)
library(xlsx)
library(AnnotationDbi)

rm(list = ls())

outputdir <- "/projects2/2024_Chiarugi/Transcriptomestic/results/unstranded/STAR_HTSeq_nosample1/ORA_padj"
cmd <- paste("mkdir -p", outputdir)
system(cmd)

###### ORA
df <- read.csv("/projects2/2024_Chiarugi/Transcriptomestic/results/unstranded/STAR_HTSeq_nosample1/results.filtered.padj.csv")

df <- df[which(df$padj < 0.05 | df$padj==0),]

genelist <- df$log2FoldChange
names(genelist) <- df$entrez

## GO classification
GOoutdir <- paste(file.path(outputdir, "GO"))
cmd <- paste("mkdir -p", GOoutdir)
system(cmd)

ontologies <- c("MF", "BP", "CC")

for (ont in ontologies){
  
  outdir_ <- paste(file.path(GOoutdir, ont))
  system(paste("mkdir -p", outdir_))
  message(paste0("Analyze...",ont))
  
  ## GO enrichment
  ego <- enrichGO(gene          = names(genelist),
                  OrgDb         = org.Hs.eg.db,
                  ont           = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE,
                  minGSSize     = 2,
                  maxGSSize     = 500)
  
  ego_table <- ego@result  
  write.table(ego_table, file.path(outdir_,"ego.tsv"), sep = "\t", quote = FALSE)
  write.csv(ego_table, file.path(outdir_,"ego.csv"), quote = FALSE)  
  #write.xlsx(ego_table, file.path(outdir_,"ego.xlsx"))
  
  if (length(which(ego@result$p.adjust < 0.05)) > 10){
    
    ##Plot
    
    #barplot
    message(paste0("Drawing barplot"))
    p <- barplot(ego, showCategory = 10, colorby = 'pvalue', font.size = 8) + ggtitle("GO Barplot")
    ggsave(file.path(outdir_,"barplot.pdf"))
    
    #dotplot show most significant enriched terms 
    message(paste0("Drawing dotplot"))
    p <-dotplot(ego, showCategory = 10) + ggtitle("GO Dotplot")
    ggsave(file.path(outdir_, "dotplot.pdf"), height = 9, width = 9)
    
    #cnetplot show genes involved in significant terms. In order to consider the potentially biological complexities in which a gene may belong to multiple annotation categories and provide information of numeric changes if available, cnetplot function extracts the complex association 
    message(paste0("Drawing cnetplot p1 for"))
    p1 <- cnetplot(ego, showCategory = 10) + ggtitle("GO Gene-concept Network plot") #category size scaled by gene number 
    ggsave(file.path(outdir_, "cnetplot.1.pdf"))
    
    ## categorySize can be scaled by 'pvalue' or 'geneNum'
    message(paste0("Drawing cnetplot p2"))
    p2 <- cnetplot(ego, categorySize="pvalue", showCategory = 10,foldChange=genelist,max.overlaps=20) + ggtitle("GO Gene-concept Network scaled by pvalue")
    ggsave(file.path(outdir_, "cnetplot.2.pdf"))
    
    message(paste0("Drawing cnetplot p3"))
    p3 <- cnetplot(ego, circular = TRUE, colorEdge = TRUE, foldChange=genelist, showCategory = 10) + ggtitle("GO Gene-concept Network circular")
    ggsave(file.path(outdir_, "cnetplot.3.pdf"))
    
    #cowplot::plot_grid(p1,p2, p3, ncol=1, labels=LETTERS[1:2])
    #ggsave(file.path(outdir_, "cnetplot.total.pdf"))
    
    #Heatmap-like functional classification
    p1 <- heatplot(ego, showCategory=5) + theme(axis.text.x = element_text(size = 4, angle = 90))
    p2 <- heatplot(ego, foldChange=genelist, showCategory=5) + theme(axis.text.x = element_text(size = 4, angle = 90))
    cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])
    ggsave(file.path(outdir_, "heatmap.pdf"), height = 10, width = 25)
    
    #Tree plot
    edox2 <- pairwise_termsim(ego)
    p1 <- treeplot(edox2, label_format = 4, showCategory = 10,nCluster=4)
    ggsave(file.path(outdir_, "clust.wardD2.pdf"), height=12, width= 15)
    p2 <- treeplot(edox2, hclust_method = "average",label_format = 4,showCategory = 10,nCluster=4)
    ggsave(file.path(outdir_, "clust.average.pdf"), height=12, width= 15)
    #aplot::plot_list(p1, p2, tag_levels='A')
    #ggsave(file.path(outdir_, "clust.pdf"))
    
    #emapplot show different genesets connected by overlapping genes. Overlapping gene sets are tend to cluster together 
    message(paste0("Drawing emapplot for"))
    p <-emapplot(edox2, showCategory = 10,layout="kk", node_scale = 2)  + ggtitle("GO Enrichment Map plot")
    ggsave(file.path(outdir_, "emapplot.pdf"))
    
    #UpSet Plot
    
    pdf(file.path(outdir_,"upset.pdf"))
    upsetplot(ego) + theme(plot.margin = unit(c(1,1,1,4), "cm"), axis.text.y = element_text(size = 15/.pt))
    dev.off()
  }    
}


## Reactome enrichment analysis
reactomepath <- paste(file.path(outputdir, "REACTOME"))
cmd <- paste("mkdir -p", reactomepath)
system(cmd)

## REACTOME
x <- enrichPathway(gene=names(genelist),pvalueCutoff=0.05, qvalueCutoff  = 0.05, readable=T, organism = "human",minGSSize = 2,maxGSSize = 500)

reactome <- as.data.frame(x@result) 
write.table(reactome,file.path(reactomepath,"reactome_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(reactome,file.path(reactomepath,"reactome_table.csv"), row.names = FALSE, quote = FALSE)
#write.xlsx(reactome,file.path(reactomepath,"reactome_table.xlsx"))

if (length(which(x@result$p.adjust < 0.05)) > 10){
  
  ##Plot
  
  #barplot
  message(paste0("Drawing barplot"))
  p <- barplot(x, showCategory = 10, colorby = 'pvalue') + ggtitle("REACTOME Barplot")
  ggsave(file.path(reactomepath,"barplot.pdf"))
  
  #dotplot show most significant enriched terms 
  message(paste0("Drawing dotplot"))
  p <- dotplot(x, showCategory = 10) + ggtitle("REACTOME Dotplot")
  ggsave(file.path(reactomepath, "dotplot.pdf"))
  
  #cnetplot show genes involved in significant terms. In order to consider the potentially biological complexities in which a gene may belong to multiple annotation categories and provide information of numeric changes if available, cnetplot function extracts the complex association 
  message(paste0("Drawing cnetplot p1 for"))
  p1 <- cnetplot(x, showCategory = 10) + ggtitle("REACTOME Gene-concept Network plot") #category size scaled by gene number 
  ggsave(file.path(reactomepath, "cnetplot.1.pdf"))
  
  ## categorySize can be scaled by 'pvalue' or 'geneNum'
  message(paste0("Drawing cnetplot p2"))
  p2 <- cnetplot(x, categorySize="pvalue", showCategory = 10,foldChange=genelist) + ggtitle("REACTOME Gene-concept Network scaled by pvalue")
  ggsave(file.path(reactomepath, "cnetplot.2.pdf"))
  
  message(paste0("Drawing cnetplot p3"))
  p3 <- cnetplot(x, circular = TRUE, colorEdge = TRUE, showCategory = 10,foldChange=genelist) + ggtitle("REACTOME Gene-concept Network circular")
  ggsave(file.path(reactomepath, "cnetplot.3.pdf"))
  
  #cowplot::plot_grid(p1,p2, p3, ncol=3, labels=LETTERS[1:2])
  #ggsave(file.path(reactomepath, "cnetplot.total.pdf"))
  
  #Heatmap-like functional classification
  p1 <- heatplot(x, showCategory=10)
  p2 <- heatplot(x, foldChange=genelist, showCategory=10) #+ theme(axis.text.x = element_text(size = 4, angle = 90))
  cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])
  ggsave(file.path(reactomepath, "heatmap.pdf"))
  
  #Tree plot
  edox2 <- pairwise_termsim(x)
  p1 <- treeplot(edox2, label_format = 4, showCategory = 50, nCluster=4)
  ggsave(file.path(reactomepath, "clust.wardD2.pdf"), height=12, width=15)
  p2 <- treeplot(edox2, hclust_method = "average",label_format = 4,nCluster=4)
  ggsave(file.path(reactomepath, "clust.average.pdf"), height=12, width=15)
  #aplot::plot_list(p1, p2, tag_levels='A')
  #ggsave(file.path(outdir_, "clust.pdf"))
  
  #emapplot show different genesets connected by overlapping genes. Overlapping gene sets are tend to cluster together 
  message(paste0("Drawing emapplot for"))
  p <- emapplot(edox2, showCategory = 20, node_scale = 2)  + ggtitle("REACTOME Enrichment Map plot")
  ggsave(file.path(reactomepath, "emapplot.pdf"))
  
  pdf(file.path(reactomepath,"upset.pdf"))
  upsetplot(x) + theme(plot.margin = unit(c(1,1,1,4), "cm"), axis.text.y = element_text(size = 15/.pt))
  dev.off()
  
}


## KEGG enrichment analysis

# install the packages
remotes::install_github("YuLab-SMU/createKEGGdb")
# import the library and create a KEGG database locally 
library(createKEGGdb)
species <-c("ath","hsa","mmu", "rno","dre","dme","cel")
createKEGGdb::create_kegg_db(species)
# You will get KEGG.db_1.0.tar.gz file in your working directory

KEGGoutdir <- paste(file.path(outputdir, "KEGG"))
cmd <- paste("mkdir -p", KEGGoutdir)
system(cmd)

hs <- search_kegg_organism("Homo sapiens", by = 'scientific_name')

kk <- enrichKEGG(gene         = names(genelist),
                 organism     = 'hsa',
                 pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.05,
                 minGSSize    = 2, 
                 maxGSSize    = 500)


kks <- setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')

#convert entrez to symbol
#symbol <- clusterProfiler::bitr(names(genelist), fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
#kk@gene2Symbol <- symbol$SYMBOL

kks_table <- kks@result
write.table(kks_table, file.path(KEGGoutdir,"KEGG.tsv"), sep = "\t", quote = FALSE)
write.csv(kks_table, file.path(KEGGoutdir,"KEGG.csv"), quote = FALSE)
write.xlsx(kks_table, file.path(KEGGoutdir,"KEGG.xlsx"))

if (length(which(kks_table$p.adjust < 0.05)) > 10){
  
  ##Plot
  
  #barplot
  message(paste0("Drawing barplot"))
  p <- barplot(kks, showCategory = 10, colorby = 'pvalue') + ggtitle("KEGG Barplot")
  ggsave(file.path(KEGGoutdir,"barplot.pdf"))
  
  #dotplot show most significant enriched terms 
  message(paste0("Drawing dotplot"))
  p <- dotplot(kks, showCategory = 10) + ggtitle("KEGG Dotplot")
  ggsave(file.path(KEGGoutdir, "dotplot.pdf"))
  
  #cnetplot show genes involved in significant terms. In order to consider the potentially biological complexities in which a gene may belong to multiple annotation categories and provide information of numeric changes if available, cnetplot function extracts the complex association 
  message(paste0("Drawing cnetplot p1 for"))
  p1 <- cnetplot(kks, showCategory = 10) + ggtitle("KEGG Gene-concept Network plot") #category size scaled by gene number 
  ggsave(file.path(KEGGoutdir, "cnetplot.1.pdf"))
  
  ## categorySize can be scaled by 'pvalue' or 'geneNum'
  message(paste0("Drawing cnetplot p2"))
  p2 <- cnetplot(kks, categorySize="pvalue", showCategory = 10,foldChange=genelist) + ggtitle("KEGG Gene-concept Network scaled by pvalue")
  ggsave(file.path(KEGGoutdir, "cnetplot.2.pdf"))
  
  message(paste0("Drawing cnetplot p3"))
  p3 <- cnetplot(kks, circular = TRUE, colorEdge = TRUE, showCategory = 10,foldChange=genelist) + ggtitle("KEGG Gene-concept Network circular")
  ggsave(file.path(KEGGoutdir, "cnetplot.3.pdf"))
  
  #cowplot::plot_grid(p1,p2, p3, ncol=3, labels=LETTERS[1:2])
  #ggsave(file.path(reactomepath, "cnetplot.total.pdf"))
  
  #Heatmap-like functional classification
  p1 <- heatplot(kks, showCategory=10)
  p2 <- heatplot(kk, foldChange=genelist, showCategory=10) #+ theme(axis.text.x = element_text(size = 4, angle = 90))
  cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])
  ggsave(file.path(KEGGoutdir, "heatmap.pdf"))
  
  #Tree plot
  edox2 <- pairwise_termsim(kks)
  p1 <- treeplot(edox2, label_format = 4, showCategory = 10,nCluster=3)
  ggsave(file.path(KEGGoutdir, "clust.wardD2.pdf"))
  p2 <- treeplot(edox2, hclust_method = "average",label_format = 3,nCluster=3)
  ggsave(file.path(KEGGoutdir, "clust.average.pdf"))
  #aplot::plot_list(p1, p2, tag_levels='A')
  #ggsave(file.path(outdir_, "clust.pdf"))
  
  #emapplot show different genesets connected by overlapping genes. Overlapping gene sets are tend to cluster together 
  message(paste0("Drawing emapplot for"))
  p <- emapplot(edox2, showCategory = 10, node_scale = 2)  + ggtitle("KEGG Enrichment Map plot")
  ggsave(file.path(KEGGoutdir, "emapplot.pdf"))
  
  #UpSet Plot
  pdf(file.path(KEGGoutdir,"upset.pdf"))
  upsetplot(kks) + theme(plot.margin = unit(c(1,1,1,4), "cm"), axis.text.y = element_text(size = 15/.pt))
  dev.off()
  
}


## Module KEGG enrichment analysis
MKEGGoutdir <- paste(file.path(outputdir, "MKEGG"))
cmd <- paste("mkdir -p", MKEGGoutdir)
system(cmd)

mkk <- enrichMKEGG(gene = names(genelist),
                   organism = 'hsa',
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   minGSSize = 2, 
                   maxGSSize = 500)                   

mkks <- setReadable(mkk, 'org.Hs.eg.db', 'ENTREZID')

#convert entrez to symbol
#symbol <- clusterProfiler::bitr(names(genelist), fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
#kk@gene2Symbol <- symbol$SYMBOL

mkks_table <- mkks@result
write.table(mkks_table, file.path(MKEGGoutdir,"MKEGG.tsv"), sep = "\t", quote = FALSE)
write.csv(mkks_table, file.path(MKEGGoutdir,"MKEGG.csv"), quote = FALSE)
write.xlsx(mkks_table, file.path(MKEGGoutdir,"MKEGG.xlsx"))

if (length(which(mkks_table$p.adjust < 0.05)) > 10){
  
  ##Plot
  
  #barplot
  message(paste0("Drawing barplot"))
  barplot(mkks, showCategory = 20, colorby = 'pvalue') + ggtitle("KEGG Barplot")
  ggsave(file.path(MKEGGoutdir,"barplot.pdf"))
  
  #dotplot show most significant enriched terms 
  message(paste0("Drawing dotplot"))
  dotplot(mkks, showCategory = 20) + ggtitle("KEGG Dotplot")
  ggsave(file.path(MKEGGoutdir, "dotplot.pdf"))
  
  #cnetplot show genes involved in significant terms. In order to consider the potentially biological complexities in which a gene may belong to multiple annotation categories and provide information of numeric changes if available, cnetplot function extracts the complex association 
  message(paste0("Drawing cnetplot p1 for"))
  p1 <- cnetplot(mkks, showCategory = 10) + ggtitle("KEGG Gene-concept Network plot") #category size scaled by gene number 
  ggsave(file.path(MKEGGoutdir, "cnetplot.1.pdf"))
  
  ## categorySize can be scaled by 'pvalue' or 'geneNum'
  message(paste0("Drawing cnetplot p2"))
  p2 <- cnetplot(mkks, categorySize="pvalue", showCategory = 10,foldChange=genelist) + ggtitle("KEGG Gene-concept Network scaled by pvalue")
  ggsave(file.path(MKEGGoutdir, "cnetplot.2.pdf"))
  
  message(paste0("Drawing cnetplot p3"))
  p3 <- cnetplot(mkks, circular = TRUE, colorEdge = TRUE, showCategory = 10,foldChange=genelist) + ggtitle("KEGG Gene-concept Network circular")
  ggsave(file.path(MKEGGoutdir, "cnetplot.3.pdf"))
  
  #cowplot::plot_grid(p1,p2, p3, ncol=3, labels=LETTERS[1:2])
  #ggsave(file.path(reactomepath, "cnetplot.total.pdf"))
  
  #Heatmap-like functional classification
  p1 <- heatplot(mkks, showCategory=10)
  p2 <- heatplot(mkks, foldChange=genelist, showCategory=10) #+ theme(axis.text.x = element_text(size = 4, angle = 90))
  cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])
  ggsave(file.path(MKEGGoutdir, "heatmap.pdf"))
  
  #Tree plot
  edox2 <- pairwise_termsim(mkks)
  p1 <- treeplot(edox2, label_format = 4, showCategory = 50,nCluster=7)
  ggsave(file.path(MKEGGoutdir, "clust.wardD2.pdf"))
  p2 <- treeplot(edox2, hclust_method = "average",label_format = 4)
  ggsave(file.path(MKEGGoutdir, "clust.average.pdf"))
  #aplot::plot_list(p1, p2, tag_levels='A')
  #ggsave(file.path(outdir_, "clust.pdf"))
  
  #emapplot show different genesets connected by overlapping genes. Overlapping gene sets are tend to cluster together 
  message(paste0("Drawing emapplot for"))
  emapplot(edox2, showCategory = 10, node_scale = 2)  + ggtitle("KEGG Enrichment Map plot")
  ggsave(file.path(MKEGGoutdir, "emapplot.pdf"))
  
  #UpSet Plot
  pdf(file.path(MKEGGoutdir,"upset.pdf"))
  upsetplot(mkks) + theme(plot.margin = unit(c(1,1,1,4), "cm"), axis.text.y = element_text(size = 15/.pt))
  dev.off()  
}


#WikiPathways analysis
wikipath <- paste(file.path(outputdir, "Wiki"))
cmd <- paste("mkdir -p", wikipath)
system(cmd)

## Wiki
wiki <- enrichWP(names(genelist), 
                 organism = "Homo sapiens", 
                 pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.05, 
                 minGSSize = 5,
                 maxGSSize = 500)
wiki <- setReadable(wiki, 'org.Hs.eg.db', 'ENTREZID')

wiki_table <- as.data.frame(wiki@result) 
write.table(wiki_table,file.path(wikipath,"wiki_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(wiki_table,file.path(wikipath,"wiki_table.csv"), row.names = FALSE, quote = FALSE)
#write.xlsx(wiki_table,file.path(wikipath,"wiki_table.xlsx"))

if (length(which(wiki_table$p.adjust < 0.05)) > 10){
  
  ##Plot
  
  #barplot
  message(paste0("Drawing barplot"))
  p <-barplot(wiki, showCategory = 20, colorby = 'pvalue') + ggtitle("Wiki Barplot")
  ggsave(file.path(wikipath,"barplot.pdf"))
  
  #dotplot show most significant enriched terms 
  message(paste0("Drawing dotplot"))
  p <-dotplot(wiki, showCategory = 20) + ggtitle("Wiki Dotplot")
  ggsave(file.path(wikipath, "dotplot.pdf"))
  
  #cnetplot show genes involved in significant terms. In order to consider the potentially biological complexities in which a gene may belong to multiple annotation categories and provide information of numeric changes if available, cnetplot function extracts the complex association 
  message(paste0("Drawing cnetplot p1 for"))
  p1 <- cnetplot(wiki, showCategory = 10) + ggtitle("Wiki Gene-concept Network plot") #category size scaled by gene number 
  ggsave(file.path(wikipath, "cnetplot.1.pdf"))
  
  ## categorySize can be scaled by 'pvalue' or 'geneNum'
  message(paste0("Drawing cnetplot p2"))
  p2 <- cnetplot(wiki, categorySize="pvalue", showCategory = 10,foldChange=genelist) + ggtitle("Wiki Gene-concept Network scaled by pvalue")
  ggsave(file.path(wikipath, "cnetplot.2.pdf"))
  
  message(paste0("Drawing cnetplot p3"))
  p3 <- cnetplot(wiki, circular = TRUE, colorEdge = TRUE, showCategory = 10,foldChange=genelist) + ggtitle("Wiki Gene-concept Network circular")
  ggsave(file.path(wikipath, "cnetplot.3.pdf"))
  
  #cowplot::plot_grid(p1,p2, p3, ncol=3, labels=LETTERS[1:2])
  #ggsave(file.path(reactomepath, "cnetplot.total.pdf"))
  
  #Heatmap-like functional classification
  p1 <- heatplot(wiki, showCategory=10)
  p2 <- heatplot(wiki, foldChange=genelist, showCategory=10)# + theme(axis.text.x = element_text(size = 4, angle = 90))
  cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])
  ggsave(file.path(wikipath, "heatmap.pdf"))
  
  #Tree plot
  edox2 <- pairwise_termsim(wiki)
  p1 <- treeplot(edox2, label_format = 4, showCategory = 50)
  ggsave(file.path(wikipath, "clust.wardD2.pdf"))
  p2 <- treeplot(edox2, hclust_method = "average",label_format = 4)
  ggsave(file.path(wikipath, "clust.average.pdf"))
  
  #emapplot show different genesets connected by overlapping genes. Overlapping gene sets are tend to cluster together 
  message(paste0("Drawing emapplot for"))
  p <- emapplot(edox2, showCategory = 10, node_scale = 2)  + ggtitle("Wiki Enrichment Map plot")
  ggsave(file.path(wikipath, "emapplot.pdf"))
  
  #UpSet Plot
  pdf(file.path(wikipath,"upset.pdf"))
  upsetplot(wiki) + theme(plot.margin = unit(c(1,1,1,4), "cm"), axis.text.y = element_text(size = 15/.pt))
  dev.off()   
} 

dftot <- read.csv("/projects2/2024_Chiarugi/Transcriptomestic/results/unstranded/STAR_HTSeq_nosample1/results.csv") 

# Disease Ontology
dopath <- paste(file.path(outputdir, "DO"))
cmd <- paste("mkdir -p", dopath)
system(cmd)

## DO
x <- enrichDO(gene         = names(genelist),
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.05, 
              readable     = T,
              minGSSize    = 2,
              maxGSSize    = 500)

do <- as.data.frame(x@result) 
write.table(do,file.path(dopath,"do.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(do,file.path(dopath,"do.csv"), row.names = FALSE, quote = FALSE)
#write.xlsx(do,file.path(dopath,"do.xlsx"), row.names = FALSE)

if (length(which(do$p.adjust < 0.05)) > 10){
  
  ##Plot
  
  #barplot
  message(paste0("Drawing barplot"))
  barplot(x, showCategory = 10, colorby = 'pvalue') + ggtitle("DO Barplot")
  ggsave(file.path(dopath,"barplot.pdf"))
  
  #dotplot show most significant enriched terms 
  message(paste0("Drawing dotplot"))
  dotplot(x, showCategory = 10) + ggtitle("DO Dotplot")
  ggsave(file.path(dopath, "dotplot.pdf"))
  
  #cnetplot show genes involved in significant terms. In order to consider the potentially biological complexities in which a gene may belong to multiple annotation categories and provide information of numeric changes if available, cnetplot function extracts the complex association 
  message(paste0("Drawing cnetplot p1 for"))
  p1 <- cnetplot(x, showCategory = 10) + ggtitle("DO Gene-concept Network plot") #category size scaled by gene number 
  ggsave(file.path(dopath, "cnetplot.1.pdf"))
  
  ## categorySize can be scaled by 'pvalue' or 'geneNum'
  message(paste0("Drawing cnetplot p2"))
  p2 <- cnetplot(x, categorySize="pvalue", showCategory = 10,foldChange=genelist) + ggtitle("DO Gene-concept Network scaled by pvalue")
  ggsave(file.path(dopath, "cnetplot.2.pdf"))
  
  message(paste0("Drawing cnetplot p3"))
  p3 <- cnetplot(x, circular = TRUE, colorEdge = TRUE, showCategory = 10,foldChange=genelist) + ggtitle("DO Gene-concept Network circular")
  ggsave(file.path(dopath, "cnetplot.3.pdf"))
  
  #cowplot::plot_grid(p1,p2, p3, ncol=3, labels=LETTERS[1:2])
  #ggsave(file.path(reactomepath, "cnetplot.total.pdf"))
  
  #Heatmap-like functional classification
  p1 <- heatplot(x, showCategory=10)
  p2 <- heatplot(x, foldChange=genelist, showCategory=10) + theme(axis.text.x = element_text(size = 4, angle = 90))
  cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])
  ggsave(file.path(dopath, "heatmap.pdf"))
  
  #Tree plot
  edox2 <- pairwise_termsim(x)
  p1 <- treeplot(edox2, label_format = 4, showCategory = 50,nCluster=7)
  ggsave(file.path(dopath, "clust.wardD2.pdf"))
  p2 <- treeplot(edox2, hclust_method = "average",label_format = 4)
  ggsave(file.path(dopath, "clust.average.pdf"))
  
  #emapplot show different genesets connected by overlapping genes. Overlapping gene sets are tend to cluster together 
  message(paste0("Drawing emapplot for"))
  emapplot(edox2, showCategory = 10, node_scale = 2)  + ggtitle("DO Enrichment Map plot")
  ggsave(file.path(dopath, "emapplot.pdf"))
  
  #UpSet Plot
  pdf(file.path(dopath,"upset.pdf"))
  upsetplot(x) + theme(plot.margin = unit(c(1,1,1,4), "cm"), axis.text.y = element_text(size = 15/.pt))
  dev.off() 
  
}


# Medical Subject Headings
library(AnnotationHub)
meshpath <- paste(file.path(outputdir, "MeSH"))
cmd <- paste("mkdir -p", meshpath)
system(cmd)

ah <- AnnotationHub(localHub=FALSE)
hsa <- query(ah, c("MeSHDb", "Homo sapiens"))
file_hsa <- hsa[[1]]
db <- MeSHDbi::MeSHDb(file_hsa)

## MeSH
x <- enrichMeSH(names(genelist), MeSHDb = db, database='gendoo', category = 'C',minGSSize = 2,maxGSSize = 500, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
x <- setReadable(x, 'org.Hs.eg.db', 'ENTREZID')

MeSH <- as.data.frame(x@result) 
write.table(MeSH,file.path(meshpath,"mesh.gendoo.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(MeSH,file.path(meshpath,"mesh.gendoo.csv"), row.names = FALSE, quote = FALSE)

if (length(which(x$p.adjust < 0.05)) > 10){
  
  ##Plot
  
  #barplot
  message(paste0("Drawing barplot"))
  p <- barplot(x, showCategory = 10, colorby = 'pvalue') + ggtitle("MeSH gendoo Barplot")
  ggsave(file.path(dopath,"barplot.pdf"))
  
  #dotplot show most significant enriched terms 
  message(paste0("Drawing dotplot"))
  dotplot(x, showCategory = 10) + ggtitle("MeSH gendoo Dotplot")
  ggsave(file.path(dopath, "dotplot.pdf"))
  
  #cnetplot show genes involved in significant terms. In order to consider the potentially biological complexities in which a gene may belong to multiple annotation categories and provide information of numeric changes if available, cnetplot function extracts the complex association 
  message(paste0("Drawing cnetplot p1 for"))
  p1 <- cnetplot(x, showCategory = 10) + ggtitle("MeSH gendoo Gene-concept Network plot") #category size scaled by gene number 
  ggsave(file.path(dopath, "cnetplot.1.pdf"))
  
  ## categorySize can be scaled by 'pvalue' or 'geneNum'
  message(paste0("Drawing cnetplot p2"))
  p2 <- cnetplot(x, categorySize="pvalue", showCategory = 10,foldChange=genelist) + ggtitle("MeSH gendoo Gene-concept Network scaled by pvalue")
  ggsave(file.path(dopath, "cnetplot.2.pdf"))
  
  message(paste0("Drawing cnetplot p3"))
  p3 <- cnetplot(x, circular = TRUE, colorEdge = TRUE, showCategory = 10,foldChange=genelist) + ggtitle("MeSH gendoo Gene-concept Network circular")
  ggsave(file.path(dopath, "cnetplot.3.pdf"))
  
  #cowplot::plot_grid(p1,p2, p3, ncol=3, labels=LETTERS[1:2])
  #ggsave(file.path(reactomepath, "cnetplot.total.pdf"))
  
  #Heatmap-like functional classification
  p1 <- heatplot(x, showCategory=5) + theme(axis.text.x = element_text(size = 3, angle = 90))
  p2 <- heatplot(x, foldChange=genelist, showCategory=5) + theme(axis.text.x = element_text(size = 3, angle = 90))
  cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])
  ggsave(file.path(dopath, "heatmap.pdf"), height = 12, width=20)
  
  #Tree plot
  edox2 <- pairwise_termsim(x)
  p1 <- treeplot(edox2, label_format = 4, showCategory = 10,nCluster=5)
  ggsave(file.path(dopath, "clust.wardD2.pdf"))
  p2 <- treeplot(edox2, hclust_method = "average",label_format = 4, showCategory = 10,nCluster=5)
  ggsave(file.path(dopath, "clust.average.pdf"))
  
  #emapplot show different genesets connected by overlapping genes. Overlapping gene sets are tend to cluster together 
  message(paste0("Drawing emapplot for"))
  p <- emapplot(edox2, showCategory = 10, node_scale = 2)  + ggtitle("MeSH gendoo Enrichment Map plot")
  ggsave(file.path(dopath, "emapplot.pdf"))
  
  #UpSet Plot
  pdf(file.path(dopath,"upset.pdf"))
  upsetplot(x) + theme(plot.margin = unit(c(1,1,1,4), "cm"), axis.text.y = element_text(size = 15/.pt))
  dev.off()   
} 

## MeSH
#x <- enrichMeSH(names(genelist), MeSHDb = db, database='gene2pubmed', category = 'C',minGSSize = 2,maxGSSize = 500)

#do <- as.data.frame(x@result) 
#write.table(x,file.path(meshpath,"mesh.gene2pubmed.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
#write.csv(x,file.path(meshpath,"mesh.gene2pubmed.csv"), row.names = FALSE, quote = FALSE)


#if (length(which(x$p.adjust < 0.05)) > 10){
  
  ##Plot
  
  #barplot
#  message(paste0("Drawing barplot"))
#  barplot(x, showCategory = 10, colorby = 'pvalue') + ggtitle("MeSH gene2pubmed Barplot")
#  ggsave(file.path(dopath,"barplot.pdf"))
  
  #dotplot show most significant enriched terms 
  #message(paste0("Drawing dotplot"))
  #dotplot(x, showCategory = 10) + ggtitle("MeSH gene2pubmed Dotplot")
  #ggsave(file.path(dopath, "dotplot.pdf"))
  
  #cnetplot show genes involved in significant terms. In order to consider the potentially biological complexities in which a gene may belong to multiple annotation categories and provide information of numeric changes if available, cnetplot function extracts the complex association 
  #message(paste0("Drawing cnetplot p1 for"))
  #p1 <- cnetplot(x, showCategory = 10) + ggtitle("MeSH gene2pubmed Gene-concept Network plot") #category size scaled by gene number 
  #ggsave(file.path(dopath, "cnetplot.1.pdf"))
  
  ## categorySize can be scaled by 'pvalue' or 'geneNum'
  #message(paste0("Drawing cnetplot p2"))
  #p2 <- cnetplot(x, categorySize="pvalue", showCategory = 10,foldChange=genelist) + ggtitle("MeSH gene2pubmed Gene-concept Network scaled by pvalue")
  #ggsave(file.path(dopath, "cnetplot.2.pdf"))
  
  #message(paste0("Drawing cnetplot p3"))
  #p3 <- cnetplot(x, circular = TRUE, colorEdge = TRUE, showCategory = 10,foldChange=genelist) + ggtitle("MeSH gene2pubmed Gene-concept Network circular")
  #ggsave(file.path(dopath, "cnetplot.3.pdf"))
  
  #cowplot::plot_grid(p1,p2, p3, ncol=3, labels=LETTERS[1:2])
  #ggsave(file.path(reactomepath, "cnetplot.total.pdf"))
  
  #Heatmap-like functional classification
  #p1 <- heatplot(x, showCategory=10)
  #p2 <- heatplot(x, foldChange=genelist, showCategory=10) #+ theme(axis.text.x = element_text(size = 4, angle = 90))
  #cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])
  #ggsave(file.path(dopath, "heatmap.pdf"))
  
  #Tree plot
  #edox2 <- pairwise_termsim(x)
  #p1 <- treeplot(edox2, label_format = 4, showCategory = 50,nCluster=7)
  #ggsave(file.path(dopath, "clust.wardD2.pdf"))
  #p2 <- treeplot(edox2, hclust_method = "average",label_format = 4)
  #ggsave(file.path(dopath, "clust.average.pdf"))
  
  #emapplot show different genesets connected by overlapping genes. Overlapping gene sets are tend to cluster together 
  #message(paste0("Drawing emapplot for"))
  #emapplot(edox2, showCategory = 10, node_scale = 2)  + ggtitle("MeSH gene2pubmed Enrichment Map plot")
  #ggsave(file.path(dopath, "emapplot.pdf"))
  
  #UpSet Plot
  #pdf(file.path(dopath,"upset.pdf"))
  #upsetplot(x) + theme(plot.margin = unit(c(1,1,1,4), "cm"), axis.text.y = element_text(size = 15/.pt))
  #dev.off()   
  
#}

#msigdbr
all_gene_sets = msigdbr(species = "Homo sapiens")
#head(all_gene_sets)

m_t2g_H <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

#Msigdbr analysis
Msigdbrpath <- paste(file.path(outputdir, "msigdbr"))
cmd <- paste("mkdir -p", Msigdbrpath)
system(cmd)

## H
H <- enricher(names(genelist), TERM2GENE = m_t2g_H, minGSSize = 2,maxGSSize = 500,pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
H <- setReadable(H, 'org.Hs.eg.db', 'ENTREZID')

Hpath <- paste(file.path(Msigdbrpath, "H"))
cmd <- paste("mkdir -p", Hpath)
system(cmd)

H_table <- as.data.frame(H@result) 
write.table(H_table,file.path(Hpath,"H_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(H_table,file.path(Hpath,"H_table.csv"), row.names = FALSE, quote = FALSE)
write.xlsx(H_table,file.path(Hpath,"H_table.xlsx"), row.names = FALSE)


if (length(which(H$p.adjust < 0.05)) > 10){
  
  ##Plot
  
  #barplot
  message(paste0("Drawing barplot"))
  p <- barplot(H, showCategory = 10, colorby = 'pvalue') + ggtitle("HALLMARK Barplot")
  ggsave(file.path(Hpath,"barplot.pdf"))
  
  #dotplot show most significant enriched terms 
  message(paste0("Drawing dotplot"))
  p <- dotplot(H, showCategory = 10) + ggtitle("HALLMARK Dotplot")
  ggsave(file.path(Hpath, "dotplot.pdf"))
  
  #cnetplot show genes involved in significant terms. In order to consider the potentially biological complexities in which a gene may belong to multiple annotation categories and provide information of numeric changes if available, cnetplot function extracts the complex association 
  message(paste0("Drawing cnetplot p1 for"))
  p1 <- cnetplot(H, showCategory = 10) + ggtitle("HALLMARK Gene-concept Network plot") #category size scaled by gene number 
  ggsave(file.path(Hpath, "cnetplot.1.pdf"))
  
  ## categorySize can be scaled by 'pvalue' or 'geneNum'
  message(paste0("Drawing cnetplot p2"))
  p2 <- cnetplot(H, categorySize="pvalue", showCategory = 10,foldChange=genelist) + ggtitle("HALLMARK Gene-concept Network scaled by pvalue")
  ggsave(file.path(Hpath, "cnetplot.2.pdf"))
  
  message(paste0("Drawing cnetplot p3"))
  p3 <- cnetplot(H, circular = TRUE, colorEdge = TRUE, showCategory = 10,foldChange=genelist) + ggtitle("HALLMARK Gene-concept Network circular")
  ggsave(file.path(Hpath, "cnetplot.3.pdf"))
  
  #cowplot::plot_grid(p1,p2, p3, ncol=3, labels=LETTERS[1:2])
  #ggsave(file.path(reactomepath, "cnetplot.total.pdf"))
  
  #Heatmap-like functional classification
  p1 <- heatplot(H, showCategory=10) + theme(axis.text.x = element_text(size = 2.5, angle = 90))
  p2 <- heatplot(H, foldChange=genelist, showCategory=10) + theme(axis.text.x = element_text(size = 2.5, angle = 90))
  cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])
  ggsave(file.path(Hpath, "heatmap.pdf"), height = 12, width=20)
  
  #Tree plot
  edox2 <- pairwise_termsim(H)
  p1 <- treeplot(edox2, label_format = 4, showCategory = 10, fontsize=3,nCluster = 3)
  ggsave(file.path(Hpath, "clust.wardD2.pdf"))
  p2 <- treeplot(edox2, hclust_method = "average",label_format = 4, showCategory = 10, fontsize=3,nCluster = 3)
  ggsave(file.path(Hpath, "clust.average.pdf"))
  
  #emapplot show different genesets connected by overlapping genes. Overlapping gene sets are tend to cluster together 
  message(paste0("Drawing emapplot for"))
  p <- emapplot(edox2, showCategory = 10, node_scale = 2)  + ggtitle("HALLMARK Enrichment Map plot")
  ggsave(file.path(Hpath, "emapplot.pdf"))
  
  #UpSet Plot
  pdf(file.path(dopath,"upset.pdf"))
  upsetplot(x) + theme(plot.margin = unit(c(1,1,1,4), "cm"), axis.text.y = element_text(size = 15/.pt))
  dev.off()   
  
} 


## Network of Cancer Gene (NCG) 
NCGpath <- paste(file.path(outputdir, "NCG"))
cmd <- paste("mkdir -p", NCGpath)
system(cmd)


x <- enrichNCG(gene=names(genelist),pvalueCutoff=0.05, qvalueCutoff=0.05, readable=T,minGSSize = 2,maxGSSize = 500)

NGC <- as.data.frame(x@result) 
write.table(NGC,file.path(NCGpath,"NGC_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(NGC,file.path(NCGpath,"NGC_table.csv"), row.names = FALSE, quote = FALSE)
write.xlsx(NGC,file.path(NCGpath,"NGC_table.xlsx"))

if (length(which(x@result$p.adjust < 0.05)) > 10){
  
  ##Plot
  
  #barplot
  message(paste0("Drawing barplot"))
  p <- barplot(x, showCategory = 20, colorby = 'pvalue') + ggtitle("NGC Barplot")
  ggsave(file.path(NCGpath,"barplot.pdf"))
  
  #dotplot show most significant enriched terms 
  message(paste0("Drawing dotplot"))
  p <- dotplot(x, showCategory = 20) + ggtitle("NGC Dotplot")
  ggsave(file.path(NCGpath, "dotplot.pdf"))
  
  #cnetplot show genes involved in significant terms. In order to consider the potentially biological complexities in which a gene may belong to multiple annotation categories and provide information of numeric changes if available, cnetplot function extracts the complex association 
  message(paste0("Drawing cnetplot p1 for"))
  p1 <- cnetplot(x, showCategory = 10) + ggtitle("NGC Gene-concept Network plot") #category size scaled by gene number 
  ggsave(file.path(NCGpath, "cnetplot.1.pdf"))
  
  ## categorySize can be scaled by 'pvalue' or 'geneNum'
  message(paste0("Drawing cnetplot p2"))
  p2 <- cnetplot(x, categorySize="pvalue", showCategory = 10,foldChange=genelist) + ggtitle("NGC Gene-concept Network scaled by pvalue")
  ggsave(file.path(NCGpath, "cnetplot.2.pdf"))
  
  message(paste0("Drawing cnetplot p3"))
  p3 <- cnetplot(x, circular = TRUE, colorEdge = TRUE, showCategory = 10,foldChange=genelist) + ggtitle("NGC Gene-concept Network circular")
  ggsave(file.path(NCGpath, "cnetplot.3.pdf"))
  
  #cowplot::plot_grid(p1,p2, p3, ncol=3, labels=LETTERS[1:2])
  #ggsave(file.path(reactomepath, "cnetplot.total.pdf"))
  
  #Heatmap-like functional classification
  p1 <- heatplot(x, showCategory=10)
  p2 <- heatplot(x, foldChange=genelist, showCategory=10) #+ theme(axis.text.x = element_text(size = 4, angle = 90))
  cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])
  ggsave(file.path(NCGpath, "heatmap.pdf"))
  
  #Tree plot
  edox2 <- pairwise_termsim(x)
  p1 <- treeplot(edox2, label_format = 4, showCategory = 50, nCluster=7)
  ggsave(file.path(NCGpath, "clust.wardD2.pdf"))
  p2 <- treeplot(edox2, hclust_method = "average",label_format = 4)
  ggsave(file.path(NCGpath, "clust.average.pdf"))
  #aplot::plot_list(p1, p2, tag_levels='A')
  #ggsave(file.path(outdir_, "clust.pdf"))
  
  #emapplot show different genesets connected by overlapping genes. Overlapping gene sets are tend to cluster together 
  message(paste0("Drawing emapplot for"))
  p <- emapplot(edox2, showCategory = 20, node_scale = 2)  + ggtitle("NGC Enrichment Map plot")
  ggsave(file.path(NCGpath, "emapplot.pdf"))
  
  #UpSet Plot
  pdf(file.path(dopath,"upset.pdf"))
  upsetplot(x) + theme(plot.margin = unit(c(1,1,1,4), "cm"), axis.text.y = element_text(size = 15/.pt))
  dev.off()   
}

#disease gene network (DGN)
DGNpath <- paste(file.path(outputdir, "DGN"))
cmd <- paste("mkdir -p", DGNpath)
system(cmd)


x <- enrichDGN(gene=names(genelist),pvalueCutoff=0.05, qvalueCutoff=0.05, readable=T,minGSSize = 2,maxGSSize = 500)

DGN <- as.data.frame(x@result) 
write.table(DGN,file.path(DGNpath,"DGN_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(DGN,file.path(DGNpath,"DGN_table.csv"), row.names = FALSE, quote = FALSE)
write.xlsx(DGN,file.path(DGNpath,"DGN_table.xlsx"))

if (length(which(x@result$p.adjust < 0.05)) > 10){
  
  ##Plot
  
  #barplot
  message(paste0("Drawing barplot"))
  p <- barplot(x, showCategory = 10, colorby = 'pvalue') + ggtitle("DGN Barplot")
  ggsave(file.path(DGNpath,"barplot.pdf"))
  
  #dotplot show most significant enriched terms 
  message(paste0("Drawing dotplot"))
  p <- dotplot(x, showCategory = 10) + ggtitle("DGN Dotplot")
  ggsave(file.path(DGNpath, "dotplot.pdf"))
  
  #cnetplot show genes involved in significant terms. In order to consider the potentially biological complexities in which a gene may belong to multiple annotation categories and provide information of numeric changes if available, cnetplot function extracts the complex association 
  message(paste0("Drawing cnetplot p1 for"))
  p1 <- cnetplot(x, showCategory = 10) + ggtitle("DGN Gene-concept Network plot") #category size scaled by gene number 
  ggsave(file.path(DGNpath, "cnetplot.1.pdf"))
  
  ## categorySize can be scaled by 'pvalue' or 'geneNum'
  message(paste0("Drawing cnetplot p2"))
  p2 <- cnetplot(x, categorySize="pvalue", showCategory = 10,foldChange=genelist) + ggtitle("DGN Gene-concept Network scaled by pvalue")
  ggsave(file.path(DGNpath, "cnetplot.2.pdf"))
  
  message(paste0("Drawing cnetplot p3"))
  p3 <- cnetplot(x, circular = TRUE, colorEdge = TRUE, showCategory = 10,foldChange=genelist) + ggtitle("DGN Gene-concept Network circular")
  ggsave(file.path(DGNpath, "cnetplot.3.pdf"))
  
  #cowplot::plot_grid(p1,p2, p3, ncol=3, labels=LETTERS[1:2])
  #ggsave(file.path(reactomepath, "cnetplot.total.pdf"))
  
  #Heatmap-like functional classification
  p1 <- heatplot(x, showCategory=10) + theme(axis.text.x = element_text(size = 2, angle = 90))
  p2 <- heatplot(x, foldChange=genelist, showCategory=10) + theme(axis.text.x = element_text(size = 2, angle = 90))
  cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])
  ggsave(file.path(DGNpath, "heatmap.pdf"), height=12, width=20)
  
  #Tree plot
  edox2 <- pairwise_termsim(x)
  p1 <- treeplot(edox2, label_format = 4, showCategory = 10, nCluster=5)
  ggsave(file.path(DGNpath, "clust.wardD2.pdf"))
  p2 <- treeplot(edox2, hclust_method = "average", showCategory = 10, nCluster=5)
  ggsave(file.path(DGNpath, "clust.average.pdf"))
  #aplot::plot_list(p1, p2, tag_levels='A')
  #ggsave(file.path(outdir_, "clust.pdf"))
  
  #emapplot show different genesets connected by overlapping genes. Overlapping gene sets are tend to cluster together 
  message(paste0("Drawing emapplot for"))
  p <- emapplot(edox2, showCategory = 10, node_scale = 2)  + ggtitle("DGN Enrichment Map plot")
  ggsave(file.path(DGNpath, "emapplot.pdf"))
  
  #UpSet Plot
  pdf(file.path(dopath,"upset.pdf"))
  upsetplot(x) + theme(plot.margin = unit(c(1,1,1,4), "cm"), axis.text.y = element_text(size = 15/.pt))
  dev.off()   
}

#C2

m_t2g_C2 <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)

## C2
C2 <- enricher(names(genelist), TERM2GENE = m_t2g_C2, minGSSize = 2,maxGSSize = 500,pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
C2 <- setReadable(C2, 'org.Hs.eg.db', 'ENTREZID')

C2path <- paste(file.path(Msigdbrpath, "C2"))
cmd <- paste("mkdir -p", C2path)
system(cmd)

C2_table <- as.data.frame(C2@result) 
write.table(C2_table,file.path(C2path,"C2_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(C2_table,file.path(C2path,"C2_table.csv"), row.names = FALSE, quote = FALSE)
write.xlsx(C2_table,file.path(C2path,"C2_table.xlsx"), row.names = FALSE)


if (length(which(C2$p.adjust < 0.05)) > 10){
  
  ##Plot
  
  #barplot
  message(paste0("Drawing barplot"))
  p <- barplot(C2, showCategory = 10, colorby = 'pvalue') + ggtitle("C2 Barplot")
  ggsave(file.path(C2path,"barplot.pdf"))
  
  #dotplot show most significant enriched terms 
  message(paste0("Drawing dotplot"))
  p <- dotplot(C2, showCategory = 10) + ggtitle("C2 Dotplot")
  ggsave(file.path(C2path, "dotplot.pdf"))
  
  #cnetplot show genes involved in significant terms. In order to consider the potentially biological complexities in which a gene may belong to multiple annotation categories and provide information of numeric changes if available, cnetplot function extracts the complex association 
  message(paste0("Drawing cnetplot p1 for"))
  p1 <- cnetplot(C2, showCategory = 10) + ggtitle("C2 Gene-concept Network plot") #category size scaled by gene number 
  ggsave(file.path(C2path, "cnetplot.1.pdf"))
  
  ## categorySize can be scaled by 'pvalue' or 'geneNum'
  message(paste0("Drawing cnetplot p2"))
  p2 <- cnetplot(C2, categorySize="pvalue", showCategory = 10,foldChange=genelist) + ggtitle("C2 Gene-concept Network scaled by pvalue")
  ggsave(file.path(C2path, "cnetplot.2.pdf"))
  
  message(paste0("Drawing cnetplot p3"))
  p3 <- cnetplot(C2, circular = TRUE, colorEdge = TRUE, showCategory = 10,foldChange=genelist) + ggtitle("C2 Gene-concept Network circular")
  ggsave(file.path(C2path, "cnetplot.3.pdf"))
  
  #cowplot::plot_grid(p1,p2, p3, ncol=3, labels=LETTERS[1:2])
  #ggsave(file.path(reactomepath, "cnetplot.total.pdf"))
  
  #Heatmap-like functional classification
  p1 <- heatplot(C2, showCategory=5) + theme(axis.text.x = element_text(size = 2.5, angle = 90))
  p2 <- heatplot(C2, foldChange=genelist, showCategory=5) + theme(axis.text.x = element_text(size = 2.5, angle = 90))
  cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])
  ggsave(file.path(C2path, "heatmap.pdf"), height = 12, width = 25)
  
  #Tree plot
  edox2 <- pairwise_termsim(C2)
  p1 <- treeplot(edox2, label_format = 4, fontsize=2.5,showCategory = 10,nCluster=4)
  ggsave(file.path(C2path, "clust.wardD2.pdf"), height=12, width=12)
  p2 <- treeplot(edox2, hclust_method = "average",label_format = 4, fontsize=2.5,showCategory = 10,nCluster=4)
  ggsave(file.path(C2path, "clust.average.pdf"), height=12, width=12)
  
  #emapplot show different genesets connected by overlapping genes. Overlapping gene sets are tend to cluster together 
  message(paste0("Drawing emapplot for"))
  p <- emapplot(edox2, showCategory = 10, node_scale = 2)  + ggtitle("C2 Enrichment Map plot")
  ggsave(file.path(C2path, "emapplot.pdf"))
  
  #UpSet Plot
  pdf(file.path(dopath,"upset.pdf"))
  upsetplot(x) + theme(plot.margin = unit(c(1,1,1,4), "cm"), axis.text.y = element_text(size = 15/.pt))
  dev.off()   
}
#C6

m_t2g_C6 <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::select(gs_name, entrez_gene)

## C6
C6 <- enricher(names(genelist), TERM2GENE = m_t2g_C6, minGSSize = 2,maxGSSize = 500,pvalueCutoff  = 0.05, qvalueCutoff  = 0.05)
C6 <- setReadable(C6, 'org.Hs.eg.db', 'ENTREZID')

C6path <- paste(file.path(Msigdbrpath, "C6"))
cmd <- paste("mkdir -p", C6path)
system(cmd)

C6_table <- as.data.frame(C6@result) 
write.table(C6_table,file.path(C6path,"C6_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(C6_table,file.path(C6path,"C6_table.csv"), row.names = FALSE, quote = FALSE)
write.xlsx(C6_table,file.path(C6path,"C6_table.xlsx"), row.names = FALSE)


if (length(which(C6$p.adjust < 0.05)) > 10){
  
  ##Plot
  
  #barplot
  message(paste0("Drawing barplot"))
  p <- barplot(C6, showCategory = 10, colorby = 'pvalue') + ggtitle("C6 Barplot")
  ggsave(file.path(C6path,"barplot.pdf"))
  
  #dotplot show most significant enriched terms 
  message(paste0("Drawing dotplot"))
  p <- dotplot(C6, showCategory = 10) + ggtitle("C6 Dotplot")
  ggsave(file.path(C6path, "dotplot.pdf"))
  
  #cnetplot show genes involved in significant terms. In order to consider the potentially biological complexities in which a gene may belong to multiple annotation categories and provide information of numeric changes if available, cnetplot function extracts the complex association 
  message(paste0("Drawing cnetplot p1 for"))
  p1 <- cnetplot(C6, showCategory = 10) + ggtitle("C6 Gene-concept Network plot") #category size scaled by gene number 
  ggsave(file.path(C6path, "cnetplot.1.pdf"))
  
  ## categorySize can be scaled by 'pvalue' or 'geneNum'
  message(paste0("Drawing cnetplot p2"))
  p2 <- cnetplot(C6, categorySize="pvalue", showCategory = 10,foldChange=genelist) + ggtitle("C6 Gene-concept Network scaled by pvalue")
  ggsave(file.path(C6path, "cnetplot.2.pdf"))
  
  message(paste0("Drawing cnetplot p3"))
  p3 <- cnetplot(C6, circular = TRUE, colorEdge = TRUE, showCategory = 10,foldChange=genelist) + ggtitle("C6 Gene-concept Network circular")
  ggsave(file.path(C6path, "cnetplot.3.pdf"))
  
  #cowplot::plot_grid(p1,p2, p3, ncol=3, labels=LETTERS[1:2])
  #ggsave(file.path(reactomepath, "cnetplot.total.pdf"))
  
  #Heatmap-like functional classification
  p1 <- heatplot(C6, showCategory=5) + theme(axis.text.x = element_text(size = 2.5, angle = 90))
  p2 <- heatplot(C6, foldChange=genelist, showCategory=5) + theme(axis.text.x = element_text(size = 2.5, angle = 90))
  cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])
  ggsave(file.path(C6path, "heatmap.pdf"), height = 12, width=20)
  
  #Tree plot
  edox2 <- pairwise_termsim(C6)
  p1 <- treeplot(edox2, label_format = 4, showCategory = 10,nCluster=4)
  ggsave(file.path(C6path, "clust.wardD2.pdf"))
  p2 <- treeplot(edox2, hclust_method = "average",label_format = 4, showCategory=10,nCluster=4)
  ggsave(file.path(C6path, "clust.average.pdf"))
  
  #emapplot show different genesets connected by overlapping genes. Overlapping gene sets are tend to cluster together 
  message(paste0("Drawing emapplot for"))
  p <- emapplot(edox2, showCategory = 10, node_scale = 2)  + ggtitle("C6 Enrichment Map plot")
  ggsave(file.path(C6path, "emapplot.pdf"))
  
  #UpSet Plot
  pdf(file.path(dopath,"upset.pdf"))
  upsetplot(x) + theme(plot.margin = unit(c(1,1,1,4), "cm"), axis.text.y = element_text(size = 15/.pt))
  dev.off()   
}