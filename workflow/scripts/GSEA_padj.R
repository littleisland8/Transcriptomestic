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
library(biomaRt)

rm(list = ls())

outputdir <- "/projects2/2024_Chiarugi/Transcriptomestic/results/unstranded/STAR_HTSeq_nosample1/GSEA_padj"
cmd <- paste("mkdir -p", outputdir)
system(cmd)

###### GSEA
df <- read.csv("/projects2/2024_Chiarugi/Transcriptomestic/results/unstranded/STAR_HTSeq_nosample1/results.filtered.padj.csv")

df <- df[which(df$padj < 0.05 | df$padj==0),]

genelist <- df$log2FoldChange
names(genelist) <- df$entrez

genelist <- genelist[order(genelist, decreasing=TRUE)]

## GO classification
GOoutdir <- paste(file.path(outputdir, "GO"))
cmd <- paste("mkdir -p", GOoutdir)
system(cmd)

#CC
cmd <- paste("mkdir -p", file.path(GOoutdir,"CC"))
system(cmd)

gseaCC <- gseGO(geneList     = genelist,
               OrgDb        = org.Hs.eg.db,
               ont          = "CC",
               minGSSize    = 2,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               verbose      = TRUE,
               by           = "DOSE",
               nPerm        = 10000,
               seed         = TRUE)

gseaCC <- setReadable(gseaCC, 'org.Hs.eg.db', 'ENTREZID')

gsea_tableCC <- gseaCC@result
write.table(gsea_tableCC,file.path(GOoutdir, "CC","gsea.table.txt"), sep = "\t", quote = FALSE)
write.csv(gsea_tableCC,file.path(GOoutdir, "CC","gsea.table.csv"), quote = FALSE)
#write.xlsx(gsea_tableCC,file.path(GOoutdir, "CC","gsea.table.xlsx"))

p <- ridgeplot(gseaCC,core_enrichment = TRUE)
ggsave(file.path(GOoutdir, "CC","ridgeplot.pdf"), height = 12, width = 10)

#p1 <- gseaplot(egoMF, geneSetID = 1, by = "runningScore", title = egoMF$Description[1])
#p2 <- gseaplot(egoMF, geneSetID = 1, by = "preranked", title = egoMF$Description[1])
#p3 <- gseaplot(egoMF, geneSetID = 1, title = egoMF$Description[1])
#cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])
#ggsave("/home/simone/mnt/part1/mouse_rnaseq/analysis/LSK/ORA/GO/MF/gsea.top.plot.pdf")

#upsetplot(egoMF)
#ggsave("/home/simone/mnt/part1/mouse_rnaseq/analysis/LSK/ORA/GO/MF/upset.gsea.pdf")

#BP
cmd <- paste("mkdir -p", file.path(GOoutdir,"BP"))
system(cmd)


gseaBP <- gseGO(geneList     = genelist,
               OrgDb        = org.Hs.eg.db,
               ont          = "BP",
               minGSSize    = 2,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               verbose      = TRUE,
               by           = "DOSE",
               nPerm        = 10000,
               seed         = TRUE)

#gseaBP <- setReadable(gseaBP, 'org.Hs.eg.db', 'ENTREZID')

#gsea_tableBP <- gseaBP@result
#write.table(gsea_tableBP,file.path(GOoutdir, "BP","gsea.table.txt"), sep = "\t", quote = FALSE)
#write.csv(gsea_tableBP,file.path(GOoutdir, "BP","gsea.table.csv"), quote = FALSE)
#write.xlsx(gsea_tableBP,file.path(GOoutdir, "BP","gsea.table.xlsx"))


#ridgeplot(gseaBP,core_enrichment = TRUE)
#ggsave(file.path(GOoutdir, "BP","ridgeplot.pdf"), height = 15, width = 12)

#MF
cmd <- paste("mkdir -p", file.path(GOoutdir,"MF"))
system(cmd)


egoMF <- gseGO(geneList     = genelist,
               OrgDb        = org.Hs.eg.db,
               ont          = "MF",
               minGSSize    = 2,
               maxGSSize    = 500,
               pvalueCutoff = 0.1,
               verbose      = FALSE,
               by           = "DOSE",
               nPerm        = 10000, 
               seed         = TRUE)


gsea_tableMF <- egoMF@result
write.table(gsea_tableMF,file.path(GOoutdir, "MF","gsea.table.txt"), sep = "\t", quote = FALSE)
write.csv(gsea_tableMF, file.path(GOoutdir, "MF","gsea.table.csv"), quote = FALSE)
write.xlsx(gsea_tableMF, file.path(GOoutdir, "MF","gsea.table.xlsx"))

p <- ridgeplot(egoMF,core_enrichment = TRUE)
ggsave(file.path(GOoutdir, "MF","ridgeplot.pdf"))

p1 <- gseaplot(egoMF, geneSetID = 1, by = "runningScore", title = egoMF$Description[1])
p2 <- gseaplot(egoMF, geneSetID = 1, by = "preranked", title = egoMF$Description[1])
p3 <- gseaplot(egoMF, geneSetID = 1, title = egoMF$Description[1])
p <- cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])
ggsave("/home/simone/mnt/part1/mouse_rnaseq/analysis/LSK/ORA/GO/MF/gsea.top.plot.pdf")

#upsetplot(egoMF)
#ggsave("/home/simone/mnt/part1/mouse_rnaseq/analysis/LSK/ORA/GO/MF/upset.gsea.pdf"

#GSEA REACTOME
reactomepath <- paste(file.path(outputdir, "REACTOME"))
cmd <- paste("mkdir -p", reactomepath)
system(cmd)

y <- gsePathway(genelist, 
                pvalueCutoff  = 0.1,
                pAdjustMethod = "BH", 
                verbose       = TRUE,
                minGSSize     = 2,
                maxGSSize     = 500,
                by            = "DOSE",
                nPerm         = 10000,
                seed          = TRUE)

#y <- setReadable(y, 'org.Hs.eg.db', 'ENTREZID')

#gsea_tableReactome <- y@result
#write.table(gsea_tableReactome,file.path(reactomepath,"gsea.table.txt"), sep = "\t", quote = FALSE)
#write.csv(gsea_tableReactome, file.path(reactomepath, "gsea.table.csv"), quote = FALSE)
#write.xlsx(gsea_tableReactome, file.path(reactomepath,"gsea.table.xlsx"))

#ridgeplot(y,core_enrichment = TRUE)
#ggsave(file.path(reactomepath,"ridgeplot.pdf"))

#GSEA KEGG
KEGGoutdir <- paste(file.path(outputdir, "KEGG"))
cmd <- paste("mkdir -p", KEGGoutdir)
system(cmd)

kk2 <- gseKEGG(geneList      = genelist,
               organism      = 'hsa',
               minGSSize     = 2,
               maxGSSize     = 500,
               pvalueCutoff  = 0.1,
               verbose       = TRUE,
               nPerm         = 10000,
               by            ="DOSE",
               seed          = TRUE)
               

#kk2 <- setReadable(kk2, 'org.Hs.eg.db', 'ENTREZID')

#gsea_tableKEGG <- kk2@result
#write.table(gsea_tableKEGG,file.path(KEGGoutdir,"gsea.table.txt"), sep = "\t", quote = FALSE)
#write.csv(gsea_tableKEGG, file.path(KEGGoutdir, "gsea.table.csv"), quote = FALSE)
#write.xlsx(gsea_tableKEGG, file.path(KEGGoutdir,"gsea.table.xlsx"))

#ridgeplot(kk2,core_enrichment = TRUE)
#ggsave(file.path(reactomepath,"ridgeplot.pdf"))

#GSEA Wiki
wikipath <- paste(file.path(outputdir, "Wiki"))
cmd <- paste("mkdir -p", wikipath)
system(cmd)

wikigsea <- gseWP(genelist, 
          organism      = "Homo sapiens",
          pvalueCutoff  = 0.1,
          minGSSize     = 2,
          maxGSSize     = 500,
          by            = "DOSE",
          nPerm         = 10000,
          seed          = TRUE,
          verbose       = TRUE)       

#wikigsea <- setReadable(wikigsea, 'org.Hs.eg.db', 'ENTREZID')
#gsea_tableWiki <- wikigsea@result
#write.table(gsea_tableWiki,file.path(wikipath,"gsea.table.txt"), sep = "\t", quote = FALSE)
#write.csv(gsea_tableWiki, file.path(wikipath, "gsea.table.csv"), quote = FALSE)
#write.xlsx(gsea_tableWiki, file.path(wikipath,"gsea.table.xlsx"))

#ridgeplot(wikigsea,core_enrichment = TRUE)
#ggsave(file.path(wikipath,"ridgeplot.pdf"))

#GSEA DO
dopath <- paste(file.path(outputdir, "DO"))
cmd <- paste("mkdir -p", dopath)
system(cmd)

y <- gseDO(genelist,
           minGSSize     = 2,
           maxGSSize     = 500,
           pvalueCutoff  = 0.1,
           pAdjustMethod = "BH",
           verbose       = TRUE,
           by            = "DOSE",
           nPerm         = 10000,
           seed          = TRUE)

#DOgsea <- setReadable(y, 'org.Hs.eg.db', 'ENTREZID')
#gsea_tableDO <- DOgsea@result
#write.table(gsea_tableDO,file.path(dopath,"gsea.table.txt"), sep = "\t", quote = FALSE)
#write.csv(gsea_tableDO, file.path(dopath, "gsea.table.csv"), quote = FALSE)
#write.xlsx(gsea_tableDO, file.path(dopath,"gsea.table.xlsx"))

#ridgeplot(DOgsea,core_enrichment = TRUE)
#ggsave(file.path(dopath,"ridgeplot.pdf"))

#GSEA HALLMARK of CANCER
Msigdbrpath <- paste(file.path(outputdir, "msigdbr"))
cmd <- paste("mkdir -p", Msigdbrpath)
system(cmd)

Hpath <- file.path(Msigdbrpath,"H")

cmd <- paste("mkdir -p", file.path(Msigdbrpath,"H"))
system(cmd)

m_t2g_H <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)


y <- GSEA(genelist, 
          TERM2GENE     = m_t2g_H,
          minGSSize     = 2,
          maxGSSize     = 500,
          pvalueCutoff  = 0.1,
          pAdjustMethod = "BH",
          verbose       = TRUE,
          by            = "DOSE",
          nPerm         = 10000,
          seed          = TRUE)

Hgsea <- setReadable(y, 'org.Hs.eg.db', 'ENTREZID')
gsea_tableH <- Hgsea@result
write.table(gsea_tableH,file.path(Hpath,"gsea.table.txt"), sep = "\t", quote = FALSE)
write.csv(gsea_tableH, file.path(Hpath, "gsea.table.csv"), quote = FALSE)
write.xlsx(gsea_tableH, file.path(Hpath,"gsea.table.xlsx"))

p <- ridgeplot(Hgsea,core_enrichment = TRUE)
ggsave(file.path(Hpath,"ridgeplot.pdf"))


#GSEA NCG
NCGpath <- paste(file.path(outputdir, "NCG"))
cmd <- paste("mkdir -p", NCGpath)
system(cmd)

ncg <- gseNCG(genelist,
              pvalueCutoff  = 0.1,
              minGSSize     = 2,
              maxGSSize     = 500,
              pAdjustMethod = "BH",
              verbose       = TRUE,
              by            = "DOSE",
              nPerm         = 10000,
              seed          = TRUE)

ncg <- setReadable(ncg, 'org.Hs.eg.db')

gsea_tablencg <- ncg@result
write.table(gsea_tablencg,file.path(NCGpath,"gsea.table.txt"), sep = "\t", quote = FALSE)
write.csv(gsea_tablencg, file.path(NCGpath, "gsea.table.csv"), quote = FALSE)
write.xlsx(gsea_tablencg, file.path(NCGpath,"gsea.table.xlsx"))

p <- ridgeplot(ncg,core_enrichment = TRUE)
ggsave(file.path(NCGpath,"ridgeplot.pdf"))


#GSEA DGN
dgnpath <- paste(file.path(outputdir, "DGN"))
cmd <- paste("mkdir -p", dgnpath)
system(cmd)

dgn <- gseDGN(genelist,
              pvalueCutoff  = 0.1,
              minGSSize     = 2,
              maxGSSize     = 500,
              pAdjustMethod = "BH",
              verbose       = TRUE,
              by            = "DOSE",
              nPerm         = 10000,
              seed          = TRUE)

dgn <- setReadable(dgn, 'org.Hs.eg.db')

#gsea_tabledgn <- dgn@result
#write.table(gsea_tabledgn,file.path(dgnpath,"gsea.table.txt"), sep = "\t", quote = FALSE)
#write.csv(gsea_tabledgn, file.path(dgnpath, "gsea.table.csv"), quote = FALSE)
#write.xlsx(gsea_tabledgn, file.path(dgnpath,"gsea.table.xlsx"))

#ridgeplot(dgn,core_enrichment = TRUE)
#ggsave(file.path(dgnpath,"ridgeplot.pdf"))

#GSEA C2
C2path <- file.path(Msigdbrpath,"C2")

cmd <- paste("mkdir -p", file.path(Msigdbrpath,"C2"))
system(cmd)

m_t2g_C2 <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)

y <- GSEA(genelist, 
          TERM2GENE     = m_t2g_C2,
          minGSSize     = 2,
          maxGSSize     = 500,
          pvalueCutoff  = 0.1,
          pAdjustMethod = "BH",
          verbose       = TRUE,
          by            = "DOSE",
          nPerm         = 10000,
          seed          = TRUE)

C2gsea <- setReadable(y, 'org.Hs.eg.db', 'ENTREZID')
gsea_tableC2 <- C2gsea@result
write.table(gsea_tableC2,file.path(C2path,"gsea.table.txt"), sep = "\t", quote = FALSE)
write.csv(gsea_tableC2, file.path(C2path, "gsea.table.csv"), quote = FALSE)
write.xlsx(gsea_tableC2, file.path(C2path,"gsea.table.xlsx"))

ridgeplot(C2gsea,core_enrichment = TRUE)
ggsave(file.path(C2path,"ridgeplot.pdf"))


#GSEA Biocarta
Biocartapath <- file.path(Msigdbrpath,"BIOCARTA")

cmd <- paste("mkdir -p", file.path(Msigdbrpath,"BIOCARTA"))
system(cmd)

m_t2g_biocarta <- msigdbr(species = "Homo sapiens", category = "C2",subcategory="CP:BIOCARTA") %>% 
  dplyr::select(gs_name, entrez_gene)

y <- GSEA(genelist, 
          TERM2GENE     = m_t2g_biocarta,
          minGSSize     = 2,
          maxGSSize     = 500,
          pvalueCutoff  = 0.1,
          pAdjustMethod = "BH",
          verbose       = TRUE,
          by            = "DOSE",
          nPerm         = 10000,
          eps           = 10,
          exponent      = 0.01,
          seed          = TRUE)

#Biocartagsea <- setReadable(y, 'org.Hs.eg.db', 'ENTREZID')
#gsea_tableBiocarta <- Biocartagsea@result
#write.table(gsea_tableBiocarta,file.path(Biocartapath,"gsea.table.txt"), sep = "\t", quote = FALSE)
#write.csv(gsea_tableBiocarta, file.path(Biocartapath, "gsea.table.csv"), quote = FALSE)
#write.xlsx(gsea_tableBiocarta, file.path(Biocartapath,"gsea.table.xlsx"))

#ridgeplot(Biocartagsea,core_enrichment = TRUE)
#ggsave(file.path(C2path,"ridgeplot.pdf"))

#GSEA PID
PIDpath <- file.path(Msigdbrpath,"PID")

cmd <- paste("mkdir -p", file.path(Msigdbrpath,"PID"))
system(cmd)

m_t2g_pid <- msigdbr(species = "Homo sapiens", category = "C2",subcategory="CP:PID") %>% 
  dplyr::select(gs_name, entrez_gene)

y <- GSEA(genelist, 
          TERM2GENE     = m_t2g_pid,
          minGSSize     = 2,
          maxGSSize     = 500,
          pvalueCutoff  = 0.1,
          pAdjustMethod = "BH",
          verbose       = TRUE,
          by            = "DOSE",
          nPerm         = 10000,
          eps           = 10,
          exponent      = 0.01,
          seed          = TRUE)

#PIDgsea <- setReadable(y, 'org.Hs.eg.db', 'ENTREZID')
#gsea_tablePID <- PIDgsea@result
#write.table(gsea_tablePID,file.path(PIDpath,"gsea.table.txt"), sep = "\t", quote = FALSE)
#write.csv(gsea_tablePID, file.path(PIDpath, "gsea.table.csv"), quote = FALSE)
#write.xlsx(gsea_tablePID, file.path(PIDpath,"gsea.table.xlsx"))

#ridgeplot(PIDgsea,core_enrichment = TRUE)
#ggsave(file.path(PIDpath,"ridgeplot.pdf"))

#GSEA C6
C6path <- file.path(Msigdbrpath,"C6")

cmd <- paste("mkdir -p", file.path(Msigdbrpath,"C6"))
system(cmd)

m_t2g_C6 <- msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::select(gs_name, entrez_gene)

y <- GSEA(genelist, 
          TERM2GENE     = m_t2g_C6,
          minGSSize     = 2,
          maxGSSize     = 500,
          pvalueCutoff  = 0.1,
          pAdjustMethod = "BH",
          verbose       = TRUE,
          by            = "DOSE",
          nPerm         = 10000,
          eps           = 10,
          exponent      = 0.01,
          seed          = TRUE)

C6gsea <- setReadable(y, 'org.Hs.eg.db', 'ENTREZID')
gsea_tableC6 <- C6gsea@result
write.table(gsea_tableC6,file.path(C6path,"gsea.table.txt"), sep = "\t", quote = FALSE)
write.csv(gsea_tableC6, file.path(C6path, "gsea.table.csv"), quote = FALSE)
write.xlsx(gsea_tableC6, file.path(C6path,"gsea.table.xlsx"))

ridgeplot(C6gsea,core_enrichment = TRUE)
#ggsave(file.path(C6path,"ridgeplot.pdf"))