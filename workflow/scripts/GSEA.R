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

outputdir <- "/projects2/2024_Chiarugi/Transcriptomestic/results/unstranded/STAR_HTSeq_nosample1/GSEA_all"
cmd <- paste("mkdir -p", outputdir)
system(cmd)

###### GSEA
df <- read.csv("/projects2/2024_Chiarugi/Transcriptomestic/results/unstranded/STAR_HTSeq_nosample1/results.csv")

genelist <- df$log2FoldChange
names(genelist) <- df$entrez

genelist <- genelist[order(genelist, decreasing=TRUE)]
genelist <- genelist[-which(is.na(genelist))]

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
               pvalueCutoff = 0.1,
               verbose      = TRUE,
               by           = "DOSE",
               eps           = 0.0001,
               nPerm         = 1000,
               exponent      = 0.01, 
               seed         = TRUE)

gseaCC <- setReadable(gseaCC, 'org.Hs.eg.db', 'ENTREZID')

gsea_tableCC <- gseaCC@result
write.table(gsea_tableCC,file.path(GOoutdir, "CC","gsea.table.txt"), sep = "\t", quote = FALSE)
write.csv(gsea_tableCC,file.path(GOoutdir, "CC","gsea.table.csv"), quote = FALSE)
#write.xlsx(gsea_tableCC,file.path(GOoutdir, "CC","gsea.table.xlsx"))

p <- ridgeplot(gseaCC,core_enrichment = TRUE, showCategory=10)
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
               pvalueCutoff = 0.1,
               verbose      = TRUE,
               by           = "DOSE",
               eps           = 0.01,
               nPerm         = 1000,
               exponent      = 0.01, 
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
               verbose      = TRUE,
               by           = "DOSE",
               eps           = 0.0001,
               nPerm         = 1000,
               exponent      = 0.01, 
               seed         = TRUE)


gsea_tableMF <- egoMF@result
write.table(gsea_tableMF,file.path(GOoutdir, "MF","gsea.table.txt"), sep = "\t", quote = FALSE)
write.csv(gsea_tableMF, file.path(GOoutdir, "MF","gsea.table.csv"), quote = FALSE)
#write.xlsx(gsea_tableMF, file.path(GOoutdir, "MF","gsea.table.xlsx"))

p <- ridgeplot(egoMF,core_enrichment = TRUE)
ggsave(file.path(GOoutdir, "MF","ridgeplot.pdf"))

p1 <- gseaplot(egoMF, geneSetID = 1, by = "runningScore", title = egoMF$Description[1])
p2 <- gseaplot(egoMF, geneSetID = 1, by = "preranked", title = egoMF$Description[1])
p3 <- gseaplot(egoMF, geneSetID = 1, title = egoMF$Description[1])
p4 <- gseaplot2(egoMF, geneSetID = 1, title = egoMF$Description[1])
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
                eps           = 0.0001,
                nPerm         = 1000,
                exponent      = 0.01,
                seed          = TRUE)

y <- setReadable(y, 'org.Hs.eg.db', 'ENTREZID')

save(y,file=file.path(reactomepath,"GSEA.rda"))

gsea_tableReactome <- y@result
write.table(gsea_tableReactome,file.path(reactomepath,"gsea.table.txt"), sep = "\t", quote = FALSE)
write.csv(gsea_tableReactome, file.path(reactomepath, "gsea.table.csv"), quote = FALSE)
#write.xlsx(gsea_tableReactome, file.path(reactomepath,"gsea.table.xlsx"))

p <- ridgeplot(y,core_enrichment = TRUE, showCategory=10)
ggsave(file.path(reactomepath,"ridgeplot.pdf"))

#RHO GTPase Effectors
selected_ <- "RHO GTPase Effectors"
p1 <- gseaplot(y, geneSetID = 1, by = "runningScore", title = selected_)
p2 <- gseaplot(y, geneSetID = 1, by = "preranked", title = selected_)
p3 <- gseaplot(y, geneSetID = 1, title = selected_)
p4 <- gseaplot2(y, geneSetID = 1, title = selected_)
#p <- cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])
pdf(file.path(reactomepath,"gsea.top.RHO.pdf"))
gseaplot(y, geneSetID =  9, title = y$Description[9])
dev.off()

#Transcriptional Regulation by TP53
selected_ <- "Transcriptional Regulation by TP53"
p4 <- gseaplot2(y, geneSetID = 1, title = selected_)
#p <- cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])
pdf(file.path(reactomepath,"gsea.top.TP53.pdf"))
gseaplot(y, geneSetID = 18, title = y$Description[18])
dev.off()

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
               eps           = 0.0001,
               nPerm         = 1000,
               exponent      = 0.01,
               seed          = TRUE,
               by            = "DOSE")
               

kk2 <- setReadable(kk2, 'org.Hs.eg.db', 'ENTREZID')

save(kk2,file=file.path(KEGGoutdir,"GSEA.rda"))

gsea_tableKEGG <- kk2@result
write.table(gsea_tableKEGG,file.path(KEGGoutdir,"gsea.table.txt"), sep = "\t", quote = FALSE)
write.csv(gsea_tableKEGG, file.path(KEGGoutdir, "gsea.table.csv"), quote = FALSE)
write.xlsx(gsea_tableKEGG, file.path(KEGGoutdir,"gsea.table.xlsx"))

ridgeplot(kk2,core_enrichment = TRUE, showCategory=10)
ggsave(file.path(KEGGoutdir,"ridgeplot.pdf"))

#Cell cycle
selected_ <- "Cell cycle"
p4 <- gseaplot2(kk2, geneSetID = 2, title = selected_)
#p <- cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])
pdf(file.path(KEGGoutdir,"gsea.top.cellcycle.pdf"))
gseaplot(y, geneSetID = 2, title = selected_)
dev.off()

#Apoptosis
selected_ <- "Apoptosis"
p4 <- gseaplot2(kk2, geneSetID = 19, title = selected_)
#p <- cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])
pdf(file.path(KEGGoutdir,"gsea.top.apoptosis.pdf"))
gseaplot(y, geneSetID = 19, title = selected_)
dev.off()

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
          nPerm         = 1000,
          eps           = 0.01,
          exponent      = 0.01,
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
           nPerm         = 1000,
           seed          = TRUE,
           eps           = 0.01,
           exponent      = 0.01,
)

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
          nPerm         = 1000,
          seed          = TRUE,
          eps           = 0.01,
          exponent      = 0.0001)

Hgsea <- setReadable(y, 'org.Hs.eg.db', 'ENTREZID')
save(Hgsea,file=file.path(Hpath,"GSEA.rda"))

gsea_tableH <- Hgsea@result
write.table(gsea_tableH,file.path(Hpath,"gsea.table.txt"), sep = "\t", quote = FALSE)
write.csv(gsea_tableH, file.path(Hpath, "gsea.table.csv"), quote = FALSE)
write.xlsx(gsea_tableH, file.path(Hpath,"gsea.table.xlsx"))

p <- ridgeplot(Hgsea,core_enrichment = TRUE, showCategory=10)
ggsave(file.path(Hpath,"ridgeplot.pdf"))


pdf(file.path(Hpath,"gsea.top.p53.pdf"))
gseaplot(Hgsea, geneSetID = 5, title = Hgsea$Description[5])
dev.off()

#HALLMARK_E2F_TARGETS
selected_ <- "HALLMARK_E2F_TARGETS"
#p4 <- gseaplot2(Hgsea, geneSetID = 1, title = Hgsea$Description[1])
#p <- cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])
pdf(file.path(Hpath,"gsea.top.E2F.pdf"))
gseaplot(Hgsea, geneSetID = 1, title = Hgsea$Description[1])
dev.off()

#HALLMARK_TNFA_SIGNALING_VIA_NFKB
selected_ <- "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
#p <- cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])
pdf(file.path(Hpath,"gsea.top.TNFA.pdf"))
gseaplot(Hgsea, geneSetID = 2, title = Hgsea$Description[2])
dev.off()

#HALLMARK_G2M_CHECKPOINT
selected_ <- "HALLMARK_G2M_CHECKPOINT"
#p <- cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])
pdf(file.path(Hpath,"gsea.top.G2M.pdf"))
gseaplot(Hgsea, geneSetID = 3, title = Hgsea$Description[3])
dev.off()

#HALLMARK_MYC_TARGETS_V1
selected_ <- "HALLMARK_MYC_TARGETS_V1"
#p <- cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])
pdf(file.path(Hpath,"gsea.top.MYC.pdf"))
gseaplot(Hgsea, geneSetID = 4, title = Hgsea$Description[4])
dev.off()

#HALLMARK_KRAS_SIGNALING_DN
selected_ <- "HALLMARK_KRAS_SIGNALING_DN"
#p <- cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])
pdf(file.path(Hpath,"gsea.top.KRAS.pdf"))
gseaplot(Hgsea, geneSetID = 9, title = Hgsea$Description[9])
dev.off()

#HALLMARK_PI3K_AKT_MTOR_SIGNALING
selected_ <- "HALLMARK_PI3K_AKT_MTOR_SIGNALING"
#p <- cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])
pdf(file.path(Hpath,"gsea.top.PI3K.pdf"))
gseaplot(Hgsea, geneSetID = 25, title = Hgsea$Description[25])
dev.off()


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
              eps           = 0.01,
              exponent      = 0.0001,
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
              eps           = 0.01,
              exponent      = 0.0001,
              nPerm         = 1000,
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
            pvalueCutoff  = 0.05,
            pAdjustMethod = "BH",
            verbose       = TRUE,
            by            = "DOSE",
            nPerm         = 1000,
            eps           = 0.00000001,
            exponent      = 0.1,
            seed          = TRUE)


C2gsea <- setReadable(y, 'org.Hs.eg.db', 'ENTREZID')
gsea_tableC2 <- C2gsea@result
write.table(gsea_tableC2,file.path(C2path,"gsea.table.txt"), sep = "\t", quote = FALSE)
write.csv(gsea_tableC2, file.path(C2path, "gsea.table.csv"), quote = FALSE)
#write.xlsx(gsea_tableC2, file.path(C2path,"gsea.table.xlsx"))

p <- ridgeplot(C2gsea,core_enrichment = TRUE, showCategory=10)
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
          pvalueCutoff  = 0.05,
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
          pvalueCutoff  = 0.05,
          pAdjustMethod = "BH",
          verbose       = TRUE,
          by            = "DOSE",
          seed          = TRUE,
          nPerm         = 1000,
          eps           = 0.00000001,
          exponent      = 0.1
)

C6gsea <- setReadable(y, 'org.Hs.eg.db', 'ENTREZID')

save(C6gsea, file=file.path(C6path,"GSEA.rda"))

gsea_tableC6 <- C6gsea@result
write.table(gsea_tableC6,file.path(C6path,"gsea.table.txt"), sep = "\t", quote = FALSE)
write.csv(gsea_tableC6, file.path(C6path, "gsea.table.csv"), quote = FALSE)
#write.xlsx(gsea_tableC6, file.path(C6path,"gsea.table.xlsx"))

p <- ridgeplot(C6gsea,core_enrichment = TRUE, showCategory=10)
ggsave(file.path(C6path,"ridgeplot.pdf"))

#P53_DN.V1_UP
selected_ <- "P53_DN.V1_UP"
pdf(file.path(C6path,"gsea.top.p53.pdf"))
gseaplot(C6gsea, geneSetID = 17, title = selected_)
dev.off()

#KRAS.600_UP.V1_DN
selected_ <- "KRAS.600_UP.V1_DN"
#p <- cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])
pdf(file.path(C6path,"gsea.top.kras.pdf"))
gseaplot(C6gsea, geneSetID = 15, title = selected_)
dev.off()

#AKT_UP.V1_DN
selected_ <- "AKT_UP.V1_DN"
#p <- cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])
pdf(file.path(C6path,"gsea.top.akt.pdf"))
gseaplot(C6gsea, geneSetID = 18, title = C6gsea$Description[18])
dev.off()

#GSEA C7
C7path <- file.path(Msigdbrpath,"C7")

cmd <- paste("mkdir -p", file.path(Msigdbrpath,"C7"))
system(cmd)

m_t2g_C7 <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_name, entrez_gene)

y <- GSEA(genelist, 
          TERM2GENE     = m_t2g_C7,
          minGSSize     = 2,
          maxGSSize     = 500,
          pvalueCutoff  = 0.05,
          pAdjustMethod = "BH",
          verbose       = TRUE,
          by            = "DOSE",
          seed          = TRUE,
          nPerm         = 1000,
          eps           = 0.00001,
          exponent      = 0.1
)

C7gsea <- setReadable(y, 'org.Hs.eg.db', 'ENTREZID')
gsea_tableC7 <- C7gsea@result
write.table(gsea_tableC7,file.path(C7path,"gsea.table.txt"), sep = "\t", quote = FALSE)
write.csv(gsea_tableC7, file.path(C7path, "gsea.table.csv"), quote = FALSE)
#write.xlsx(gsea_tableC6, file.path(C6path,"gsea.table.xlsx"))

p <- ridgeplot(C6gsea,core_enrichment = TRUE, showCategory=10)
ggsave(file.path(C6path,"ridgeplot.pdf"))