#### ###  #### ###
##LungTx Compath##
##ALL stat analyses and figures##
#### ###  #### ### 20191201

#libraries 
############# 
library("survival")
library("survminer")
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)
library("ggsci")
library(dplyr)
library(ComplexHeatmap)
library(ggcorrplot)
library(corrplot)
library(gridExtra)
library(RColorBrewer)
library(reshape2)
library(ggstatsplot)
library(circlize)
#############

#Functions
#############
colNA <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
ggcorrplotModified <- function (corr, method = c("square", "circle"), type = c("full",
                                                                               "lower", "upper"), ggtheme = ggplot2::theme_minimal, title = "",
                                show.legend = TRUE, legend.title = "Corr", show.diag = FALSE,
                                colors = c("blue", "white", "red"), outline.color = "gray",
                                hc.order = FALSE, hc.method = "complete", lab = FALSE, lab_col = "black",
                                lab_size = 4, p.mat = NULL, sig.level = 0.05, insig = c("pch",
                                                                                        "blank"), pch = 4, pch.col = "black", pch.cex = 5, tl.cex = 12,
                                tl.col = "black", tl.srt = 45, lab.notsig="")
{
  type <- match.arg(type)
  method <- match.arg(method)
  insig <- match.arg(insig)
  if (!is.matrix(corr) & !is.data.frame(corr))
    stop("Need a matrix or data frame!")
  corr <- as.matrix(corr)
  if (hc.order) {
    ord <- ggcorrplot:::.hc_cormat_order(corr)
    corr <- corr[ord, ord]
    if (!is.null(p.mat))
      p.mat <- p.mat[ord, ord]
  }
  if (type == "lower") {
    corr <- ggcorrplot:::.get_lower_tri(corr, show.diag)
    p.mat <- ggcorrplot:::.get_lower_tri(p.mat, show.diag)
  }
  else if (type == "upper") {
    corr <- ggcorrplot:::.get_upper_tri(corr, show.diag)
    p.mat <- ggcorrplot:::.get_upper_tri(p.mat, show.diag)
  }
  corr <- reshape2::melt(corr, na.rm = TRUE)
  corr$pvalue <- rep(NA, nrow(corr))
  corr$signif <- rep(NA, nrow(corr))
  if (!is.null(p.mat)) {
    p.mat <- reshape2::melt(p.mat, na.rm = TRUE)
    corr$coef <- corr$value
    corr$pvalue <- p.mat$value
    corr$signif <- as.numeric(p.mat$value <= sig.level)
    p.mat <- subset(p.mat, p.mat$value > sig.level)
    if (insig == "blank")
      corr$value <- corr$value * corr$signif
  }
  corr$abs_corr <- abs(corr$value) * 10
  p <- ggplot2::ggplot(corr, ggplot2::aes_string("Var1", "Var2",
                                                 fill = "value"))
  if (method == "square")
    p <- p + ggplot2::geom_tile(color = outline.color)
  else if (method == "circle") {
    p <- p + ggplot2::geom_point(color = outline.color, shape = 21,
                                 ggplot2::aes_string(size = "abs_corr")) + ggplot2::scale_size(range = c(4,
                                                                                                         10)) + ggplot2::guides(size = FALSE)
  }
  p <- p + ggplot2::scale_fill_gradient2(low = colors[1], high = colors[3],
                                         mid = colors[2], midpoint = 0, limit = c(-1, 1), space = "Lab",
                                         name = legend.title) + ggtheme() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = tl.srt,
                                                                                                                               vjust = 1, size = tl.cex, hjust = 1), axis.text.y = ggplot2::element_text(size = tl.cex)) +
    ggplot2::coord_fixed()
  label <- as.character(round(corr[, "value"], 2))
  label[label=="0"] <- lab.notsig
  if (lab)
    p <- p + ggplot2::geom_text(ggplot2::aes_string("Var1",
                                                    "Var2", label = "label"), color = lab_col, size = lab_size)
  if (!is.null(p.mat) & insig == "pch") {
    p <- p + ggplot2::geom_point(data = p.mat, ggplot2::aes_string("Var1",
                                                                   "Var2"), shape = pch, size = pch.cex, color = pch.col)
  }
  if (title != "")
    p <- p + ggplot2::ggtitle(title)
  if (!show.legend)
    p <- p + ggplot2::theme(legend.position = "none")
  p <- p + ggcorrplot:::.no_panel()
  p
}
#############

#Data
#############

#TRACERx and TCGA data
load("X:/Dropbox (ICR)/yuanlab/Manuscripts/Lung Tx1/revision1/190501/data/tx100_compath.RData")

#1
#regional: tx100 multiregion cohort (n=275 for 85 patients). 
#pathology_LymEstimate: pathology lym cell fraction scores for the discordant regions (n=42). 
#ASCAT and VAF: DNA-based tumor purity. cd8.score.danaher: Danaher et al RNA-seq signature of CD8+. orig_immune_cluster: Rosenthal et al RNA-seq based immune classification. 

#2
#diagnsotic: tx100 diagnostic cohort from both HE, IHC slides (n=100 for 100 patients). 
#Stromal.TIL: pathology scores on the first 100 slides. 

#3
#both regional and disgnostic scores combined together with key genomic scores, summarised in:
#tx: combined tx100 regional and diagnsotic cohort at patient level (n=85, patients with a regional sample).

#4 TCGA dataset
#tcga: contains NSCLC patient IDs, lym%, fd, immuneClass (n=933). 
#tcgaluad_immune: spatial TILs for LUAD patients with RNA immune signatures (n=464, after removing NAs).

#5
#reg_tils_immune: deconvulted immune cells onto the HE TIL classes for the same registered sections.
#reg_eval: manual evaluation of HE-IHC registered sections.  

#6
#rosenthal_immune_signatures.RData: RNA-seq Danaher et al immune signatures for 142 tracerx regions

#LATTICe-A data
load("X:/Dropbox (ICR)/yuanlab/Manuscripts/Lung Tx1/revision1/190501/data/latticea_compath.RData")
#latticea: manual estimations for fraction of lym and path. TILs, and corresponding automated scores (n=80 patients LATTICe-A cohort)
#latticea_bioai: list contains 3 DFs for TTF1, CD45, SMA. Rows correspond to image patches from cores and cols show the IHC/HE cell counts by deep learning

#latticea_970_finalCohort.RData: entire LATTICe-A cohort (n=970 patients), all patient level and regional data

#############

#### ###  #### ### #### ###  #### ###
#Figures
#### ###  #### ### #### ###  #### ###

#1e-f and S2 g: biological AI corrs 
#############
ggscatter(latticea_bioai$TTF1, x = "tumour_norm", y = "ttf1_norm", add = "reg.line", xlab = "H&E-based cancer cell / 100um", 
              ylab = "TTF1+ cell / 100um",
              color = "#00ff00",
              conf.int = TRUE, size = 2,
              add.params = list(color = "grey50", fill = "azure3"), 
              cor.coef = TRUE, cor.method = "spearman")

ggscatter(latticea_bioai$CD3, x = "lym_norm", y = "cd3_norm", add = "reg.line", xlab = "H&E-based lymphocyte cell / 100um", 
              ylab = "CD45+ cell / 100um",
              color = "#0000ff",
              conf.int = TRUE, size = 2,
              add.params = list(color = "grey50", fill = "azure3"), 
              cor.coef = TRUE, cor.method = "spearman") 

ggscatter(latticea_bioai$SMA, x = "stromal_norm", y = "sma_norm", add = "reg.line", xlab = "H&E-based stromal cell / 100um", 
              ylab = "SMA cell / 100um",
              color = "goldenrod",
              conf.int = TRUE, size = 2,
              add.params = list(color = "grey50", fill = "azure3"), 
              cor.coef = TRUE, cor.method = "spearman") 
#############

#2a
#############
r2=regional
r2$Histology[r2$Histology=="Other"] <- "Z"
r2 = r2[order(r2$lymphocytes_per, decreasing = F),]
r2 = r2[order(r2$Histology, decreasing = F),]
pDotplot = ggdotchart(r2, x = "PublicationID", y = "lymphocytes_per",
           color = "immuneClass_2",                                
           group = "PublicationID",
           sort.by.groups = TRUE,           
           palette = c("blue2", "red2", "#f4d5d5"), 
           #sorting = "descending",                       
           #rotate = TRUE,                               
           #add = "segments",                             
           add.params = list(color = "lightgray", size = 0.8),
           dot.size = 3,                                 
           #label = min(regionalT$lymphocytes_per),
           #y.text.col = TRUE,                            
           ggtheme = theme_pubr())
#############

#2c-d: Danaher immune signatures 
#############
load("X:/Dropbox (ICR)/yuanlab/Manuscripts/Lung Tx1/revision2/data/rosenthal_immune_signatures.RData")
r = tx_rosenthal
Heatmap(scale(r[, c(3)]), name = "comb", row_split = r$immuneClass_2, row_gap = unit(2, "mm"), border = F, 
        col = colorRamp2(c(-4.2, 0, 6), c("#4DBBD5B2", "white", "red2")), cluster_columns = FALSE,
        show_row_dend = F, show_column_dend = F, row_labels = F, show_row_names = F, cluster_row_slices = F, show_heatmap_legend = F,
        column_names_rot = 45)+
  Heatmap(scale(r[, c(5:18)]), name = "comb", row_split = r$immuneClass_2, row_gap = unit(2, "mm"), border = F, 
          col = colorRamp2(c(-4.2, 0, 6), c("#4DBBD5B2", "white", "red2")), column_title = "RNA-seq immune cell signatures", cluster_columns = FALSE,
          show_row_dend = F, show_column_dend = F, row_labels = F, show_row_names = F, cluster_row_slices = F, show_heatmap_legend = F,
          column_names_rot = 45)

#hot-cold comparisons: 
rC = r[! r$immuneClass_2 =="intermediate" ,]
rC$sample = paste0(rC$PublicationID, "_", rC$region)
rC = rC[order(rC$immuneClass_2) ,]
rC$Name = rep("tx_all")
rC = rC[, c(19:20, 4, 5:18)]
df2<-melt(rC,id.var=c("sample","immuneClass_2", "Name"))

p <- ggplot(df2, aes(variable, value,fill=immuneClass_2), color = c("blue2", "red2"))
p + geom_boxplot()+scale_fill_manual(values = c("blue2", "red2"))+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Immune signature") + ylab("Estimated infiltration (Danaher et al)") 

comList <- list( c("cold", "hot"))
p1=ggboxplot(rC, x = "immuneClass_2", 
             y = "cd8", color = "immuneClass_2", ylab = "CD8",
             add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())
p2=ggboxplot(rC, x = "immuneClass_2", 
             y = "cd4", color = "immuneClass_2", ylab = "CD4",
             add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())
p3=ggboxplot(rC, x = "immuneClass_2", 
             y = "bcell", color = "immuneClass_2", ylab = "B cell",
             add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())
p4=ggboxplot(rC, x = "immuneClass_2", 
             y = "cd45", color = "immuneClass_2", ylab = "CD45",
             add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())
p5=ggboxplot(rC, x = "immuneClass_2", 
             y = "cyto", color = "immuneClass_2", ylab = "Cyto",
             add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())
p6=ggboxplot(rC, x = "immuneClass_2", 
             y = "dend", color = "immuneClass_2", ylab = "Dendritic",
             add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())
p7=ggboxplot(rC, x = "immuneClass_2", 
             y = "cd8.exhausted", color = "immuneClass_2", ylab = "CD8 exhausted",
             add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())
p8=ggboxplot(rC, x = "immuneClass_2", 
             y = "macrophage", color = "immuneClass_2", ylab = "Macrophage",
             add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())
p9=ggboxplot(rC, x = "immuneClass_2", 
             y = "neutrophil", color = "immuneClass_2", ylab = "Neutrophil",
             add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())
p10=ggboxplot(rC, x = "immuneClass_2", 
              y = "nkcd56dim", color = "immuneClass_2", ylab = "NK CD56",
              add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())
p11=ggboxplot(rC, x = "immuneClass_2", 
              y = "nk", color = "immuneClass_2", ylab = "NK",
              add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())
p12=ggboxplot(rC, x = "immuneClass_2", 
              y = "tcells", color = "immuneClass_2", ylab = "T cells",
              add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())
p13=ggboxplot(rC, x = "immuneClass_2", 
              y = "th1", color = "immuneClass_2", ylab = "TH1",
              add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())
p14=ggboxplot(rC, x = "immuneClass_2", 
              y = "treg", color = "immuneClass_2", ylab = "Treg",
              add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())
p = c(5.8e-5, 1.5e-6, 9.4e-9, 3.4e-5, 0.00014, 0.018, 1.4e-5, 
      0.007, 0.0029, 6.6e-6, 4e-4, 2.9e-6, 3.3e-6, 9.1e-7)
p.adjust(p, method = "BH")
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, ncol=7)


#############

#2b-c and S5a: phenotype-based pairwise genomic dist
#############
#mutations summary from Jamal-Hanjani, NEJM 2017
#summarised data files saved seprately in csvs

#tx_genomicImmune_pairdist.csv
ggboxplot(tx_genomicImmune_pairdist[tx_genomicImmune_pairdist$Histology=="LUAD",], x = "immuneClass2", y = "immuneDist", 
              color = "immuneClass2", palette = c("blue", "red2"),
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))
ggboxplot(tx_genomicImmune_pairdist[tx_genomicImmune_pairdist$Histology=="LUAD",], x = "immuneClass2", y = "geneticDist", 
              color = "immuneClass2", palette = c("blue", "red2"), 
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))

ggboxplot(tx_genomicImmune_pairdist[tx_genomicImmune_pairdist$Histology=="LUSC",], x = "immuneClass2", y = "immuneDist", 
              color = "immuneClass2", palette = c("blue", "red2"), 
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))
ggboxplot(tx_genomicImmune_pairdist[tx_genomicImmune_pairdist$Histology=="LUSC",], x = "immuneClass2", y = "geneticDist", 
              color = "immuneClass2", palette = c("blue", "red2"), 
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))


#distance using phy. tree for pairs in LUAD tumors

#txluad_domclone_pairdist.csv
comList <- list( c("hot_hot", "cold_cold"))
ggboxplot(phy, x = "immuneClass2", y = "phyDist2", color = "immuneClass2", palette = c("red2", "blue"), ylab = "Dist. of dominant clones to the MRCA (phylogenetics)",
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())
ggboxplot(phy, x = "immuneClass2", y = "geneticDist", color = "immuneClass2", palette = c("red2", "blue"), ylab = "Genomic distance",
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())
ggboxplot(phy, x = "immuneClass2", y = "immuneDist", color = "immuneClass2", palette = c("red2", "blue"), ylab = "Immune distance",
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())

#############

#3b-d, S6: FD
#############

#FD vs immune class in tx and tcga
comList = list(c("cold", "hot"))
ggboxplot(regional[regional$Histology=="LUAD" & ! regional$immuneClass_2=="intermediate",], 
          x = "immuneClass_2", y = "fd", color = "immuneClass_2", ylab = "Fractal dimension",
          add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))
ggboxplot(regional[regional$Histology=="LUSC"& ! regional$immuneClass_2=="intermediate",], 
          x = "immuneClass_2", y = "fd", color = "immuneClass_2", ylab = "Fractal dimension",
          add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))
ggboxplot(tcga[tcga$disease_code=="LUAD" & ! tcga$immuneClass_2=="intermediate",], x = "immuneClass_2", y = "fd", 
          color = "immuneClass_2",ylab = "Fractal dimension",
          add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))
ggboxplot(tcga[tcga$disease_code=="LUSC" & ! tcga$immuneClass_2=="intermediate",], x = "immuneClass_2", y = "fd", 
          color = "immuneClass_2",ylab = "Fractal dimension",
          add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))

#now FD for all vs immmune class
comList = list(c("cold", "hot"))
ggboxplot(regional[! regional$immuneClass_2=="intermediate",], 
          x = "immuneClass_2", y = "fd", color = "immuneClass_2", ylab = "Fractal dimension",
          add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))

#now test stromal% instead of FD in both classes
ggboxplot(regional[! regional$immuneClass_2=="intermediate",], 
          x = "immuneClass_2", y = "fibroblasts_per", color = "immuneClass_2", ylab = "Stromal%",
          add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))
ggboxplot(regional[regional$Histology=="LUAD" & ! regional$immuneClass_2=="intermediate",], 
          x = "immuneClass_2", y = "fibroblasts_per", color = "immuneClass_2", ylab = "Stromal%",
          add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))
ggboxplot(regional[regional$Histology=="LUSC"& ! regional$immuneClass_2=="intermediate",], 
          x = "immuneClass_2", y = "fibroblasts_per", color = "immuneClass_2", ylab = "Stromal%",
          add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))


#FD vs LOHHLA (these are just the boxplots, p-vals are computed using multivariate model then corrected, see below)
comList <- list( c("Loss", "No loss"))
ggboxplot(tx[tx$Histology=="LUSC",], x = "any.HLA.loss", y = "fd_max", 
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+ 
  theme(text = element_text(size=18))
#no sign. in LUSC

#LUAD
ggboxplot(tx[tx$Histology=="LUAD",], x = "any.HLA.loss", y = "fd_max", 
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+ 
  theme(text = element_text(size=18))
ggboxplot(tx[tx$Histology=="LUAD",], x = "any.HLA.A.loss", y = "fd_max", 
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+ 
  theme(text = element_text(size=18))
ggboxplot(tx[tx$Histology=="LUAD",], x = "any.HLA.B.loss", y = "fd_max", 
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+ 
  theme(text = element_text(size=18))
ggboxplot(tx[tx$Histology=="LUAD",], x = "any.HLA.C.loss", y = "fd_max",
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+ 
  theme(text = element_text(size=18))

#check independence to ClonalNeo
fit <- glm(fd_max ~ factor(any.HLA.loss) + ClonalNeo  , data=tx[tx$Histology=="LUAD",])
summary(fit)
#repeat the above for any.HLA.A.loss, any.HLA.B.loss, any.HLA.C.loss,
p=c(0.00249, 0.00948, 0.0968, 0.0331)
p.adjust(p, method = "BH")

#3d: evasion mechanisims in LUAD
txH2LUAD = tx[tx$Histology=="LUAD",]
txH2LUAD$FDorHLALOH <- NULL
txH2LUAD$FDorHLALOH[txH2LUAD$FD == "High" | txH2LUAD$any.HLA.loss == "Loss"] <- "Yes"
txH2LUAD$FDorHLALOH[is.na(txH2LUAD$FDorHLALOH)] <- "No"

#p-value using a Fisher test
fisher.test(table(txH2LUAD$FDorHLALOH, txH2LUAD$coldP3))

#plot 3d
txH2LUAD$FDorHLALOH <- NULL
txH2LUAD$FDorHLALOH[txH2LUAD$FD == "High" & txH2LUAD$any.HLA.loss == "Loss"] <- "Both"
txH2LUAD$FDorHLALOH[txH2LUAD$FD == "High" & txH2LUAD$any.HLA.loss == "No loss"] <- "FD"
txH2LUAD$FDorHLALOH[txH2LUAD$FD == "Low" & txH2LUAD$any.HLA.loss == "Loss"] <- "HLALOH"
txH2LUAD$FDorHLALOH[txH2LUAD$FD == "Low" & txH2LUAD$any.HLA.loss == "No loss"] <- "None"
ggstatsplot::ggbarstats(
  data = txH2LUAD,
  main = FDorHLALOH,
  condition = coldP3,
  bf.message = TRUE,
  sampling.plan = "indepMulti",
  title = "Lung adenocarcinoma (n=49)",
  xlab = "Immune cold regions",
  perc.k = 1,
  x.axis.orientation = "slant",
  #ggtheme = hrbrthemes::theme_modern_rc(),
  ggstatsplot.layer = FALSE,
  #ggplot.component = ggplot2::theme(axis.text.x = ggplot2::element_text(face = "italic")),
  palette = "Set2",
  messages = FALSE,
  data.label = "both",
  ggtheme = ggplot2::theme_classic())+
  scale_fill_manual(values = c("#ADB7BA","#FF0000", "#ACFF2F",  "#C8AA1F"))

#FD distribution 
ggdensity(regional, x = "fd", color = "Histology", palette = c("#98AED6", "#682A45", "black"),
          add = "median", rug = TRUE) +
  theme(legend.position="none")

#correlations between FD and stromal%, tumor%
ggscatter(regional, x = "fibroblasts_per", y = "fd", show.legend.text = F, 
                xlab = "Stromal%", ylab = "Fractal dimension", font.label = 4,
                shape = 19, size = 2,
                color = "Histology", palette = c("#98AED6", "#682A45", "black"),
                add.params = list(color = "grey50", fill = "azure3"))+ 
  stat_cor(aes(color = Histology), label.x = 50, method="spearman", show.legend = F)+
  theme(legend.position="none")
ggscatter(regional, x = "fibroblasts_per", y = "fd", cor.coef = T)

ggscatter(regional, x = "tumour_per", y = "fd", show.legend.text = F, 
                xlab = "Tumour%", ylab = "Fractal dimension", font.label = 4,
                shape = 19, size = 2,
                color = "Histology", palette = c("#98AED6", "#682A45", "black"),
                add.params = list(color = "grey50", fill = "azure3"))+ 
  stat_cor(aes(color = Histology), label.x = 10, method="spearman", show.legend = F)+
  theme(legend.position="none")
ggscatter(regional, x = "tumour_per", y = "fd", cor.coef = T)

#average distance from a stromal cell to nearest cancer cell (each row is a region)
#read the csv stromal_tumor_cellAvgDist.csv
ggplot(distD, aes(minDistMean_um)) +
  geom_density(color="darkblue", fill="pink")+xlim(0,300) + theme_classic()

#############

#4b, d - S7 a,c,e: ATL~ClonalNeo
#############
diagnosticLUAD = diagnostic[diagnostic$histology_group=="Adenocarcinoma",]
diagnosticLUAD$ClonalNeoMedian <- NULL
diagnosticLUAD$ClonalNeoMedian[diagnosticLUAD$ClonalNeo>=quantile(diagnosticLUAD$ClonalNeo)[3]] <-">=Median"
diagnosticLUAD$ClonalNeoMedian[diagnosticLUAD$ClonalNeo<quantile(diagnosticLUAD$ClonalNeo)[3]] <-"<Median"
comList <- list( c(">=Median", "<Median"))
ggboxplot(diagnosticLUAD, x = "ClonalNeoMedian", y = "ATL_fibroRatio", 
          xlab = "Clonal neoantigens", ylab = "Adjacent-tumor lymphocytes/stroma", 
          color = "ClonalNeoMedian", palette = c("red2", "blue"), #title = "LUAD",
          add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))

#deconvulted immune subsets in ATL vs DTL vs ITL (converted as percentages)
comList = list(c("ATL", "DTL"), c("DTL", "ITL"), c("ATL", "ITL"))
ggboxplot(reg_tils_immune, x = "TILClass", y = "cd8foxp3_Ratio", color = "TILClass", palette = c("blue", "#8291F7", "black"),
          add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox", paired = T)+
  theme(legend.position="") + theme(text = element_text(size=18))

ggboxplot(reg_tils_immune, x = "TILClass", y = "cd8", color = "TILClass", palette = c("blue", "#8291F7", "black"),
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox", paired = T)+
  theme(legend.position="") + theme(text = element_text(size=18))

ggboxplot(reg_tils_immune, x = "TILClass", y = "cd4", color = "TILClass", palette = c("blue", "#8291F7", "black"),
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox", paired = T)+
  theme(legend.position="") + theme(text = element_text(size=18))


#immune features vs NeoAgs with p-val correction
toCorrLUAD <- select(tx[tx$Histology=="LUAD",], lymphocytes_per, ATL_fibroRatio, ITL_tumorRatio, cd4_per, cd8_per, cd8foxp3Ratio,
                     ClonalNeo, SubclonalNeo, Regions.LOHHLAloss.A, Regions.LOHHLAloss.B, Regions.LOHHLAloss.C)
#no HLA event (replace NAs with 0)
toCorrLUAD[is.na(toCorrLUAD)] = 0
p.rLUAD <- cor_pmat(toCorrLUAD, method=c("spearman"), use="complete.obs")
p.rLUAD <- matrix(p.adjust(p.rLUAD, method = "BH", n = length(p.rLUAD)), nrow = 11, ncol = 11) ##
corRLUAD <- round(cor(toCorrLUAD, method = c("spearman"), use="complete.obs"),2)
ggcorrplotModified(corRLUAD[7:11, 1:6], lab = T, sig.level=0.05, insig ="blank", lab.notsig="", 
                   colors = c("#4DBBD5B2","grey98","#DC0000B2"), legend.title = "Spearman Corr",
                   p.mat = p.rLUAD[7:11, 1:6])

#S7 HE-IHC reg. eval. 
ggdensity(regEval, x = "distance", color = "purple", xlab = "HE-IHC registration distance (um)",
          add = "mean", rug = TRUE)
#############

#S7d - TCGA corrs: ATL vs RNA gene signatures 
#############
p.r <- cor_pmat(tcgaluad_immune[, -c(1)], method=c("spearman"), use="complete.obs")
p.r <- matrix(p.adjust(p.r, method = "BH", n = length(p.r)), nrow = 11, ncol = 11) ##
corR <- round(cor(tcgaluad_immune[, -c(1)], method = c("spearman"), use="complete.obs"),2)
ggcorrplotModified(corR[4:11, 1:3], lab = T, sig.level=0.05, insig ="blank", lab.notsig="", 
                   colors = c("#4DBBD5B2","grey98","#DC0000B2"), legend.title = "Spearman Corr",
                   p.mat = p.r[4:11, 1:3], title = "TCGA LUAD (n=464)")
#############

#S2 a-d: DL validation corrs
#############
#diagnostic lym% and ATL vs IHC cd8 and stromal TILs
ggscatter(diagnostic, x = "Stromal.TIL", y = "ATL_fibroRatio", show.legend.text = F, 
                xlab = "Stromal TILs", ylab = "ATL/Stroma", font.label = 4,
                shape = 19, size = 2,
                color = "Histology", palette = c("#98AED6", "#682A45", "black"),
                add.params = list(color = "grey50", fill = "azure3"))+ 
  stat_cor(aes(color = Histology), label.x = 15, method="spearman", show.legend = F)+
  theme(legend.position="none")
ggscatter(diagnostic, x = "Stromal.TIL", y = "ATL_fibroRatio", cor.coef = T) 

ggscatter(diagnostic, x = "lymphocytes_per", y = "cd8_per", show.legend.text = F, 
              xlab = "H&E Lymphocytes%", ylab = "IHC CD8%", font.label = 4,
              shape = 19, size = 2,
              color = "Histology", palette = c("#98AED6", "#682A45", "black"),
              add.params = list(color = "grey50", fill = "azure3"))+ 
  stat_cor(aes(color = Histology), label.x = 5, method="spearman", show.legend = F)+
  theme(legend.position="none")
ggscatter(diagnostic, x = "lymphocytes_per", y = "cd8_per", cor.coef = T)

#regional DNA tumor purity 
ggscatter(regional, x = "ASCAT.purity", y = "tumour_per", show.legend.text = F, 
                xlab = "ASCAT Purity", ylab = "Tumour%", font.label = 4,
                shape = 19, size = 2,
                color = "Histology", palette = c("#98AED6", "#682A45", "black"),
                add.params = list(color = "grey50", fill = "azure3"))+ 
  stat_cor(aes(color = Histology), label.x = 0.75, method="spearman", show.legend = F)+
  theme(legend.position="none")
ggscatter(regional, x = "ASCAT.purity", y = "tumour_per", cor.coef = T)

ggscatter(regional, x = "VAF.purity", y = "tumour_per", show.legend.text = F, 
                xlab = "VAF Purity", ylab = "Tumour%", font.label = 4,
                shape = 19, size = 2,
                color = "Histology", palette = c("#98AED6", "#682A45", "black"),
                add.params = list(color = "grey50", fill = "azure3"))+ 
  stat_cor(aes(color = Histology), label.x = 0.75, method="spearman", show.legend = F)+
  theme(legend.position="none")
ggscatter(regional, x = "VAF.purity", y = "tumour_per", cor.coef = T)

#Danaher cd8 signature from RNA
ggscatter(regional, x = "cd8.score.danaher", y = "lymphocytes_per", show.legend.text = F, 
                xlab = "CD8 score (Danaher et. al. signature)", ylab = "Lymphocytes%", font.label = 4,
                shape = 19, size = 2,
                color = "Histology", palette = c("#98AED6", "#682A45", "black"),
                add.params = list(color = "grey50", fill = "azure3"))+ 
  stat_cor(aes(color = Histology), label.x = 1, method="spearman", show.legend = F)+
  theme(legend.position="none")
ggscatter(regional, x = "cd8.score.danaher", y = "lymphocytes_per", cor.coef = T)
#############

#S2 e: LATTICe-A validation corrs
#############
ggscatter(lattice, x = "pathology_LymEstimate", y = "lymphocytes_per", add = "reg.line", 
          xlab = "Lymphocyte cell fraction (scored by pathologist)", ylab = "Lymphocyte%", font.label = 4,
          color = "#98AED6",shape = 19,
          conf.int = TRUE, size = 2,
          add.params = list(color = "grey50", fill = "azure3"), 
          cor.coef = TRUE, cor.method = "spearman") 
ggscatter(latticeTIL, x = "Stromal.TIL", y = "ATL_fibroRatio", add = "reg.line", 
          xlab = "Stromal TILs (scored by pathologist)", ylab = "ATL/Stroma", font.label = 4,
          color = "#98AED6",shape = 19,
          conf.int = TRUE, size = 2,
          add.params = list(color = "grey50", fill = "azure3"), 
          cor.coef = TRUE, cor.method = "spearman") 

#############

#S3 a-e: regional lym dist
#############
ggdensity(tx[! tx$Histology=="Other",], x = "lymphocytes_per_sd", color = "Histology", palette = c("#98AED6", "#682A45"),
          add = "mean", rug = TRUE)
ggdensity(tx[! tx$Histology=="Other",], x = "lymphocytes_per_sd", color = "pathologyTNM", 
          add = "mean", rug = TRUE)
ggdensity(tx[! tx$Histology=="Other",], x = "lymphocytes_per_mean", color = "pathologyTNM", 
          add = "mean", rug = TRUE)

ggstatsplot::ggbetweenstats(
  data = tx[! tx$Histology=="Other",], type = "np",
  x = pathologyTNM2, xlab = "Pathological stage",
  y = lymphocytes_per_sd, ylab = "Lymphocyte% (SD)", 
  pairwise.comparisons = T,      
  pairwise.annotation = "p.value",  
  pairwise.display = "significant",
  p.adjust.method = "BH",   
  bf.message = TRUE,   
  conf.level = 0.99,   
  ggplot.component = list(   
    ggplot2::scale_y_continuous(sec.axis = ggplot2::dup_axis())),
  k = 3,
  messages = T,
  ggtheme = ggplot2::theme_classic(),
  mean.color = c("#f3776e", "#b59f31", "#32b24b", "#1dbdc2", "#9268ad"))+
  scale_color_manual(values = c("#f3776e", "#b59f31", "#32b24b", "#1dbdc2", "#9268ad"))
ggstatsplot::ggbetweenstats(
  data = tx[! tx$Histology=="Other",],  type = "np",
  x = pathologyTNM2, xlab = "Pathological stage",
  y = lymphocytes_per_mean, ylab = "Lymphocyte% (mean)", 
  pairwise.comparisons = T,      
  pairwise.annotation = "p.value",  
  pairwise.display = "significant",
  p.adjust.method = "BH",   
  bf.message = TRUE, 
  conf.level = 0.99,          
  ggplot.component = list(    
    ggplot2::scale_y_continuous(sec.axis = ggplot2::dup_axis())),
  k = 3,
  messages = T,
  ggtheme = ggplot2::theme_classic(),
  mean.color = c("#f3776e", "#b59f31", "#32b24b", "#1dbdc2", "#9268ad"))+
  scale_color_manual(values = c("#f3776e", "#b59f31", "#32b24b", "#1dbdc2", "#9268ad"))
#############

#S4: RNA vs deep learning 
#############
r = regional[order(regional$Histology, decreasing = F),]
r$orig_immune_cluster[is.na(r$orig_immune_cluster)] <- "NA"
ggdotchart(r, x = "PublicationID", y = "lymphocytes_per",
           color = "orig_immune_cluster",                                # Color by groups.
           group = "Histology",
           sort.by.groups = TRUE,           
           palette = c("red2", "blue2", "gray45"), # Custom color palette
           #sorting = "descending",                       # Sort value in descending order
           #rotate = TRUE,                                # Rotate vertically
           #add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 0.8),
           dot.size = 3,                                 # Large dot size
           #label = min(regionalT$lymphocytes_per),
           #y.text.col = TRUE,                            # Color y text by groups
           ggtheme = theme_pubr())

#test overlap
r = colNA(regional, "orig_immune_cluster")
r = r[! r$immuneClass_2=="intermediate",]
fisher.test(table(r$immuneClass_2, r$orig_immune_cluster))

#comparisons
r = regional[! regional$immuneClass_2=="intermediate",]
comList <- list( c("hot", "cold"))
r <- r[order(r$immuneClass_2),] 
ggboxplot(r, x = "immuneClass_2", y = "lymphocytes_per", color = "immuneClass_2",
              add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",comparisons = comList, method = "wilcox")+ 
  theme(text = element_text(size=18))
ggboxplot(r, x = "immuneClass_2", y = "ITL_tumorRatio", color = "immuneClass_2",
              add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",comparisons = comList, method = "wilcox")+ 
  theme(text = element_text(size=18))
ggboxplot(r, x = "immuneClass_2", y = "ATL_fibroRatio", color = "immuneClass_2",
              add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",comparisons = comList, method = "wilcox")+ 
  theme(text = element_text(size=18))
ggboxplot(r, x = "immuneClass_2", y = "ASCAT.purity", color = "immuneClass_2",
              add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",comparisons = comList, method = "wilcox")+ 
  theme(text = element_text(size=18))

r <- colNA(regional, "orig_immune_cluster")
comList <- list( c("low", "high"))
r <- r[order(r$orig_immune_cluster, decreasing = T),] 
ggboxplot(r, x = "orig_immune_cluster", y = "lymphocytes_per", color = "orig_immune_cluster",
          add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",comparisons = comList, method = "wilcox")+ 
  theme(text = element_text(size=18))
ggboxplot(r, x = "orig_immune_cluster", y = "ITL_tumorRatio", color = "orig_immune_cluster",
          add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",comparisons = comList, method = "wilcox")+ 
  theme(text = element_text(size=18))
ggboxplot(r, x = "orig_immune_cluster", y = "ATL_fibroRatio", color = "orig_immune_cluster",
          add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",comparisons = comList, method = "wilcox")+ 
  theme(text = element_text(size=18))
ggboxplot(r, x = "orig_immune_cluster", y = "ASCAT.purity", color = "orig_immune_cluster",
          add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",comparisons = comList, method = "wilcox")+ 
  theme(text = element_text(size=18))

#immune hotspot analysis
r = colNA(regional, "orig_immune_cluster")
r = r[! r$immuneClass_2=="intermediate",]
r$discHC = NULL
r$discHC[r$immuneClass_2=="hot" & r$orig_immune_cluster == "high"] = "TRUE"
r$discHC[r$immuneClass_2=="cold" & r$orig_immune_cluster == "low"] = "TRUE"
r$discHC[is.na(r$discHC)] = "FALSE"
comList <- list( c("TRUE", "FALSE"))
ggboxplot(r, x = "discHC", y = "fraction_immune_hotspot", color = "discHC", palette = c("#f4abc2", "#97f2a4"), ylab = "Fraction of immune hotspot",
          add = "dotplot", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())

#against path. lym fraction for overall discordant regions
ggscatter(regional, x = "lymphocytes_per", y = "pathology_LymEstimate", add = "reg.line",
          xlab = "Lymphocyte%", ylab = "Lymphocyte cell fraction (Pathologist 2)",
          color = "black",
          conf.int = TRUE, size = 2,
          add.params = list(color = "grey50", fill = "azure3"), 
          cor.coef = TRUE, cor.method = "spearman")+ theme(legend.position="none")+xlim(10, 60)+ylim(5, 75)+ theme(text = element_text(size=14))

#############

#3d-f and S5: survival using the number of cold regions in TRACERx
#############
txH2LUADLUSC <- tx[ which(tx$Histology=="LUAD" 
                            | tx$Histology == "LUSC"), ]
txH2LUADLUSC$Histology = as.factor(as.character(txH2LUADLUSC$Histology))
txH2LUADLUSC$sex = as.factor(as.character(txH2LUADLUSC$sex))
txH2LUADLUSC$coldP3 = as.factor(as.character(txH2LUADLUSC$coldP3))
txH2LUADLUSC$pathologyTNM2 = as.factor(as.character(txH2LUADLUSC$pathologyTNM2))
#with clonal neo local quartile for the available 79 luad lusc patients: 
#quantile(txH2LUADLUSC$ClonalNeo[txH2LUADLUSC$Histology=="LUSC"])
txH2LUADLUSC$ClonalNeoUQ[txH2LUADLUSC$Histology=="LUAD" & txH2LUADLUSC$ClonalNeo >= 189]<- "UQ.high"
txH2LUADLUSC$ClonalNeoUQ[txH2LUADLUSC$Histology=="LUAD" & txH2LUADLUSC$ClonalNeo < 189]<- "cUQ.low"
txH2LUADLUSC$ClonalNeoUQ[txH2LUADLUSC$Histology=="LUSC" & txH2LUADLUSC$ClonalNeo >= 206.25]<- "UQ.high"
txH2LUADLUSC$ClonalNeoUQ[txH2LUADLUSC$Histology=="LUSC" & txH2LUADLUSC$ClonalNeo < 206.25]<- "cUQ.low"
txH2LUADLUSC$ClonalNeoUQ <- as.factor(as.character(txH2LUADLUSC$ClonalNeoUQ))

#first test, without the genomic score
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~coldP3+age+sex+pack_years_calculated+Histology+pathologyTNM2+Adjuvant.therapy, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)

#KM
fit <- survfit(Surv(DFS_time_days, DFS_censor_variable)~rCold<=1, data = txH2LUADLUSC) 
ggsurvplot(fit, data = txH2LUADLUSC, conf.int = FALSE,
           pval = TRUE, pval.size = 5, pval.coord = c(0.2, 0.1),
           linetype = "solid",
           #surv.median.line = "hv",
           legend = "none", legend.title = "", title = "LUAD+LUSC",
           legend.labs = c("High", "Low"),
           surv.plot.height = 0.7, palette = c("red2", "blue"),
           risk.table = TRUE,
           risk.table.col = "black", break.time.by = 200,
           tables.height = 0.25,
           tables.theme = theme_cleantable(),
           tables.y.text = TRUE, risk.table.title = "Number at Risk",
           tables.x.text = "", xlim = c(0, 1400),
           xlab = "Days to Death or Recurrence", ylab = "Disease-free Survival")
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~coldP3+age+sex+pack_years_calculated+Histology+pathologyTNM2+Adjuvant.therapy+ClonalNeoUQ, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)
#also reported in text, univariate as cont.:
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~rCold, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)
#or at the median again:
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~coldP3, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)

#including total number of regions
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~coldP3+age+sex+pack_years_calculated+Histology+pathologyTNM2+Adjuvant.therapy+ClonalNeoUQ+nTotalRegions, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)

#as continuous, keeping it simple 
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~rCold+age+sex+pack_years_calculated+Histology+pathologyTNM2+Adjuvant.therapy, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)

#as continuous, testing independce against other scores
#correcting for Histology (by adding +Histology) in each fit doesn't make a big difference
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~rCold+subclonalCNA, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~rCold+mean_ASCAT.purity, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~rCold+tumour_per_mean, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)

#the original test including other histology patients (only 6 more), cant do this with ClonalNeo... 
fit <- survfit(Surv(DFS_time_days, DFS_censor_variable)~rCold<=1, data = tx) 
ggsurvplot(fit, data = tx, conf.int = FALSE,
           pval = TRUE, pval.size = 5, pval.coord = c(0.2, 0.1),
           linetype = "solid",
           #surv.median.line = "hv",
           legend = "none", legend.title = "", title = "ALL",
           legend.labs = c("High", "Low"),
           surv.plot.height = 0.7, palette = c("red2", "blue"),
           risk.table = TRUE,
           risk.table.col = "black", break.time.by = 200,
           tables.height = 0.25,
           tables.theme = theme_cleantable(),
           tables.y.text = TRUE, risk.table.title = "Number at Risk",
           tables.x.text = "", xlim = c(0, 1400),
           xlab = "Days to Death or Recurrence", ylab = "Disease-free Survival")
txS = tx
txS$Histology = as.factor(as.character(txS$Histology))
txS$sex = as.factor(as.character(txS$sex))
txS$coldP3 = as.factor(as.character(txS$coldP3))
txS$pathologyTNM2 = as.factor(as.character(txS$pathologyTNM2))
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~coldP3+age+sex+pack_years_calculated+Histology+pathologyTNM2+Adjuvant.therapy, data = txS)
ggforest(fit, data=txS)

#lastly HR comparison panel f
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~rCold, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~coldP3, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~coldProp, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~rHot, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~lymphocytes_per, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~lymphocytes_per_min, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~lymphocytes_per_mean, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~lymphocytes_per_sd, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~cd8_per, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~cd4_per, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~cd8foxp3_Ratio, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~cd4foxp3_Ratio, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)

summaryPs <- data.frame(
  Measure = c("9N cold regions", "99>1 cold regions", "99 N cold regions/total",
              "8N hot regions",
              "7Diagnostic lymphocyte%", "6Regional lymphocyte% (min)", "5Regional lymphocyte% (mean)",
              "55Regional lymphocyte% (SD)",
              "4CD8+%", "3CD4+FOXP3-%", "2CD8+%/CD4+FOXP3+%", "1CD4+FOXP3-%/CD4+FOXP3+%"),
  PVal = c(0.003, 0.005, 0.098,
           0.795,
           0.05, 0.058, 0.125,
           0.855,
           0.631, 0.517, 0.713, 0.804),
  lower95 = c(1.1, 1.3, 0.87,
              0.75,
              0.86, 0.94, 0.94,
              0.93,
              0.89, 0.91, 0.7, 0.86),
  upper95 = c(1.8, 5.3, 5.3,
              1.2,
              1, 1,1,
              1.1,
              1.1, 1, 1.3, 1.2), 
  hr = c(1.4, 2.7, 2.2,
         0.97,
         0.93, 0.97, 0.97,
         1,
         0.98, 0.98, 0.95, 1))
summaryPs$PVal_log10 <- -log10(summaryPs$PVal)
summaryPs$PvalAdjust <- p.adjust(summaryPs$PVal, method = "BH")
summaryPs$PvalAdjust_log10 <- -log10(summaryPs$PvalAdjust)
ggplot(data=summaryPs,aes(x=hr,y=Measure))+
  geom_point(aes(size=PVal_log10), colour="#CB6728", fill="#CB6728", shape=21)+
  geom_errorbarh(aes(xmin=lower95,xmax=upper95), height=0)+
  geom_vline(xintercept=1,linetype="longdash", colour="grey70")+
  scale_size_continuous(breaks=c(5000,10000,15000))+theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank()) + labs(x = "Hazard ratio", y = "Immune feature")
ggplot(data=summaryPs,aes(x=hr,y=Measure))+
  geom_point(aes(size=PvalAdjust_log10), colour="#CB6728", fill="#CB6728", shape=21)+
  geom_errorbarh(aes(xmin=lower95,xmax=upper95), height=0)+
  geom_vline(xintercept=1,linetype="longdash", colour="grey70")+
  scale_size_continuous(breaks=c(5000,10000,15000))+theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank()) + labs(x = "Hazard ratio", y = "Immune feature")
#############

#3e-g and S5: survival using the number of cold regions in LATTICe-A (validation)
#############
load(("X:/Dropbox (ICR)/yuanlab/Manuscripts/Lung Tx1/revision2/data/latticea_970_finalCohort.RData"))
tmp <- lattice
tmp$coldP3 = as.factor(as.character(tmp$coldP3))
tmp$gender = as.factor(as.character(tmp$gender))
tmp$Adjuvant.therapy = as.factor(as.character(tmp$Adjuvant.therapy))
tmp$stage = as.factor(as.character(tmp$stage))
tmplatt=tmp
levels(tmplatt$stage)[levels(tmplatt$stage)=="IIIA"] <- "III"
levels(tmplatt$stage)[levels(tmplatt$stage)=="IIIB"] <- "III"
tmplatt = tmplatt[!is.na(tmplatt$age_surgery) & !is.na(tmplatt$gender) &
                    !is.na(tmplatt$pack_years) & !is.na(tmplatt$stage) &
                    !is.na(tmplatt$Adjuvant.therapy) & !is.na(tmplatt$time_to_rfs) &
                    !is.na(tmplatt$rfs_status),]
tmp_model <- coxph(Surv(time_to_rfs, rfs_status == '1') ~ coldP3 + age_surgery + gender + pack_years + stage + Adjuvant.therapy, data=tmplatt)
ggforest(tmp_model, data = NULL, main = "", fontsize = 0.7, refLabel = "", noDigits = 2)
tmp_model <- coxph(Surv(time_to_rfs, rfs_status == '1') ~ coldP3 + age_surgery + gender + pack_years + stage + Adjuvant.therapy + nTotalRegions, data=tmplatt)
ggforest(tmp_model, data = NULL, main = "", fontsize = 0.7, refLabel = "", noDigits = 2)
tmplatt=tmp
tmplatt = tmplatt[!is.na(tmplatt$age_surgery) & !is.na(tmplatt$gender) &
                    !is.na(tmplatt$stage) &
                    !is.na(tmplatt$Adjuvant.therapy) & !is.na(tmplatt$time_to_rfs) &
                    !is.na(tmplatt$rfs_status),]
tmp_model <- coxph(Surv(time_to_rfs, rfs_status == '1') ~ coldP3 + age_surgery + gender + pack_years + stage + Adjuvant.therapy + nTotalRegions, data=tmplatt)
ggforest(tmp_model, data = NULL, main = "", fontsize = 0.7, refLabel = "", noDigits = 2)

#KM
fit <- survfit(Surv(time_to_rfs, rfs_status)~cold<=1, data = lattice)
ggsurvplot(fit, data = lattice, conf.int = FALSE,
           pval = TRUE, pval.size = 5, pval.coord = c(0.2, 0.1),
           linetype = "solid",
           #surv.median.line = "hv",
           legend = "none", legend.title = "", #title = "N cold regions (median) <= 1",
           legend.labs = c("High", "Low"),
           surv.plot.height = 0.7, palette = c("#81391b", "#ef5123"),
           risk.table = TRUE,
           risk.table.col = "black", break.time.by = 1000,
           tables.height = 0.25,
           tables.theme = theme_cleantable(),
           tables.y.text = TRUE, risk.table.title = "Number at Risk",
           tables.x.text = "", xlim = c(0, 6345),
           xlab = "Days to Death or Recurrence", ylab = "Disease-free Survival")

#head to head, image-immune features: 
# lt$immuneClass[lt$lymphocytes_per >=median(lt$lymphocytes_per)] <- "hot"
# lt$immuneClass[lt$lymphocytes_per <median(lt$lymphocytes_per)] <- "cold"
# lt$immuneClass_2 = lt$immuneClass
# lt$immuneClass_2[lt$lymphocytes_per >= (median(lt$lymphocytes_per)-sd(lt$lymphocytes_per)/4) &
#                    lt$lymphocytes_per <= (median(lt$lymphocytes_per)+sd(lt$lymphocytes_per)/4)] <- "intermediate"

rU2 = lt %>%
  group_by(ACA_ID) %>%
  summarise(lymphocytes_per_mean = mean(lymphocytes_per),
            lymphocytes_per_min = min(lymphocytes_per), 
            lymphocytes_per_sd = sd(lymphocytes_per))
lattice_lym = NULL
lattice_lym = merge(lattice, rU2)
fit <- coxph(Surv(time_to_rfs, rfs_status)~cold, data = lattice_lym)
ggforest(fit, data=lattice_lym)
fit <- coxph(Surv(time_to_rfs, rfs_status)~coldP3, data = lattice_lym)
ggforest(fit, data=lattice_lym)
lattice_lym$coldProp = lattice_lym$cold/lattice_lym$nTotalRegions
fit <- coxph(Surv(time_to_rfs, rfs_status)~cold+coldProp, data = lattice_lym)
ggforest(fit, data=lattice_lym)
fit <- coxph(Surv(time_to_rfs, rfs_status)~cold+hot, data = lattice_lym)
ggforest(fit, data=lattice_lym)
fit <- coxph(Surv(time_to_rfs, rfs_status)~cold+lymphocytes_per_min, data = lattice_lym)
ggforest(fit, data=lattice_lym)
fit <- coxph(Surv(time_to_rfs, rfs_status)~cold+lymphocytes_per_mean, data = lattice_lym)
ggforest(fit, data=lattice_lym)
fit <- coxph(Surv(time_to_rfs, rfs_status)~cold+lymphocytes_per_sd, data = lattice_lym)
ggforest(fit, data=lattice_lym)
#############

#S3 LATTICe-A: lym_sd/mean vd stage
#############
p1 = ggstatsplot::ggbetweenstats(
  data = lattice_lym, type = "np",
  x = stage, xlab = "Pathological stage",
  y = lymphocytes_per_sd, ylab = "Lymphocyte% (SD)", 
  pairwise.comparisons = T,      # display significant pairwise comparisons
  pairwise.annotation = "p.value",  # how do you want to annotate the pairwise comparisons
  pairwise.display = "significant",
  p.adjust.method = "BH",   # method for adjusting p-values for multiple comparisons
  bf.message = TRUE,                # display Bayes Factor in favor of the null hypothesis
  conf.level = 0.99,                # changing confidence level to 99%
  k = 3,
  messages = T,
  ggtheme = ggplot2::theme_classic(),
  mean.color = c("#f3776e", "#b59f31", "#32b24b", "#1dbdc2", "#9268ad"))+
  scale_color_manual(values = c("#f3776e", "#b59f31", "#32b24b", "#1dbdc2", "#9268ad"))
p2 = ggstatsplot::ggbetweenstats(
  data = lattice_lym, type = "np",
  x = stage, xlab = "Pathological stage",
  y = lymphocytes_per_mean, ylab = "Lymphocyte% (mean)", 
  pairwise.comparisons = T,      # display significant pairwise comparisons
  pairwise.annotation = "p.value",  # how do you want to annotate the pairwise comparisons
  pairwise.display = "significant",
  p.adjust.method = "BH",   # method for adjusting p-values for multiple comparisons
  bf.message = TRUE,                # display Bayes Factor in favor of the null hypothesis
  conf.level = 0.99,                # changing confidence level to 99%
  k = 3,
  messages = T,
  ggtheme = ggplot2::theme_classic(),
  mean.color = c("#f3776e", "#b59f31", "#32b24b", "#1dbdc2", "#9268ad"))+
  scale_color_manual(values = c("#f3776e", "#b59f31", "#32b24b", "#1dbdc2", "#9268ad"))

lym_clean = colNA(lattice_lym, "stage")

p4 = ggdensity(lym_clean, x = "lymphocytes_per_sd", color = "#98AED6", xlab = "Lymphocyte% (regional SD)",
               add = "median", rug = TRUE)+
  scale_color_manual(values = c("#f3776e", "#b59f31", "#32b24b", "#1dbdc2", "#9268ad"))+ theme(legend.position = "none")
p5 = ggdensity(lym_clean, x = "lymphocytes_per_sd", color = "stage", xlab = "Lymphocyte% (regional SD)",
               add = "median", rug = TRUE)+
  scale_color_manual(values = c("#f3776e", "#b59f31", "#32b24b", "#1dbdc2", "#9268ad"))+ theme(legend.position = "none")
p6 = ggdensity(lym_clean, x = "lymphocytes_per_mean", color = "stage", xlab = "Lymphocyte% (regional mean)",
               add = "median", rug = TRUE)+
  scale_color_manual(values = c("#f3776e", "#b59f31", "#32b24b", "#1dbdc2", "#9268ad"))+ theme(legend.position = "none")

#############

#S6: immune phenotype classification 
#############
#thresholds visualisation 
ggdensity(regional, x = "lymphocytes_per", 
          add = "median", rug = T)+ 
  geom_vline(xintercept = 14.55154, linetype="dotted", 
             color = "gray", size=0.9)+ 
  geom_vline(xintercept = 38.02652, linetype="dotted", 
             color = "gray", size=0.9)+ 
  geom_vline(xintercept = 20.42028, linetype="dotted", 
             color = "gray", size=0.9)+ 
  geom_vline(xintercept = 32.15778, linetype="dotted", 
             color = "gray", size=0.9)+ 
  geom_vline(xintercept = 23.35466, linetype="dotted", 
             color = "gray", size=0.9)+ 
  geom_vline(xintercept = 29.2234, linetype="dotted", 
             color = "gray", size=0.9)+ 
  geom_vline(xintercept = 22.37653, linetype="dotted", 
             color = "gray", size=0.9)+ 
  geom_vline(xintercept = 30.20153, linetype="dotted", 
             color = "gray", size=0.9)+ 
  geom_vline(xintercept = 24.33278, linetype="dotted", 
             color = "gray", size=0.9)+ 
  geom_vline(xintercept = 28.24528, linetype="dotted", 
             color = "gray", size=0.9)

#test survival in each threshold
txH2LUADLUSC <- tx[ which(tx$Histology=="LUAD" 
                          | tx$Histology == "LUSC"), ]
txH2LUADLUSC$Histology = as.factor(as.character(txH2LUADLUSC$Histology))
txH2LUADLUSC$sex = as.factor(as.character(txH2LUADLUSC$sex))
txH2LUADLUSC$coldP3 = as.factor(as.character(txH2LUADLUSC$coldP3))
txH2LUADLUSC$pathologyTNM2 = as.factor(as.character(txH2LUADLUSC$pathologyTNM2))
#with clonal neo local quartile for the available 79 luad lusc patients: 
#quantile(txH2LUADLUSC$ClonalNeo[txH2LUADLUSC$Histology=="LUSC"])
txH2LUADLUSC$ClonalNeoUQ[txH2LUADLUSC$Histology=="LUAD" & txH2LUADLUSC$ClonalNeo >= 189]<- "UQ.high"
txH2LUADLUSC$ClonalNeoUQ[txH2LUADLUSC$Histology=="LUAD" & txH2LUADLUSC$ClonalNeo < 189]<- "cUQ.low"
txH2LUADLUSC$ClonalNeoUQ[txH2LUADLUSC$Histology=="LUSC" & txH2LUADLUSC$ClonalNeo >= 206.25]<- "UQ.high"
txH2LUADLUSC$ClonalNeoUQ[txH2LUADLUSC$Histology=="LUSC" & txH2LUADLUSC$ClonalNeo < 206.25]<- "cUQ.low"
txH2LUADLUSC$ClonalNeoUQ <- as.factor(as.character(txH2LUADLUSC$ClonalNeoUQ))

reg = regional[, c("PublicationID", "region", "lymphocytes_per")]
reg$imC[reg$lymphocytes_per >=median(reg$lymphocytes_per)] <- "hot"
reg$imC[reg$lymphocytes_per <median(reg$lymphocytes_per)] <- "cold"
reg$imC1 = reg$imC
reg$imC1[reg$lymphocytes_per >= (median(reg$lymphocytes_per)-sd(reg$lymphocytes_per)/1) &
           reg$lymphocytes_per <= (median(reg$lymphocytes_per)+sd(reg$lymphocytes_per)/1)] <- "intermediate"
reg$imC2 = reg$imC
reg$imC2[reg$lymphocytes_per >= (median(reg$lymphocytes_per)-sd(reg$lymphocytes_per)/2) &
           reg$lymphocytes_per <= (median(reg$lymphocytes_per)+sd(reg$lymphocytes_per)/2)] <- "intermediate"
reg$imC3 = reg$imC
reg$imC3[reg$lymphocytes_per >= (median(reg$lymphocytes_per)-sd(reg$lymphocytes_per)/3) &
           reg$lymphocytes_per <= (median(reg$lymphocytes_per)+sd(reg$lymphocytes_per)/3)] <- "intermediate"
reg$imC6 = reg$imC
reg$imC6[reg$lymphocytes_per >= (median(reg$lymphocytes_per)-sd(reg$lymphocytes_per)/6) &
           reg$lymphocytes_per <= (median(reg$lymphocytes_per)+sd(reg$lymphocytes_per)/6)] <- "intermediate"

regS0 <- reg %>% group_by(PublicationID, imC) %>%
  dplyr::summarise(count=n())
regS0 = dcast(regS0, PublicationID ~ imC)
names(regS0) <- c("PublicationID", "rCold", "rHot")
regS0[is.na(regS0)] <- 0
tx0 = merge(txH2LUADLUSC, regS0)

regS1 <- reg %>% group_by(PublicationID, imC1) %>%
  dplyr::summarise(count=n())
regS1 = dcast(regS1, PublicationID ~ imC1)
names(regS1) <- c("PublicationID", "rCold", "rHot", "rIntermediate")
regS1[is.na(regS1)] <- 0
tx1 = merge(txH2LUADLUSC, regS1)

regS2 <- reg %>% group_by(PublicationID, imC2) %>%
  dplyr::summarise(count=n())
regS2 = dcast(regS2, PublicationID ~ imC2)
names(regS2) <- c("PublicationID", "rCold", "rHot", "rIntermediate")
regS2[is.na(regS2)] <- 0
tx2 = merge(txH2LUADLUSC, regS2)

regS3 <- reg %>% group_by(PublicationID, imC3) %>%
  dplyr::summarise(count=n())
regS3 = dcast(regS3, PublicationID ~ imC3)
names(regS3) <- c("PublicationID", "rCold", "rHot", "rIntermediate")
regS3[is.na(regS3)] <- 0
tx3 = merge(txH2LUADLUSC, regS3)

regS6 <- reg %>% group_by(PublicationID, imC6) %>%
  dplyr::summarise(count=n())
regS6 = dcast(regS6, PublicationID ~ imC6)
names(regS6) <- c("PublicationID", "rCold", "rHot", "rIntermediate")
regS6[is.na(regS6)] <- 0
tx6 = merge(txH2LUADLUSC, regS6)

tx0$coldP[tx0$rCold>median(tx0$rCold)] <- "med.High"
tx0$coldP[tx0$rCold<=median(tx0$rCold)] <- "Low"

tx1$coldP[tx1$rCold>median(tx1$rCold)] <- "med.High"
tx1$coldP[tx1$rCold<=median(tx1$rCold)] <- "Low"

tx2$coldP[tx2$rCold>median(tx2$rCold)] <- "med.High"
tx2$coldP[tx2$rCold<=median(tx2$rCold)] <- "Low"

tx3$coldP[tx3$rCold>median(tx3$rCold)] <- "med.High"
tx3$coldP[tx3$rCold<=median(tx3$rCold)] <- "Low"

tx6$coldP[tx6$rCold>median(tx6$rCold)] <- "med.High"
tx6$coldP[tx6$rCold<=median(tx6$rCold)] <- "Low"

tx0$coldP = as.factor(as.character(tx0$coldP))
tx1$coldP = as.factor(as.character(tx1$coldP))
tx2$coldP = as.factor(as.character(tx2$coldP))
tx3$coldP = as.factor(as.character(tx3$coldP))
tx6$coldP = as.factor(as.character(tx6$coldP))

fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~coldP+age+sex+pack_years_calculated+Histology+pathologyTNM2+Adjuvant.therapy+ClonalNeoUQ, data = tx0)
ggforest(fit, data=tx0, main = "Threshold: No intermediate zone")
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~coldP+age+sex+pack_years_calculated+Histology+pathologyTNM2+Adjuvant.therapy+ClonalNeoUQ, data = tx1)
ggforest(fit, data=tx1, main = "Threshold: Entire SD")
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~coldP+age+sex+pack_years_calculated+Histology+pathologyTNM2+Adjuvant.therapy+ClonalNeoUQ, data = tx2)
ggforest(fit, data=tx2, main = "Threshold: SD/2")
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~coldP+age+sex+pack_years_calculated+Histology+pathologyTNM2+Adjuvant.therapy+ClonalNeoUQ, data = tx3)
ggforest(fit, data=tx3, main = "Threshold: SD/3")
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~coldP+age+sex+pack_years_calculated+Histology+pathologyTNM2+Adjuvant.therapy+ClonalNeoUQ, data = tx6)
ggforest(fit, data=tx6, main = "Threshold: SD/6")

#cd8 danaher differnce in hot vs cold regions using all thresholds
comList <- list( c("cold", "hot"))
ggboxplot(regional[! regional$immuneClass_2=="intermediate" ,], x = "immuneClass_2", y = "cd8.score.danaher", color = "immuneClass_2", 
          add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))

#repeat for imC1, 2, ...6
reg2 = merge(reg, regional[, c("PublicationID", "region", "cd8.score.danaher"), ], by = c("PublicationID", "region"), all.x = T)
ggboxplot(reg2[! reg2$imC=="intermediate" ,], x = "imC", y = "cd8.score.danaher", color = "imC", 
          add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))

#distribution of lym% in tx and tcga for histology subtypes
ggdensity(regional, x = "lymphocytes_per", color = "Histology", palette = c("#98AED6", "#682A45", "black"),
          add = "median", rug = TRUE) +
  theme(legend.position="none")
ggdensity(tcga, x = "lym_per", color = "disease_code", palette = c("#98AED6", "#682A45"),
          add = "median", rug = TRUE) +
  theme(legend.position="none")

comList <- list( c("LUAD", "LUSC"))
ggboxplot(regional, x = "Histology", y = "lymphocytes_per", color = "Histology", title = "TRACERx (n=275)", 
          add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+
  scale_color_manual(values=c("#98AED6", "#682A45", "black"))
ggboxplot(tcga, x = "disease_code", y = "lym_per", color = "disease_code", title = "TCGA (n=939)", 
          add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+
  scale_color_manual(values=c("#98AED6", "#682A45"))

#############

#S8: summary stats heatmap
#############
txS8=tx
txS8$Histology[txS8$Histology=="Other"] <- "Z"
txS8 = txS8[order(txS8$lymphocytes_per_min, decreasing = F),]
txS8 = txS8[order(txS8$Histology, decreasing = F),]
txS8$DFS_yes[txS8$DFS_censor_variable==1] <- txS8$DFS_time_days
txS8$DFS_no[txS8$DFS_censor_variable==0] <- txS8$DFS_time_days

#dotplot same as Fig 2
pDotplot

#now the heatmap annotations for all immune scores
Heatmap(txS8$Histology, name = "Histology", width = unit(2.85, "mm"), 
        cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
        col = c("#98AED6", "#682A45", "black")) +
  Heatmap(txS8$FD, name = "FD", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = c("greenyellow", "hotpink"))+
  Heatmap(txS8$cd8_per, name = "CD8", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = colorRamp2(c(0, 17), c("#F2E2E2", "#8E0000")))+
  Heatmap(txS8$cd4_per, name = "CD4", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = colorRamp2(c(0, 28), c("#E3F0E2", "#067C00")))+
  Heatmap(txS8$foxp3_per, name = "FOXP3", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = colorRamp2(c(0, 21), c("#F9E2F2", "#CC008E")))+
  Heatmap(txS8$ATL_fibroRatio, name = "ATLF", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = colorRamp2(c(0, 1.6), c("#EDE2EF", "#650077")))+
  Heatmap(txS8$ClonalNeo, name = "cNEO", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          colorRamp2(c(0, 800), 
                     c("#E3E2EF", "#030070")))+
  Heatmap(txS8$any.HLA.loss, name = "HLA", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = c("red1", "#4DBBD5B2", "gray45"))+
  Heatmap(txS8$pathologyTNM, name = "Stage", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = c("#7FC97F" ,"#BEAED4" ,"#FDC086", "#FFFF99" ,"#386CB0", "#F0027F"))+
  Heatmap(txS8$pack_years_calculated, name = "Pack years", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = colorRamp2(c(0, 140), 
                           c("#E2E6E7", "#00252D")))+
  Heatmap(txS8$DFS_yes, name = "DFS_Yes", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = colorRamp2(c(0, 1364.0, NA), 
                           c("#FBF0E2", "#DD7800", "gray45")))+
  Heatmap(txS8$DFS_no, name = "DFS_No", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = colorRamp2(c(54.00, 1364.00, NA ), 
                           c("#E2EAE2", "#024402", "gray45")))
#############

#KA 20191201