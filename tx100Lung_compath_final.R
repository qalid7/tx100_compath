#### ###  #### ###
###Geospatial immune variability illuminates differential evolution of lung adenocarcinoma##
##ALL stat analyses and figures##
#### ###  #### ### 20200310


#R version notes
############# 
#mostly developed on R 3.5.1, but also tested on R 3.6.3 on 20200329

#for latest R if you struggle with ComplexHeatmap, do this: 
#R version 3.6.3
#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
#############

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
library(psych)
#############

#Functions
#############
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
setwd("/Users/kjabbar/Dropbox (ICR)/yuanlab/Manuscripts/Lung Tx1/code_data/")
load("./data/data_10.RData")
##
#TRACERx 100 cohort

#1# - regional: tx100 multiregion cohort (n=275 for 85 patients). includes main dig. path. scores e.g. lym$ and immune class. 
#ASCAT and VAF: DNA-based tumor purity. cd8.score.danaher: Danaher et al RNA-seq signature of CD8+. orig_immune_cluster: Rosenthal et al RNA-seq based immune classification. 
#pathology TILs: pathology TIL scoring using the Salgado et al method (n=260 regions).

#2# - diagnsotic: tx100 diagnostic cohort from both HE, IHC slides (n=100 for 100 patients). 
#Stromal.TIL: pathology scores on the first 100 slides. 

#3# - tx: combined tx100 regional and diagnsotic cohort at patient level (n=85, patients with a regional sample).
#both regional and disgnostic scores combined together with key genomic scores.

#4# - tx_surv15May2018: original tx100 survival data from 15 May 2018 - just kept here for the record. 
#5# - txclin_supp: provided to make the Supp Extended Data Fig 1 heatmap. 

##Other items
#6# - tx_genomic_pairdist: tx phenotype-based pairwise genomic dist. 
#7# - tx_genomic_pairdist_extended: tx phenotype-based pairwise genomic dist, ALL regions - to validate class. schemes. 
#8# - regional_immune_schemes: 4 more region-specific immune schemes. 

#9# - txluad_domClon_pairdist: tx phenotype-based pairwise phy. analysis. 
#10# - stromal_tumor_cellAvgDist: tx regional min and max dist between a tumor and a stromal cell for all 275 tumor regions.
#11# - reg_tils_immune: deconvulted immune cells onto the HE TIL classes for the same registered sections.
#12# - reg_eval: manual evaluation of HE-IHC registered sections.  
##

##
#LATTICe-A cohort

#13# - lattice_80: manual estimations for fraction of lym and path. TILs, and corresponding automated scores (n=80 patients LATTICe-A cohort).
#14# - latticea_bioai: list contains 3 DFs for TTF1, CD45, SMA. Rows correspond to image patches from cores and cols show the IHC/HE cell counts by deep learning.
#15# - lattice: entire LATTICe-A cohort (n=970 patients), all summarised patient level data.
#16# - lt: entire LATTICe-A immune data for all 4,324 H&E samples.
##


## data from Rosenthal et al 2019 (Nature): 
#17# - tx_rosenthal_RNA_danaher: RNA-seq Danaher et al immune signatures for 142 tracerx regions with histology/deep learning-based immune class.
#18# - tx_rosenthal_RNA_clusters: the original RNA-seq clusters from Rosetnhal et al 2019. Provided seprately for S5b. 

#############

#### ###  #### ### #### ###  #### ###
#Figures
#### ###  #### ### #### ###  #### ###

#1g-h and S2 g: biological-AI correlations  
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

#2a: regions specific immune variability in tx
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

#and the heatmap
txS8=tx
txS8$Histology[txS8$Histology=="Other"] <- "Z"
txS8 = txS8[order(txS8$lymphocytes_per_min, decreasing = F),]
txS8 = txS8[order(txS8$Histology, decreasing = F),]
txS8$DFS_yes[txS8$DFS_censor_variable==1] <- txS8$DFS_time_days
txS8$DFS_no[txS8$DFS_censor_variable==0] <- txS8$DFS_time_days

Heatmap(txS8$cd8_per, name = "CD8", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = colorRamp2(c(0, 17), c("#F2E2E2", "#8E0000")))+
  Heatmap(txS8$cd4_per, name = "CD4", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = colorRamp2(c(0, 28), c("#E3F0E2", "#067C00")))+
  Heatmap(txS8$foxp3_per, name = "FOXP3", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = colorRamp2(c(0, 21), c("#F9E2F2", "#CC008E")))
#############

#2c-d: deep learning vs Danaher immune signatures from Rosenthal et al 2019
#############
r = tx_rosenthal_RNA_danaher
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

#below for the stats. sign. for above merged panel 
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

#p-value correction to be annotated on the figure
p = c(5.8e-5, 1.5e-6, 9.4e-9, 3.4e-5, 0.00014, 0.018, 1.4e-5, 
      0.007, 0.0029, 6.6e-6, 4e-4, 2.9e-6, 3.3e-6, 9.1e-7)
p.adjust(p, method = "BH")
#grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, ncol=7)

#############

#3a-b and S5a: phenotype-based pairwise genomic dist and phy. analysis
#############
#mutations summary from Jamal-Hanjani, NEJM 2017
#then pairwise dist. for subclonal mutations for every hot-hot, cold-cold pair

comList <- list( c("hot_hot", "cold_cold"))
ggboxplot(tx_genomic_pairdist[tx_genomic_pairdist$Histology=="LUAD",], x = "immuneClass2", y = "geneticDist", 
              color = "immuneClass2", palette = c("blue", "red2"), 
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))

ggboxplot(tx_genomic_pairdist[tx_genomic_pairdist$Histology=="LUSC",], x = "immuneClass2", y = "geneticDist", 
              color = "immuneClass2", palette = c("blue", "red2"), 
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))


#distance using phy. tree for pairs in LUAD tumors
#this is checked manually from the file: all.pyclone_trees.20161220.RData, using these data: 
#trees$LTX0050$manual_tree
#trees$LTX0050$manual_mean_pyclone_ccf
#trees$LTX0050$manual_edgelength

comList <- list( c("hot_hot", "cold_cold"))
ggboxplot(txluad_domClon_pairdist, x = "immuneClass2", y = "phyDist2", color = "immuneClass2", palette = c("red2", "blue"), ylab = "Dist. of dominant clones to the MRCA (phylogenetics)",
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())

#reported in the caption, phyDist is when the dominant clone closest 
#to the most recent common ancestor of each tree was considered
ggboxplot(txluad_domClon_pairdist, x = "immuneClass2", y = "phyDist", color = "immuneClass2", palette = c("red2", "blue"), ylab = "Dist. of dominant clones to the MRCA (phylogenetics)",
          add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())

#not in the ms but just to confirm method a (from MJH Supp mutations data) and method b (Nicky's trees) are the same!
ggboxplot(txluad_domClon_pairdist, x = "immuneClass2", y = "geneticDist", color = "immuneClass2", palette = c("red2", "blue"), ylab = "Genomic distance",
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())

#############

#3d,f and S5c, d, e, f, g, h: survival analysis using the number of cold regions in TRACERx
#############

#many tests in this section, the ones are in the final paper are clearly labeled to the panel
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

#Fig 3d#
#KM
fit <- survfit(Surv(DFS_time_days, DFS_censor_variable)~rCold<=1, data = txH2LUADLUSC) 
ggsurvplot(fit, data = txH2LUADLUSC, conf.int = FALSE,
           pval = TRUE, pval.size = 5, pval.coord = c(0.2, 0.1),
           linetype = "solid",
           #surv.median.line = "hv",
           legend = "none", legend.title = "", title = "LUAD+LUSC",
           legend.labs = c("High", "Low"),
           surv.plot.height = 0.7, palette = c("#81391b", "#ef5123"),
           risk.table = TRUE,
           risk.table.col = "black", break.time.by = 200,
           tables.height = 0.25,
           tables.theme = theme_cleantable(),
           tables.y.text = TRUE, risk.table.title = "Number at Risk",
           tables.x.text = "", xlim = c(0, 1400),
           xlab = "Days to Death or Recurrence", ylab = "Disease-free Survival")
#Fig 3f#
#including total number of regions
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~coldP3+age+sex+pack_years_calculated+Histology+pathologyTNM2+Adjuvant.therapy+ClonalNeoUQ+nTotalRegions, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)

#Fig S5c#
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~coldP3+age+sex+pack_years_calculated+Histology+pathologyTNM2+Adjuvant.therapy+ClonalNeoUQ, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)
#also "alluded to" in text, univariate as cont.:
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~rCold, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)
#or at the median again:
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~coldP3, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)

#as continuous, keeping it simple with main clin data
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~rCold+age+sex+pack_years_calculated+Histology+pathologyTNM2+Adjuvant.therapy, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)

#Fig S5f#

#as continuous, testing independce against other scores
#correcting for Histology (by adding +Histology) in each fit doesn't make a big difference
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~rCold+subclonalCNA, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~rCold+mean_ASCAT.purity, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~rCold+tumour_per_mean, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)

#Fig S5g#
#the original test including other histology patients (only 6 more), cant do this with ClonalNeo UQ 
fit <- survfit(Surv(DFS_time_days, DFS_censor_variable)~rCold<=1, data = tx) 
ggsurvplot(fit, data = tx, conf.int = FALSE,
           pval = TRUE, pval.size = 5, pval.coord = c(0.2, 0.1),
           linetype = "solid",
           #surv.median.line = "hv",
           legend = "none", legend.title = "", title = "ALL",
           legend.labs = c("High", "Low"),
           surv.plot.height = 0.7, palette = c("#81391b", "#ef5123"),
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

#Fig S5e#
#lastly HR comparison panel
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

txH2LUADLUSC$cd8foxp3_Ratio = txH2LUADLUSC$cd8_per/txH2LUADLUSC$foxp3_per
txH2LUADLUSC$cd4foxp3_Ratio = txH2LUADLUSC$cd4_per/txH2LUADLUSC$foxp3_per
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
#adjust p-vals
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

#S5b: survival analysis using the number of immune cold region but with the RNA-seq-based classification (Rosenthal et al)
#############
#PS we don't use the overall patient level immune classification in Rosenthal et al 2019
#because that involves "salvaged" immune status from pathology TILs
#this is strict to N cold regions, using RNA clusters, not histology 

#why a seprate RNA df? because there might be some regions in "regional" that have RNA, but no histology
#we want the max number of tumors here which is 64 according to Rosenthal et al 2019

r = tx_rosenthal_RNA_clusters[, c(1, 2)] %>% 
  group_by(PublicationID, orig_immune_cluster) %>%
  tally()
txRNA = dcast(r, PublicationID ~ orig_immune_cluster)
txRNA[is.na(txRNA)] <- 0
names(txRNA)[2] <- "rRNA_high"
names(txRNA)[3] <- "rRNA_low"

nat = merge(txRNA, tx_surv15May2018)
median(nat$rRNA_low)

nat$coldRNAp3[nat$rRNA_low>1] <- ">1"
nat$coldRNAp3[nat$rRNA_low<=1] <- "<=1"
nat$coldRNAp3 = as.factor(as.character(nat$coldRNAp3))

fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~rRNA_low, data = nat)
u1 = ggforest(fit, data=nat)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~coldRNAp3, data = nat)
u3 = ggforest(fit, data=nat)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~rRNA_high, data = nat)
u2 = ggforest(fit, data=nat)

grid.arrange(u1, u3, u2, ncol=1)
#############

#3e,g and S5c, d, e: survival using the number of cold regions in LATTICe-A (validation cohort)
#############
tmp <- lattice
tmp$coldP3 = as.factor(as.character(tmp$coldP3))
tmp$gender = as.factor(as.character(tmp$gender))
tmp$Adjuvant.therapy = as.factor(as.character(tmp$Adjuvant.therapy))
tmp$stage = as.factor(as.character(tmp$stage))
levels(tmp$stage)[levels(tmp$stage)=="IIIA"] <- "III"
levels(tmp$stage)[levels(tmp$stage)=="IIIB"] <- "III"
tmplatt=tmp
tmplatt = tmplatt[!is.na(tmplatt$age_surgery) & !is.na(tmplatt$gender) &
                    !is.na(tmplatt$pack_years) & !is.na(tmplatt$stage) &
                    !is.na(tmplatt$Adjuvant.therapy) & !is.na(tmplatt$time_to_rfs) &
                    !is.na(tmplatt$rfs_status),]
#a simple test with the r cold on its own
tmp_model <- coxph(Surv(time_to_rfs, rfs_status == '1') ~ coldP3 + age_surgery + gender + pack_years + stage + Adjuvant.therapy, data=tmplatt)
ggforest(tmp_model, data = NULL, main = "", fontsize = 0.7, refLabel = "", noDigits = 2)

#Fig 3e#
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
surv_pvalue(fit)

#Fig 3g#
tmp_model <- coxph(Surv(time_to_rfs, rfs_status == '1') ~ coldP3 + age_surgery + gender + pack_years + stage + Adjuvant.therapy + nTotalRegions, data=tmplatt)
ggforest(tmp_model, data = NULL, main = "", fontsize = 0.7, refLabel = "", noDigits = 2)

#Fig 3g# - caption info: n=827
tmp_model <- coxph(Surv(time_to_rfs, rfs_status == '1') ~ coldP3 + age_surgery + gender + pack_years + stage + Adjuvant.therapy + nTotalRegions, data=tmp)
ggforest(tmp_model, data = NULL, main = "", fontsize = 0.7, refLabel = "", noDigits = 2)

#Fig S5c#
tmp_model <- coxph(Surv(time_to_rfs, rfs_status == '1') ~ coldP3 + age_surgery + gender + pack_years + stage + Adjuvant.therapy , data=tmplatt)
ggforest(tmp_model, data = NULL, main = "", fontsize = 0.7, refLabel = "", noDigits = 2)


#Fig S5d# with tumor size instead of total regions
tmp2 = tmplatt[!is.na(tmplatt$largest_tumour_size),]
tmp_model <- coxph(Surv(time_to_rfs, rfs_status == '1') ~ coldP3 + age_surgery + gender + pack_years + stage + Adjuvant.therapy + largest_tumour_size, data=tmp2)
ggforest(tmp_model, data = NULL, main = "", fontsize = 0.7, refLabel = "", noDigits = 2)

#Fig S5d# caption info: n=815
tmp3 = tmp
tmp3 = tmp[!is.na(tmp$age_surgery) & !is.na(tmp$gender) & !is.na(tmp$largest_tumour_size) &
                    !is.na(tmp$stage) &
                   !is.na(tmp$Adjuvant.therapy) & !is.na(tmp$time_to_rfs) &
                   !is.na(tmp$rfs_status),] 
tmp_model <- coxph(Surv(time_to_rfs, rfs_status == '1') ~ coldP3 + age_surgery + gender + pack_years + stage + Adjuvant.therapy + largest_tumour_size, data=tmp3)
ggforest(tmp_model, data = NULL, main = "", fontsize = 0.7, refLabel = "", noDigits = 2)


#Fig S5e#
#head to head, image-immune features: 
rU2 = lt %>%
  group_by(ACA_ID) %>%
  summarise(lymphocytes_per_mean = mean(lymphocytes_per),
            lymphocytes_per_min = min(lymphocytes_per), 
            lymphocytes_per_sd = sd(lymphocytes_per))
lattice_lym = NULL
lattice_lym = merge(lattice, rU2)
fit <- coxph(Surv(time_to_rfs, rfs_status)~cold, data = lattice_lym)
ggforest(fit, data=lattice_lym)
summary(fit)
fit <- coxph(Surv(time_to_rfs, rfs_status)~coldP3, data = lattice_lym)
ggforest(fit, data=lattice_lym)
summary(fit)
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

summaryPs <- data.frame(
  Measure = c("9N cold regions", "99>1 cold regions", "99 N cold regions/total","8N hot regions",
              "6Regional lymphocyte% (min)", "5Regional lymphocyte% (mean)", "55Regional lymphocyte% (SD)"),
  PVal = c(4.83e-07, 2.18e-05, 0.54, 0.233,
           0.459, 0.647, 0.732),
  lower95 = c(1.1, 1.2, 0.97, 0.98, 
              0.97, 0.97, 0.97),
  upper95 = c(1.2, 1.7, 1.2, 1.1, 
              1, 1, 1), 
  hr = c(1.1, 1.5, 0.99, 1, 
         0.99, 0.99, 1))
summaryPs$PVal_log10 <- -log10(summaryPs$PVal)
summaryPs$PvalAdjust <- p.adjust(summaryPs$PVal, method = "BH")
summaryPs$PvalAdjust_log10 <- -log10(summaryPs$PVal)

#change HR for the number of cold regions/total number of regions in LATTICe-A is adjusted from
#0.88[0.57-1.3] to 0.95[0.95-1.2].  
#not adjusted
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
#adjusted p-vals
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

#S6: validation of immune phenotype classification over 4 new schemes 
#############

#below numbers derived simply from:
reg = regional[, c("PublicationID", "region", "lymphocytes_per")]

median(reg$lymphocytes_per)+sd(reg$lymphocytes_per)/4 # or /2, /3, /6
median(reg$lymphocytes_per)-sd(reg$lymphocytes_per)/4 # or /2, /3, /6

#schemes visualisation 
ggdensity(regional, x = "lymphocytes_per", 
          add = "median", rug = T)+ 
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
txH2LUADLUSC = txH2LUADLUSC[, -c(14:16)]
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

reg$imC[reg$lymphocytes_per >=median(reg$lymphocytes_per)] <- "hot"
reg$imC[reg$lymphocytes_per <median(reg$lymphocytes_per)] <- "cold"

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
tx0 = merge(txH2LUADLUSC, regS0, by = "PublicationID")

regS2 <- reg %>% group_by(PublicationID, imC2) %>%
  dplyr::summarise(count=n())
regS2 = dcast(regS2, PublicationID ~ imC2)
names(regS2) <- c("PublicationID", "rCold", "rHot", "rIntermediate")
regS2[is.na(regS2)] <- 0
tx2 = merge(txH2LUADLUSC, regS2, by = "PublicationID")

regS3 <- reg %>% group_by(PublicationID, imC3) %>%
  dplyr::summarise(count=n())
regS3 = dcast(regS3, PublicationID ~ imC3)
names(regS3) <- c("PublicationID", "rCold", "rHot", "rIntermediate")
regS3[is.na(regS3)] <- 0
tx3 = merge(txH2LUADLUSC, regS3, by = "PublicationID")

regS6 <- reg %>% group_by(PublicationID, imC6) %>%
  dplyr::summarise(count=n())
regS6 = dcast(regS6, PublicationID ~ imC6)
names(regS6) <- c("PublicationID", "rCold", "rHot", "rIntermediate")
regS6[is.na(regS6)] <- 0
tx6 = merge(txH2LUADLUSC, regS6, by = "PublicationID")

tx0$coldP[tx0$rCold>median(tx0$rCold)] <- "med.High"
tx0$coldP[tx0$rCold<=median(tx0$rCold)] <- "Low"

tx2$coldP[tx2$rCold>median(tx2$rCold)] <- "med.High"
tx2$coldP[tx2$rCold<=median(tx2$rCold)] <- "Low"

tx3$coldP[tx3$rCold>median(tx3$rCold)] <- "med.High"
tx3$coldP[tx3$rCold<=median(tx3$rCold)] <- "Low"

tx6$coldP[tx6$rCold>median(tx6$rCold)] <- "med.High"
tx6$coldP[tx6$rCold<=median(tx6$rCold)] <- "Low"

tx0$coldP = as.factor(as.character(tx0$coldP))
tx2$coldP = as.factor(as.character(tx2$coldP))
tx3$coldP = as.factor(as.character(tx3$coldP))
tx6$coldP = as.factor(as.character(tx6$coldP))

fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~coldP+age+sex+pack_years_calculated+Histology+pathologyTNM2+Adjuvant.therapy+ClonalNeoUQ, data = tx0)
ggforest(fit, data=tx0, main = "Threshold: No intermediate zone")
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

reg2 = merge(reg, regional[, c("PublicationID", "region", "cd8.score.danaher"), ], by = c("PublicationID", "region"), all.x = T)
ggboxplot(reg2[! reg2$imC=="intermediate" ,], x = "imC", y = "cd8.score.danaher", color = "imC", 
          add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))
ggboxplot(reg2[! reg2$imC2=="intermediate" ,], x = "imC2", y = "cd8.score.danaher", color = "imC2", 
          add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))
ggboxplot(reg2[! reg2$imC3=="intermediate" ,], x = "imC3", y = "cd8.score.danaher", color = "imC3", 
          add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))
ggboxplot(reg2[! reg2$imC6=="intermediate" ,], x = "imC6", y = "cd8.score.danaher", color = "imC6", 
          add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))

##NOTE
#using R ver 6.3.3 the p-value for the SD/2 vs cd8.danaher was 1.5e-06 instead 3.7e-06
##

#now using the the same new 4 schemes, we test the pairwise genetic dist. 
#we changed the immune classification 4 times, will the gentics dist in LUAD remain sign.? 

#here we've the 4 classifications/schemes ready in a df: regional_immune_schemes
#just done exactly like the above code for survival tests
#imC: no intermediate zone
#imC2: SD/2 for the intermediate zone
#imC3: SD/3 for the intermediate zone
#imC6: SD/6 for the intermediate zone
#and all possible pairwise genetic dist. calculated in: tx_genomic_pairdist_extended
#2 dfs together we can reproduce S6 genetic dist. tests

#this immune_scheme_sub
#no immune class, we will add the 4 schemes in manually

tx_genomic_pairdist_extended$immune_scheme_sub <- NA

#start
for (i in 1:nrow(tx_genomic_pairdist_extended)){
  x=substring(tx_genomic_pairdist_extended$rsCRUK[[i]], 1,11)
  y=paste0(substring(tx_genomic_pairdist_extended$rsCRUK[[i]], 1,9), substring(tx_genomic_pairdist_extended$rsCRUK[[i]], 12,13))
  
  tx_genomic_pairdist_extended$immune_scheme_sub[[i]] <- paste0(regional_immune_schemes$imC[regional_immune_schemes$regionCRUK %in% x], 
                                                                "_", regional_immune_schemes$imC[regional_immune_schemes$regionCRUK %in% y])
} ###change here imC, imC2, imC3, imC6###
tx_genomic_pairdist_extended_ch = tx_genomic_pairdist_extended[tx_genomic_pairdist_extended$immune_scheme_sub=="hot_hot" | tx_genomic_pairdist_extended$immune_scheme_sub=="cold_cold" ,]
tx_genomic_pairdist_extended_ch$PublicationID = substring(tx_genomic_pairdist_extended_ch$rsCRUK, 1,8)
tx_genomic_pairdist_extended_ch <- tx_genomic_pairdist_extended_ch[order(tx_genomic_pairdist_extended_ch$immune_scheme_sub),]
comList <- list( c("hot_hot", "cold_cold"))
p1<-ggboxplot(tx_genomic_pairdist_extended_ch[tx_genomic_pairdist_extended_ch$histology_group=="Adenocarcinoma",], x = "immune_scheme_sub", y = "geneticDist", 
              color = "immune_scheme_sub", palette = c("blue", "red2"), title = "LUAD",
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))

p2<-ggboxplot(tx_genomic_pairdist_extended_ch[tx_genomic_pairdist_extended_ch$histology_group=="Squamous cell carcinoma",], x = "immune_scheme_sub", y = "geneticDist", 
              color = "immune_scheme_sub", palette = c("blue", "red2"), title = "LUSC",
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))

grid.arrange(p1,p2, nrow=1)

#check n
#table(tx_genomic_pairdist_extended_ch$histology_group, tx_genomic_pairdist_extended_ch$immune_scheme_sub)

tx_genomic_pairdist_extended_ch = NULL
### end of figure

#panel d: 
#distribution of lym% in tx and LATTICe-A for histology subtypes
ggdensity(regional, x = "lymphocytes_per", color = "Histology", palette = c("#98AED6", "#682A45", "black"),
          add = "mean", rug = TRUE) +
  theme(legend.position="none")
ggdensity(lt, x = "lymphocytes_per", color = "#98AED6", 
          add = "mean", rug = TRUE) +
  theme(legend.position="none")

comList <- list( c("LUAD", "LUSC"))
ggboxplot(regional, x = "Histology", y = "lymphocytes_per", color = "Histology", title = "TRACERx (n=275)", 
          add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+
  scale_color_manual(values=c("#98AED6", "#682A45", "black"))


#############

#4b-c and S7: cancer-stromal cell-cell fractal dimension analysis
#############

#FD vs immune class in tx
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

#effect size reported in main text
r = regional[! regional$immuneClass_2=="intermediate",]
r$imc = as.factor(as.character(r$immuneClass_2))
cohen.d(r$fd, r$imc)
cohen.d(r$fibroblasts_per, r$imc)

#FD vs LOHHLA - at region level
r1 = regional[!is.na(regional$LOHHLAloss.A) | !is.na(regional$LOHHLAloss.B) |
                      !is.na(regional$LOHHLAloss.C),]
r1$LOHHLAloss_any = NULL
r1$LOHHLAloss_any[r1$LOHHLAloss.A == T | r1$LOHHLAloss.B == T 
                        | r1$LOHHLAloss.C == T ] <-"Loss"
r1$LOHHLAloss_any[is.na(r1$LOHHLAloss_any)] <-"No loss"

comList <- list( c("Loss", "No loss"))
p1 = ggboxplot(r1[r1$Histology=="LUAD",], 
               x = "LOHHLAloss_any", xlab = "LOH: any HLA type", 
               y = "fd", ylab = "Fractal dimension",
               add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+ 
  theme(text = element_text(size=18))

r = regional[!is.na(regional$LOHHLAloss.A), ]
r$LOHHLAloss.A_L[r$LOHHLAloss.A == T] <-"Loss"
r$LOHHLAloss.A_L[r$LOHHLAloss.A == F] <-"No loss"
p2 = ggboxplot(r[r$Histology=="LUAD",], 
               x = "LOHHLAloss.A_L", xlab = "LOH: HLA type A", 
               y = "fd", ylab = "Fractal dimension",
               add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+ 
  theme(text = element_text(size=18))

r = regional[!is.na(regional$LOHHLAloss.B), ]
r$LOHHLAloss.B_L[r$LOHHLAloss.B == T] <-"Loss"
r$LOHHLAloss.B_L[r$LOHHLAloss.B == F] <-"No loss"
p3 = ggboxplot(r[r$Histology=="LUAD",], 
               x = "LOHHLAloss.B_L", xlab = "LOH: HLA type B", 
               y = "fd", ylab = "Fractal dimension",
               add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+ 
  theme(text = element_text(size=18))

r = regional[!is.na(regional$LOHHLAloss.C), ]
r$LOHHLAloss.C_L[r$LOHHLAloss.C == T] <-"Loss"
r$LOHHLAloss.C_L[r$LOHHLAloss.C == F] <-"No loss"
p4 = ggboxplot(r[r$Histology=="LUAD",], 
               x = "LOHHLAloss.C_L", xlab = "LOH: HLA type C", 
               y = "fd", ylab = "Fractal dimension",
               add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+ 
  theme(text = element_text(size=18))

p5 = ggboxplot(r1[r1$Histology=="LUSC",], 
               x = "LOHHLAloss_any", xlab = "LOH: any HLA type", 
               y = "fd", ylab = "Fractal dimension",
               add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+ 
  theme(text = element_text(size=18))

grid.arrange(p1, p2, p3, p4, p5, ncol=5)

#p-value correction as shown in the final figure
#correction LUAD region level analysis
p=c(0.0043, 0.074, 0.0019, 0.0029)
p.adjust(p, method = "BH")


#now FD vs LOHHLA at tumor level, using fd max to represent a patient
#(these are just the boxplots, p-vals are computed using multivariate model then corrected, see below)
txf = tx[!is.na(tx$any.HLA.loss), ]
comList <- list( c("Loss", "No loss"))
ggboxplot(txf[txf$Histology=="LUSC",], x = "any.HLA.loss", y = "fd_max", 
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+ 
  theme(text = element_text(size=18))
#no sign. in LUSC
fit <- glm(fd_max ~ factor(any.HLA.loss) + ClonalNeo  , data=txf[txf$Histology=="LUSC",])
summary(fit)

#LUAD
#again only vizualisation 
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
fit <- glm(fd_max ~ factor(any.HLA.loss) + ClonalNeo  , data=txf[txf$Histology=="LUAD",])
summary(fit)
#repeat the above for any.HLA.A.loss, any.HLA.B.loss, any.HLA.C.loss,
#p-value correction as shown in the final figure
p=c(0.00249, 0.00948, 0.0968, 0.0331)
p.adjust(p, method = "BH")


#FD distribution per histology 
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


#distribution of the min dist between a stromal cell to a cancer cell
#each df is a tumor region
#each row is a cell, we take the mean of each region to plot the distribution
stromal_tumor_cellAvgDist_summary = data.frame(t(sapply(stromal_tumor_cellAvgDist,colMeans)))

#average distance from a stromal cell to nearest cancer cell (each row is now a region)
ggplot(stromal_tumor_cellAvgDist_summary, aes(min_dist_f_to_t)) +
  geom_density(color="darkblue", fill="pink")+xlim(0,300) + theme_classic()

#############

#4e,g and S8: ATL/stromal correlation with ClonalNeo
#############
diagnosticLUAD = diagnostic[diagnostic$histology_group=="Adenocarcinoma",]
diagnosticLUAD$ClonalNeoMedian <- NULL
diagnosticLUAD$ClonalNeoMedian[diagnosticLUAD$ClonalNeo>=quantile(diagnosticLUAD$ClonalNeo)[3]] <-">=Median"
diagnosticLUAD$ClonalNeoMedian[diagnosticLUAD$ClonalNeo<quantile(diagnosticLUAD$ClonalNeo)[3]] <-"<Median"
comList <- list( c(">=Median", "<Median"))
ggboxplot(diagnosticLUAD, x = "ClonalNeoMedian", y = "ATL_fibroRatio", 
          xlab = "Clonal neoantigens", ylab = "Adjacent-tumor lymphocytes/stroma", 
          color = "ClonalNeoMedian", palette = c("blue", "red2"), #title = "LUAD",
          add = "dotplot", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))
##NOTE
#using R ver 6.3.3 the p-value was 0.0068 instead of 0.0074
##

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

ggboxplot(reg_tils_immune, x = "TILClass", y = "foxp3", color = "TILClass", palette = c("blue", "#8291F7", "black"),
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox", paired = T)+
  theme(legend.position="") + theme(text = element_text(size=18))


#immune features vs NeoAgs with p-val correction
toCorrLUAD <- select(tx[tx$Histology=="LUAD",], lymphocytes_per, ATL_fibroRatio, ITL_tumorRatio, cd4_per, cd8_per, cd8foxp3Ratio,
                     ClonalNeo, SubclonalNeo, Regions.LOHHLAloss.A, Regions.LOHHLAloss.B, Regions.LOHHLAloss.C)

#no confirmed HLA status (i.e. NAs) will be ignored to be safe

p.rLUAD <- cor_pmat(toCorrLUAD, method=c("spearman"), use="complete.obs")
p.rLUAD <- matrix(p.adjust(p.rLUAD, method = "BH", n = length(p.rLUAD)), nrow = 11, ncol = 11) ##
corRLUAD <- round(cor(toCorrLUAD, method = c("spearman"), use="complete.obs"),2)
ggcorrplotModified(corRLUAD[7:11, 1:6], lab = T, sig.level=0.05, insig ="blank", lab.notsig="", 
                   colors = c("#4DBBD5B2","grey98","#DC0000B2"), legend.title = "Spearman Corr",
                   p.mat = p.rLUAD[7:11, 1:6])

#S8 HE-IHC reg. eval. 
ggdensity(reg_eval, x = "distance", color = "purple", xlab = "HE-IHC registration distance (um)",
          add = "mean", rug = TRUE)
#############

#S2 a-d: correlations of deep learning scores vs DNA, RNA, pathology estimations in tx
#############
#diagnostic lym% and ATL vs IHC cd8 and stromal TILs

#the green annotation is computed using the overall correlation for every plot
#replacing the black (Other histology) in all corr plots

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

#S2 e: correlations of deep learning scores in LATTICe-A
#############
ggscatter(lattice_80, x = "pathology_LymEstimate", y = "lymphocytes_per", add = "reg.line", 
          xlab = "Lymphocyte cell fraction (scored by pathologist)", ylab = "Lymphocyte%", font.label = 4,
          color = "#98AED6",shape = 19,
          conf.int = TRUE, size = 2,
          add.params = list(color = "grey50", fill = "azure3"), 
          cor.coef = TRUE, cor.method = "spearman") 
ggscatter(lattice_80, x = "Stromal.TIL", y = "ATL_fibroRatio", add = "reg.line", 
          xlab = "Stromal TILs (scored by pathologist)", ylab = "ATL/Stroma", font.label = 4,
          color = "#98AED6",shape = 19,
          conf.int = TRUE, size = 2,
          add.params = list(color = "grey50", fill = "azure3"), 
          cor.coef = TRUE, cor.method = "spearman") 

#############

#S3: regional lym% distribution against stage in tx and LATTICe-A
#############
p1 = ggdensity(tx[! tx$Histology=="Other",], x = "lymphocytes_per_sd", color = "Histology", palette = c("#98AED6", "#682A45"), xlab = "Lymphocyte% (regional SD)",
          add = "median", rug = TRUE)+theme(legend.position = "none")
p2 = ggdensity(tx[! tx$Histology=="Other",], x = "lymphocytes_per_sd", color = "pathologyTNM2", xlab = "Lymphocyte% (regional SD)",
          add = "median", rug = TRUE)+
  scale_color_manual(values = c("#f3776e", "#b59f31", "#32b24b", "#1dbdc2", "#9268ad"))+ theme(legend.position = "none")
p3 = ggdensity(tx[! tx$Histology=="Other",], x = "lymphocytes_per_mean", color = "pathologyTNM2", xlab = "Lymphocyte% (regional mean)",
          add = "median", rug = TRUE)+
  scale_color_manual(values = c("#f3776e", "#b59f31", "#32b24b", "#1dbdc2", "#9268ad"))+ theme(legend.position = "none")

lattice$stage2 = lattice$stage
lattice$stage2[lattice$stage =="IIIA"] <- "III"
lattice$stage2[lattice$stage =="IIIB"] <- "III"

rU2 = lt %>% select(ACA_ID, lymphocytes_per) %>%
  dplyr::group_by(ACA_ID) %>%
  dplyr::summarise(lymphocytes_per_sd = sd(lymphocytes_per),
                   lymphocytes_per_mean = mean(lymphocytes_per),
                   lymphocytes_per_min = min(lymphocytes_per))
lattice_lym = merge(lattice, rU2, all.x = T)
lym_clean = lattice_lym[! is.na(lattice_lym$stage) ,]

p4 = ggdensity(lym_clean, x = "lymphocytes_per_sd", color = "#98AED6", xlab = "Lymphocyte% (regional SD)",
          add = "median", rug = TRUE)+theme(legend.position = "none")
p5 = ggdensity(lym_clean, x = "lymphocytes_per_sd", color = "stage2", xlab = "Lymphocyte% (regional SD)",
          add = "median", rug = TRUE)+
  scale_color_manual(values = c("#f3776e", "#b59f31", "#32b24b", "#1dbdc2", "#9268ad"))+ theme(legend.position = "none")
p6 = ggdensity(lym_clean, x = "lymphocytes_per_mean", color = "stage2", xlab = "Lymphocyte% (regional mean)",
          add = "median", rug = TRUE)+
  scale_color_manual(values = c("#f3776e", "#b59f31", "#32b24b", "#1dbdc2", "#9268ad"))+ theme(legend.position = "none")

grid.arrange(p1, p2, p3, p4, p5, p6, ncol=3)

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

ggstatsplot::ggbetweenstats(
  data = lattice_lym, type = "np",
  x = stage2, xlab = "Pathological stage",
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
ggstatsplot::ggbetweenstats(
  data = lattice_lym, type = "np",
  x = stage2, xlab = "Pathological stage",
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
#############

#S4: concordance of deep learning based immune classification with RNA/pathology TILs (from Rosenthal et al 2019)
#############

#against path TILs for all available regions
comList <- list( c("cold", "hot"))
ggboxplot(regional[! regional$immuneClass_2 == "intermediate" ,], x = "immuneClass_2", y = "salgado_TILs", 
          color = "immuneClass_2", palette = c("red2", "blue"), ylab = "Pathology TILs",
          add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())

#test overlap
r = regional[! is.na(regional$orig_immune_cluster) ,]
r = r[! r$immuneClass_2=="intermediate",]
fisher.test(table(r$immuneClass_2, r$orig_immune_cluster))
table(r$immuneClass_2, r$orig_immune_cluster)

#immune hotspot analysis
r = regional[! is.na(regional$orig_immune_cluster) ,]
r = r[! r$immuneClass_2=="intermediate",]
r$discHC = NULL
r$discHC[r$immuneClass_2=="hot" & r$orig_immune_cluster == "high"] = "TRUE"
r$discHC[r$immuneClass_2=="cold" & r$orig_immune_cluster == "low"] = "TRUE"
r$discHC[is.na(r$discHC)] = "FALSE"
comList <- list( c("TRUE", "FALSE"))
ggboxplot(r, x = "discHC", y = "fraction_immune_hotspot", color = "discHC", palette = c( "#97f2a4", "#f4abc2"), ylab = "Fraction of immune hotspot",
          add = "dotplot", border = "white")+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())

#comparisons
r = regional[! regional$immuneClass_2=="intermediate",]
comList <- list( c("hot", "cold"))
r <- r[order(r$immuneClass_2),] 
p1 = ggboxplot(r, x = "immuneClass_2", y = "lymphocytes_per", color = "immuneClass_2", ylab = "Lymphocyte%", 
              add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",comparisons = comList, method = "wilcox")+ 
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())
p2 = ggboxplot(r, x = "immuneClass_2", y = "ITL_tumorRatio", color = "immuneClass_2", ylab = "Intra-tumor lym/tumor",
              add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",comparisons = comList, method = "wilcox")+ 
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())
p3 = ggboxplot(r, x = "immuneClass_2", y = "ATL_fibroRatio", color = "immuneClass_2", ylab = "Adjacent-tumor lym/stroma",
              add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",comparisons = comList, method = "wilcox")+ 
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())
p4 = ggboxplot(r, x = "immuneClass_2", y = "ASCAT.purity", color = "immuneClass_2",
              add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",comparisons = comList, method = "wilcox")+ 
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())

r = regional[! is.na(regional$orig_immune_cluster) ,]
comList <- list( c("low", "high"))
r <- r[order(r$orig_immune_cluster, decreasing = T),] 
p5 = ggboxplot(r, x = "orig_immune_cluster", y = "lymphocytes_per", color = "orig_immune_cluster", ylab = "Lymphocyte%",
          add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",comparisons = comList, method = "wilcox")+ 
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())
p6 = ggboxplot(r, x = "orig_immune_cluster", y = "ITL_tumorRatio", color = "orig_immune_cluster", ylab = "Intra-tumor lym/tumor",
          add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",comparisons = comList, method = "wilcox")+ 
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())
p7 = ggboxplot(r, x = "orig_immune_cluster", y = "ATL_fibroRatio", color = "orig_immune_cluster", ylab = "Adjacent-tumor lym/stroma",
          add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",comparisons = comList, method = "wilcox")+ 
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())
p8 = ggboxplot(r, x = "orig_immune_cluster", y = "ASCAT.purity", color = "orig_immune_cluster",
          add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",comparisons = comList, method = "wilcox")+ 
  theme(legend.position="") + theme(text = element_text(size=18))+theme(axis.title.x=element_blank())

#############

#S1: patient char. in tx, and tx-LATTICe-A demographic/cohort comparison
#############

#txclin_supp is used, which is tx~diagnostic (to plot patients without region specific sampels)
#kept seprately to avoi any confusions! 

#TRACERx 100 patient characteristics
Heatmap(txclin_supp$Histology, name = "Histology", width = unit(2.85, "mm"), 
        cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
        col = c("#98AED6", "#682A45", "black")) +
  Heatmap(txclin_supp$pathologyTNM, name = "Stage", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = c("yellow2", "chartreuse1", "cadetblue1", "burlywood4", "orchid", "firebrick"))+
  Heatmap(txclin_supp$pack_years_calculated, name = "Pack years", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = colorRamp2(c(0, 20, 40, 60, 80, 100, 120, 140), 
                           c("#568E9B", "#427884", "#386A75", "#395760", "#4F4F4F", "#3D3D3D", "#333333", "black")))+
  Heatmap(txclin_supp$DFS_censor_variable, name = "Event", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = c("lightblue", "midnightblue"))+
  Heatmap(txclin_supp$nTotalRegions, name = "N regions", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = c("#FADDA2", "#EEBF88", "#E3A371", "#D9885C", "#CF7049", "#C65738", "#BC3C29"))+
  Heatmap(txclin_supp$nRNARegions, name = "N RNA regions", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = c("#FCF9CE", "#FADDA2", "#EEBF88", "#E3A371", "#D9885C", "#CF7049", "#C65738", "gray45"))+
  Heatmap(txclin_supp$nTRegions, name = "N hist. regions", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = c("#FCF9CE", "#FADDA2", "#EEBF88", "#E3A371", "#D9885C", "#CF7049", "#C65738", "#BC3C29", "gray45"))+
  Heatmap(txclin_supp$rCold, name = "N cold regions", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = c("#71c7ec", "#1ebbd7", "#189ad3", "#107dac", "#005073", "#033547","gray45"))


#compare cohorts
p1 = ggdensity(tx, x = "age",
               add = "median", rug = TRUE,
               color = "sex", fill = "sex", xlab = "Age",
               palette = c("#F39B7FB2", "#4DBBD5B2"))
p2 = ggdensity(tx, x = "pack_years_calculated", xlab = "Pack years", 
               add = "median", rug = T)

stagetx = tx %>% select(PublicationID, pathologyTNM2) %>% group_by(pathologyTNM2) %>% tally()
stagetx$per = stagetx$n / sum(stagetx$n)

p3 =ggplot(stagetx, aes(x="", y=per, fill=pathologyTNM2)) + geom_bar(stat="identity", width=1)+ 
  coord_polar("y", start=0) + geom_text(aes(label = paste0(round(per*100), "%")), position = position_stack(vjust = 0.5))+ 
  scale_fill_manual(values=c("#f3776e", "#b59f31", "#32b24b", "#1dbdc2", "#9268ad"))+ 
  labs(x = NULL, y = NULL, fill = NULL)+ theme_classic() + theme(axis.line = element_blank(),
                                                                 axis.text = element_blank(),
                                                                 axis.ticks = element_blank(),
                                                                 plot.title = element_text(hjust = 0.5, color = "#666666"))

tmp2=lattice
tmp2$stage[tmp2$stage=="IIIA"] <- "III"
tmp2$stage[tmp2$stage=="IIIB"] <- "III"
p4 = ggdensity(tmp2, x = "age_surgery",
               add = "median", rug = TRUE,
               color = "gender", fill = "gender", xlab = "Age",
               palette = c("#F39B7FB2", "#4DBBD5B2"))
p5 = ggdensity(tmp2, x = "pack_years", xlab = "Pack years", 
               add = "median", rug = T)

stageLattice = tmp2[! is.na(tmp2$stage) ,] %>% select(ACA_ID, stage) %>% group_by(stage) %>% tally()
stageLattice$per = stageLattice$n / sum(stageLattice$n)

# Create a basic bar
p6 = ggplot(stageLattice, aes(x="", y=per, fill=stage)) + geom_bar(stat="identity", width=1)+ 
  coord_polar("y", start=0) + geom_text(aes(label = paste0(round(per*100), "%")), position = position_stack(vjust = 0.5))+ 
  scale_fill_manual(values=c("#f3776e", "#b59f31", "#32b24b", "#1dbdc2", "#9268ad"))+ 
  labs(x = NULL, y = NULL, fill = NULL)+ theme_classic() + theme(axis.line = element_blank(),
                                                                 axis.text = element_blank(),
                                                                 axis.ticks = element_blank(),
                                                                 plot.title = element_text(hjust = 0.5, color = "#666666"))

grid.arrange(p1, p2, p3, p4, p5, p6, nrow=2)


#############

#S9: summary stats heatmap
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

#KA 20200310