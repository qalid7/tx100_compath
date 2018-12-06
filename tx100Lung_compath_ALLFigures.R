#### ###  #### ###

##LungTx Compath##
##ALL stat analyses and figures##

#### ###  #### ###

setwd('X:/Dropbox (ICR)/yuanlab/Projects/lung/tracerx')
path = file.path("X:/Dropbox (ICR)/yuanlab/Projects/lung/tracerx/tracerx100/results/LungTx1Figures/")

source('./code/Kfun.R')
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

#LOAD - this is just an example, a similar data summary stats file is loaded serveal times 
load("./tracerx100/results/trait/trait5_08Oct2018.RData")
#mutations summary from Jamal-Hanjani, NEJM 2017
mut <- read.csv("./tracerx100/data/tracerx100tableS3Mutations.csv", stringsAsFactors = F)

#definitions:
#diagnostic: summary file for all diangostic immune scores for tx100 patients, 100 patients
#regional: summary file for all regional immune scores for tx100 patients, 275 tumor regions, 85 patients
#tx: summary regional+diangostic (85 patients)
#txH2: summary of tx for luad and lusc patients only (79 patients)
#genetic: majority of tx100 genetic and genmoic summary data from Jamal-Hanjani NEJM 2017, McGranahan Cell 2017 
#and Rosenthal co-submitted to Nature 2018 
#also sometimes we need a seprate df for LUAD, LUSC, etc

##Correlations with Spearman:
#all correlations are made for each histology subtype (luad, lusc or other) within the ggscatter function
#followed by a whole cohort corr (highlighted in green in all corr figures)

##Test of difference amnong groups:
#we use ggboxplot's and set comparison list to show exact p val within the figure 
#e.g. c("hot", "cold), then assign T test as the method
#sometimes where apropriate, we allow the paired option 

#most of the sections below re-load the needed data, it needs to be "melted" or "re-shaped" often
#hence, clear the vars and re-load


##fig S3, deep learning validation using pathology TILs, DNA purity, RNA cd8, IHC~HE
########################################################################################################################

diagnostic$Histology[diagnostic$histology_group== "Adenocarcinoma"] <- "LUAD"
diagnostic$Histology[diagnostic$histology_group== "Squamous cell carcinoma"] <- "LUSC"
diagnostic$Histology[is.na(diagnostic$Histology)] <- "Other"

#load pathology Stromal TILs by Roberto and co
diagnostic$ATL_fibroRatio = diagnostic$ATL/diagnostic$fibroblasts
pathTILs <- read.csv(file="X:/Dropbox (ICR)/TRACERx-lung/tracerx100/tx100.roberto.scores.csv", stringsAsFactors = F)
pathTILs$PATIENT_SCREENING_ID = sub("(.{3})(.*)", "\\10\\2", pathTILs$PATIENT_SCREENING_ID) #add 0 in LTX
diagnosticpathTILs <- merge(diagnostic, pathTILs, by = "PATIENT_SCREENING_ID", all.y = T)

p1 <- ggscatter(diagnosticpathTILs, x = "Stromal.TIL", y = "ATL_fibroRatio", show.legend.text = F, 
                xlab = "Stromal TILs", ylab = "ATL/Stroma", font.label = 4,
                shape = 19, size = 2,
                color = "Histology", palette = c("#98AED6", "#682A45", "black"),
                add.params = list(color = "grey50", fill = "azure3"))+ 
  stat_cor(aes(color = Histology), label.x = 15, method="spearman", show.legend = F)+
  theme(legend.position="none")
ggscatter(diagnosticpathTILs, x = "Stromal.TIL", y = "ATL_fibroRatio", cor.coef = T)

p2<-ggscatter(diagnostic, x = "lymphocytes_per", y = "cd8_per", show.legend.text = F, 
          xlab = "H&E Lymphocytes%", ylab = "IHC CD8%", font.label = 4,
          shape = 19, size = 2,
          color = "Histology", palette = c("#98AED6", "#682A45", "black"),
          add.params = list(color = "grey50", fill = "azure3"))+ 
  stat_cor(aes(color = Histology), label.x = 5, method="spearman", show.legend = F)+
  theme(legend.position="none")
ggscatter(diagnostic, x = "lymphocytes_per", y = "cd8_per", cor.coef = T)


#regional DNA purity 
regionalDNA = colNA(regional, "ASCAT.purity")
regionalDNA$Histology[regionalDNA$Histology== "Invasive adenocarcinoma"] <- "LUAD"
regionalDNA$Histology[regionalDNA$Histology== "Squamous cell carcinoma"] <- "LUSC"
regionalDNA$Histology[regionalDNA$Histology== "Large cell carcinoma"] <- "Other"
regionalDNA$Histology[regionalDNA$Histology== "Large Cell Neuroendocrine"] <- "Other"
regionalDNA$Histology[regionalDNA$Histology== "Carcinosarcoma"] <- "Other"
regionalDNA$Histology[regionalDNA$Histology== "Adenosquamous carcinoma"] <- "Other"

p3 <- ggscatter(regionalDNA, x = "ASCAT.purity", y = "tumour_per", show.legend.text = F, 
                xlab = "ASCAT Purity", ylab = "Tumour%", font.label = 4,
                shape = 19, size = 2,
                color = "Histology", palette = c("#98AED6", "#682A45", "black"),
                add.params = list(color = "grey50", fill = "azure3"))+ 
  stat_cor(aes(color = Histology), label.x = 0.75, method="spearman", show.legend = F)+
  theme(legend.position="none")
ggscatter(regionalDNA, x = "ASCAT.purity", y = "tumour_per", cor.coef = T)

p4 <- ggscatter(regionalDNA, x = "VAF.purity", y = "tumour_per", show.legend.text = F, 
                xlab = "VAF Purity", ylab = "Tumour%", font.label = 4,
                shape = 19, size = 2,
                color = "Histology", palette = c("#98AED6", "#682A45", "black"),
                add.params = list(color = "grey50", fill = "azure3"))+ 
  stat_cor(aes(color = Histology), label.x = 0.75, method="spearman", show.legend = F)+
  theme(legend.position="none")
ggscatter(regionalDNA, x = "VAF.purity", y = "tumour_per", cor.coef = T)


#Danaher cd8 signature from RNA
load("X:/Dropbox (ICR)/yuanlab/Projects/lung/tracerx/tracerx100/data/cd8_and_lohhla.RData")
colnames(cd8.scores)[which(names(cd8.scores) == "RegionID")] <- "file_name"
cd8.scores$file_name = gsub(":", "_T1", cd8.scores$file_name)
regional$Histology[regional$Histology== "Invasive adenocarcinoma"] <- "LUAD"
regional$Histology[regional$Histology== "Squamous cell carcinoma"] <- "LUSC"
regional$Histology[regional$Histology== "Large cell carcinoma"] <- "Other"
regional$Histology[regional$Histology== "Large Cell Neuroendocrine"] <- "Other"
regional$Histology[regional$Histology== "Carcinosarcoma"] <- "Other"
regional$Histology[regional$Histology== "Adenosquamous carcinoma"] <- "Other"
rCD8 = merge(cd8.scores, regional, all.x = T)
rCD8 = colNA(rCD8, "cd8.score.danaher")
rCD8 = colNA(rCD8, "lymphocytes_per")

p5 <- ggscatter(rCD8, x = "cd8.score.danaher", y = "lymphocytes_per", show.legend.text = F, 
                xlab = "CD8 score (Danaher et. al. signature)", ylab = "Lymphocytes%", font.label = 4,
                shape = 19, size = 2,
                color = "Histology", palette = c("#98AED6", "#682A45", "black"),
                add.params = list(color = "grey50", fill = "azure3"))+ 
  stat_cor(aes(color = Histology), label.x = 1, method="spearman", show.legend = F)+
  theme(legend.position="none")
ggscatter(rCD8, x = "cd8.score.danaher", y = "lymphocytes_per", cor.coef = T)

pdf(file=paste0(path, "/figS3.pdf"), width = 10 , height = 16)
grid.arrange(p1,p2,p3,p4,p5,  ncol=2)
dev.off()
########################################################################################################################

##fig S4, biological validation with TTF1 (independent cohort)
########################################################################################################################
IHCDir<-"R:/tracerx/HE_IHC_John/AEC3/TTF1/Misc/IHC-annotations-20180530/rect_labels-180718/rectlabels/3_HE.ndpi"
CELLS<-list.files(path = IHCDir,pattern=glob2rx("Da*.csv"))
CELLS <- CELLS[!CELLS %in% "Da171.csv"]
options(stringsAsFactors = FALSE)
TTF1_list<- vector("list", length(CELLS))
names(TTF1_list) <- CELLS
IHCDir2<-"./tracerx100/results/tx1_validation/csv/TTF1-20180530/detected_points/Manualp_3_TTF1.tif"
for (iCELL in CELLS){
  limit <- read.csv(paste0(IHCDir, '/',iCELL))
  Datmp2 <- read.csv(paste0(IHCDir2, '/',iCELL))
  print(iCELL)
  
  Datmp2 = Datmp2[Datmp2$V2 >= limit$xf1 & Datmp2$V2 <= limit$xf2 & 
                    Datmp2$V3 >= limit$yf1 & Datmp2$V3 <= limit$yf2 ,]
  TTF1_list[[iCELL]] <- nrow(Datmp2[Datmp2$V1=="None" ,])
}

IHCDir<-"R:/tracerx/HE_IHC_John/AEC3/TTF1/Misc/IHC-annotations-20180530/rect_labels-180718/rectlabels/3_HE.ndpi"
CELLS<-list.files(path = IHCDir,pattern=glob2rx("Da*.csv"))
CELLS <- CELLS[!CELLS %in% "Da171.csv"]
options(stringsAsFactors = FALSE)
HETTF1_list<- vector("list", length(CELLS))
names(HETTF1_list) <- CELLS
IHCDir2<-"./tracerx100/results/tx1_validation/csv/HETTF1-20180327/detected_points/3_HE.ndpi"
for (iCELL in CELLS){
  limit <- read.csv(paste0(IHCDir, '/',iCELL))
  Datmp2 <- read.csv(paste0(IHCDir2, '/',iCELL))
  print(iCELL)
  
  Datmp2 = Datmp2[Datmp2$V2 >= limit$xf1 & Datmp2$V2 <= limit$xf2 & 
                    Datmp2$V3 >= limit$yf1 & Datmp2$V3 <= limit$yf2 ,]
  HETTF1_list[[iCELL]] <- nrow(Datmp2[Datmp2$V1=="t" ,])
}

TTF1_list=data.frame(t(do.call(cbind.data.frame, TTF1_list)))
TTF1_list$Da = rownames(TTF1_list)
names(TTF1_list) = c("cellCount_ttf1", "tile")
HETTF1_list=data.frame(t(do.call(cbind.data.frame, HETTF1_list)))
HETTF1_list$Da = rownames(HETTF1_list)
names(HETTF1_list) = c("cellCount_tumour", "tile")
TTF1 = merge(TTF1_list, HETTF1_list)

tissTTF1 = read.csv("X:/Dropbox (ICR)/TheTemp/SMACD3TIFF1/Manualp_3_TTF1.tif.csv", stringsAsFactors = F)
tissHETTF1 = read.csv("X:/Dropbox (ICR)/TheTemp/SMACD3TIFF1/3_HE.ndpi.csv", stringsAsFactors = F)

tissTTF1$tile = paste0(tissTTF1$tile, ".csv")
tissHETTF1$tile = paste0(tissHETTF1$tile, ".csv")

names(tissTTF1)[2:4] = paste0(names(tissTTF1)[2:4], "_ihc")
names(tissHETTF1)[2:4] = paste0(names(tissHETTF1)[2:4], "_he")

TTF1 = merge(TTF1, tissTTF1)
TTF1 = merge(TTF1, tissHETTF1)

TTF1 = TTF1[TTF1$tissper_he>0.05 | TTF1$tissper_ihc>0.05 ,]

#normalise the tissue area by 16 (SS1 image is x16 lower), 0.3225*1000 is for 20x nanozoomer
TTF1$ttf1_norm = TTF1$cellCount_ttf1/TTF1$tissuearea_grayImage_ihc*16*0.3225*1000
TTF1$tumour_norm = TTF1$cellCount_tumour/TTF1$tissuearea_grayImage_he*16*0.3225*1000

pdf(file=paste0(path, "/countTissuareaCorrs_tiles_bw_updated.pdf"), width = 5, height = 5)
ggscatter(TTF1, x = "tumour_norm", y = "ttf1_norm", add = "reg.line", 
          color = "green",
          conf.int = TRUE, size = 2,
          add.params = list(color = "grey50", fill = "azure3"), 
          cor.coef = TRUE, cor.method = "spearman") 
dev.off()
########################################################################################################################

##fig 2, immune classification and survival for number of immune cold regions + fig S7
########################################################################################################################
#make the same dotplot for all 275 regions, 85 patients, this is used in figs 2 and 6
load("./tracerx100/results/trait/trait5_08Oct2018.RData")
tx$Histology[tx$histology_group == "Squamous cell carcinoma"] <- "Squamous-cell carcinoma"
tx$Histology[tx$histology_group == "Adenocarcinoma"] <- "Adenocarcinoma"
tx$Histology[tx$histology_group=="Pleomorphic carcinoma"] <- "Z"
tx$Histology[tx$histology_group=="Adenosquamous carcinoma"] <- "Z"
tx$Histology[tx$histology_group=="Large Cell Neuroendocrine Carcinoma"] <- "Z"
tx$Histology[tx$histology_group=="sarcomatoid carcinoma of pleomorphic type arising from adenocarcinoma"] <- "Z"
regional$Histology[regional$Histology == "Squamous cell carcinoma"] <- "Squamous-cell carcinoma"
regional$Histology[regional$Histology == "Invasive adenocarcinoma"] <- "Adenocarcinoma"
regional$Histology[regional$Histology=="Carcinosarcoma"] <- "Z"
regional$Histology[regional$Histology=="Adenosquamous carcinoma"] <- "Z"
regional$Histology[regional$Histology=="Large Cell Neuroendocrine"] <- "Z"
regional = regional[order(regional$lymphocytes_per, decreasing = F),]
regional = regional[order(regional$Histology, decreasing = F),]

ggdotchart(regional, x = "PublicationID", y = "lymphocytes_per",
           color = "immuneClass_2",                                # Color by groups.
           group = "PublicationID",
           sort.by.groups = TRUE,           
           palette = c("blue2", "red2", "#f4d5d5"), # Custom color palette
           #sorting = "descending",                       # Sort value in descending order
           #rotate = TRUE,                                # Rotate vertically
           #add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 0.8),
           dot.size = 3,                                 # Large dot size
           #label = min(regionalT$lymphocytes_per),
           #y.text.col = TRUE,                            # Color y text by groups
           ggtheme = theme_pubr())

#now the survival analysis using the number of cold regions
load("./tracerx100/results/trait/trait5_15May2018.RData")
diagnosticLUAD <- diagnostic[diagnostic$histology_group=="Adenocarcinoma",]
diagnosticLUSC <- diagnostic[diagnostic$histology_group=="Squamous cell carcinoma",]
geneticLUAD <- genetic[genetic$histology_group=="Adenocarcinoma",]
geneticLUSC <- genetic[genetic$histology_group=="Squamous cell carcinoma",]
diagnosticIHC = colNA(diagnostic, "cd8_per")
diagnosticIHC$cd8foxp3Ratio = diagnosticIHC$cd8_per/diagnosticIHC$foxp3_per
diagnosticIHCLUAD <- diagnosticIHC[diagnosticIHC$histology_group=="Adenocarcinoma",]
diagnosticIHCLUSC <- diagnosticIHC[diagnosticIHC$histology_group=="Squamous cell carcinoma",]
tx$coldP = as.factor(as.character(tx$coldP))
tx$coldP3 = as.factor(as.character(tx$coldP3))
tx$pathologyTNM = as.factor(as.character(tx$pathologyTNM))
txH2 = tx
txH2$Histology[txH2$histology_group=="Adenocarcinoma" ] <- "LUAD"
txH2$Histology[txH2$histology_group=="Squamous cell carcinoma" ] <- "LUSC"
txH2$Histology[is.na(txH2$Histology)] <- "Other"
txH2$Histology = as.factor(as.character(txH2$Histology))
txH2$sex = as.factor(as.character(txH2$sex))
txH2$pathologyTNM2=txH2$pathologyTNM
levels(txH2$pathologyTNM2)[levels(txH2$pathologyTNM2)=="IIIA"] <- "III"
levels(txH2$pathologyTNM2)[levels(txH2$pathologyTNM2)=="IIIB"] <- "III"
txH2$AdjuvantTherapy = txH2$adjuvant_treatment_given
txH2$AdjuvantTherapyC[txH2$adjuvant_treatment_given == ""]<- 1
txH2$AdjuvantTherapyC[txH2$adjuvant_treatment_given == "Platinum chemo"]<- 0
txH2$AdjuvantTherapyC[txH2$adjuvant_treatment_given == "Radiotherapy"]<- 0
txH2$AdjuvantTherapyC[txH2$adjuvant_treatment_given == "Platinum chemo/radiotherapy"]<- 0
txH2LUADLUSC <- txH2[ which(txH2$Histology=="LUAD" 
                            | txH2$Histology == "LUSC"), ]
txH2LUADLUSC$Histology <- factor(txH2LUADLUSC$Histology)

#include number of tumor regions/samples, this is from the DNA paper
num = readRDS("./tracerx100/data/number_regions.RDS")
num$REGTrialNo = rownames(num)
names(num) <- c("nTotalRegions", "nRNARegions", "PATIENT_SCREENING_ID")
num$PATIENT_SCREENING_ID = sub("(.{3})(.*)", "\\10\\2", num$PATIENT_SCREENING_ID)
txH2LUADLUSC = merge(txH2LUADLUSC, num, all.x = T)
txH2 = merge(txH2, num, all.x = T)

#with clonal neo but local quartile for the available 79 luad lusc patients: quantile(txH2LUADLUSC$ClonalNeo[txH2LUADLUSC$Histology=="LUSC"])
txH2LUADLUSC$ClonalNeoUQ[txH2LUADLUSC$Histology=="LUAD" & txH2LUADLUSC$ClonalNeo >= 189]<- "UQ.high"
txH2LUADLUSC$ClonalNeoUQ[txH2LUADLUSC$Histology=="LUAD" & txH2LUADLUSC$ClonalNeo < 189]<- "cUQ.low"
txH2LUADLUSC$ClonalNeoUQ[txH2LUADLUSC$Histology=="LUSC" & txH2LUADLUSC$ClonalNeo >= 206.25]<- "UQ.high"
txH2LUADLUSC$ClonalNeoUQ[txH2LUADLUSC$Histology=="LUSC" & txH2LUADLUSC$ClonalNeo < 206.25]<- "cUQ.low"
txH2LUADLUSC$ClonalNeoUQ <- as.factor(as.character(txH2LUADLUSC$ClonalNeoUQ))

#first test, without the genomic score
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~coldP3+age+sex+pack_years_calculated+Histology+pathologyTNM2+Adjuvant.therapy, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)

#MAIN fig 2 TEST
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

#supp test S7a - which includes total number of regions
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~coldP3+age+sex+pack_years_calculated+Histology+pathologyTNM2+Adjuvant.therapy+ClonalNeoUQ+nTotalRegions, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)

#supp test S7b - as continuous, keeping it simple 
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~rCold+age+sex+pack_years_calculated+Histology+pathologyTNM2+Adjuvant.therapy, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)

#supp test S7c - as continuous, testing independce against other scores
regional = read.csv("./tracerx100/results/trait/100xRegional_TGene.csv", stringsAsFactors = F )
ascat <- regional %>% group_by(PublicationID) %>%
  summarise(mean_ASCAT.purity=mean(ASCAT.purity))
txH2LUADLUSC <- merge(txH2LUADLUSC, ascat, all.x=T)

#correcting for Histology (by adding +Histology) in each fit doesn't make a big difference
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~rCold+subclonalCNA, data = txH2LUADLUSC)
u1<-ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~rCold+mean_ASCAT.purity, data = txH2LUADLUSC)
u2<-ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~rCold+tumour_per_mean, data = txH2LUADLUSC)
u3<-ggforest(fit, data=txH2LUADLUSC)
pdf(file=paste0(path, "/forestDFS_rCold_validationS2.pdf"), width = 18 , height = 8)
grid.arrange(u1, u2, u3, nrow=3)
dev.off()

#supp test S7d - the original test including other histology patients (only 6 more), cant do this with ClonalNeo... 
fit <- survfit(Surv(DFS_time_days, DFS_censor_variable)~rCold<=1, data = txH2) 
ggsurvplot(fit, data = txH2, conf.int = FALSE,
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
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~coldP3+age+sex+pack_years_calculated+Histology+pathologyTNM2+Adjuvant.therapy, data = txH2)
ggforest(fit, data=txH2)

#ratio cold/total
txH2LUADLUSC$coldtoAll_ratio = txH2LUADLUSC$rCold/txH2LUADLUSC$nTotalRegions
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~coldtoAll_ratio, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~coldtoAll_ratio+age+sex+pack_years_calculated+Histology+pathologyTNM2+Adjuvant.therapy, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)

########################################################################################################################

##fig 2, genomic and immune distances of paired hot-hot and cold-cold regions
########################################################################################################################
load("./tracerx100/results/trait/trait5_08Oct2018.RData")
regional <- colNA(regional, "PublicationID")
regional <- subset(regional, select = c(PublicationID, region, lymphocytes_per))
regional$region <- substring(regional$region, 3)
regional$regionCRUK <- with(regional, paste0(PublicationID, "_",region))
regional <- regional[order(regional$regionCRUK),]
library(plyr)
singleR <- count(regional, "PublicationID")
removeCRUK <- singleR$PublicationID[singleR$freq==1]
regional <- regional[!(regional$PublicationID%in%removeCRUK),]
idsToGO <- unique(regional$PublicationID)
regional2 <- subset(regional, select = -c(regionCRUK))

lymL <- vector("list", length(idsToGO))
names(lymL) <- idsToGO
for (i in 1:length(idsToGO)){
  lymL[[i]] = regional2[regional2$PublicationID %in% idsToGO[i],]
  lymL[[i]] = subset(lymL[[i]], select = -c(PublicationID))
  lymL[[i]] = t(lymL[[i]])
  colnames(lymL[[i]]) <- lymL[[i]][c(1),]
  lymL[[i]] <- lymL[[i]][-c(1),]
}
for (i in 1:length(lymL)){
  lymL[[i]] <- dist(lymL[[i]])
}
for (i in 1:length(lymL)){
  lymL[[i]] <- as.matrix(lymL[[i]])
  lymL[[i]] <- melt(lymL[[i]])[melt(upper.tri(lymL[[i]]))$value,]
  names(lymL[[i]]) <- c("r1", "r2", "immuneDist")
  lymL[[i]]$rsCRUK <- paste0(names(lymL)[i], "_",lymL[[i]]$r1, lymL[[i]]$r2)
}

#the same binary matrix needed for all mutations to compute the distance
#for each patient in a list, we make the cols are regions and rows as muts (1=mutated, 0=not)
#then simply take the Ecu. dist. for every two pairs

mut <- read.csv("./tracerx100/data/tracerx100tableS3Mutations.csv", stringsAsFactors = F)
colnames(mut)[which(names(mut) == "SampleID")] <- "PublicationID"
mut <- subset(mut, select = c(MutationID, PublicationID, Hugo_Symbol, RegionSum))
mutL <- vector("list", length(idsToGO))
names(mutL) <- idsToGO
for (i in 1:length(idsToGO)){
  mutL[[i]] = mut[mut$PublicationID %in% idsToGO[i],]
  mutL[[i]] <- mutL[[i]] %>%
    separate(RegionSum, c(regional$region[regional$PublicationID == idsToGO[i]]), ";")
}
for (i in 1:length(mutL)){
  for (f in 1:(length(mutL[[i]][, -c(1:3)]))){
    names(mutL[[i]])[f+3] <- strsplit(mutL[[i]][1,f+3], ":")[[1]][1]
  }
}
for (i in 1:length(mutL)){
  mutL[[i]][, -c(1:3)] <- gsub(".*:", "", as.matrix(mutL[[i]][, -c(1:3)]))
}
for (i in 1:length(mutL)){
  mutL[[i]] <- mutL[[i]][!is.na(names(mutL[[i]]))]
}
for (i in 1:length(mutL)){
  mutL[[i]][, -c(1:3)] <- as.numeric(gsub("/.*", "", as.matrix(mutL[[i]][, -c(1:3)])))
}
for (i in 1:length(mutL)){
  mutL[[i]][, -c(1:3)][mutL[[i]][, -c(1:3)]>0] <- 1
}

mutLBin <- mutL
for (i in 1:length(mutLBin)){
  mutLBin[[i]] <- mutLBin[[i]][, -c(1:3)] 
}
for (i in 1:length(mutLBin)){
  mutLBin[[i]] <- dist(t(mutLBin[[i]]))
}
for (i in 1:length(mutLBin)){
  mutLBin[[i]] <- as.matrix(mutLBin[[i]])
  mutLBin[[i]] <- melt(mutLBin[[i]])[melt(upper.tri(mutLBin[[i]]))$value,]
  names(mutLBin[[i]]) <- c("r1", "r2", "geneticDist")
  mutLBin[[i]]$rsCRUK <- paste0(names(mutLBin)[i], "_",mutLBin[[i]]$r1, mutLBin[[i]]$r2)
}

lyms <- do.call(rbind.data.frame, lymL)
lyms <- subset(lyms, select = c(rsCRUK, immuneDist))

mutations <- do.call(rbind.data.frame, mutLBin)
mutations <- subset(mutations, select = c(rsCRUK, geneticDist))

lyms_mut <- merge(lyms, mutations, by = "rsCRUK")

lyms_mut$PublicationID = substring(lyms_mut$rsCRUK, 1,8)
lyms_mut <- merge(lyms_mut, diagnostic[, c(3, 48)])

lyms_mut$Histology[! (lyms_mut$histology_group== "Adenocarcinoma" | 
                        lyms_mut$histology_group== "Squamous cell carcinoma") ] <- "Other"
lyms_mut$Histology[lyms_mut$histology_group== "Adenocarcinoma" ] <- "LUAD"
lyms_mut$Histology[lyms_mut$histology_group== "Squamous cell carcinoma" ] <- "LUSC"
detach("package:plyr", unload=TRUE)

#now we subset for hothot and coldcold
#hothot vs coldcold
r=regional
load("./tracerx100/results/trait/trait5_08Oct2018.RData")
regional$region <- substring(regional$region, 3)
r = merge(r, regional[, c(4, 5, 65)], by = c("PublicationID","region"), all.x = T)

lym_mutPh = lyms_mut
lym_mutPh$immuneClass2 <- NA
for (i in 1:nrow(lym_mutPh)){
  x=substring(lym_mutPh$rsCRUK[[i]], 1,11)
  y=paste0(substring(lym_mutPh$rsCRUK[[i]], 1,9), substring(lym_mutPh$rsCRUK[[i]], 12,13))
  
  lym_mutPh$immuneClass2[[i]] <- paste0(r$immuneClass_2[r$regionCRUK %in% x], "_", r$immuneClass_2[r$regionCRUK %in% y])
}

lym_mutPh_ch = lym_mutPh[lym_mutPh$immuneClass2=="hot_hot" | lym_mutPh$immuneClass2=="cold_cold" ,]
lym_mutPh_ch$PublicationID = substring(lym_mutPh_ch$rsCRUK, 1,8)

lym_mutPh_ch <- merge(lym_mutPh_ch, diagnostic[, c(3, 48)])
lym_mutPh_ch <- lym_mutPh_ch[order(lym_mutPh_ch$immuneClass2),]
comList <- list( c("hot_hot", "cold_cold"))
p1<-ggboxplot(lym_mutPh_ch[lym_mutPh_ch$histology_group=="Adenocarcinoma",], x = "immuneClass2", y = "immuneDist", 
              color = "immuneClass2", palette = c("blue", "red2"), title = "LUAD",
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "t.test")+
  theme(legend.position="") + theme(text = element_text(size=18))
p2<-ggboxplot(lym_mutPh_ch[lym_mutPh_ch$histology_group=="Adenocarcinoma",], x = "immuneClass2", y = "geneticDist", 
              color = "immuneClass2", palette = c("blue", "red2"), title = "LUAD",
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "t.test")+
  theme(legend.position="") + theme(text = element_text(size=18))

p3<-ggboxplot(lym_mutPh_ch[lym_mutPh_ch$histology_group=="Squamous cell carcinoma",], x = "immuneClass2", y = "immuneDist", 
              color = "immuneClass2", palette = c("blue", "red2"), title = "LUSC",
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "t.test")+
  theme(legend.position="") + theme(text = element_text(size=18))
p4<-ggboxplot(lym_mutPh_ch[lym_mutPh_ch$histology_group=="Squamous cell carcinoma",], x = "immuneClass2", y = "geneticDist", 
              color = "immuneClass2", palette = c("blue", "red2"), title = "LUSC",
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "t.test")+
  theme(legend.position="") + theme(text = element_text(size=18))
pdf(file=paste0(path, "/dist_hothotVscoldcold_histology.pdf"), width = 10 , height = 10)
grid.arrange(p1,p2,p3,p4,ncol=2)
dev.off()
########################################################################################################################

##fig S5, concordance of deep learning immune classification vs RNA
########################################################################################################################
regional = read.csv("./tracerx100/results/trait/100xRegional_TGene.csv", stringsAsFactors = F )
regional$Histology[regional$Histology == "Squamous cell carcinoma"] <- "Squamous-cell carcinoma"
regional$Histology[regional$Histology == "Invasive adenocarcinoma"] <- "Adenocarcinoma"
regional$Histology[regional$Histology=="Carcinosarcoma"] <- "Z"
regional$Histology[regional$Histology=="Adenosquamous carcinoma"] <- "Z"
regional$Histology[regional$Histology=="Large Cell Neuroendocrine"] <- "Z"
regional = regional[order(regional$lymphocytes_per, decreasing = F),]
regional = regional[order(regional$Histology, decreasing = F),]
RNA =readRDS("./tracerx100/data/RNA/20180926-immune-clusters-rnaseq.RDS")
RNA$sample = gsub(":", "_", rownames(RNA))
r = regional
r = merge(r, RNA, all.x=T)
r = r[order(r$Histology, decreasing = F),]
r$orig_immune_cluster[is.na(r$orig_immune_cluster)] <- "NA"
pdf(file=paste0(path, "/RNAClustervsLym_dotPlot.pdf"), width = 20 , height = 8)
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
           ggtheme = theme_pubr())+
  geom_hline(yintercept=29.2, linetype="dashed", color = "tomato")+
  geom_hline(yintercept=23.2, linetype="dashed", color = "cadetblue1")

dev.off()

#boxplots for key immune scores for deep learning based and RNA based immune classifications
regional = read.csv("./tracerx100/results/trait/100xRegional_TGene.csv", stringsAsFactors = F )
regional <- regional[regional$immuneClass_2=="hot" | regional$immuneClass_2=="cold",] 
comList <- list( c("hot", "cold"))
regional <- regional[order(regional$immuneClass_2),] 
p1<-ggboxplot(regional, x = "immuneClass_2", y = "lymphocytes_per", color = "immuneClass_2",
              add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",comparisons = comList, method = "t.test")+ 
  theme(text = element_text(size=18))
p2<-ggboxplot(regional, x = "immuneClass_2", y = "ITL_tumorRatio", color = "immuneClass_2",
              add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",comparisons = comList, method = "t.test")+ 
  theme(text = element_text(size=18))
p3<-ggboxplot(regional, x = "immuneClass_2", y = "ATL_fibroRatio", color = "immuneClass_2",
              add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",comparisons = comList, method = "t.test")+ 
  theme(text = element_text(size=18))
p4<-ggboxplot(regional, x = "immuneClass_2", y = "ASCAT.purity", color = "immuneClass_2",
              add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",comparisons = comList, method = "t.test")+ 
  theme(text = element_text(size=18))
regional = read.csv("./tracerx100/results/trait/100xRegional_TGene.csv", stringsAsFactors = F )
regional <- colNA(regional, "RNAcluster")
comList <- list( c("low", "high"))
regional <- regional[order(regional$RNAcluster, decreasing = T),] 
p5<-ggboxplot(regional, x = "RNAcluster", y = "lymphocytes_per", color = "RNAcluster",
              add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",comparisons = comList, method = "t.test")+ 
  theme(text = element_text(size=18))
p6<-ggboxplot(regional, x = "RNAcluster", y = "ITL_tumorRatio", color = "RNAcluster",
              add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",comparisons = comList, method = "t.test")+ 
  theme(text = element_text(size=18))
p7<-ggboxplot(regional, x = "RNAcluster", y = "ATL_fibroRatio", color = "RNAcluster",
              add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",comparisons = comList, method = "t.test")+ 
  theme(text = element_text(size=18))
p8<-ggboxplot(regional, x = "RNAcluster", y = "ASCAT.purity", color = "RNAcluster",
              add = "jitter", border = "white", palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",comparisons = comList, method = "t.test")+ 
  theme(text = element_text(size=18))

pdf(file=paste0(path, "/RNAvsDeppL_class.pdf"), height = 9, width = 18)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8,  ncol=4)
dev.off()

#scores provided by email for CRUK0065
a = r[r$PublicationID=="CRUK0065",]
a$pathology_TILs = c(30, 5, 80, 10, 70, 60, NA)
a=colNA(a, "pathology_TILs")
ggscatter(a, x = "lymphocytes_per", y = "pathology_TILs",
          font.label = 12, 
          shape = 19, size = 2, label = "region",  cor.coef = F, cor.method = "spearman",
          color = "orig_immune_cluster", palette = c("red2", "blue2", "gray45"))

#JLQ seprately scored lym infilt. for 47 regions
John = read.csv("./tracerx100/data/John_lymScores/RNAImaging_concordance_181019.csv", stringsAsFactors = F)
John = merge(John, r, all.x = T)

names(John)[names(John) == "John"] <- "pathology_LymEstimate"
ggscatter(John, x = "lymphocytes_per", y = "pathology_LymEstimate", add = "reg.line", 
          xlab = "Lymphocyte%", ylab = "Pathology lym estimate (%)",
          font.label = 12, cor.coef = T, cor.method = "spearman", 
          shape = 19, size = 2)+ theme(legend.position="none")+xlim(10, 60)+ylim(5, 75)

########################################################################################################################

##fig S6, distribution of lym% over stage, histology 
########################################################################################################################
load("./tracerx100/results/trait/trait5_15May2018.RData")
diagnosticLUAD <- diagnostic[diagnostic$histology_group=="Adenocarcinoma",]
diagnosticLUSC <- diagnostic[diagnostic$histology_group=="Squamous cell carcinoma",]
geneticLUAD <- genetic[genetic$histology_group=="Adenocarcinoma",]
geneticLUSC <- genetic[genetic$histology_group=="Squamous cell carcinoma",]
diagnosticIHC = colNA(diagnostic, "cd8_per")
diagnosticIHC$cd8foxp3Ratio = diagnosticIHC$cd8_per/diagnosticIHC$foxp3_per
diagnosticIHCLUAD <- diagnosticIHC[diagnosticIHC$histology_group=="Adenocarcinoma",]
diagnosticIHCLUSC <- diagnosticIHC[diagnosticIHC$histology_group=="Squamous cell carcinoma",]
tx$coldP = as.factor(as.character(tx$coldP))
tx$coldP3 = as.factor(as.character(tx$coldP3))
tx$pathologyTNM = as.factor(as.character(tx$pathologyTNM))
txH2 = tx
txH2$Histology[txH2$histology_group=="Adenocarcinoma" ] <- "LUAD"
txH2$Histology[txH2$histology_group=="Squamous cell carcinoma" ] <- "LUSC"
txH2$Histology[is.na(txH2$Histology)] <- "Other"
txH2$Histology = as.factor(as.character(txH2$Histology))
txH2$sex = as.factor(as.character(txH2$sex))
txH2$pathologyTNM2=txH2$pathologyTNM
levels(txH2$pathologyTNM2)[levels(txH2$pathologyTNM2)=="IIIA"] <- "III"
levels(txH2$pathologyTNM2)[levels(txH2$pathologyTNM2)=="IIIB"] <- "III"
txH2$AdjuvantTherapy = txH2$adjuvant_treatment_given
txH2$AdjuvantTherapyC[txH2$adjuvant_treatment_given == ""]<- 1
txH2$AdjuvantTherapyC[txH2$adjuvant_treatment_given == "Platinum chemo"]<- 0
txH2$AdjuvantTherapyC[txH2$adjuvant_treatment_given == "Radiotherapy"]<- 0
txH2$AdjuvantTherapyC[txH2$adjuvant_treatment_given == "Platinum chemo/radiotherapy"]<- 0
txH2LUADLUSC <- txH2[ which(txH2$Histology=="LUAD" 
                            | txH2$Histology == "LUSC"), ]
txH2LUADLUSC$Histology <- factor(txH2LUADLUSC$Histology)

r = regional[, c(6, 12)]
rU = uniqueLTX_survPrepare2(r)
rU = rU[, c(1:2, 4, 6:8)]

txH2LUADLUSC = merge(txH2LUADLUSC, rU, all.x = T)
txH2 = merge(txH2, rU, all.x = T)


p1<-ggboxplot(txH2LUADLUSC, x = "pathologyTNM2", y = "lymphocytes_per_sd", color = "pathologyTNM2",title = "LUADLUSC 79",
              add = "jitter", border = "white")+
  stat_compare_means(label = "p.signif",method = "t.test", paired = F, ref.group = ".all.", hide.ns = F)+
  theme(legend.position="") + theme(text = element_text(size=18)) 
p2<-ggboxplot(txH2LUADLUSC, x = "pathologyTNM2", y = "lymphocytes_per_mean", color = "pathologyTNM2",title = "LUADLUSC 79",
              add = "jitter", border = "white")+
  stat_compare_means(label = "p.signif",method = "t.test", paired = F, ref.group = ".all.", hide.ns = F)+
  theme(legend.position="") + theme(text = element_text(size=18)) 

pdf(file="./tracerx100/results/Regional/lymLandscpae_SD_density_Stage.pdf", width = 10 , height = 6)
ggdensity(txH2LUADLUSC, x = "lymphocytes_per_sd", color = "Histology", palette = c("#98AED6", "#682A45"),
             add = "mean", rug = TRUE)
ggdensity(txH2LUADLUSC, x = "lymphocytes_per_sd", color = "pathologyTNM", 
             add = "mean", rug = TRUE)
ggdensity(txH2LUADLUSC, x = "lymphocytes_per_mean", color = "pathologyTNM", 
             add = "mean", rug = TRUE)
grid.arrange(p1,p2, nrow=2)
dev.off()
########################################################################################################################

##fig S8, corr of lym% in regional vs diagnsotic samples + fig S9b lym% normal vs lym% tumor min
########################################################################################################################
diagnostic = read.csv("./tracerx100/results/trait/100xDiagnostic.csv", stringsAsFactors = F)
tx = read.csv("./tracerx100/results/trait/100xCombined_2.csv", stringsAsFactors = F)
regional = read.csv("./tracerx100/results/trait/100xRegional_TGene.csv", stringsAsFactors = F)
r <- regional[, c(4, 12)]
r <- r %>% group_by(PublicationID) %>%
  summarize(regional_min = min(lymphocytes_per),
            regional_max = max(lymphocytes_per),
            regional_mean = mean(lymphocytes_per))
diL <- diagnostic[, c(3, 37, 46:48)]
colnames(diL)[which(names(diL) == "lymphocytes_per")] <- "diagnostic"
l <- merge(diL, r, all.x = T)
l <- colNA(l, "regional_min")
load("./tracerx100/results/trait/trait5_15May2018.RData")
l <- merge(l, tx[, c(2,86)], all.x=T)
l$Histology[l$histology_group== "Adenocarcinoma"] <- "LUAD"
l$Histology[l$histology_group== "Squamous cell carcinoma"] <- "LUSC"
l$Histology[l$histology_group== "sarcomatoid carcinoma of pleomorphic type arising from adenocarcinoma"] <- "Other"
l$Histology[l$histology_group== "Large Cell Neuroendocrine Carcinoma"] <- "Other"
l$Histology[l$histology_group== "Pleomorphic carcinoma"] <- "Other"
l$Histology[l$histology_group== "Adenosquamous carcinoma"] <- "Other"
p1 <- ggscatter(l, x = "diagnostic", y = "regional_mean", show.legend.text = F, 
                xlab = "Diagnostic lymphocytes%", ylab = "Regional mean lymphocytes%", font.label = 4,
                shape = 19, size = 2,
                color = "Histology", palette = c("#98AED6", "#682A45", "black"),
                add.params = list(color = "grey50", fill = "azure3"))+ 
  stat_cor(aes(color = Histology), label.x = 25, method="spearman", show.legend = F)+
  theme(legend.position="none")+
  geom_errorbar(aes(ymin = regional_min, ymax = regional_max), 
                alpha = 0.25, width=0.2, size=1, color="gray20")
ggscatter(l, x = "diagnostic", y = "regional_mean", cor.coef = T)
p2 <- ggscatter(l, x = "diagnostic", y = "regional_min", show.legend.text = F, 
                xlab = "Diagnostic lymphocytes%", ylab = "Regional min lymphocytes%", font.label = 4,
                shape = 19, size = 2,
                color = "Histology", palette = c("#98AED6", "#682A45", "black"),
                add.params = list(color = "grey50", fill = "azure3"))+ 
  stat_cor(aes(color = Histology), label.x = 25, method="spearman", show.legend = F)+
  theme(legend.position="none")
ggscatter(l, x = "diagnostic", y = "regional_min", cor.coef = T)
p3 <- ggscatter(l, x = "cd8_per", y = "regional_min", show.legend.text = F, 
                xlab = "CD8%", ylab = "Regional min lymphocytes%", font.label = 4,
                shape = 19, size = 2,
                color = "Histology", palette = c("#98AED6", "#682A45", "black"),
                add.params = list(color = "grey50", fill = "azure3"))+ 
  stat_cor(aes(color = Histology), label.x = 15, method="spearman", show.legend = F)+
  theme(legend.position="none")
ggscatter(l, x = "cd8_per", y = "regional_min", cor.coef = T)

## Normal

diagnostic = read.csv("./tracerx100/results/trait/100xDiagnostic.csv", stringsAsFactors = F)
tx = read.csv("./tracerx100/results/trait/100xCombined_2.csv", stringsAsFactors = F)
regional = read.csv("./tracerx100/results/trait/100xRegional.csv", stringsAsFactors = F)
regional <- colNA(regional, "PublicationID")
N <- regional[regional$region=="N",]
N <- N[, c(1, 10)]
colnames(N)[which(names(N) == "lymphocytes_per")] <- "normal"
regional = read.csv("./tracerx100/results/trait/100xRegional.csv", stringsAsFactors = F)
regional <- colNA(regional, "PublicationID")
LN <- regional[regional$region=="LN1" | regional$region=="LN2",]
LN <- LN[, c(1, 10)]
colnames(LN)[which(names(LN) == "lymphocytes_per")] <- "lymph_node"
regional = read.csv("./tracerx100/results/trait/100xRegional_TGene.csv", stringsAsFactors = F)
regional <- regional[regional$immuneClass_2=="hot" | regional$immuneClass_2=="cold",] 
rP <- regional[, c(1, 3:4, 11, 64)]
rP <- merge(rP, N, all.x = T)
rP <- merge(rP, LN, all.x = T)

LN.t <- merge(LN, rP, all.x = T)
N.t <- merge(N, rP, all.x = T)
LN.tt <- merge(LN, tx[, c(1, 3:77)], all.x = T)
LN.tt <- merge(LN.tt, r, all.x = T)
LN.tt <- LN.tt[, c(1, 2, 10, 78:80)]
LN.tt <- merge(LN.tt, N, all.x = T)
colnames(LN.tt)[which(names(LN.tt) == "lymphocytes_per")] <- "diagnostic"
LN.tt <- reshape::melt(LN.tt) 
colnames(LN.tt)[which(names(LN.tt) == "variable")] <- "site"
colnames(LN.tt)[which(names(LN.tt) == "value")] <- "infiltration"
N.tt <- merge(N, tx[, c(1, 3:77)], all.x = T)
N.tt <- merge(N.tt, r, all.x = T)
N.tt$Histology[N.tt$Histology== "Invasive adenocarcinoma"] <- "LUAD"
N.tt$Histology[N.tt$Histology== "Squamous cell carcinoma"] <- "LUSC"
N.tt$Histology[N.tt$Histology== "Large cell carcinoma"] <- "Other"
N.tt$Histology[N.tt$Histology== "Large Cell Neuroendocrine"] <- "Other"
N.tt$Histology[N.tt$Histology== "Carcinosarcoma"] <- "Other"
N.tt$Histology[N.tt$Histology== "Adenosquamous carcinoma"] <- "Other"
p4 <- ggscatter(N.tt, x = "normal", y = "regional_min", show.legend.text = F, 
                xlab = "Normal lymphocytes%", ylab = "Regional min lymphocytes%", font.label = 4,
                shape = 19, size = 2,
                color = "Histology", palette = c("#98AED6", "#682A45", "black"),
                add.params = list(color = "grey50", fill = "azure3"))+ 
  stat_cor(aes(color = Histology), label.x = 10, method="spearman", show.legend = F)+
  theme(legend.position="none")
ggscatter(N.tt, x = "normal", y = "regional_min", cor.coef = T)

pdf(file=paste0(path, "/figSLym.pdf"), width = 10 , height = 16)
grid.arrange(p1,p2,p3,p4,  ncol=2)
dev.off()

########################################################################################################################

##fig 3, lym% least vs most survival in Tx100 and TCGA
########################################################################################################################

#need to load TCGA:
luad <- read.csv(file = "/Users/kjabbar/Dropbox (ICR)/yuanlab/Projects/lung/tcga/dev/immuneScoringPaper/LUADSummary_AllSlides.csv", stringsAsFactors = F)
lusc <- read.csv(file = "/Users/kjabbar/Dropbox (ICR)/yuanlab/Projects/lung/tcga/dev/immuneScoringPaper/LUSCSummary_AllSlides.csv", stringsAsFactors = F)
luad$PATIENT_SCREENING_ID <- substr(luad$file_name, 1, 12)
lusc$PATIENT_SCREENING_ID <- substr(lusc$file_name, 1, 12)
luad = luad[, c(19, 17)]
lusc = lusc[, c(19, 17)]
load("/Users/kjabbar/Dropbox (ICR)/yuanlab/Projects/lung/tcga/trait/TCGAmin.RData")
luad = luad[luad$PATIENT_SCREENING_ID %in% LUADmin$bcr_patient_barcode,]
lusc = lusc[lusc$PATIENT_SCREENING_ID %in% LUSCmin$bcr_patient_barcode,]
luad = uniqueLTX_survPrepare2(luad)
lusc = uniqueLTX_survPrepare2(lusc)
names(luad)[names(luad) == "PATIENT_SCREENING_ID"] <- "bcr_patient_barcode"
names(lusc)[names(lusc) == "PATIENT_SCREENING_ID"] <- "bcr_patient_barcode"
LUSCmean = read.csv("/Users/kjabbar/Dropbox (ICR)/yuanlab/Projects/lung/tcga/dev/immuneScoringPaper/results/LUSCmean.csv")
LUADmean = read.csv("/Users/kjabbar/Dropbox (ICR)/yuanlab/Projects/lung/tcga/dev/immuneScoringPaper/results/LUADmean.csv")
LUSCclin = LUSCmean[, c(8:10, 13, 22:26)]
LUADclin = LUADmean[, c(8:10, 13, 22:26)]
luad= merge(luad, LUADclin, all.x=T)
lusc = merge(lusc, LUSCclin, all.x = T)
luad$Histology = rep("LUAD")
lusc$Histology = rep("LUSC")
tcga <- rbind(luad, lusc)
tcga$Histology = as.factor(as.character(tcga$Histology))
fit <- coxph(Surv(time, event)~lym_per_min+Histology, data = tcga)
u1=ggforest(fit, data=tcga)
fit <- coxph(Surv(time, event)~lym_per_mean+Histology, data = tcga)
u2=ggforest(fit, data=tcga)
fit <- coxph(Surv(time, event)~lym_per_max+Histology, data = tcga)
u3=ggforest(fit, data=tcga)
fit <- coxph(Surv(time, event)~lym_per_median+Histology, data = tcga)
u4=ggforest(fit, data=tcga)
fit <- coxph(Surv(time, event)~lym_per_dist+Histology, data = tcga)
u5=ggforest(fit, data=tcga)
fit <- coxph(Surv(time, event)~lym_per_sd+Histology, data = tcga)
u6=ggforest(fit, data=tcga)
#now a simple summary df from the above tests
summaryPsTCGA <- data.frame(
  Measure = c("6Minimum", "5Mean", "4Maximum",
              "3Median", "2Distance", "1SD"),
  PVal = c(0.035, 0.055, 0.149, 0.058, 0.817, 0.726),
  lower95 = c(0.0063, 0.011, 0.043, 0.011, 0.15, 0.083),
  upper95 = c(0.83, 1.1, 1.6, 1.1, 2, 2), #10.9, 35.4
  hr = c(0.073, 0.11, 0.26, 0.11, 1.3, 1.7))
summaryPsTCGA$PVal_log10 <- -log10(summaryPsTCGA$PVal)

#now the same for tx
load("./tracerx100/results/trait/trait5_15May2018.RData")
diagnosticLUAD <- diagnostic[diagnostic$histology_group=="Adenocarcinoma",]
diagnosticLUSC <- diagnostic[diagnostic$histology_group=="Squamous cell carcinoma",]
geneticLUAD <- genetic[genetic$histology_group=="Adenocarcinoma",]
geneticLUSC <- genetic[genetic$histology_group=="Squamous cell carcinoma",]
diagnosticIHC = colNA(diagnostic, "cd8_per")
diagnosticIHC$cd8foxp3Ratio = diagnosticIHC$cd8_per/diagnosticIHC$foxp3_per
diagnosticIHCLUAD <- diagnosticIHC[diagnosticIHC$histology_group=="Adenocarcinoma",]
diagnosticIHCLUSC <- diagnosticIHC[diagnosticIHC$histology_group=="Squamous cell carcinoma",]
tx$coldP = as.factor(as.character(tx$coldP))
tx$coldP3 = as.factor(as.character(tx$coldP3))
tx$pathologyTNM = as.factor(as.character(tx$pathologyTNM))
txH2 = tx
txH2$Histology[txH2$histology_group=="Adenocarcinoma" ] <- "LUAD"
txH2$Histology[txH2$histology_group=="Squamous cell carcinoma" ] <- "LUSC"
txH2$Histology[is.na(txH2$Histology)] <- "Other"
txH2$Histology = as.factor(as.character(txH2$Histology))
txH2$sex = as.factor(as.character(txH2$sex))
txH2$pathologyTNM2=txH2$pathologyTNM
levels(txH2$pathologyTNM2)[levels(txH2$pathologyTNM2)=="IIIA"] <- "III"
levels(txH2$pathologyTNM2)[levels(txH2$pathologyTNM2)=="IIIB"] <- "III"
txH2$AdjuvantTherapy = txH2$adjuvant_treatment_given
txH2$AdjuvantTherapyC[txH2$adjuvant_treatment_given == ""]<- 1
txH2$AdjuvantTherapyC[txH2$adjuvant_treatment_given == "Platinum chemo"]<- 0
txH2$AdjuvantTherapyC[txH2$adjuvant_treatment_given == "Radiotherapy"]<- 0
txH2$AdjuvantTherapyC[txH2$adjuvant_treatment_given == "Platinum chemo/radiotherapy"]<- 0
txH2LUADLUSC <- txH2[ which(txH2$Histology=="LUAD" 
                            | txH2$Histology == "LUSC"), ]
txH2LUADLUSC$Histology <- factor(txH2LUADLUSC$Histology)
r = regional[, c(6, 12)]
rU = uniqueLTX_survPrepare2(r)
rU = rU[, c(1:2, 4, 6:8)]
txH2LUADLUSC = merge(txH2LUADLUSC, rU, all.x = T)
txH2 = merge(txH2, rU, all.x = T)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~lymphocytes_per_min+Histology, data = txH2LUADLUSC)
u1<-ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~lymphocytes_per+Histology, data = txH2LUADLUSC)
u2<-ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~lymphocytes_per_mean+Histology, data = txH2LUADLUSC)
u3<-ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~lymphocytes_per_max+Histology, data = txH2LUADLUSC)
u4<-ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~lymphocytes_per_median+Histology, data = txH2LUADLUSC)
u5<-ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~lymphocytes_per_dist+Histology, data = txH2LUADLUSC)
u6<-ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~lymphocytes_per_sd+Histology, data = txH2LUADLUSC)
u7<-ggforest(fit, data=txH2LUADLUSC)
#add more scores from IHC diagnostic
txH2LUADLUSC$cd8foxp3Ratio = txH2LUADLUSC$cd8_per/txH2LUADLUSC$foxp3_per
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~cd8foxp3Ratio+Histology, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~cd8_per+Histology, data = txH2LUADLUSC)
ggforest(fit, data=txH2LUADLUSC)
#another simple summary df
summaryPs <- data.frame(
  Measure = c("6Minimum", "5Mean", "4Maximum",
              "3Median", "2Distance", "1SD", "0Diagnostic", "00CD8%", "000CD8/FOXP3"),
  PVal = c(0.062, 0.124, 0.472, 0.119, 0.129, 0.995, 0.048, 0.6, 0.62),
  lower95 = c(0.94, 0.94, 0.96, 0.94, 0.99, 0.92, 0.86, 0.88, 0.67),
  upper95 = c(1, 1, 1, 1, 1.1, 1.1, 1, 1.1, 1.3), #10.9, 35.4
  hr = c(0.97, 0.97, 0.99, 0.97, 1.03, 1, 0.93, 0.97, 0.92))
summaryPs$PVal_log10 <- -log10(summaryPs$PVal)

#forest plots for the two
pdf(file=paste0(path, "/NatForest_lym_txtcga.pdf"), height = 5, width = 5)
ggplot(data=summaryPsTCGA,aes(x=hr,y=Measure))+
  geom_point(aes(size=PVal_log10), colour="#567C3B", fill="#567C3B",shape=21)+
  geom_errorbarh(aes(xmin=lower95,xmax=upper95), height=0)+
  geom_vline(xintercept=1,linetype="dashed")+
  scale_size_continuous(breaks=c(5000,10000,15000))+theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank()) + labs(x = "Hazard ratio", y = "Lymphocytic score")
ggplot(data=summaryPs,aes(x=hr,y=Measure))+
  geom_point(aes(size=PVal_log10), colour="#CB6728", fill="#CB6728", shape=21)+
  geom_errorbarh(aes(xmin=lower95,xmax=upper95), height=0)+
  geom_vline(xintercept=1,linetype="dashed")+
  scale_size_continuous(breaks=c(5000,10000,15000))+theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.x = element_blank()) + labs(x = "Hazard ratio", y = "Lymphocytic score")

dev.off()
########################################################################################################################

##fig 3, tumor-normal divergence and survival analysis
########################################################################################################################
load("./tracerx100/results/trait/trait5_08Oct2018.RData")
txH2=merge(txH2, normal, all.x = T)
txH2LUAD=merge(txH2LUAD, normal, all.x = T)
txH2LUSC=merge(txH2LUSC, normal, all.x = T)
txH2 = colNA(txH2, "normal")
txH2$normalDiff = txH2$regional_mean-txH2$normal
txH2$LQ.normalMeanDiff[txH2$normalDiff <= quantile(txH2$normalDiff)[[2]]] <- "bLow"
txH2$LQ.normalMeanDiff[txH2$normalDiff > quantile(txH2$normalDiff)[[2]]] <- "High"
txH2$LQ.normalMeanDiff = as.factor(as.character(txH2$LQ.normalMeanDiff))
txH2LUAD = txH2[txH2$Histology=="LUAD" ,]
txH2LUSC = txH2[txH2$Histology=="LUSC" ,]

#paired boxplots for tumor normal dive for luad and lusc
rLUAD = txH2LUAD[, c(1, 118, 115)]
rLUAD = melt(rLUAD, variable.name = "Tissue", id.vars = "PublicationID")

rLUSC = txH2LUSC[, c(1, 118, 115)]
rLUSC = melt(rLUSC, variable.name = "Tissue", id.vars = "PublicationID")

rLUAD = merge(rLUAD, txH2LUAD[, c(1, 123)])
rLUSC = merge(rLUSC, txH2LUSC[, c(1, 123)])

p1<-ggpaired(rLUAD, x = "Tissue", y = "value", title = "LUAD",
             color = "Tissue", line.color = "gray", line.size = 0.4,
             palette = "npg")+
  stat_compare_means(paired = TRUE, method = "t.test")
p2<-ggpaired(rLUSC, x = "Tissue", y = "value", title = "LUSC",
             color = "Tissue", line.color = "gray", line.size = 0.4,
             palette = "npg")+
  stat_compare_means(paired = TRUE, method = "t.test")

p3<-ggplot(rLUAD, aes(x = Tissue, y = value)) +
  geom_boxplot(aes(fill = Tissue), alpha = 0.2, col = "grey") +
  geom_point() +
  geom_line(aes(group = PublicationID, col = LQ.normalMeanDiff)) +
  scale_colour_manual(values = c("blue", "red2")) + theme_classic() +theme(legend.position="none")
p4<-ggplot(rLUSC, aes(x = Tissue, y = value)) +
  geom_boxplot(aes(fill = Tissue), alpha = 0.2, col = "grey") +
  geom_point() +
  geom_line(aes(group = PublicationID, col = LQ.normalMeanDiff)) +
  scale_colour_manual(values = c("blue", "red2")) + theme_classic() +theme(legend.position="none")

pdf(file=paste0(path, "/paired_normalToMean.pdf"), width = 8 , height = 6)
grid.arrange(p1, p2, ncol=2)
dev.off()
pdf(file=paste0(path, "/paired_normalToMean_colored.pdf"), width = 8 , height = 6)
grid.arrange(p3,p4, ncol=2)
dev.off()

#KM curves in 3d
splots <- list()
fit <- survfit(Surv(DFS_time_days, DFS_censor_variable)~
                 normalDiff<=quantile(txH2LUAD$normalDiff)[[2]], data = txH2LUAD)
splots[[1]]<-ggsurvplot(fit, data = txH2LUAD, conf.int = FALSE,
                        pval = TRUE, pval.size = 5, pval.coord = c(0.2, 0.1),
                        linetype = "solid",
                        #surv.median.line = "hv",
                        legend = "none", legend.title = "", title = "Tx LUAD",
                        legend.labs = c("High", "Low"),
                        surv.plot.height = 0.7, palette = c("red2", "blue"),
                        risk.table = TRUE,
                        risk.table.col = "black", break.time.by = 200,
                        tables.height = 0.25,
                        tables.theme = theme_cleantable(),
                        tables.y.text = TRUE, risk.table.title = "Number at Risk",
                        tables.x.text = "", xlim = c(0, 1400),
                        xlab = "Days to Death or Recurrence", ylab = "Disease-free Survival")
fit <- survfit(Surv(DFS_time_days, DFS_censor_variable)~
                 normalDiff<=quantile(txH2LUSC$normalDiff)[[2]], data = txH2LUSC)
splots[[2]]<-ggsurvplot(fit, data = txH2LUSC, conf.int = FALSE,
                        pval = TRUE, pval.size = 5, pval.coord = c(0.2, 0.1),
                        linetype = "solid",
                        #surv.median.line = "hv",
                        legend = "none", legend.title = "", title = "Tx LUSC",
                        legend.labs = c("High", "Low"),
                        surv.plot.height = 0.7, palette = c("red2", "blue"),
                        risk.table = TRUE,
                        risk.table.col = "black", break.time.by = 200,
                        tables.height = 0.25,
                        tables.theme = theme_cleantable(),
                        tables.y.text = TRUE, risk.table.title = "Number at Risk",
                        tables.x.text = "", xlim = c(0, 1400),
                        xlab = "Days to Death or Recurrence", ylab = "Disease-free Survival")
pdf(file=paste0(path, "/KMDFS_rMeanMinusNormal_lowerQ.pdf"), width = 12, height = 8)
arrange_ggsurvplots(splots, print = TRUE,ncol = 2, nrow = 1)
dev.off()

#MAIN TEST fig 3e
#with clonal neo but local quartile: quantile(txH2$ClonalNeo[txH2$Histology=="LUAD"])
txH2$ClonalNeoUQ[txH2$Histology=="LUAD" & txH2$ClonalNeo >= 203.25]<- "UQ.high"
txH2$ClonalNeoUQ[txH2$Histology=="LUAD" & txH2$ClonalNeo < 203.25]<- "cUQ.low"
txH2$ClonalNeoUQ[txH2$Histology=="LUSC" & txH2$ClonalNeo >= 225]<- "UQ.high"
txH2$ClonalNeoUQ[txH2$Histology=="LUSC" & txH2$ClonalNeo < 225]<- "cUQ.low"
txH2$ClonalNeoUQ <- as.factor(as.character(txH2$ClonalNeoUQ))
fit <- coxph(Surv(DFS_time_days, DFS_censor_variable)~LQ.normalMeanDiff+age+sex+pack_years_calculated+Histology+pathologyTNM2+Adjuvant.therapy+ClonalNeoUQ, data = txH2)
ggforest(fit, data=txH2)





#side test which is reported in text for number of cold regions vs tumor normal dive

#does this mean the high rcold patients all have lower tumour lym compared to their normals?
comList <- list( c("<=1", ">1"))
p1<-ggboxplot(txH2, x = "coldP3", y = "normalDiff", color = "coldP3", 
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "t.test")+
  theme(legend.position="") + theme(text = element_text(size=18))
comList <- list( c("High", "Low"))
p2<-ggboxplot(txH2, x = "LQ.normalMeanDiff", y = "rCold", color = "LQ.normalMeanDiff", 
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "t.test")+
  theme(legend.position="") + theme(text = element_text(size=18))
pdf(file=paste0(path, "/rColdVSLQnormalMeanDiff.pdf"), width = 8 , height = 6)
grid.arrange(p1,p2,ncol=2)
dev.off()

########################################################################################################################

##fig S9a, levels of lym infiltration in normal lung vs tumor
########################################################################################################################
load("./tracerx100/results/trait/trait5_08Oct2018.RData")
normal = colNA(normal, "normal")
normal = merge(normal, txH2[, c(2, 106)])
normal = normal[, -c(6:7)]
normal = na.omit(normal)
r = melt(normal)
r = r[! r$variable=="diagnostic",]
normalLUAD = normal[normal$Histology=="LUAD",]
normalLUSC = normal[normal$Histology=="LUSC",]

#doing all reported paired t tests manually 
t.test(normalLUSC$normal, normalLUSC$regional_min, paired = TRUE, alternative = "two.sided")
t.test(normalLUAD$normal, normalLUAD$regional_min, paired = TRUE, alternative = "two.sided")
#etc for the others

comList <- list( c("normal", "regional_min"), c("normal", "regional_max"), c("normal", "regional_mean"))
ggboxplot(r[r$Histology=="LUAD",], x = "variable", y = "value", color = "variable", title = "LUAD",
          add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "t.test")+
  theme(legend.position="") + theme(text = element_text(size=18)) +scale_color_d3()

ggboxplot(r[r$Histology=="LUSC",], x = "variable", y = "value", color = "variable", title = "LUSC",
          add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "t.test")+
  theme(legend.position="") + theme(text = element_text(size=18)) +scale_color_d3()
#NOTE: the statistical tests for these groups are done manually just like reported above, it works more accurately
#because we're comparing normal to all others
########################################################################################################################

##fig 4 + fig S10a, spatial immune subsets in tx LUAD  
########################################################################################################################
tx = read.csv("./tracerx100/results/trait/100xCombined_5.csv", stringsAsFactors = F)
diagnostic = read.csv("./tracerx100/results/trait/100xDiagnostic.csv", stringsAsFactors = F)
diagnostic$cd8foxp3Ratio = diagnostic$cd8_per/diagnostic$foxp3_per
diagnostic$ATL_fibroRatio <- diagnostic$ATL/diagnostic$fibroblasts

dATL = diagnostic[diagnostic$Histology == "Invasive adenocarcinoma",]
dATL = dATL[order(dATL$ATL_fibroRatio, decreasing = T),]
mydist <- function(row){
  dists <- (row[["x"]] - df2$x)^2 + (row[["y"]]- df2$y)^2
  return(cbind(df2[which.min(dists),], distance = min(dists)))
}
IHCDir<-'./tracerx100/results/Diagnostic/ATLNeo/LUAD_markedATL/classification/csv/cellPos/'
TILDir<-'./tracerx100/results/Diagnostic/ATLNeo/LUAD_TILsCellPos/'
CELLS<-dir(path = IHCDir,pattern = '*.csv')
options(stringsAsFactors = FALSE)

LUADATL_list <- vector("list", length(CELLS))
names(LUADATL_list) <- CELLS
i=0

for (iCELL in CELLS){
  IHC <- read.csv(paste0(IHCDir, iCELL))
  IHC = IHC[, c(2:3, 1)]
  
  print(iCELL)
  TIL <- read.csv(paste0(TILDir, iCELL))
  
  # TIL2 = TIL[TIL$x<=max(IHC$x),]
  # TIL2 = TIL2[TIL2$x>=min(IHC$x),]
  # TIL2 = TIL2[TIL2$y<=max(IHC$y),]
  # TIL2 = TIL2[TIL2$y>=min(IHC$y),]
  
  #cd8-foxp3 regulation
  cd8 = IHC[IHC$class=="cd8" ,]
  foxp3 = IHC[IHC$class =="p3" ,]
  #take out CD8s from IHC
  IHC = IHC[! IHC$class=="cd8",]
  
  df2 = foxp3
  df1 = cd8
  cd8R <- cbind(df1, do.call(rbind, lapply(1:nrow(df1), function(x) mydist(df1[x,])))) 
  names(cd8R) <- c("x", "y", "c1", "xfoxp3", "yfoxp3", "c2", "distance")
  # cd8R$class[cd8R$distance>=median(cd8R$distance)] <- "cd8"
  # cd8R$class[cd8R$distance<median(cd8R$distance)] <- "cd8Reg"
  cd8R$class[cd8R$distance>=5] <- "cd8"
  cd8R$class[cd8R$distance<5] <- "cd8Reg"
  cd8R=cd8R[, c(1,2,8)]
  
  #25 is tooo much! change to 20um 
  #16*0.3225*4 = 20.64, 16*0.3225*5 = 25.8
  
  #bring back CD8s regulated into IHC
  IHC = rbind(IHC, cd8R)
  
  #
  df2 = TIL
  df1 = IHC #smaller df is df2
  dR <- cbind(df1, do.call(rbind, lapply(1:nrow(df1), function(x) mydist(df1[x,])))) 
  names(dR) <- c("xihc", "yihc", "immuneClass", "xtil", "ytil", "TILClass", "distance")
  i=i+1
  LUADATL_list[[i]] = dR
}

#distance < 20 (approximate pixels radius)
LUADATL_list2 = LUADATL_list
LUADATL_list2[[6]] <- NULL
LUADATL_list2[[8]] <- NULL
LUADATL_bind <- bind_rows(LUADATL_list2, .id = "file_name")
LUADATL_bind = LUADATL_bind[LUADATL_bind$distance<20 ,]

#statitical sign. 
p1 = t(table(LUADATL_list2[[1]]$immuneClass, LUADATL_list2[[1]]$TILClass))
p2 = t(table(LUADATL_list2[[2]]$immuneClass, LUADATL_list2[[2]]$TILClass))
p3 = t(table(LUADATL_list2[[3]]$immuneClass, LUADATL_list2[[3]]$TILClass))
p4 = t(table(LUADATL_list2[[4]]$immuneClass, LUADATL_list2[[4]]$TILClass))
p5 = t(table(LUADATL_list2[[5]]$immuneClass, LUADATL_list2[[5]]$TILClass))
p6 = t(table(LUADATL_list2[[6]]$immuneClass, LUADATL_list2[[6]]$TILClass))
p7 = t(table(LUADATL_list2[[7]]$immuneClass, LUADATL_list2[[7]]$TILClass))
pSummin = rbind(p1,p2,p3,p4,p5,p6,p7)
pSumm = data.frame(rbind(p1,p2,p3,p4,p5,p6,p7))
pSumm$TILClass = rep(c("ATL", "DTL", "ITL"))

load("./tracerx100/results/Diagnostic/ATLNeo/R_datafiles/exec_2.RData")
#stats per images - split the main final lists
IHCDir<-'./tracerx100/results/Diagnostic/ATLNeo/Map/'
CELLS<-dir(path = IHCDir,pattern = '*.csv')
options(stringsAsFactors = FALSE)
Da_list <- vector("list", length(CELLS))
names(Da_list) <- CELLS
i=0
for (iCELL in CELLS){
  print(iCELL)
  i = i+1
  Da_list[[i]] = read.csv(paste0(IHCDir, iCELL))
}

#split the LUADATL_list based on the tiles 
das = LUADATL_list
for (i in 1:length(das)){
  
  das[[i]]$tile <- NULL
  for (k in 1:nrow(das[[i]])){
    x = das[[i]]$xtil[[k]]
    y = das[[i]]$ytil[[k]]
    
    ada = Da_list[[i]]$Daname[x >= Da_list[[i]]$start_w &
                                x <= Da_list[[i]]$end_w &
                                y >= Da_list[[i]]$start_h &
                                y <= Da_list[[i]]$end_h]
    
    das[[i]]$tile[[k]] <- ada
    
  }
}

for (i in 1:length(das)){
  
  das[[i]]$daSample = paste0(names(das)[[i]], '_',das[[i]]$tile)
  
}

#distance < 20
das_bind <- bind_rows(das, .id = "file_name")
das_bind = das_bind[das_bind$distance<20 ,]
das_bind = das_bind[, c(4,7,10)]

das_bindP = das_bind %>% group_by(daSample, TILClass, immuneClass) %>%
  tally()

DaSamples = unique(das_bind$daSample)
das_bind_list <- vector("list", length(DaSamples))
names(das_bind_list) <- DaSamples

for (i in 1:length(DaSamples)){
  tmDaS = das_bind[das_bind$daSample == DaSamples[i],]
  tmDaS = t(table(tmDaS$immuneClass, tmDaS$TILClass))
  das_bind_list[[i]] <- tmDaS
}

pSummDas = data.frame(rbind(das_bind_list[[1]], das_bind_list[[2]], das_bind_list[[3]], das_bind_list[[4]], das_bind_list[[5]],
                            das_bind_list[[6]], das_bind_list[[7]], das_bind_list[[8]], das_bind_list[[9]], das_bind_list[[10]],
                            das_bind_list[[11]], das_bind_list[[12]], das_bind_list[[13]], das_bind_list[[14]], das_bind_list[[15]],
                            das_bind_list[[16]], das_bind_list[[17]], das_bind_list[[18]], das_bind_list[[19]], das_bind_list[[20]]))
pSummDas = pSummDas/rowSums(pSummDas)*100
pSummDas$TILClass = rep(c("ATL", "DTL", "ITL"))
pSummDas$cd8foxp3_Ratio = pSummDas$cd8/pSummDas$p3
pSummDas$cd8foxp3_Ratio[pSummDas$p3==0]<-0
comList = list(c("ATL", "DTL"), c("DTL", "ITL"), c("ATL", "ITL"))
p1<-ggboxplot(pSummDas, x = "TILClass", y = "cd8", color = "TILClass", palette = c("blue", "#8291F7", "black"),
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "t.test", paired = T)+
  theme(legend.position="") + theme(text = element_text(size=18))
p2<-ggboxplot(pSummDas, x = "TILClass", y = "cd8foxp3_Ratio", color = "TILClass", palette = c("blue", "#8291F7", "black"),
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "t.test", paired = T)+
  theme(legend.position="") + theme(text = element_text(size=18))
p3<-ggboxplot(pSummDas, x = "TILClass", y = "cd4", color = "TILClass", palette = c("blue", "#8291F7", "black"),
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "t.test", paired = T)+
  theme(legend.position="") + theme(text = element_text(size=18))
p4<-ggboxplot(pSummDas, x = "TILClass", y = "cd8Reg", color = "TILClass", palette = c("blue", "#8291F7", "black"),
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "t.test", paired = T)+
  theme(legend.position="") + theme(text = element_text(size=18))
p5<-ggboxplot(pSummDas, x = "TILClass", y = "p3", color = "TILClass", palette = c("blue", "#8291F7", "black"),
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "t.test", paired = T)+
  theme(legend.position="") + theme(text = element_text(size=18))

pdf(file=paste0(path, "/DaTiles_immuneVsTILs_percentage.pdf"), height = 6, width = 16)
grid.arrange(p1, p2, p3, p4, p5, ncol=5)
dev.off()
########################################################################################################################

##fig 10b, TCGA LUAD corrs of RNAseq signatures vs spatial TILs 
########################################################################################################################

#load TCGA TILs data
load("./results/TILs/LUAD_TILClass.rdata")
tilLUAD = data.frame(mat.l)
tilLUAD$file_name <- rownames(tilLUAD)
luad <- read.csv(file = "./dev/immuneScoringPaper/LUADSummary_AllSlides.csv", stringsAsFactors = F)
luad$PATIENT_SCREENING_ID <- substr(luad$file_name, 1, 12)
luad = luad[, c(3, 12:14, 16:18, 19)]
luad = merge(luad, tilLUAD, all.x=T)
names(luad)[names(luad) == "PATIENT_SCREENING_ID"] <- "bcr_patient_barcode"
luad = luad[, -c(1)]

luadM = luad %>% 
  group_by(bcr_patient_barcode) %>%
  summarise_all("mean")

LUADmean = read.csv("./dev/immuneScoringPaper/results/LUADmean.csv")
LUADclin = LUADmean[, c(8:10, 13, 22:26)]
luadM= merge(luadM, LUADclin, all.y=T)

luadM$ATL_fibroRatio = luadM$Adjacent/luadM$stromal
luadM$ITLR = luadM$Intra/luadM$tumour

#summary of all immune gene signatures for TCGA luad samples
load("./data/geneSignatures.RData")
LUADimmune <- Reduce(function(x, y) merge(x, y, by ="bcr_patient_barcode", all.x=TRUE), list(luadM, estLUAD2, BoLi2, Davoli))
toCorr <- LUADimmune[, c(8:10, 28, 31, 40:42)]
p.r <- cor_pmat(toCorr, method=c("spearman"), use="complete.obs")
#p val correction with BH
p.r <- matrix(p.adjust(p.r, method = "BH", n = length(p.r)), nrow = 8, ncol = 8) ##
corR <- round(cor(toCorr, method = c("spearman"), use="complete.obs"),2)

pdf(file=paste0(path, "/TILsvsRNASignatures_BH_.pdf"), width = 8, height = 6)
ggcorrplotModified(corR[4:8, 1:3], lab = T, sig.level=0.05, insig ="blank", lab.notsig="", 
                   colors = c("#4DBBD5B2","grey98","#DC0000B2"), legend.title = "Spearman Corr",
                   p.mat = p.r[4:8, 1:3], title = "TCGA LUAD (n=464)")
dev.off()
########################################################################################################################

##fig 5 + fig S11, fractal dim analysis + corrs
########################################################################################################################

#
tx = read.csv("./tracerx100/results/trait/100xCombined_4.csv", stringsAsFactors = F)
dFrac <- read.csv("./tracerx100/results/tumorFibro_frac/diagnostic/fracBoxcount-171130.csv", stringsAsFactors = F)
diagnostic = read.csv("./tracerx100/results/trait/100xDiagnostic.csv", stringsAsFactors = F)
dFrac$file_name <- substr(dFrac$file_name,1,nchar(dFrac$file_name)-12)
diagnostic <- merge(diagnostic, dFrac, all.x = T)
dF <- diagnostic[, c(3, 60:66)]
tx <- merge(tx, dF, all.x = T)

regional = read.csv("./tracerx100/results/trait/100xRegional_TGene.csv", stringsAsFactors = F )
rFrac <- read.csv("./tracerx100/results/tumorFibro_frac/regional/controlBoxSize/fracBoxcount49-180105.csv", 
                  stringsAsFactors = F)
colnames(rFrac)[which(names(rFrac) == "FileName")] <- "file_name"
colnames(rFrac)[which(names(rFrac) == "std")] <- "fd_std"
regional <- merge(regional, rFrac, all.x = T)
regional <- regional[regional$immuneClass_2=="hot" | regional$immuneClass_2=="cold",] 
comList <- list( c("hot", "cold"))
regional <- regional[order(regional$immuneClass_2),] 
rLUAD <- regional[regional$Histology=="Invasive adenocarcinoma",]
rLUSC <- regional[regional$Histology=="Squamous cell carcinoma",]

#5c boxes
p1<-ggboxplot(regional, x = "immuneClass_2", y = "fd", color = "immuneClass_2", title = "ALL, n=219",
              add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))
p2<-ggboxplot(rLUAD, x = "immuneClass_2", y = "fd", color = "immuneClass_2", title = "LUAD, n=113", 
              add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))
p3<-ggboxplot(rLUSC, x = "immuneClass_2", y = "fd", color = "immuneClass_2", title = "LUSC, n=84", 
              add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))
p4<-ggboxplot(regional, x = "immuneClass_2", y = "fibroblasts_per", color = "immuneClass_2",
              add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))
p5<-ggboxplot(rLUAD, x = "immuneClass_2", y = "fibroblasts_per", color = "immuneClass_2",
              add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))
p6<-ggboxplot(rLUSC, x = "immuneClass_2", y = "fibroblasts_per", color = "immuneClass_2",
              add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(comparisons = comList, method = "wilcox")+
  theme(legend.position="") + theme(text = element_text(size=18))

pdf(file=paste0(path, "/ver3_frac_hist_controlBoxSize_hotcold.pdf"), height = 12, width = 16)
grid.arrange(p1, p2, p3, 
             p4, p5, p6, ncol=3)
dev.off()

#distribution of dist between a stromal to a tumor
con = 0.3225
#20x  : 0.3225 um/pixel or 3.10 pixels/um
distD = data.frame()
Dir <- './tracerx100/results/tumorFibro_frac/regional/dist'
dirtmp <- dir(Dir)
fflist <- grep('_dist.csv', dirtmp)
for (i in 1:length(dirtmp[fflist])){
  dist <- read.csv(paste0('./tracerx100/results/tumorFibro_frac/regional/dist/',dirtmp[fflist][i]))
  dist = dist[dist$min_dist_f_to_t>12,]
  dist$minDist_f_to_t_um = con*dist$min_dist_f_to_t
  dist$maxDist_f_to_t_um = con*dist$max_dist_f_to_t
  minM = mean(dist$minDist_f_to_t_um)
  maxM = mean(dist$maxDist_f_to_t_um)
  df <- data.frame(minDistMean_um = minM, maxDistMean_um = maxM)
  distD <- rbind(distD,df)
}
pdf(file=paste0(path, "/meanMinDist_ftot_allRs.pdf"))
ggplot(distD, aes(minDistMean_um)) +
  geom_density(color="darkblue", fill="pink") + theme_classic()
ggplot(distD, aes(minDistMean_um)) +
  geom_density(color="darkblue", fill="pink")+xlim(0,300) + theme_classic()
dev.off()

#corrs in S11b
p1 <- ggscatter(regional, x = "fibroblasts_per", y = "fd", show.legend.text = F, 
                xlab = "Fibroblasts%", ylab = "Fractal dimension", font.label = 4,
                shape = 19, size = 2,
                color = "Histology", palette = c("#98AED6", "#682A45", "black"),
                add.params = list(color = "grey50", fill = "azure3"))+ 
  stat_cor(aes(color = Histology), label.x = 50, method="spearman", show.legend = F)+
  theme(legend.position="none")
ggscatter(regional, x = "fibroblasts_per", y = "fd", cor.coef = T)
p2 <- ggscatter(regional, x = "tumour_per", y = "fd", show.legend.text = F, 
                xlab = "Tumour%", ylab = "Fractal dimension", font.label = 4,
                shape = 19, size = 2,
                color = "Histology", palette = c("#98AED6", "#682A45", "black"),
                add.params = list(color = "grey50", fill = "azure3"))+ 
  stat_cor(aes(color = Histology), label.x = 10, method="spearman", show.legend = F)+
  theme(legend.position="none")
ggscatter(regional, x = "tumour_per", y = "fd", cor.coef = T)
pdf(file=paste0(path, "/fig4Corrs.pdf"), width = 10 , height = 16)
grid.arrange(p1,p2,  ncol=2)
dev.off()

#repeat almost the same for TCGA samples
#fig S11c
luad <- read.csv(file = "X:/Dropbox (ICR)/yuanlab/Projects/lung/tcga/dev/immuneScoringPaper/LUADSummary_AllSlides.csv", stringsAsFactors = F)
lusc <- read.csv(file = "X:/Dropbox (ICR)/yuanlab/Projects/lung/tcga/dev/immuneScoringPaper/LUSCSummary_AllSlides.csv", stringsAsFactors = F)
load("X:/Dropbox (ICR)/yuanlab/Projects/lung/tcga/trait/TCGAmin.RData")
LUSCmean = read.csv("X:/Dropbox (ICR)/yuanlab/Projects/lung/tcga/dev/immuneScoringPaper/results/LUSCmean.csv")
LUADmean = read.csv("X:/Dropbox (ICR)/yuanlab/Projects/lung/tcga/dev/immuneScoringPaper/results/LUADmean.csv")
clinLUAD = LUADmean[, c(8:10, 13:26)]
clinLUSC = LUSCmean[, c(8:10, 13:26)]
fracLUAD = read.csv("X:/Dropbox (ICR)/yuanlab/Projects/lung/tcga/results/tumorFibro_frac/fracBoxcount49_LUAD.csv")
fracLUSC= read.csv("X:/Dropbox (ICR)/yuanlab/Projects/lung/tcga/results/tumorFibro_frac/fracBoxcount49_LUSC.csv")
fracLUAD$bcr_patient_barcode = substr(fracLUAD$file_name ,1, 12)
fracLUSC$bcr_patient_barcode = substr(fracLUSC$file_name ,1, 12)

fracLUADsum = fracLUAD[, c(2,4)] %>% group_by(bcr_patient_barcode) %>%
  summarise_all(funs(fd_min = min, 
                     fd_median = median, 
                     fd_max = max,
                     fd_mean = mean))
fracLUSCsum = fracLUSC[, c(2,4)] %>% group_by(bcr_patient_barcode) %>%
  summarise_all(funs(fd_min = min, 
                     fd_median = median, 
                     fd_max = max,
                     fd_mean = mean))
fracLUSCsum2 = colNA(fracLUSCsum, "fd_max")
fracLUADsum = merge(clinLUAD, fracLUADsum, all.x = T)
fracLUSCsum = merge(clinLUSC, fracLUSCsum, all.x = T)
tcgaBoth = rbind(fracLUADsum, fracLUSCsum2)
tcgaBoth1 = rbind(LUADmean, LUSCmean)

trait = read.csv(file = "X:/Dropbox (ICR)/yuanlab/Projects/lung/tcga/trait/PanTCGA_immuneClassV3.csv", stringsAsFactors = F)
tcgaBoth = merge(tcgaBoth, trait, all.x = T)
tcgaBoth = merge(tcgaBoth, tcgaBoth1[, c(8,2)], all.x = T)
tcgaBoth = tcgaBoth[! tcgaBoth$immuneClass_2=="intermediate",]

p1<-ggboxplot(tcgaBoth[tcgaBoth$disease_code=="LUAD",], x = "immuneClass_2", y = "fd_max", color = "immuneClass_2", title = "LUAD", 
              add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",method = "t.test", ref.group = ".all.", hide.ns = T)+
  theme(legend.position="") + theme(text = element_text(size=18))

p2<-ggboxplot(tcgaBoth[tcgaBoth$disease_code=="LUSC",], x = "immuneClass_2", y = "fd_max", color = "immuneClass_2", title = "LUSC", 
              add = "jitter", border = "white",palette = c("blue2", "red2"))+
  stat_compare_means(label = "p.signif",method = "t.test", ref.group = ".all.", hide.ns = T)+
  theme(legend.position="") + theme(text = element_text(size=18))

pdf(file=paste0(path, "/ver3_frac_hist_controlBoxSize.pdf"), height = 20, width = 12)
grid.arrange(p1, p2, ncol=2)
dev.off()
########################################################################################################################

##fig 5d + fig S12, fractal vs LOH for each HLA type
########################################################################################################################
tx = read.csv("./tracerx100/results/trait/100xCombined_5.csv", stringsAsFactors = F)
tx$Histology[tx$Histology == "Squamous cell carcinoma"] <- "Squamous-cell carcinoma"
tx$Histology[tx$Histology == "Invasive adenocarcinoma"] <- "Adenocarcinoma"
tx = colNA(tx, "any.HLA.loss")
txLUAD <- tx[tx$Histology == "Adenocarcinoma" , ]
txLUSC <- tx[tx$Histology == "Squamous-cell carcinoma", ]

txLUAD = txLUAD[order(txLUAD$any.HLA.B.loss),]
txLUAD = txLUAD[order(txLUAD$any.HLA.C.loss),]

txLUSC = txLUSC[order(txLUSC$any.HLA.B.loss),]
txLUSC = txLUSC[order(txLUSC$any.HLA.C.loss),]

comList <- list( c("Loss", "No loss"))
p1<-ggboxplot(txLUSC, x = "any.HLA.loss", y = "fd_max", title = "LUSC (n=29)", 
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "t.test")+ 
  theme(text = element_text(size=18))
p2<-ggboxplot(txLUSC, x = "any.HLA.A.loss", y = "fd_max", title = "LUSC (n=29)",
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "t.test")+ 
  theme(text = element_text(size=18))
p3<-ggboxplot(txLUSC, x = "any.HLA.B.loss", y = "fd_max", title = "LUSC (n=29)",
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "t.test")+ 
  theme(text = element_text(size=18))
p4<-ggboxplot(txLUSC, x = "any.HLA.C.loss", y = "fd_max", title = "LUSC (n=29)",
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "t.test")+ 
  theme(text = element_text(size=18))
p5<-ggboxplot(txLUAD, x = "any.HLA.loss", y = "fd_max", title = "LUAD (n=48)",
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "t.test")+ 
  theme(text = element_text(size=18))
p6<-ggboxplot(txLUAD, x = "any.HLA.A.loss", y = "fd_max", title = "LUAD (n=48)",
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "t.test")+ 
  theme(text = element_text(size=18))
p7<-ggboxplot(txLUAD, x = "any.HLA.B.loss", y = "fd_max", title = "LUAD (n=48)",
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "t.test")+ 
  theme(text = element_text(size=18))
p8<-ggboxplot(txLUAD, x = "any.HLA.C.loss", y = "fd_max", title = "LUAD (n=48)",
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "t.test")+ 
  theme(text = element_text(size=18))
pdf(file=paste0(path, "/patientHLA_LUADLUSC_fd_max.pdf"), height = 10, width = 16)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8,  ncol=4)
dev.off()
########################################################################################################################

##fig S13, distribution of lym% between luad and lusc tx and tcga tumors
########################################################################################################################
lymTCGA = read.csv("X:/Dropbox (ICR)/yuanlab/Projects/lung/tcga/trait/PanTCGA_immuneClassV3.csv", stringsAsFactors = F)
comList <- list( c("LUAD", "LUSC"))
lymTCGA = lymTCGA[lymTCGA$disease_code=="LUAD" | lymTCGA$disease_code=="LUSC" ,]
lymTCGA$lym_per = lymTCGA$lym_per*100
lymTCGA = lymTCGA[order(lymTCGA$disease_code),]
p=ggdensity(lymTCGA, x = "lym_per", color = "disease_code", palette = c("#98AED6", "#682A45"),
            add = "mean", rug = TRUE)
p2<-ggboxplot(lymTCGA, x = "disease_code", y = "lym_per", color = "disease_code", title = "TCGA (n=939)", 
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "t.test")+
  theme(legend.position="") + theme(text = element_text(size=18))+
  scale_color_manual(values=c("#98AED6", "#682A45", "orange", "#00AFBB"))

regional = read.csv("./tracerx100/results/trait/100xRegional_TGene.csv", stringsAsFactors = F )
regional$Histology[regional$Histology == "Squamous cell carcinoma"] <- "LUSC"
regional$Histology[regional$Histology == "Invasive adenocarcinoma"] <- "LUAD"
regional$Histology[regional$Histology=="Carcinosarcoma"] <- "Other"
regional$Histology[regional$Histology=="Adenosquamous carcinoma"] <- "Other"
regional$Histology[regional$Histology=="Large Cell Neuroendocrine"] <- "Other"
regional = regional[order(regional$Histology),]
p3=ggdensity(regional, x = "lymphocytes_per", color = "Histology", palette = c("#98AED6", "#682A45", "black"),
             add = "mean", rug = TRUE)
p4<-ggboxplot(regional, x = "Histology", y = "lymphocytes_per", color = "Histology", title = "TRACERx (n=275)", 
              add = "jitter", border = "white")+
  stat_compare_means(comparisons = comList, method = "t.test")+
  theme(legend.position="") + theme(text = element_text(size=18))+
  scale_color_manual(values=c("#98AED6", "#682A45", "black"))

pdf(file=("./results/Pan_lymInfiltration_luadlusc.pdf"), height = 10, width = 10)
grid.arrange(p4, p3, p2, p, ncol=2)
dev.off()
########################################################################################################################

##fig 6a, summary statistics heatmap
########################################################################################################################
load("./tracerx100/results/trait/trait5_08Oct2018.RData")
tx$Histology[tx$histology_group == "Squamous cell carcinoma"] <- "Squamous-cell carcinoma"
tx$Histology[tx$histology_group == "Adenocarcinoma"] <- "Adenocarcinoma"
tx$Histology[tx$histology_group=="Pleomorphic carcinoma"] <- "Z"
tx$Histology[tx$histology_group=="Adenosquamous carcinoma"] <- "Z"
tx$Histology[tx$histology_group=="Large Cell Neuroendocrine Carcinoma"] <- "Z"
tx$Histology[tx$histology_group=="sarcomatoid carcinoma of pleomorphic type arising from adenocarcinoma"] <- "Z"
regional$Histology[regional$Histology == "Squamous cell carcinoma"] <- "Squamous-cell carcinoma"
regional$Histology[regional$Histology == "Invasive adenocarcinoma"] <- "Adenocarcinoma"
regional$Histology[regional$Histology=="Carcinosarcoma"] <- "Z"
regional$Histology[regional$Histology=="Adenosquamous carcinoma"] <- "Z"
regional$Histology[regional$Histology=="Large Cell Neuroendocrine"] <- "Z"
regional = regional[order(regional$lymphocytes_per, decreasing = F),]
regional = regional[order(regional$Histology, decreasing = F),]
#normal
rub=merge(tx, normal, all.x = T)
rub = colNA(rub, "normal")
rub$normalDiff = rub$regional_mean-rub$normal
rub$LQ.normalMeanDiff[rub$normalDiff <= quantile(rub$normalDiff)[[2]]] <- "Low"
rub$LQ.normalMeanDiff[rub$normalDiff > quantile(rub$normalDiff)[[2]]] <- "High"
rub$LQ.normalMeanDiff = as.factor(as.character(rub$LQ.normalMeanDiff))
tx = merge(tx, rub[, c(1,115)], all.x = T)
#set
tx = tx[order(tx$lymphocytes_per_min, decreasing = F),]
tx = tx[order(tx$Histology, decreasing = F),]
tx$DFS_yes[tx$DFS_censor_variable==1] <- tx$DFS_time_days
tx$DFS_no[tx$DFS_censor_variable==0] <- tx$DFS_time_days
p1<-ggdotchart(regional, x = "PublicationID", y = "lymphocytes_per",
               color = "immuneClass_2",                                # Color by groups.
               group = "PublicationID",
               sort.by.groups = TRUE,           
               palette = c("blue2", "red2", "#f4d5d5"), # Custom color palette
               #sorting = "descending",                       # Sort value in descending order
               #rotate = TRUE,                                # Rotate vertically
               #add = "segments",                             # Add segments from y = 0 to dots
               add.params = list(color = "lightgray", size = 0.8),
               dot.size = 3,                                 # Large dot size
               #label = min(regionalT$lymphocytes_per),
               #y.text.col = TRUE,                            # Color y text by groups
               ggtheme = theme_pubr())

#now the heatmap annotations for all immune scores
pdf(file=paste0(path, "/mainHeatmap_ADDrows.pdf"), width = 9 , height = 21.5)

Heatmap(tx$Histology, name = "Histology", width = unit(2.85, "mm"), 
        cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
        col = c("#98AED6", "#682A45", "black")) +
  Heatmap(tx$FD, name = "FD", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = c("greenyellow", "hotpink"))+
  Heatmap(tx$LQ.normalMeanDiff, name = "TNDive", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = c("#C40013", "#0080FF", "gray45"))+
  Heatmap(tx$cd8_per, name = "CD8", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = colorRamp2(c(0, 17), c("#F2E2E2", "#8E0000")))+
  Heatmap(tx$cd4_per, name = "CD4", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = colorRamp2(c(0, 28), c("#E3F0E2", "#067C00")))+
  Heatmap(tx$foxp3_per, name = "FOXP3", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = colorRamp2(c(0, 21), c("#F9E2F2", "#CC008E")))+
  Heatmap(tx$ATL_fibroRatio, name = "ATLF", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = colorRamp2(c(0, 1.6), c("#EDE2EF", "#650077")))+
  
  
  Heatmap(tx$ClonalNeo, name = "cNEO", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          colorRamp2(c(0, 800), 
                     c("#E3E2EF", "#030070")))+
  Heatmap(tx$SubclonalNeo, name = "subNEO", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          colorRamp2(c(0, 400), 
                     c("#E3E2EF", "#030070")))+
  Heatmap(tx$any.HLA.loss, name = "HLA", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = c("red1", "#4DBBD5B2", "gray45"))+
  
  
  Heatmap(tx$pathologyTNM, name = "Stage", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = c("#7FC97F" ,"#BEAED4" ,"#FDC086", "#FFFF99" ,"#386CB0", "#F0027F"))+
  Heatmap(tx$pack_years_calculated, name = "Pack years", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = colorRamp2(c(0, 140), 
                           c("#E2E6E7", "#00252D")))+
  Heatmap(tx$DFS_yes, name = "DFS_Yes", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = colorRamp2(c(0, 1364.0, NA), 
                           c("#FBF0E2", "#DD7800", "gray45")))+
  Heatmap(tx$DFS_no, name = "DFS_No", width = unit(2.85, "mm"), 
          cluster_rows = F, cluster_columns = F, gap = unit(2, "mm"),
          col = colorRamp2(c(54.00, 1364.00, NA ), 
                           c("#E2EAE2", "#024402", "gray45")))
dev.off()
########################################################################################################################

diagnosticLUAD <- diagnostic[diagnostic$histology_group=="Adenocarcinoma",]
diagnosticLUSC <- diagnostic[diagnostic$histology_group=="Squamous cell carcinoma",]
geneticLUAD <- genetic[genetic$histology_group=="Adenocarcinoma",]
geneticLUSC <- genetic[genetic$histology_group=="Squamous cell carcinoma",]
diagnosticIHC = colNA(diagnostic, "cd8_per")
diagnosticIHC$cd8foxp3Ratio = diagnosticIHC$cd8_per/diagnosticIHC$foxp3_per
diagnosticIHCLUAD <- diagnosticIHC[diagnosticIHC$histology_group=="Adenocarcinoma",]
diagnosticIHCLUSC <- diagnosticIHC[diagnosticIHC$histology_group=="Squamous cell carcinoma",]
tx$coldP = as.factor(as.character(tx$coldP))
tx$coldP3 = as.factor(as.character(tx$coldP3))
tx$pathologyTNM = as.factor(as.character(tx$pathologyTNM))
txH2 = tx
txH2$Histology[txH2$histology_group=="Adenocarcinoma" ] <- "LUAD"
txH2$Histology[txH2$histology_group=="Squamous cell carcinoma" ] <- "LUSC"
txH2$Histology[is.na(txH2$Histology)] <- "Other"
txH2$Histology = as.factor(as.character(txH2$Histology))
txH2$sex = as.factor(as.character(txH2$sex))
txH2$pathologyTNM2=txH2$pathologyTNM
levels(txH2$pathologyTNM2)[levels(txH2$pathologyTNM2)=="IIIA"] <- "III"
levels(txH2$pathologyTNM2)[levels(txH2$pathologyTNM2)=="IIIB"] <- "III"
txH2$AdjuvantTherapy = txH2$adjuvant_treatment_given
txH2$AdjuvantTherapyC[txH2$adjuvant_treatment_given == ""]<- 1
txH2$AdjuvantTherapyC[txH2$adjuvant_treatment_given == "Platinum chemo"]<- 0
txH2$AdjuvantTherapyC[txH2$adjuvant_treatment_given == "Radiotherapy"]<- 0
txH2$AdjuvantTherapyC[txH2$adjuvant_treatment_given == "Platinum chemo/radiotherapy"]<- 0

txH2$CD8FOXP3Ratio <- txH2$cd8_per/txH2$foxp3_per
txH2$CD4FOXP3Ratio <- txH2$cd4_per/txH2$foxp3_per

txH2LUADLUSC <- txH2[ which(txH2$Histology=="LUAD" 
                            | txH2$Histology == "LUSC"), ]
txH2LUADLUSC$Histology <- factor(txH2LUADLUSC$Histology)

num = readRDS("./tracerx100/data/number_regions.RDS")
num$REGTrialNo = rownames(num)
names(num) <- c("nTotalRegions", "nRNARegions", "PATIENT_SCREENING_ID")
num$PATIENT_SCREENING_ID = sub("(.{3})(.*)", "\\10\\2", num$PATIENT_SCREENING_ID)
txH2LUADLUSC = merge(txH2LUADLUSC, num, all.x = T)
txH2 = merge(txH2, num, all.x = T)
