library(knitr)
library(limma)
library(minfi)
library(RColorBrewer)
library(stringr)
library(wateRmelon)
library(future)
library(gtools)
library(matrixStats)
library(data.table)
library(MASS)
library(sandwich)
library(lmtest)
library(parallel)
library(dplyr)
library(tidyverse)
library(compare)
library("tidyr")

rm(list=ls())
beta2 <- readRDS("/vol/projects/BIIM/2000HIV/AnalysisDV/Discovery_MethylationDataBased/1qc6/2000HIV.Mvalue.rds")
beta3 <- beta2[ ,-c(318,325,328)]

dis <- readRDS("/vol/projects/BIIM/2000HIV/AnalysisDV/Discovery_MethylationDataBased/2Pheno/Phe_dis_clean_PCs.rds")
ids <- dis$ID
phe <- read.csv("/vol/projects/BIIM/2000HIV/Phenotype/Phenotype_2000HIV_all_01.tsv", sep="\t")
colnames(phe)[1] = "ID"
phe <- dplyr::filter(phe, phe$ID%in%ids)

phe <- phe[, c("ID", "Cohort", "SEX_BIRTH", "AGE", "BMI_BASELINE", "DNAm_SamplePlate",
               "CMV_IgG_Serology", "CMV_IgG_IU.mL", "ETHNICITY", "season_sin", "season_cos")]
#phe <- dplyr::filter(phe, phe$ETHNICITY=="White")

cellCount <- dis %>%
    dplyr::select("ID", "NK", "Mono", "Neu", "CD4T", "CD8T", "Bcell", "PC1", "PC2", "PC3", "PC4", "PC5")

phe <- merge(cellCount, phe, by="ID")

pheno5 <- dplyr::filter(phe, !is.na(phe$CMV_IgG_Serology))
idx <- intersect(colnames(beta3), pheno5$ID)
beta3 <- beta3[, idx]
idx2 <- match(colnames(beta3), pheno5$ID)
pheno6 <- pheno5[idx2, ]

compare(colnames(beta3), pheno6$ID)
cat(dim(pheno6), dim(beta3))

EWAS_pheno=data.frame()
for (i in 1:length(pheno6$ID)) {
    id = pheno6$ID[i]
    #cat(i, id, "|")
    EWAS_pheno[i, "ID_idat"] = id
    INR = dplyr::filter(pheno6, pheno6$ID==id)[, "CMV_IgG_Serology"] 
    if (is.na(INR)) {
        EWAS_pheno[i, "Control"] = NA
    } else {
        INR = as.numeric(INR)
        EWAS_pheno[i, "Control"] = INR
    }
}

EWAS_pheno$Control <- as.factor(EWAS_pheno$Control)
ix <- (which(is.na(EWAS_pheno$Control)))
if (length(ix)>0) {
    beta3 <- beta3[ ,-ix]
    EWAS_pheno <- EWAS_pheno[-ix, ]
    pheno6 <- pheno6[-ix, ]
}

cat("Adter removing the NA individuals. Do the Pheno5 and EWAS_Pheno have the same ID_idat: ")
compare(pheno6$ID, EWAS_pheno$ID_idat)
cat("Adter removing the NA individuals. Do the beta3's colnames are the same as ID_idat in pheno5 and EWAS_pheno: ")
compare(colnames(beta3), EWAS_pheno$ID_idat)
cat("Adter removing the NA individuals. The dimentions of EWAS_pheno, beta3 and pheno5: \n")
cat(dim(EWAS_pheno), dim(beta3), dim(pheno6))

removeOutliers<-function(probes){
  require(matrixStats)
  if(nrow(probes) < ncol(probes)) warning("expecting probes are rows (long dataset)")
  rowIQR <- rowIQRs(probes, na.rm = T)
  row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = T)
  maskL <- probes < row2575[,1] - 3 * rowIQR 
  maskU <- probes > row2575[,2] + 3 * rowIQR 
  initial_NAs<-rowSums(is.na(probes))
  probes[maskL] <- NA
  removed_lower <- rowSums(is.na(probes))-initial_NAs
  probes[maskU] <- NA
  removed_upper <- rowSums(is.na(probes))-removed_lower-initial_NAs
  N_for_probe<-rowSums(!is.na(probes))
  Log<-data.frame(initial_NAs,removed_lower,removed_upper,N_for_probe)
  return(list(probes, Log))
}

#Remove outliers from METH (methylation data where probes are rows and samples are columns)
system.time(OutlierResults <- removeOutliers(beta3))  
METH.2 <- OutlierResults[[1]]
Log <- OutlierResults[[2]]
rm(beta2)
rm(beta3)
save(Log,file = "Outlier_log_All.Rdata") #save log

mVals.T0 = t(METH.2)
pheno.T0 = pheno6
TrImmRes.T0 = EWAS_pheno

genes <- unique(colnames(mVals.T0))#test remove [1:100]
df = mVals.T0 
df = df %>% as.data.frame()

cor_test=function(gene)
{ 
  sub <- df[, gene] %>% data.frame()  #df=normalized gene count matrix
  colnames(sub) <- 'methylation'
  test.df <- cbind(pheno.T0, sub)
  cyto1 <- TrImmRes.T0[,2]
  bad <- as.numeric(rep(NA, 4))
  names(bad) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  result <- bad  
  tryCatch(
      {
        ML = rlm(methylation ~ cyto1 + AGE + SEX_BIRTH + 
                 BMI_BASELINE + SamplePlate + season_sin + season_cos +
                 CD8T + CD4T + NK + Bcell + Mono + Neu + 
                 PC1 + PC2 + PC3 + PC4 + PC5,
                 data=test.df, maxit=200)
        cf <- try(coeftest(ML, vcov=vcovHC(ML, type="HC0")))
        x <<- x+1
        a <- sum(x:length(P_value_Index))
        #cat(x, ":",length(P_value_Index),".\n",sep="")
        if (class(cf)=="try-error") {
          bad <- as.numeric(rep(NA, 4))
          names(bad)<- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
          result <- bad
        }
        else{
          result <- cf[2, c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]
        }  
      },
      error=function(error_message) {
        message("And below is the error message from R:")
        message(error_message)
        return(result)
      }
    )
  return(result)
}

pheno.T0$SEX_BIRTH = as.factor(pheno.T0$SEX_BIRTH)
pheno.T0$SamplePlate = droplevels(as.factor(pheno.T0$DNAm_SamplePlate))
pheno.T0$AGE = as.numeric(pheno.T0$AGE)
pheno.T0$BMI_BASELINE = as.numeric(pheno.T0$BMI_BASELINE)
# pheno.T0$HIV_DURATION = as.numeric(pheno.T0$HIV_DURATION)
# pheno.T0$RISK_BEHAV_MSM = as.factor(pheno.T0$RISK_BEHAV_MSM)
# pheno.T0$RISK_BEHAV_Heterosexual = as.factor(pheno.T0$RISK_BEHAV_Heterosexual)
# pheno.T0$RISK_BEHAV_Unknown = as.factor(pheno.T0$RISK_BEHAV_Unknown)

L = list()
P_value_Index=genes
cor_info <- data.frame()
x <- 0
Result=vapply(P_value_Index, cor_test, numeric(4))
cor_info=as.data.frame(t(Result))
colnames(cor_info)=c("Estimate", "StdError", "z-score", "pval")
cor_info[,5]=rownames(cor_info)#cor_info[,4]=genes
colnames(cor_info)=c("Estimate", "StdError", "z-score", "pval", "CpGsite")
cor_info$FDR=p.adjust(cor_info$pval, method = "fdr")
saveRDS(cor_info,"EWAS_outputs/corEWAS_dis_CMV_season_BMI_SP_cor_rlm.rds",compress="xz")



