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
library(foreach)
library(doParallel)
library(dplyr)
library(tidyverse)
library(compare)
library("tidyr")

rm(list=ls())
beta2 <- readRDS("/vol/projects/BIIM/2000HIV/AnalysisDV/Validation_MethylationDataBased/1qc6/2000HIV.Mvalue.rds")
beta3 <- beta2

val <- readRDS("/vol/projects/BIIM/2000HIV/AnalysisDV/Validation_MethylationDataBased/2Pheno/Phe_Rep_Clean_GenotypePCs.rds")
ids <- val$ID
phe <- read.csv("/vol/projects/BIIM/2000HIV/Phenotype/Phenotype_2000HIV_all_01.tsv", sep="\t")
colnames(phe)[1] = "ID"
phe <- dplyr::filter(phe, phe$ID%in%ids)

phe <- phe[, c("ID", "Cohort", "SEX_BIRTH", "AGE", "BMI_BASELINE", "DNAm_SamplePlate",
               "CMV_IgG_Serology", "CMV_IgG_IU.mL", "ETHNICITY", "season_sin", "season_cos")]


cellCount <- val %>%
    dplyr::select("ID", "NK", "Mono", "Neu", "CD4T", "CD8T", "Bcell", "PC1", "PC2", "PC3", "PC4", "PC5")
head(cellCount)

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
  rowIQR <- rowIQRs(probes, na.rm = T, useNames = TRUE)
  row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = T, useNames = TRUE)
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
save(Log,file = "Outlier_log_rep.Rdata") #save log

mVals.T0 = t(METH.2)
pheno.T0 = pheno6
TrImmRes.T0 = EWAS_pheno

genes <- unique(colnames(mVals.T0))#test remove [1:100]
df = mVals.T0 
df = df %>% as.data.frame()
#gene = genes

cor_test <- function(gene, df, pheno.T0, TrImmRes.T0) {
  sub <- df[, gene, drop = FALSE]
  colnames(sub) <- 'methylation'
  test.df <- cbind(pheno.T0, sub)
  cyto1 <- TrImmRes.T0[, 2]
  
  tryCatch({
    ML <- rlm(methylation ~ cyto1 + AGE + SEX_BIRTH, data = test.df, maxit = 200)
    cf <- try(coeftest(ML, vcov = vcovHC(ML, type = "HC0")), silent = TRUE)
    
    if (class(cf) == "try-error") {
      result <- c(CpGsite = gene, Estimate = NA, `Std. Error` = NA, `z value` = NA, `Pr(>|z|)` = NA)
    } else {
      result <- c(CpGsite = gene, cf[2, c("Estimate", "Std. Error", "z value", "Pr(>|z|)")])
    }
  }, error = function(error_message) {
    message("Error in processing gene: ", gene)
    message(error_message)
    result <- c(CpGsite = gene, Estimate = NA, `Std. Error` = NA, `z value` = NA, `Pr(>|z|)` = NA)
  })
  
  return(result)
}

pheno.T0$SEX_BIRTH = as.factor(pheno.T0$SEX_BIRTH)
pheno.T0$SamplePlate = droplevels(as.factor(pheno.T0$DNAm_SamplePlate))
pheno.T0$AGE = as.numeric(pheno.T0$AGE)
pheno.T0$BMI_BASELINE = as.numeric(pheno.T0$BMI_BASELINE)


registerDoParallel(8)
Result <- foreach(gene=genes, .combine = rbind) %dopar% {
  tryCatch({
              cor_test(gene, df, pheno.T0, TrImmRes.T0)
           }, error = function(e) {
                message("Error on iteration ", i, ": ", e)
                NULL  # or return some default value
            })
}

colnames(Result)=c("CpGsite", "Estimate", "StdError", "z_score", "pval")
Result <- as.data.frame(Result)
Result$Estimate <- as.numeric(Result$Estimate)
Result$pval <- as.numeric(Result$pval)
Result$StdError <- as.numeric(Result$StdError)
Result$z_score <- as.numeric(Result$z_score)
Result$FDR=p.adjust(Result$pval, method = "fdr")

saveRDS(Result,"../EWAS_outputs/corEWAS_rep_CMV_AGE_SEX_rlm_doParallel.rds",compress="xz")