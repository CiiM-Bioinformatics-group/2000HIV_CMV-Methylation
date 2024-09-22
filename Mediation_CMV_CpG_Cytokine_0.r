library(tidyverse)
library(data.table)
library(mediation)
library(lmtest)
library(MASS)
rm(list=ls())
load("../temp_data/ValidatedCMVAssoDMSmValue_forCytokine.RData")
df$DNAm_SamplePlate <- droplevels(as.factor(df$DNAm_SamplePlate))
inormal <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
file0="CMV_validatedDMS_validateIRTCytokine_mediationTest0.tsv"
cytokines <- c("pbmc_7d_il10_pha","pbmc_7d_il22_spneu", "pbmc_7d_il22_mtb","pbmc_24h_il1ra_il1a", "pbmc_7d_il22_saureus", "pbmc_24h_il1b_cmv", "pbmc_24h_il1ra_cmv", "pbmc_24h_il8_cmv", "pbmc_24h_mcp1_cmv")
options(warn=-1)
cat("CpG\tCytokine\tES_ACME\tpval_ACME\tES_ADE\tpval_ADE\tES_ProMediated\tpval_ProMediated\tES_total\tpval_total\n", file=file0, sep="", append=F)
for (cpg in cpgs[1:250]) {
for (cytokine in cytokines){
        cat(cpg, "	", cytokine, file=file0, sep="", append=T)
        dat = df[, c(cpg, cytokine, "CMV_IgG_Serology", "AGE", "SEX_BIRTH", "BMI_BASELINE", "DNAm_SamplePlate", "Institute.Abbreviation", "season_sin", "season_cos")]
        colnames(dat) = c("Betavalue", "cytokine", "CMV_IgG_Serology", "AGE", "SEX_BIRTH", "BMI_BASELINE", "Sample_Plate", "Center", "season_sin", "season_cos")
        dat = na.omit(dat)
        dat$cytokine <- inormal(dat$cytokine)
        med_fit <- try(rlm(Betavalue ~ CMV_IgG_Serology + AGE + SEX_BIRTH + BMI_BASELINE +  Sample_Plate + Center + season_sin + season_cos, data=dat, maxit=200))
        out_fit <- try(rlm(cytokine ~ Betavalue + CMV_IgG_Serology + AGE + SEX_BIRTH + BMI_BASELINE + Sample_Plate + Center + season_sin + season_cos, data=dat, maxit=200))
        med_out <- try(mediate(med_fit, out_fit, treat="CMV_IgG_Serology",  mediator = "Betavalue", robustSE = T, sims=1000))


        if (class(med_fit)[1]=="try-error" | class(out_fit)[1]=="try-error" | class(med_out)=="try-error") {
            ES_ACME <- NA # ACME Estimate
            pval_ACME <- NA # ACME p-value
            ES_ADE <- NA # ADE Estimate 
            pval_ADE <- NA # ADE p-value
            ES_PropMediated <- NA # prop.mediated Estimate
            pval_PropMediated <- NA # prop.mediated p-value
            ES_Total <- NA # Total Effect Estimate
            pval_Total <- NA # Total Effect p-value
        } else{
            cf <- try(summary(med_out))
            ES_ACME <- cf$d0 # ACME Estimate
            pval_ACME <- cf$d0.p # ACME p-value
            ES_ADE <- cf$z0 # ADE Estimate 
            pval_ADE <- cf$z0.p # ADE p-value
            ES_PropMediated <- cf$n0 # prop.mediated Estimate
            pval_PropMediated <- cf$n0.p # prop.mediated p-value
            ES_Total <- cf$tau.coef # Total Effect Estimate
            pval_Total <- cf$tau.p # Total Effect p-value
        }
        cat("\t", file=file0, sep="", append=T)
        cat(ES_ACME, "\t", pval_ACME, "\t", ES_ADE, "\t", pval_ADE, "\t", ES_PropMediated, "\t", pval_PropMediated, "\t", ES_Total, "\t", pval_Total, "\n",
            file=file0, sep="", append=T)
    }
    
}