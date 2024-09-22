# library(BBMRIomics)
# library(bacon)
library(limma)
library(R.utils)
library(edgeR)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(FDb.InfiniumMethylation.hg19)
library(igraph)
library(ggrepel)
library(tidyverse)

rm(list=ls())
genomeFile <- "/vol/projects/BIIM/2000HIV/AnalysisDV/Discovery_MethylationDataBased/4EWAS_AllEthnicity/Genome_annotation/gencode.v45lift37.basic.annotation.gff3.gz"
cpgBedFile <- "/vol/projects/BIIM/2000HIV/AnalysisDV/Discovery_MethylationDataBased/4EWAS_AllEthnicity/CMV_rlm/Enrichment_analysis/BedFile_GreatWeb/validated_CMV_14536DMS.bed"
load("/vol/projects/BIIM/2000HIV/AnalysisDV/Discovery_MethylationDataBased/4EWAS_AllEthnicity/CMV_rlm/temp_data/CMVAssoDMSmValue_forRNAseq.RData")
GENEs <- genes # genes was used as a variable in this script. So assign it to another one
data$Sample_Plate <- droplevels(as.factor(data$DNAm_SamplePlate))
data$SEX_BIRTH <- as.factor(data$SEX_BIRTH)
data$Center <- as.factor(data$`Institute.Abbreviation`)

# annotated genome data
Genome <- read.table(genomeFile)
colnames(Genome) <- c("CHR", "V2", "annotation", "start", "end", "V6", "strand", "V8", "description")


# the position of the DMS
df_DMS <- read.table(cpgBedFile, header = FALSE, sep="\t") %>% 
            setNames(c("CHR", "start", "end", "CpGsite"))
df_DMS <- dplyr::filter(df_DMS, CpGsite%in%cpgs) # only check the validated DMS


for (i in 1:500){ #dim(df_DMS)[1]
    cpgChr <- df_DMS[i, "CHR"]
    cpgStart <- df_DMS[i, "start"]
    cpgEnd <- df_DMS[i, "end"]
    cpg <- df_DMS[i, "CpGsite"]
    upstream <- cpgStart - 250000
    downstream <- cpgEnd + 250000
    arr <- dplyr::filter(Genome, CHR==cpgChr & start>upstream & end<downstream)
    genes = list()
    for (description in arr$description){
        if (grepl("gene_name=", description)){
            gene <- str_split_1(description, "gene_name=")[2]
            gene <- str_split_1(gene, ";")[1] # get the gene name in the description column
            if (gene %in% GENEs){
            genes <- append(genes, gene)
            }
        }
    }
    genes <- unique(genes) %>% unlist()
    
    if (length(genes)>0) {
    
        dat = data[, c(GENEs, cpg, c("BMI_BASELINE", #"SMOKING", 
                            "AGE", "SEX_BIRTH", "Sample_Plate", "PC1", "PC2", "PC3", "PC4", "BMI_BASELINE",
                            "PC5", "CD8T", "CD4T", "Bcell", "Mono", "NK", "Neu", 
                            "season_sin", "season_cos", "Center", "ID"))]
        rownames(dat) <- dat$ID
        dat <- na.omit(dat)
        vars <- dat[, c(cpg, "BMI_BASELINE", #"SMOKING",
                            "AGE", "SEX_BIRTH", "Sample_Plate", "PC1", "PC2", "PC3", "PC4", "BMI_BASELINE",
                            "PC5", "CD8T", "CD4T", "Bcell", "Mono", "NK", "Neu", 
                            "season_sin", "season_cos", "Center")]
        design <- model.matrix(~ ., vars)
        counts. <- t(dat[, GENEs])
        counts. <- DGEList(counts.)
        counts. <- counts.[rowSums(counts.$counts > 0) > 0.8 * ncol(counts.), ] #filter: min 80% samples at least 1 count
        counts. <- calcNormFactors(counts.)
        counts. <- voom(counts., design)
        fit <- lmFit(counts., design)
        fit <- eBayes(fit)
        results <- limma::topTable(fit, coef=2, n=Inf)
        sig <-  dplyr::filter(results, results$`adj.P.Val`< 0.05)
        if (dim(sig)[1]!=0){
            sig$CpG <- cpg
            sig$Gene <- rownames(sig)
            rownames(sig) <- paste0(sig$CpG, "_", sig$Gene)
            if (i==1){
                total_result <- sig
            } else{
                total_result <- rbind(total_result, sig)
            }
        }
        results <- results[genes, ]
        results$CpG <- cpg
        results$Gene <- rownames(results)
        rownames(results) <- paste0(results$CpG, "_", results$Gene)
        if (i==1){
            cis_result <- results
        } else{
            cis_result <- rbind(cis_result, results)
        }
    } 
}

#Total_result$`adj.P.Val` <- p.adjust(Total_result$`P.Value`, method = "fdr")
write.table(total_result, file = "eQTM_CMV_2478validatedDMS_Limma.tsv00", row.names = FALSE, sep="\t")
write.table(cis_result, file = "cis-eQTM_CMV_2478validatedDMS_Limma.tsv00", row.names = FALSE, sep="\t")


