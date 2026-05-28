path <- "C:/Users/pvgam/Documents/CASB185" #change as needed!

setwd(path)
install.packages("tidyverse")
install.packages("dplyr")
install.packages("janitor")
install.packages("edgeR")

library(tidyverse)
library(limma)
library(dplyr)
library(janitor)

install.packages("BiocManager")
BiocManager::install("edgeR")
library(edgeR)

excDF <- read.delim(
  "merged.tpm.Exc.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

inhDF <- read.delim(
  "merged.tpm.Inh (1).txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)


metaData <- read.delim(
  "metadata_snRNAseq_SEA-AD_DLPFC_with_mapping.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

controlMeta <- metaData %>%
  filter(`Consensus clinical diagnosis` == "Control")

ALZMeta <- metaData %>% 
  filter(`Consensus clinical diagnosis` != "Control")

control_samples <- controlMeta$individualID
alz_samples <- ALZMeta$individualID

processDF <- function(DF,controlNames,ALZNames){
  
  gene_names <- DF$gene_name
  expr <- DF[, -1]
  expr <- log10(expr + 1)
  keep <- rowSums(expr > 1) >= 3
  expr <- expr[keep, ]
  expr <- cbind(
    gene_name = gene_names[keep],
    expr
  )
  
  return(expr)
}

CC_exc <- excDF %>% select(gene_name, any_of(control_samples))
CC_exc <- processDF(CC_exc,control_samples,alz_samples) 
ALZ_exc <- excDF %>% select(gene_name, any_of(alz_samples)) 
ALZ_exc <- processDF(ALZ_exc,control_samples,alz_samples) 
CC_inh <- inhDF %>% select(gene_name, any_of(control_samples))
CC_inh <- processDF(CC_inh,control_samples,alz_samples) 
ALZ_inh <- inhDF %>% select(gene_name, any_of(alz_samples)) 
ALZ_inh <- processDF(ALZ_inh,control_samples,alz_samples)



run_DEG <- function(
    df1,
    df2,
    group1_name,
    group2_name,
    min_expr = 1,
    min_samples = 3,
    logFC_cutoff = 1,
    adjP_cutoff = 0.05
) {
  
  
  #Limma uniquely requires gene names to be a "row name" instead of a column, this function ensures that, as well as that the entire matrix is numerical
  prep_matrix <- function(df) {
    
    df <- as.data.frame(df)
    rownames(df) <- df$gene_name
    df$gene_name <- NULL
    mat <- as.matrix(df)
    storage.mode(mat) <- "numeric"
    return(mat)
  }
  
  #All results will be the change from MAT 2 --> MAT 1
  mat1 <- prep_matrix(df1)
  mat2 <- prep_matrix(df2)
  
  
  #Some genes were filtered out due to minimal expression in one or both datasets, so just align the genes in each matrix here!
  common_genes <- intersect(
    rownames(mat1),
    rownames(mat2)
  )
  mat1 <- mat1[
    common_genes,
    ,drop = FALSE]
  mat2 <- mat2[
    common_genes,
    , drop = FALSE]
  
  #After sucessfull filtering, combine them here!
  COMBINED <- cbind(mat1, mat2)
  

  #Remove very low expression genes
  keep <- rowSums(COMBINED > min_expr) >= min_samples
  COMBINED <- COMBINED[keep, ]
  
  #Remove housekeeping genes with no variance
  gene_var <- apply(COMBINED, 1, var)
  COMBINED <- COMBINED[gene_var > 0, ]
  
  # Format matrix for Limma, here we denote each sample as belonging to group identity 1 or 2. The beauty of Limma is that it can perform high-fidelity DEG across multiple samples in a group
  group <- factor(c(
    rep(group1_name, ncol(mat1)),
    rep(group2_name, ncol(mat2))
  ))
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  
  # Run the FIT model
  fit <- lmFit(COMBINED, design)
  
  #This is the formula for Fold Change
  contrast_formula <- paste0(
    group2_name,
    " - ",
    group1_name
  )
  
  #New Matrix is created w/ fold change/DEG contrasts
  contrast.matrix <- makeContrasts(
    contrasts = contrast_formula,
    levels = design
  )
  
  
  #This is where the statistical testing occurs and we get the relevant p-adjusted values and parameters of linear fit (Bayesian/Bejamini-Hochberg)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
 #Format Results Table
  results <- topTable(
    fit2,
    number = Inf,
    adjust.method = "BH"
  )
  
 #Select for DEGs at alpha threshold 
  sig_genes <- results %>%
    filter(
      adj.P.Val < adjP_cutoff,
      abs(logFC) > logFC_cutoff
    )
  
  #Volcano Plot!
  volcanoplot(
    fit2,
    coef = 1,
    highlight = 20,
    main = paste(group2_name, "vs", group1_name)
  )
  
  hist(
    results$adj.P.Val,
    breaks = 50,
    main = "Adjusted P-Value Distribution",
    xlab = "Adjusted P-Value",
    col = "skyblue",
    border = "black"
  )  
  #Return in case we want to use these further outside the method. 
  
  return(list(
    all_results = results,
    significant_genes = sig_genes,
    fit = fit2,
    combined_matrix = COMBINED,
    design_matrix = design
  ))
}

results1 <- run_DEG(ALZ_exc,ALZ_inh,"Excitatory_Neuron_ALZ","Inhibitory_Neuron_ALZ")
results2 <- run_DEG(CC_exc, CC_inh, "Excitatory_Neuron_CC","Inhibitory_Neuron_CC")

