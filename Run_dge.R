#Function for Differential Peak regions
RunDGE <- function(dge, dfcondition){
  # Setting the design as factor
  dfcondition[["celltyp"]] <- as.factor(dfcondition[["celltyp"]])
  
  # Making empty dataframes to do save Padj, FoldChange and Pvals data for each 
  #Pairwaise comparison
  dfQval <- data.frame()
  dfFC <- data.frame()
  dfPval <- data.frame()
  
  # Calculating dge across all the pairwise combination
  for(i in seq_along(dfPairwiseCond)){
    comparision <- print(paste0(dfPairwiseCond[1,i], "vs", dfPairwiseCond[2,i]))
    
    #Save the results of DGE as dataframe for downstream analysis
    res <- results(dge, contrast = c("celltyp", dfPairwiseCond[1,i],dfPairwiseCond[2,i]), format = "DataFrame")
    # Add results to the respective data frames
    dfQval[rownames(res), paste0(comparision, "_padj")] <- res$padj
    dfFC[rownames(res), paste0(comparision, "_fc")] <- res$log2FoldChange
    dfPval[rownames(res), paste0(comparision, "_pval")] <- res$pvalue
  }
  # Return df
  return(list(qvals = dfQval, fc = dfFC, pvals = dfPval))
}