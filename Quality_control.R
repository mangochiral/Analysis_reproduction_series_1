#Function to run DESeq
RunDeseq <- function(df,dfcondition){
  dfcondition[["celltyp"]] <- as.factor(dfcondition[["celltyp"]]) # set the T cell type as factor
  dfcondition[["time"]] <- as.factor(dfcondition[["time"]]) # set the time as factor
  #Running DESeq
  dds <- DESeqDataSetFromMatrix(countData = df, colData = dfcondition, ~celltyp+time)
  deRNA <- DESeq(dds)
  return(deRNA)
}

# Plotting PCA function
RunPCA <- function(dfnorm,dfcondition){
  # To run PCA only regions or genes which contributed 80% of the variance between design is used
  dfNormFilter <- varFilter(as.matrix(dfnorm), var.cutoff = 0.2, filterByQuantile = T)
  
  # Data is scaled and centered for plotting PCA
  dfNormPCA <- prcomp(t(dfNormFilter), scale. = T, center = T)
  
  # The variance contribution is calculated for labelling contribution of variance by each PC
  percentVar <- ((dfNormPCA$sdev)^2 / sum(dfNormPCA$sdev^2)) * 100
  
  #Plotting PCA of PC1 and PC2
  p1 <- ggplot(dfNormPCA$x, aes(x= dfNormPCA$x[,1], y = dfNormPCA$x[, 2], colour = dfcondition$celltyp))+
    geom_point(size = 5)+
    xlab( paste0("PC1: ", round(percentVar[1]), "%"))+
    ylab(paste0("PC2: ", round(percentVar[2]), "%"))+
    scale_color_manual(name = "", values = c("Naive" = "grey",
                                             "Tmem" = "lightblue",
                                             "Tex" = "salmon"))+
    theme_classic()
  return(p1)
}