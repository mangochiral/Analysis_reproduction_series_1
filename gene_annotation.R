GeneAnnot <- function(dfDge){
  bed <- dfDge|> select("ID")
  data <- str_split_fixed(bed[["ID"]], "-", 3)
  bed$chr <- data[,1]
  bed$start <- data[,2]
  bed$end <- data[,3]
  bed <- bed |> select(-c(1))
  #Making GRange object with Gene Ontology ID and  
  res = great(toGRanges(bed), "GO:BP", "txdb:mm10")# for mouse
  
  # Getting gene names associated with chromosome location id
  dfRegionsGenes <- as.data.frame(getRegionGeneAssociations(res))
  
  # Pasting the chromosome start end and id to join with Differential peaks dataframe
  bed.associated_genes <- dfRegionsGenes |> select(1:3, 6) |> 
    mutate(ID = paste0(seqnames,"-",start,"-",end)) |> 
    select(-c(1:3)) 
  rownames(bed.associated_genes) <- NULL
  return(bed.associated_genes)
}

# Color Scheme

ColorScheme <- function(joined_df){
  joined_df[["color"]] <- case_when(joined_df[["RNA_FC"]] > 0~ "salmon", 
                                    joined_df[["RNA_FC"]] < 0~"lightblue",
                                    joined_df[["RNA_FC"]] == 0~ "white")
  return(joined_df)
}

# Scatter Plot
plotCorrelation <- function(joined_df, col1, col2, color){
  p2 <- ggplot(joined_df, aes(x = .data[[col1]], y = .data[[col2]])) +
    geom_point(aes(color = .data[[color]]), alpha = 0.6, size = 1) +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    labs(
      x = as.character(col1),
      y = as.character(col2),
      title = expression(T[EX]~"versus"~T[MEM])
    ) +
    theme_minimal() +
    scale_color_manual(values = c("salmon", "lightblue", "white")) +
    theme(legend.position = "none")
  return(p2)
} 
