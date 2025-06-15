
# Loading the packages required for preprocessing
suppressPackageStartupMessages({library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
  library(tidyverse)
  library(DESeq2)
  library(openxlsx)
  library(cowplot)
  library(readr)
  library(edgeR)
  library(rhdf5)
  library(stringr)
  library(sva)
  library(vsn)
  library(genefilter)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(rGREAT)
  library(org.Mm.eg.db)})

path <- "C:/practice_bulk/"
dfCounts.H3K4me <- read_tsv(paste0(path, "GSE285245_RawRC_H3K4me3_cleanID.txt"))
dfCounts.H3K9me3 <- read_tsv(paste0(path, "GSE285245_RawRC_H3K9me3_cleanID.txt"))
dfCounts.H3K27ac <- read_tsv(paste0(path, "GSE285245_RawRC_H3K27ac_cleanID.txt"))
dfCounts.H3K27me3 <- read_tsv(paste0(path, "GSE285245_RawRC_H3K27me3_cleanID.txt"))

# Bedfiles
bedData.H3K4me <- dfCounts.H3K4me[,2:6]
bedData.H3K9me3 <- dfCounts.H3K9me3[,2:6]
bedData.H3K27ac <- dfCounts.H3K27ac[,2:6]
bedData.H3K27me3 <- dfCounts.H3K27me3[,2:6]


#Function of  rownames of counts matrix
SetRow <- function(df){
  df <- df |> 
    mutate(ID = paste0(Chr,"-",Start,"-",End)) |> 
    select(-c(1:6)) |> 
    column_to_rownames("ID")
  return(df)
}

# Function to make condition df

SetCondition <- function(df){
  dfCondition <- as.data.frame(colnames(df))
  colnames(dfCondition) <- "Samples"
  data <- str_split_fixed(dfCondition$Samples, "_", 4)
  dfCondition$Group <- data[,1]
  dfCondition$time <- data[,3]
  dfCondition <- dfCondition |> mutate(celltyp = case_when(Group == "Arm" ~ "Tmem",
                                                                           Group == "Cl13" ~ "Tex",
                                                                           .default = Group))
  return(dfCondition)
}

#make the rownames as chromosome loci
dfCounts.H3K4me <- SetRow(dfCounts.H3K4me)
dfConditions.H3K4me <- SetCondition(dfCounts.H3K4me)

dfCounts.H3K9me3 <- SetRow(dfCounts.H3K9me3)
dfConditions.H3K9me3 <- SetCondition(dfCounts.H3K9me3)

dfCounts.H3K27ac <- SetRow(dfCounts.H3K27ac)
dfConditions.H3K27ac <- SetCondition(dfCounts.H3K27ac)

dfCounts.H3K27me3 <- SetRow(dfCounts.H3K27me3)
dfConditions.H3K27me3 <- SetCondition(dfCounts.H3K27ac)


dfConditions.H3K4me <- as.data.frame(colnames(dfCounts.H3K4me))
colnames(dfConditions.H3K4me) <- "Samples"
data <- str_split_fixed(dfConditions.H3K4me$Samples, "_", 4)
dfConditions.H3K4me$Group <- data[,1]
dfConditions.H3K4me$time <- data[,3]
dfConditions.H3K4me <- dfConditions.H3K4me |> mutate(celltyp = case_when(Group == "Arm" ~ "Tmem",
                                                           Group == "Cl13" ~ "Tex",
                                                           .default = Group))


#Function to run DESeq
RunDeseq <- function(df,dfcondition){
  dfcondition[["celltyp"]] <- as.factor(dfcondition[["celltyp"]])
  dfcondition[["time"]] <- as.factor(dfcondition[["time"]])
  dds <- DESeqDataSetFromMatrix(countData = df, colData = dfcondition, ~celltyp+time)
  deRNA <- DESeq(dds)
  return(deRNA)
}

deRNA.H3K4me <- RunDeseq(dfCounts.H3K4me, dfConditions.H3K4me)
deRNA.H3K9me3 <- RunDeseq(dfCounts.H3K9me3, dfConditions.H3K9me3)
deRNA.H3K27ac <- RunDeseq(dfCounts.H3K27ac,dfConditions.H3K27ac)
deRNA.H3K27me3 <- RunDeseq(dfCounts.H3K27me3 , dfConditions.H3K27me3)

# Normalized data
dfNormalizedCounts <- as.data.frame(counts(deRNA, normalized = T))

# Plotting PCA function
RunPCA <- function(dfnorm,dfcondition){
  dfNormFilter <- varFilter(as.matrix(dfnorm), var.cutoff = 0.2, filterByQuantile = T)
  dfNormPCA <- prcomp(t(dfNormFilter), scale. = T, center = T)
  percentVar <- ((dfNormPCA$sdev)^2 / sum(dfNormPCA$sdev^2)) * 100
  p1 <- ggplot(dfNormPCA$x, aes(x= dfNormPCA$x[,1], y = dfNormPCA$x[, 2], colour = dfcondition$celltyp))+
    geom_point(size = 4)+
    xlab( paste0("PC1: ", round(percentVar[1]), "%"))+
    ylab(paste0("PC2: ", round(percentVar[2]), "%"))+
    scale_color_manual(name = "", values = c("Naive" = "grey",
                                             "Tmem" = "lightblue",
                                             "Tex" = "salmon"))+
    theme_classic()
  return(p1)
}

# Normalized data
dfNormalizedCounts.H3K4me <- as.data.frame(counts(deRNA.H3K4me, normalized = T))
p.H3K4me <- RunPCA(dfNormalizedCounts.H3K4me, dfConditions.H3K4me)+ggtitle("H3K4me")

dfNormalizedCounts.H3K9me3<- as.data.frame(counts(deRNA.H3K9me3, normalized = T))
p.H3K9me3 <- RunPCA(dfNormalizedCounts.H3K9me3, dfConditions.H3K9me3)+ggtitle("H3K9me3")

dfNormalizedCounts.H3K27ac <- as.data.frame(counts(deRNA.H3K27ac, normalized = T))
p.H3K27ac <- RunPCA(dfNormalizedCounts.H3K27ac, dfConditions.H3K27ac)+ggtitle("H3K27ac")

dfNormalizedCounts.H3K27me3 <- as.data.frame(counts(deRNA.H3K27me3, normalized = T))
p.H3K27me3 <- RunPCA(dfNormalizedCounts.H3K27me3, dfConditions.H3K27me3)+ggtitle("H3K27me3")

plot_grid(ncol = 2, p.H3K4me,p.H3K9me3,p.H3K27ac, p.H3K27me3)


# Taking in combinations of two samples to create pairwise samples
dfPairwiseCond <- as.data.frame(combn(levels(as.factor(dfConditions.H3K27ac$celltyp)), 2)) 


#Function for Differential Peak regions
RunDPE <- function(deRNA, dfcondition){
  dfcondition[["celltyp"]] <- as.factor(dfcondition[["celltyp"]])
  dfQval <- data.frame()
  dfFC <- data.frame()
  dfPval <- data.frame()
  for(i in seq_along(dfPairwiseCond)){
    comparision <- print(paste0(dfPairwiseCond[1,i], "vs", dfPairwiseCond[2,i]))
    res <- results(deRNA, contrast = c("celltyp", dfPairwiseCond[1,i],dfPairwiseCond[2,i]), format = "DataFrame")
    # Add results to the respective data frames
    dfQval[rownames(res), paste0(comparision, "_padj")] <- res$padj
    dfFC[rownames(res), paste0(comparision, "_fc")] <- res$log2FoldChange
    dfPval[rownames(res), paste0(comparision, "_pval")] <- res$pvalue
  }
  # Return df
  return(list(qvals = dfQval, fc = dfFC, pvals = dfPval))
}

DPE.H3K4me <- RunDPE(deRNA.H3K4me, dfConditions.H3K4me)
DPE.H3K9me3 <- RunDPE(deRNA.H3K9me3, dfConditions.H3K9me3)
# Setting padj NA values to 1 to make not significant based on FDR correction (padj)
dfQvals.H3K9me3 <- DPE.H3K9me3$qvals
dfQvals.H3K9me3[is.na(dfQvals.H3K9me3)] <- 1

matQValsPerComparison.H3K9me3 <- as.matrix(dfQvals.H3K9me3)
matFCPerComparison.H3K9me3 <- as.matrix(DPE.H3K9me3$fc)
# Converting padj value df to matrix
# matBelowQ <- as.matrix(dfQvals)
matBelowQ.H3K9me3 <- matQValsPerComparison.H3K9me3 < 0.05
matAboveFC.H3K9me3 <- abs(matFCPerComparison.H3K9me3) > 1.5

# Multiply padj matrix to logFC dataframe
matPeakType.H3K9me3 <- matBelowQ.H3K9me3 * matAboveFC.H3K9me3*DPE.H3K9me3$fc

# Renaming conditions and Significance  
groupLabels <- sapply(dfPairwiseCond, function(x) {
  paste0(x[[1]], "vs", x[[2]])})
colnames(matPeakType.H3K9me3) <- paste0(groupLabels, "_Sig")

# Marking 1 for upregulated -1 for downregulated
matPeakType.H3K9me3[matPeakType.H3K9me3 > 0] <- 1
matPeakType.H3K9me3[matPeakType.H3K9me3 < 0] <- (-1)

all_comp_signi.H3K9me3 <- matPeakType.H3K9me3 %>% 
  # select(contains("Sig")) %>%
  pivot_longer(cols = everything(), names_to = "comparison", values_to = "sig")%>% 
  filter(sig != 0)

# Compiling the number of picks
all_comp_signi.H3K9me3 <- all_comp_signi.H3K9me3 |>
  group_by(comparison, sig) |>
  summarise(count = n()) |> mutate(com_peaks = sig*count) |> select(c(1,4)) 

DPE.H3K27ac <- RunDPE(deRNA.H3K27ac, dfConditions.H3K27ac)
DPE.H3K27me3 <- RunDPE(deRNA.H3K27me3, dfConditions.H3K27me3)


dfCounts <- read_tsv(paste0(path, "GSE285248_RawRC_RNA.txt"))
bedData <- dfCounts[,1:6]

dfCounts <- dfCounts[!duplicated(dfCounts$Gene.Name),]
bedData <- dfCounts[,1:6]

dfCounts <- dfCounts |> 
  select(c(5,7:ncol(dfCounts))) |> 
  column_to_rownames("Gene.Name")

dfConditions <- SetCondition(dfCounts)

# Filtering genes greater variance cutoff 20% percentile
dfNormFilter <- varFilter(as.matrix(dfNormalizedCounts), var.cutoff = 0.2, filterByQuantile = T)


dfNormPCA <- prcomp(t(dfNormFilter), scale. = T, center = T)

# Calculate the Percentage of Variance Explained by PCs
percentVar <- ((dfNormPCA$sdev)^2 / sum(dfNormPCA$sdev^2)) * 100


dfOriginal <- read_tsv(paste0(path, "GSE285248_NormRC_filter_wanno.txt"))

dfConditions$celltyp <- as.factor(dfConditions$celltyp)
dfConditions$time <- as.factor(dfConditions$time)

dfCounts <- subset(dfCounts, subset = apply(dfCounts, 1,max) > 5)

dds <- DESeqDataSetFromMatrix(countData = dfCounts,
                              colData = dfConditions,
                              design = ~ celltyp + time)
deRNA <- DESeq(dds)

# Normalized data
dfNormalizedCounts <- as.data.frame(counts(deRNA, normalized = T))

DGE.rna <- RunDPE(deRNA, dfConditions)
dge_sig.rna <- as.data.frame(cbind(DGE.rna$qvals$TexvsTmem_padj,DGE.rna$fc$TexvsTmem_fc))
colnames(dge_sig.rna) <- c("RNA_padj", "RNA_FC")
rownames(dge_sig.rna) <- rownames(DGE.rna$qvals)
dge_sig.rna <- subset(dge_sig.rna, RNA_padj < 0.05 & abs(RNA_FC) > 1.5)
dge_sig.rna <- dge_sig.rna |> rownames_to_column("gene")


dge_sig.H3K27ac <- as.data.frame(cbind(DPE.H3K27ac$qvals$TexvsTmem_padj,DPE.H3K27ac$fc$TexvsTmem_fc))
colnames(dge_sig.H3K27ac) <- c("H3K27ac_padj", "H3K27ac_FC")
rownames(dge_sig.H3K27ac) <- rownames(DPE.H3K27ac$qvals)
dge_sig.H3K27ac <- subset(dge_sig.H3K27ac, H3K27ac_padj < 0.05 & abs(H3K27ac_FC) > 1.5)
dge_sig.H3K27ac <- dge_sig.H3K27ac |> rownames_to_column("ID")
data <- str_split_fixed(dge_sig.H3K27ac$ID, "-", 3)
dge_sig.H3K27ac$chr <- data[,1]
dge_sig.H3K27ac$start <- data[,2]
dge_sig.H3K27ac$end <- data[,3]
dge_sig.H3K27ac <- cbind(dge_sig.H3K27ac[,4:6], dge_sig.H3K27ac[,2:3])
dge_sig.H3K27ac.bed <- dge_sig.H3K27ac[,1:3]

library(ChIPpeakAnno)

#Making GRange object with Gene Ontology ID and  
res = great(toGRanges(dge_sig.H3K27ac.bed), "GO:BP", "txdb:mm10")# for mouse

# Getting gene names associated with chromosome location id
dfRegionsGenes <- as.data.frame(getRegionGeneAssociations(res))

associated_genes <- dfRegionsGenes |> select(1:3, 6) |> 
  mutate(ID = paste0(seqnames,"-",start,"-",end)) |> 
  select(-c(1:3)) 
rownames(associated_genes) <- NULL
associated_genes <- associated_genes |> 
  column_to_rownames("ID")
dge_sig.H3K27ac <- cbind(dfRegionsGenes$annotated_genes, dge_sig.H3K27ac)

joined_df <- inner_join(dge_sig.H3K27ac, associated_genes, by = "ID")
joined_df <- unnest(joined_df, annotated_genes)
colnames(joined_df)[4] <- "gene"

joined_df.rna.H3K27ac <- inner_join(joined_df, dge_sig.rna, by = "gene")
joined_df.rna.H3K27ac$RNA_FC[is.na(joined_df.rna.H3K27ac$RNA_FC)] <- 0
# Color based on RNA direction (optional, similar to red/blue in your plot)
joined_df.rna.H3K27ac$color <- ifelse(df$RNA_log2FC > 0, "red", "blue")

# Color based on RNA direction (optional, similar to red/blue in your plot)
joined_df.rna.H3K27ac$color <- case_when(joined_df.rna.H3K27ac$RNA_FC > 0~ "salmon", 
                                         joined_df.rna.H3K27ac$RNA_FC < 0~"lightblue",
                                         joined_df.rna.H3K27ac$RNA_FC == 0~ "white")
p1 <- ggplot(joined_df.rna.H3K27ac, aes(x = H3K27ac_FC, y = RNA_FC)) +
  geom_point(aes(color = color), alpha = 0.6, size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    x = "H3K27ac log2FC",
    y = "RNA log2FC",
    title = expression(T[EX]~"versus"~T[MEM])
  ) +
  theme_minimal() +
  scale_color_manual(values = c("salmon", "lightblue", "white")) +
  theme(legend.position = "none")
p1


ggplot(joined_df.rna.H3K27ac, aes(x = H3K27ac_FC, y = RNA_FC, color = color)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_manual(values = c("salmon", "lightblue", "white")) +
  labs(
    x = "H3K27ac log2FC",
    y = "RNA log2FC",
    title = expression(T[EX]~"versus"~T[MEM])
  ) +
  annotate("text", x = min(joined_df.rna.H3K27ac$H3K27ac_FC), y = max(joined_df.rna.H3K27ac$RNA_FC), 
           label = paste("R =", cor(joined_df.rna.H3K27ac$H3K27ac_FC, joined_df.rna.H3K27ac$RNA_FC), ", p < 0.001")