
#Function of  rownames of counts matrix
SetRow <- function(df){
  df <- df |> 
    mutate(ID = paste0(Chr,"-",Start,"-",End)) |> # paste chromosome location and save in ID column
    select(-c(1:6)) |>  # remove columns not required for DGE
    column_to_rownames("ID") # set ID as rownames
  return(df)
}

# Function to make coldata df
SetCondition <- function(df){
  dfCondition <- as.data.frame(colnames(df)) # Get the design from colname names
  colnames(dfCondition) <- "Samples"
  data <- str_split_fixed(dfCondition$Samples, "_", 4) # Split the names to get design names
  dfCondition$Group <- data[,1]
  dfCondition$time <- data[,3]
  dfCondition <- dfCondition |> mutate(celltyp = case_when(Group == "Arm" ~ "Tmem",
                                                           Group == "Cl13" ~ "Tex",
                                                           .default = Group))
  return(dfCondition)
}

#Function to run DESeq
RunDeseq <- function(df,dfcondition){
  dfcondition[["celltyp"]] <- as.factor(dfcondition[["celltyp"]]) # set the T cell type as factor
  dfcondition[["time"]] <- as.factor(dfcondition[["time"]]) # set the time as factor
  dds <- DESeqDataSetFromMatrix(countData = df, colData = dfcondition, ~celltyp+time)
  deRNA <- DESeq(dds)
  return(deRNA)
}