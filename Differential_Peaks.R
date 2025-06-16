# Function for counts of Differentail peaks
CountDPE <- function(list_data, comparison_name){
  # Extract the qval columns and fold change colums from the dataframe inside a list
  df <- as.data.frame(cbind(list_data[["qvals"]]["TexvsTmem_padj"],list_data[["fc"]]["TexvsTmem_fc"]))
  #Set the rownames to qvals dataframe
  rownames(df) <- rownames(list_data[["qvals"]])
  #Set colnames
  colnames(df) <- c("padj","fc")
  
  # Extract significant picks based upon FDR and absolute Fold change cut according to the paper
  df.sig <- df |> 
    filter(padj < 0.05 & abs(fc) > 1.5) |> 
    # Set the direction positives for left negatives for rightside of the comparison
    mutate(!!comparison_name := if_else(padj*fc < 0 , -1, 1)) |> 
    select(3) # Keep only direction table
  
  # Pivot longer for setting the comparison name
  all_comp_signi <- df.sig |>
    pivot_longer(cols = everything(), names_to = "comparison", values_to = "sig") |>
    group_by(comparison, sig) |>
    summarise(count = n())  |>  # Count the number of peak in each direction
    mutate(com_peaks = sig*count) |> # Setting the direction of the counts
    select(c(1,4)) # Only keep the comparison and peaks counts column with direct
}
