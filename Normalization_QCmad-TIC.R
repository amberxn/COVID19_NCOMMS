library(xlsx)
options(stringsAsFactors = F, warn = -1) 

# Prepare NA-replaced Metabolomics Raw Peak Area Data ---------------------------

# NA-replaced peak area
df <- read.csv("5-worksheet after filling gaps.csv", sep = ",", header = T) # SuppData 2

# how many batch do you have? 
total_batch_num <- 2

# mean level of each metabolite in all QC samples: df.QC$all_mean
if (T) {
  df.QC <- df[, grep(pattern="QC", colnames(df))]
  rownames(df.QC) <- df$Metabolite
  colnames(df.QC) <- gsub("[.]", "_", colnames(df.QC))
  df.QC$all_mean <- apply(df.QC, 1, mean)
}

# Start Normalization ---------------------------

df_final.S <- as.data.frame(df[, 1])
df_final.QC <- as.data.frame(df[, 1])
for (i in 1:total_batch_num) {
  # batch name
  batch_name <- paste0("Batch", sprintf("%02d", i))
  
  # sample sequencing
  sam_seq <- read.xlsx("untargeted.samseq.xlsx", # not provided
                       sheetIndex = i, header = F, startRow = 1, as.data.frame = T)
  # QC index in each batch
  QC_index <- grep("QC_", sam_seq$X1)
  batch_QC_num <- length(QC_index)
  
  for (j in 1:batch_QC_num) {
    # adjacent QC
    s <- ifelse(j!= 1 && QC_index[j] - QC_index[j - 1] > 1, T, F)
    
    if(s){
      # adjacent QC in each group，j-1, j 
      QC1 <- paste0(batch_name, "_QC_", sprintf("%02d", j - 1))
      QC2 <- paste0(batch_name, "_QC_", sprintf("%02d", j))
      group_QC <- df.QC[, which(colnames(df.QC) == QC1 | colnames(df.QC) == QC2)]
      
      # mean level of adjacent QC in each group and factor
      group_QC_mean <- apply(group_QC, 1, mean)
      group_QC_factor <- as.data.frame(df.QC$all_mean / group_QC_mean)
      
      # samples in each group
      group_j <- sam_seq$X1[c((QC_index[j - 1] + 1):(QC_index[j] - 1))]
      group_j <- gsub("Sars_Covid_2020-06-10_", paste(batch_name, "_", sep = ""), group_j)
      group_index <- match(group_j, colnames(df))
      
      if(complete.cases(group_index)){
        # sample QCmad Normalization
        df_S.QCm_normalized <- df[,group_index] * group_QC_factor[, 1]
        df_QC.QCm_normalized <- group_QC*group_QC_factor[,1]
        
        # Sample TIC Normalization
        df_S.TIC_normalized <- apply(df_S.QCm_normalized, 2, function(x) x / sum(x))
        df_QC.TIC_normalized <- apply(df_QC.QCm_normalized, 2, function(x) x / sum(x))
        
        # merge all batches
        df_final.S <- cbind(df_final.S, df_S.TIC_normalized, stringsAsFactors = F)
        df_final.QC <- cbind(df_final.QC, df_QC.TIC_normalized, stringsAsFactors = F)
      }
    }
  }
}
colnames(df_final.S)[1] <- "Metabolite"
colnames(df_final.QC)[1] <- "Metabolite"
df_final <- cbind(df_final.S[, -1], df_final.QC[, -1])

# save
write.csv(df_final, "COVID-19.untargeted.QCmad-TIC.csv", 
          quote = F, row.names = T)
