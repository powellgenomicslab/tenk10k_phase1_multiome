library(tidyverse)
library(data.table)
library(glue)
library(qvalue)

chr_num <- ""

folder <- glue("/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Results_merged/chr{chr_num}")
output <- glue("/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Results_merged/AllResults/chr{chr_num}_ACAT.csv")

df_list <- list()

print(glue("Processing chromosome {chr_num}"))
file_list <- list.files(glue(folder), pattern = "_ACAT$", full.names = TRUE)

pb <- txtProgressBar(0, length(file_list), style = 3)
for (i in seq_along(file_list)) {
  df_list[[basename(file_list[i])]] <- fread(file_list[i])
  setTxtProgressBar(pb, i)
}
close(pb)

df <- rbindlist(df_list)
df$qvalue <- qvalue(df$ACAT_p)$qvalues

fwrite(df, output)
