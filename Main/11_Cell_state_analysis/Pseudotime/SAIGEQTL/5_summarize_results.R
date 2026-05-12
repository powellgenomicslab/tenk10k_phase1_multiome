
library(data.table)
library(qvalue)
library(glue)

# for (chr in 1:22) {
#   print(glue("Processing chromosome {chr}"))
#   df <- fread(glue("/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Results_merged/AllResults/chr{chr}.csv"))
#   df$qval <- qvalue(df$pval_ge)$qvalues
#   df$lfdr <- qvalue(df$pval_ge)$lfdr
#   df <- df[lfdr < 0.05]
#   if (chr == 1) {
#     all_df <- df
#   } else {
#     all_df <- rbind(all_df, df)
#   }
# }
# fwrite(all_df, "/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Results_merged/AllResults/All_chr_snp_sig.csv")


for (chr in 1:22) {
  print(glue("Processing chromosome {chr}"))
  df <- fread(glue("/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Results_merged/AllResults/chr{chr}.csv"))
  df$qval <- qvalue(df$pval_ge)$qvalues
  df$lfdr <- qvalue(df$pval_ge)$lfdr
  # df <- df[lfdr < 0.05]
  if (chr == 1) {
    all_df <- df
  } else {
    all_df <- rbind(all_df, df)
  }
}
fwrite(all_df, "/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Results_merged/AllResults/All_chr_snp.csv")