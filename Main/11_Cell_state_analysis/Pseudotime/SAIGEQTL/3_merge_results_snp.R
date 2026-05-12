library(tidyverse)
library(data.table)
library(glue)
library(qvalue)

chr_num <- ""

folder <- glue("/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Results_merged/chr{chr_num}")
output <- glue("/g/data/ei56/jf1058/TenK10K/Multiome/SAIGEQTL/Dynamic/Results_merged/AllResults/chr{chr_num}.csv")

print(glue("Processing chromosome {chr_num}"))
file_list <- list.files(folder, pattern = "^chr[^/]+_\\d+_\\d+$", full.names = TRUE)
file_list <- file_list[!grepl("_ACAT$", file_list)]

df <- rbindlist(lapply(file_list, function(f) {
  if (dir.exists(f)) return(NULL)
  peak <- basename(f)
  d <- fread(f)
  d[, Peak := peak]
  
  # Calculate q-value per peak
  tryCatch({
    qobj <- qvalue(p = d$p.value)
    d[, qvalue := qobj$qvalues]
  }, error = function(e) {
    # qvalue can fail if p-values are poorly behaved (e.g., all ~1)
    warning(paste("qvalue failed for peak:", peak, "-", e$message))
    d[, qvalue := NA_real_]
  })
  
  d
}), use.names = TRUE, fill = TRUE)

setcolorder(df, c("Peak", setdiff(names(df), "Peak")))

fwrite(df, output)

