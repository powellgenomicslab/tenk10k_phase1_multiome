# Load libraries
library(tidyverse)
library(data.table)
library(glue)

d <- fread(PATH)
d <- d[!is.na(d$GWAS_loci) & d$GWAS_loci != "", ]
n_loci <- nrow(d)

steps <- list(
  list(lab = "eQTLGen",                  col = "#9A5B8A", column = "Colocalized_with_eQTLGen_eQTL"),
  list(lab = "Add by TenK eQTL 100 kb",  col = "#8DB87F", column = "Colocalized_with_100kb_eQTL"),
  list(lab = "Add by TenK eQTL 1 Mb",    col = "#5E8FA0", column = "Colocalized_with_eQTL"),
  list(lab = "Add by TenK caQTL",        col = "#E5A07F", column = "Colocalized_with_caQTL"),
  list(lab = "Add by TenK eQTL mv",      col = "#8CBFB0", column = "Colocalized_with_multi_variant_coloc_eQTL"),
  list(lab = "Add by TenK caQTL mv",     col = "#F3D0B5", column = "Colocalized_with_multi_variant_coloc_caQTL"),
  list(lab = "Add by GTEx eQTL",         col = "#87B9DC", column = "Colocalized_with_GTEx_eQTL"),
  list(lab = "Add by IBDverse eQTL",     col = "#CD6B67", column = "Colocalized_with_IBDverse_eQTL")
)
missing <- setdiff(sapply(steps, `[[`, "column"), names(d))
if (length(missing)) stop("Missing columns in input: ", paste(missing, collapse = ", "))

assigned <- rep(FALSE, n_loci)
cnt <- integer(0)
for (st in steps) {
  hit <- as_lgl(d[[st$column]]) & !assigned
  cnt <- c(cnt, sum(hit))
  assigned <- assigned | hit
}
labs <- c(sapply(steps, `[[`, "lab"), "None")
cols <- c(sapply(steps, `[[`, "col"), "#C9CDCE")
cnt  <- c(cnt, sum(!assigned))
frac <- cnt / n_loci
stopifnot(sum(cnt) == n_loci)

summary_df <- data.frame(category = labs, n_loci = cnt, fraction = round(frac, 4))
print(summary_df, row.names = FALSE)
cat(sprintf("Total loci: %d; colocalizing with any molQTL: %.1f%%\n",
            n_loci, 100 * (1 - tail(frac, 1))))
write.csv(summary_df, paste0(out_pre, "_summary.csv"), row.names = FALSE)

# ---- draw ------------------------------------------------------------------
draw <- function() {
  layout(matrix(1:2, 1), widths = c(3.2, 1.2))
  # panel 1: bar
  par(mar = c(1, 1, 3, 0.5), family = "sans", xpd = NA)
  plot(NA, xlim = c(-0.04, 1.04), ylim = c(-0.6, 0.4), axes = FALSE, xlab = "", ylab = "",
       xaxs = "i", yaxs = "i")
  edges <- c(0, cumsum(frac)); keep <- frac > 0
  rect(edges[-length(edges)][keep], -0.35, edges[-1][keep], -0.05, col = cols[keep], border = NA)
  ay <- 0.05
  segments(0, ay, 1, ay, lwd = 1.5)
  segments(seq(0, 1, 0.25), ay, seq(0, 1, 0.25), ay + 0.07, lwd = 1.5)
  segments(seq(0.125, 0.875, 0.25), ay, seq(0.125, 0.875, 0.25), ay + 0.035, lwd = 1)
  text(seq(0, 1, 0.25), ay + 0.17, c("0", "0.25", "0.50", "0.75", "1.00"), cex = 1.4)
  mtext("IBD GWAS clumped loci colocalizing with molQTLs (%)", side = 3, line = 0.3, cex = 1.2, adj = 0.5)
  # panel 2: legend (top to bottom = reverse of bar order, as in the published key)
  par(mar = c(1, 0.5, 3, 1))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  text(0, 1.04, "GWAS loci type", adj = c(0, 0), cex = 1.1)
  n <- length(labs); ys <- seq(0.88, 0.06, length.out = n)
  rect(0, ys - 0.04, 0.1, ys + 0.04, col = rev(cols), border = NA)
  text(0.15, ys, sprintf("%s (n=%d)", rev(labs), rev(cnt)), adj = c(0, 0.5), cex = 1.0)
}

png(paste0(out_pre, ".png"), width = 3200, height = 1000, res = 300); draw(); invisible(dev.off())
