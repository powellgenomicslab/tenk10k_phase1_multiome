# Inputs: one locus list per file (chr:pos:ref:alt) in this folder; loci are matched on chr:pos
in_dir  <- "."
out_pre <- "Figure_3D"
read_loci <- function(f) unique(sub("^([^:]+:[^:]+):.*", "\\1", readLines(file.path(in_dir, f))))

gwas <- read_loci("GWAS_clump_loci.txt")
n_loci <- length(gwas)

# Ordered steps; each locus is assigned to the first step it colocalizes in
steps <- list(
  list(lab = "eQTLGen",                  col = "#9A5B8A", file = "eQTL_eQTLGen.txt"),
  list(lab = "Add by TenK eQTL 100 kb",  col = "#8DB87F", file = "eQTL_100kb_single_variants.txt"),
  list(lab = "Add by TenK eQTL 1 Mb",    col = "#5E8FA0", file = "eQTL_1Mb_single_variants.txt"),
  list(lab = "Add by TenK caQTL",        col = "#E5A07F", file = "caQTL_1Mb_single_variants.txt"),
  list(lab = "Add by TenK eQTL mv",      col = "#8CBFB0", file = "eQTL_1Mb_multi_variants.txt"),
  list(lab = "Add by TenK caQTL mv",     col = "#F3D0B5", file = "caQTL_1Mb_multi_variants.txt"),
  list(lab = "Add by GTEx eQTL",         col = "#87B9DC", file = "eQTL_GTEx.txt"),
  list(lab = "Add by IBDverse eQTL",     col = "#CD6B67", file = "eQTL_IBDVerse.txt")
)

assigned <- rep(FALSE, n_loci)
cnt <- integer(0)
for (st in steps) {
  hit <- gwas %in% read_loci(st$file) & !assigned
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
  layout(matrix(1:2, 1), widths = c(3.2, 1.8))
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
  text(0.15, ys, sprintf("%s (n=%d)", rev(labs), rev(cnt)), adj = c(0, 0.5), cex = 0.9)
}

png(paste0(out_pre, ".png"), width = 3200, height = 1000, res = 300); draw(); invisible(dev.off())
