# Snapshot fields of BM-CF-BACARENA for the microenvironment figures of S12 (make_microenv.py), from the saved state, no rerun.
# usage (in the run directory): Rscript extract_bacarena_fields.R <out prefix>
# Writes
#   <prefix>_media.csv.gz  step_h, x, y, glucose, H2, acetate, methane: concentration at every site (mM) at every hourly step
#   <prefix>_cells.csv.gz  step_h, x, y, species, biomass_pg: every individual at every hourly step
# Site lookup as BacArena's getSublb and extract_bacarena_fig6_qual.R: matrix(diffmat, nrow = m, ncol = n)[y, x], i.e. vector
# index (x - 1) * m + y; mM = amount / conv with conv = 1e12 * 0.01 * scale (as in extract_bacarena_fig6_qual.R).
.libPaths(c("/home/ohayoun/benchmark_tools/R_libs", .libPaths()))
suppressMessages(library(BacArena))
pre <- commandArgs(trailingOnly = TRUE)[1]
s <- readRDS("bacarena_fig6_80x80_closed_nopool.rds")
SP <- c("c_beijerinckii", "m_barkeri")
MET <- c(glucose = "EX_glc_D_e", H2 = "EX_h2_e", acetate = "EX_ac_e", methane = "EX_ch4_e")
conv <- 1e12 * 0.01 * s@scale
grid <- expand.grid(y = 1:s@m, x = 1:s@n)          # column-major: row index (x - 1) * m + y
med <- list(); cel <- list()
for (i in seq_along(s@medlist)) {
  md <- s@medlist[[i]]
  df <- data.frame(step_h = i - 1, x = grid$x, y = grid$y)
  for (k in names(MET)) df[[k]] <- if (is.null(md[[MET[[k]]]])) NA_real_ else as.numeric(md[[MET[[k]]]]) / conv
  med[[i]] <- df
  p <- s@simlist[[i]]
  cel[[i]] <- data.frame(step_h = i - 1, x = p$x, y = p$y, species = SP[p$type], biomass_pg = p$biomass)
}
write.csv(do.call(rbind, med), gzfile(paste0(pre, "_media.csv.gz")), row.names = FALSE)
write.csv(do.call(rbind, cel), gzfile(paste0(pre, "_cells.csv.gz")), row.names = FALSE)
cat("steps", length(s@medlist), "sites", s@n, "x", s@m, "\n")
