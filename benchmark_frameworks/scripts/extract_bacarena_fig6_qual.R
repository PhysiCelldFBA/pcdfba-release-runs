# Qualitative CF-REF criteria inputs from the saved BacArena v2 run (bacarena_fig6_80x80_closed_nopool.rds), no rerun.
# usage (in the run directory): Rscript extract_bacarena_cf_qual.R
# Writes
#   bacarena_fig6_closed_nopool_exchange_signs.csv  step_h, species, metabolite, flux: species-summed exchange flux (mfluxlist) at 10/20/30 h
#   bacarena_fig6_closed_nopool_mb_local.csv        step_h, species, biomass_pg, n, H2_mM, acetate_mM: per hourly step, total biomass of
#                                   the species and the biomass-weighted mean concentration at its individuals' sites
# Site lookup as BacArena's getSublb: matrix(diffmat, nrow = m, ncol = n)[y, x], i.e. index (x - 1) * m + y.
# Concentrations converted to mM with the factor of run_bacarena_cfref_48x48_v2.R (conv = 1e12 * 0.01 * scale).
.libPaths(c("/home/ohayoun/benchmark_tools/R_libs", .libPaths()))
suppressMessages(library(BacArena))
s <- readRDS("bacarena_fig6_80x80_closed_nopool.rds")
SP <- c("c_beijerinckii", "m_barkeri")            # orgdat type 1, 2 (order of addOrg)
EX <- c(glucose = "EX_glc_D_e", acetate = "EX_ac_e", H2 = "EX_h2_e", CO2 = "EX_co2_e", methane = "EX_ch4_e")
conv <- 1e12 * 0.01 * s@scale
stopifnot(length(s@mfluxlist[[1]]) == 2)

ex <- list()
for (h in c(10, 20, 30, 40, 50)) {
  f <- s@mfluxlist[[h + 1]]
  for (k in 1:2) {
    v <- f[[k]]
    for (nm in names(EX)) ex[[length(ex) + 1]] <- data.frame(step_h = h, species = SP[k], metabolite = nm,
                                                            flux = if (EX[[nm]] %in% names(v)) unname(v[EX[[nm]]]) else NA)
  }
}
write.csv(do.call(rbind, ex), "bacarena_fig6_closed_nopool_exchange_signs.csv", row.names = FALSE)

loc <- list()
for (i in seq_along(s@simlist)) {
  p <- s@simlist[[i]]; md <- s@medlist[[i]]
  for (k in 1:2) {
    q <- p[p$type == k, ]
    idx <- (q$x - 1) * s@m + q$y
    w <- q$biomass
    loc[[length(loc) + 1]] <- data.frame(step_h = i - 1, species = SP[k], biomass_pg = sum(w), n = nrow(q),
      H2_mM = if (nrow(q)) sum(w * md[["EX_h2_e"]][idx]) / sum(w) / conv else NA,
      acetate_mM = if (nrow(q)) sum(w * md[["EX_ac_e"]][idx]) / sum(w) / conv else NA)
  }
}
write.csv(do.call(rbind, loc), "bacarena_fig6_closed_nopool_mb_local.csv", row.names = FALSE)
cat("done\n")
