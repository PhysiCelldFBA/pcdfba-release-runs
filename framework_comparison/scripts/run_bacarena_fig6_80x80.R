# v2 (2026-11-19, env BOUNDARY=const, default): the boundary of the PhysiCelldFBA run with BacArena's own
#   constant-boundary diffusion model (pde "ConstBoundDiff2d", edge value boundS, in the same per-site amount as the
#   medium): glucose, CO2 1 mM; NH4, Pi, SO4 5 mM; H2, acetate, methane 0 mM. Diffusion coefficients of the PhysiCelldFBA
#   configuration (um2/min x 6e-7 = cm2/h) instead of BacArena's default 0.02412 cm2/h for every substance.
#   BOUNDARY=closed: closed arena (BacArena's batch mode), same diffusion coefficients. INIT_POOL: initial H2/acetate (mM).
#   OUT_SFX names the outputs.
# v2 also restores the SBML exchange bounds (env SBML_BOUNDS=1, default): Bac() sets the lower bound of every exchange
#   to -1000, so in v1 the untracked nutrients supplied at 1000 mM had no uptake limit (M. barkeri used H2S, SBML
#   bound -10 mmol/gDW/h, as electron donor at ~80 mmol/gDW/h, releasing SO3) and C. beijerinckii could take up H2,
#   acetate and CO2 (Vmax 0 in PhysiCelldFBA, SBML lower bound 0). Now: substrates with kinetics open (limited by the
#   Michaelis-Menten rate, as in PhysiCelldFBA); every other exchange at its SBML lower bound.
# Consortium of manuscript Fig. 6 (PhysiCelldFBA run CF3-NC50_s50) in BacArena 1.8.2, 2026-11-19.
# From run_bacarena_fig6_80x80.R; changes, all from the CF3-NC50_s50 configuration:
#   box 200 um -> 80 x 80 sites of 2.5 um; 50 h; Km(H2) of M. barkeri 0.01 mM; H2 and acetate start at 1 mM;
#   models without the FVA caps: the bounds of *_trimmed_nocaps.xml (same reactions and stoichiometry as the trimmed
#   models, checked with cobra) applied to the trimmed modelorgs from nocaps_bounds_{cb,mb}.csv, before the
#   exchange and biomass bounds below. Unchanged: closed arena (BacArena has no fixed-concentration boundary, so the
#   [v1 only] 0 mM sink of H2 and acetate at the edge not represented), random mixed seeding block at the centre (seed 50),
#   unit mM in BacArena's convention, untracked open exchanges at 1000 mM.

Sys.setenv(R_LIBS_USER="/home/ohayoun/benchmark_tools/R_libs")
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(library(BacArena))
suppressMessages(library(sybil))

t0 <- Sys.time()

cb_model <- readRDS("cb_modelorg.rds")
mb_model <- readRDS("mb_modelorg.rds")
apply_nocaps <- function(m, f) {
  b <- read.csv(f, stringsAsFactors = FALSE)
  stopifnot(all(b$react %in% react_id(m)))
  changeBounds(m, react = b$react, lb = b$lb, ub = b$ub)
}
cb_model <- apply_nocaps(cb_model, "nocaps_bounds_cb.csv")
mb_model <- apply_nocaps(mb_model, "nocaps_bounds_mb.csv")

ex_glc <- "EX_glc_D_e"; ex_h2 <- "EX_h2_e"; ex_nh4 <- "EX_nh4_e"; ex_pi <- "EX_pi_e"
ex_so4 <- "EX_so4_e"; ex_ac <- "EX_ac_e"; ex_co2 <- "EX_co2_e"; ex_ch4 <- "EX_ch4_e"

SBML_BOUNDS <- Sys.getenv("SBML_BOUNDS", "1") == "1"
cb_kin <- c(ex_glc, ex_nh4, ex_pi, ex_so4); mb_kin <- c(ex_co2, ex_h2, ex_nh4, ex_pi, ex_ac)
cb_open <- if (SBML_BOUNDS) cb_kin else c(ex_glc, ex_h2, ex_nh4, ex_pi, ex_so4, ex_ac, ex_co2)
cb_model <- changeBounds(cb_model, react = cb_open, lb = rep(-1000, length(cb_open)), ub = rep(1000, length(cb_open)))
cb_model <- changeBounds(cb_model, react = "biomass", lb = 0, ub = 0.4)

mb_open <- if (SBML_BOUNDS) mb_kin else c(ex_co2, ex_h2, ex_ch4, ex_nh4, ex_pi, ex_ac)
mb_model <- changeBounds(mb_model, react = mb_open, lb = rep(-1000, length(mb_open)), ub = rep(1000, length(mb_open)))
mb_model <- changeBounds(mb_model, react = "Mb_biomass_65", lb = 0, ub = 0.12)

cb_bac <- Bac(cb_model, deathrate = 0)
mb_bac <- Bac(mb_model, deathrate = 0)
restore_lb <- function(bac, model) {
  ids <- react_id(model); ex <- intersect(names(bac@lbnd), ids[grepl("^EX_", ids)])
  bac@lbnd[ex] <- lowbnd(model)[match(ex, ids)]; bac
}
if (SBML_BOUNDS) { cb_bac <- restore_lb(cb_bac, cb_model); mb_bac <- restore_lb(mb_bac, mb_model) }
cat("lbnd EX_h2s_e (Mb):", mb_bac@lbnd["EX_h2s_e"], " EX_h2_e (Cb):", cb_bac@lbnd["EX_h2_e"], " EX_ch4_e (Mb):", mb_bac@lbnd["EX_ch4_e"],
    " EX_h2_e (Mb):", mb_bac@lbnd["EX_h2_e"], "\n")

cb_bac <- setKinetics(cb_bac, exchangeR = ex_glc, Km = 1, vmax = 7.5)
cb_bac <- setKinetics(cb_bac, exchangeR = ex_nh4, Km = 0.001, vmax = 15.0)
cb_bac <- setKinetics(cb_bac, exchangeR = ex_pi,  Km = 0.001, vmax = 15.0)
cb_bac <- setKinetics(cb_bac, exchangeR = ex_so4, Km = 0.001, vmax = 15.0)

mb_bac <- setKinetics(mb_bac, exchangeR = ex_co2, Km = 0.0001, vmax = 100.0)
mb_bac <- setKinetics(mb_bac, exchangeR = ex_h2,  Km = 0.01,   vmax = 20.0)
mb_bac <- setKinetics(mb_bac, exchangeR = ex_nh4, Km = 0.0001, vmax = 10.0)
mb_bac <- setKinetics(mb_bac, exchangeR = ex_pi,  Km = 0.0001, vmax = 10.0)
mb_bac <- setKinetics(mb_bac, exchangeR = ex_ac,  Km = 0.0001, vmax = 3.0)

# PROPERLY-SIZED grid: 48x48 sites at 2.5um each = 120x120um, matching CF-REF's own
# box exactly at BacArena's minimum feasible site size (sqrt(cellarea)=2.1um). Our
# earlier 20x20 grid covered only 50x50um of physical space at this site size --
# a self-imposed undersizing, not something CF-REF or BacArena's design required.
N_BOXES <- 80
site_width_cm <- 2.5e-4
arena <- Arena(n = N_BOXES, m = N_BOXES, Lx = N_BOXES * site_width_cm, Ly = N_BOXES * site_width_cm)

set.seed(50)
n_cb <- 62; n_mb <- 49
block_r <- ceiling(sqrt(n_cb + n_mb))
center <- N_BOXES %/% 2
coords <- expand.grid(x = (center - block_r %/% 2):(center + block_r %/% 2),
                       y = (center - block_r %/% 2):(center + block_r %/% 2))
coords <- coords[sample(nrow(coords)), ]
cb_coords <- coords[1:n_cb, ]
mb_coords <- coords[(n_cb + 1):(n_cb + n_mb), ]
arena <- addOrg(arena, cb_bac, amount = n_cb, x = cb_coords$x, y = cb_coords$y)
arena <- addOrg(arena, mb_bac, amount = n_mb, x = mb_coords$x, y = mb_coords$y)

cat("n organisms placed:", nrow(arena@orgdat), " total sites:", N_BOXES*N_BOXES, "\n")

BOUNDARY <- Sys.getenv("BOUNDARY", "const")
# initial and boundary concentrations (mM) and diffusion coefficients (um2/min) of CF3-NC50_s50
POOL <- as.numeric(Sys.getenv("INIT_POOL", "1"))   # initial H2 and acetate (mM); 0 = no initial pool
init  <- c(EX_glc_D_e = 1, EX_co2_e = 1, EX_nh4_e = 5, EX_pi_e = 5, EX_so4_e = 5, EX_h2_e = POOL, EX_ac_e = POOL, EX_ch4_e = 0)
bound <- c(EX_glc_D_e = 1, EX_co2_e = 1, EX_nh4_e = 5, EX_pi_e = 5, EX_so4_e = 5, EX_h2_e = 0, EX_ac_e = 0, EX_ch4_e = 0)
D_um2min <- c(EX_glc_D_e = 40000, EX_co2_e = 115000, EX_nh4_e = 118000, EX_pi_e = 53000, EX_so4_e = 64000,
              EX_h2_e = 270000, EX_ac_e = 72000, EX_ch4_e = 90000)
conv0 <- 1e12 * 0.01 * arena@scale   # per-site amount per BacArena-mM (as in the extraction)
for (r in names(init)) {
  arena <- addSubs(arena, smax = init[[r]], mediac = r, unit = "mM", difspeed = D_um2min[[r]] * 6e-7,
                   pde = if (BOUNDARY == "const") "ConstBoundDiff2d" else "Diff2d")
  if (BOUNDARY == "const") arena@media[[r]]@boundS <- bound[[r]] * conv0
}
cat("boundary:", BOUNDARY, " H2 difspeed:", arena@media[["EX_h2_e"]]@difspeed, " glc boundS:", arena@media[["EX_glc_D_e"]]@boundS,
    " glc smax-site:", max(arena@media[["EX_glc_D_e"]]@diffmat), "\n")

tracked <- c(ex_glc, ex_co2, ex_nh4, ex_pi, ex_so4, ex_h2, ex_ac, ex_ch4)
open_ex <- function(m) { ids <- react_id(m); ex <- ids[grepl("^EX_", ids)]; ex[lowbnd(m)[match(ex, ids)] < 0] }
extra <- setdiff(union(open_ex(cb_model), open_ex(mb_model)), tracked)
cat("extra medium at 1000 mM:", extra, "\n")
for (r in extra) arena <- addSubs(arena, smax = 1000.0, mediac = r, unit = "mM")

MAXCYCLES <- as.integer(Sys.getenv("CF_MAXCYCLES", "50"))
sim1 <- simEnv(arena, time = MAXCYCLES)
saveRDS(sim1, paste0("bacarena_fig6_80x80", Sys.getenv("OUT_SFX", "_cb"), ".rds"))
t_run <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
conv <- 1e12 * 0.01 * sim1@scale   # amount per site per BacArena-mM
rows <- list()
for (i in seq_along(sim1@simlist)) { p <- sim1@simlist[[i]]; md <- sim1@medlist[[i]]
  for (ty in 1:2) { q <- p[p$type == ty, ]
    rows[[length(rows) + 1]] <- data.frame(step_h = i - 1, species = c("c_beijerinckii", "m_barkeri")[ty],
      n_cells = nrow(q), biomass_pg = sum(q$biomass),
      glucose_mM = mean(md[[ex_glc]]) / conv, H2_mM = mean(md[[ex_h2]]) / conv,
      acetate_mM = mean(md[[ex_ac]]) / conv, methane_mM = mean(md[[ex_ch4]]) / conv) } }
ts <- do.call(rbind, rows)
write.csv(ts, paste0("bacarena_fig6_80x80", Sys.getenv("OUT_SFX", "_cb"), "_timeseries.csv"), row.names = FALSE)
writeLines(sprintf('{"elapsed_s": %.3f, "threads": 1, "host": "%s"}', t_run, Sys.info()[["nodename"]]), paste0("bacarena_fig6_80x80", Sys.getenv("OUT_SFX", "_cb"), "_timing.json"))
print(ts[ts$step_h %in% c(0, 10, 20, 30, 40, 50), ])
cat("TIMING", t_run, "\n")
