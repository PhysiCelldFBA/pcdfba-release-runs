# ACE-REF in BacArena 1.8.2, v6 (2026-11-19): v5 (medium and Km scaled to the PhysiCelldFBA box, F_SCALE) with the
#   model bounds that BacArena does not convert to a step converted as it converts the exchange bounds
#   (constrain(): exchange lb x dryweight/cellweight_mean x tstep; the LP is solved per step):
#   - ATP maintenance lower bound 8.39 x TSTEP (BacArena leaves it at 8.39 per step; dryweight/cellweight_mean taken as 1);
#   - biomass upper bound 0.86 x TSTEP (growExp grows the biomass by the flux of one step);
#   - Vmax passed as Vmax / TSTEP, because constrain() multiplies the Michaelis-Menten bound by tstep twice.
#   Env TSTEP (default 0.05 h), PULSE_H 3.5. Untouched otherwise.
#   Env GRID_N (default 9): lattice of GRID_N x GRID_N sites over the same 36 um; at 17 the site (2.1 um, 4.5 um2) holds
#   one cell (bacnum 1; consume() removes uptake x bacnum). Env DEATHRATE (default: BacArena's 0.21 1/h) overrides the
#   shrinking of individuals with no growth (diagnostic).
# v5: (2026-11-19): v4 with the medium scaled to the PhysiCelldFBA box (env F_SCALE, default
#   10 cm / 36 um = 2777.8): BacArena's mM refers to a 10 cm column per site, the PhysiCelldFBA box is 36 um deep, so the
#   initial glucose and O2, the O2 pulse and the Km of glucose, O2 and acetate (converted by BacArena with the same
#   column factor) are divided by F_SCALE; Vmax unchanged; outputs reported in PhysiCelldFBA-equivalent mM (x F_SCALE).
#   Env TSTEP (h, default 1) sets BacArena's step and PULSE_H (default 3 at TSTEP 1, else 3.5) the O2 pulse;
#   ATPM_LB overrides the ATP maintenance lower bound (diagnostic only). F_SCALE=1 reproduces v4. The untracked inorganic exchanges stay at 1000 mM (non-limiting).
# v4 (2026-11-19): as v2, plus the SBML exchange bounds restored after Bac() (Bac() sets every
# exchange lower bound to -1000, so in v2 the cells could take up again the formate, ethanol, lactate, succinate, ... they
# secrete; Ec_core lower bound 0, as in PhysiCelldFBA). Glucose, O2 and acetate keep their kinetics.
# v2 notes (2026-11-18), changes from job 667a90ce:
#  - only the inorganic exchanges left open in PhysiCelldFBA (h, h2o, nh4, pi) are added at 1000 mM;
#    adding every other exchange at 1000 mM turned the glucose field into NA after the first step;
#  - concentrations converted with the factor read from t = 0 (amount per mM), not a hard-coded value.
.libPaths(c("/home/ohayoun/benchmark_tools/R_libs", .libPaths()))
suppressMessages({library(BacArena); library(sybil); library(jsonlite)})
F <- as.numeric(Sys.getenv("F_SCALE", as.character(0.1 / 36e-6)))
cat("F_SCALE", F, "\n")
t0 <- Sys.time()
data(Ec_core); model <- Ec_core; model0 <- Ec_core; ids <- react_id(model)
ex_glc <- "EX_glc(e)"; ex_o2 <- "EX_o2(e)"; ex_ac <- "EX_ac(e)"; ex_co2 <- "EX_co2(e)"
model <- changeBounds(model, react=c(ex_glc,ex_o2,ex_ac,ex_co2), lb=rep(-1000,4), ub=rep(1000,4))
bio_id <- ids[which(obj_coef(model) != 0)]
model <- changeBounds(model, react=bio_id, lb=0, ub=0.86)
TSTEP <- as.numeric(Sys.getenv("TSTEP", "0.05"))
model <- changeBounds(model, "ATPM", lb = lowbnd(model0)[ids == "ATPM"] * TSTEP)
model <- changeBounds(model, react=bio_id, lb=0, ub=0.86 * TSTEP)
cat("ATPM lb", lowbnd(model)[ids == "ATPM"], " biomass ub", uppbnd(model)[ids == bio_id], "\n")
bac <- Bac(model)
kin <- c(ex_glc, ex_o2, ex_ac)
exr <- setdiff(intersect(names(bac@lbnd), ids[grepl("^EX_", ids)]), kin)
bac@lbnd[exr] <- lowbnd(model0)[match(exr, ids)]
cat("lbnd EX_for(e):", bac@lbnd["EX_for(e)"], " EX_co2(e):", bac@lbnd["EX_co2(e)"], " EX_glc(e):", bac@lbnd["EX_glc(e)"], "\n")
if (Sys.getenv("DEATHRATE") != "") bac@deathrate <- as.numeric(Sys.getenv("DEATHRATE"))
cat("deathrate", bac@deathrate, "\n")
bac <- setKinetics(bac, exchangeR=ex_glc, Km=0.02/F,  vmax=8.0/TSTEP)
bac <- setKinetics(bac, exchangeR=ex_o2,  Km=0.002/F, vmax=18.0/TSTEP)
bac <- setKinetics(bac, exchangeR=ex_ac,  Km=1.0/F,   vmax=5.0/TSTEP)
PULSE_H <- as.numeric(Sys.getenv("PULSE_H", "3.5"))
n1 <- round(PULSE_H / TSTEP); n2 <- round((10 - PULSE_H) / TSTEP)
GRID_N <- as.integer(Sys.getenv("GRID_N", "9"))
arena <- Arena(n=GRID_N, m=GRID_N, Lx=36e-4, Ly=36e-4, tstep=TSTEP)
cat("TSTEP", TSTEP, "PULSE_H", PULSE_H, "bacnum (cells per site):", round(arena@scale / (bac@cellarea * 10^(-8))), "\n")
set.seed(50); arena <- addOrg(arena, bac, amount=20)
arena <- addSubs(arena, smax=8.0/F, mediac=ex_glc, unit="mM")
arena <- addSubs(arena, smax=5.0/F, mediac=ex_o2,  unit="mM")
arena <- addSubs(arena, smax=0.0, mediac=ex_ac,  unit="mM")
arena <- addSubs(arena, smax=0.0, mediac=ex_co2, unit="mM")
for (r in c("EX_h(e)","EX_h2o(e)","EX_nh4(e)","EX_pi(e)")) arena <- addSubs(arena, smax=1000, mediac=r, unit="mM")
# leg 1: 0-3 h; O2 pulse at 3 h (BacArena steps are whole hours at tstep = 1; PhysiCelldFBA pulses at 210 min)
sim1 <- suppressWarnings(simEnv(arena, time=n1))
arena2 <- getArena(sim1, time=n1)
arena2 <- suppressWarnings(changeSub(arena2, smax=5.0/F, mediac=ex_o2, unit="mM"))
sim2 <- suppressWarnings(simEnv(arena2, time=n2))
t_total <- as.numeric(difftime(Sys.time(), t0, units="secs"))
conv <- mean(sim1@medlist[[1]][[ex_glc]]) / 8.0   # per PhysiCelldFBA-equivalent mM (glucose 8 mM <-> 8/F BacArena-mM)
cat("glucose per site at t0 (fmol):", mean(sim1@medlist[[1]][[ex_glc]]), "\n")
ser <- function(sim, off){ do.call(rbind, lapply(seq_along(sim@medlist), function(i){ md <- sim@medlist[[i]]; p <- sim@simlist[[i]]
  data.frame(time_h=off+(i-1)*TSTEP, n_cells=nrow(p), biomass_pg=sum(p$biomass),
    glucose_mM=mean(md[[ex_glc]])/conv, oxygen_mM=mean(md[[ex_o2]])/conv, acetate_mM=mean(md[[ex_ac]])/conv, co2_mM=mean(md[[ex_co2]])/conv)}))}
full <- rbind(ser(sim1,0), ser(sim2,PULSE_H)[-1,]); full$time_min <- full$time_h*60
write.csv(full, paste0("bacarena_aceref_v6", Sys.getenv("OUT_SFX"), "_timeseries.csv"), row.names=FALSE)
write(toJSON(list(elapsed_s=t_total, threads=1, host=Sys.info()[["nodename"]], bacarena=as.character(packageVersion("BacArena"))), auto_unbox=TRUE), paste0("bacarena_aceref_v6", Sys.getenv("OUT_SFX"), "_timing.json"))
print(full); cat("TIMING", t_total, "\n")
