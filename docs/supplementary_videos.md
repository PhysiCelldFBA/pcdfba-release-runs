# Supplementary Videos S1–S4

The four supplementary videos of the revised manuscript are drawn from the same runs as the main figures. The runs are in the outputs record doi:10.5281/zenodo.24093301 (one archive per use case; each run folder has every full-data snapshot).

| Video | Content | Run | Archive | Command (from the use-case folder) |
|---|---|---|---|---|
| S1 | Acetate exchange of each cell, one quadrant of the colony, 0–10 h, every 20 min, 4 fps | `COL-REF_s50` (Fig. 4) | `colony_runs.tar.zst` | `python scripts/make_movie.py runs/COL-REF_s50 Video_S1.mp4 4` |
| S2 | Growth rate of each cell (Fig. 5A), 0–72 h, every hour, 8 fps | `TUM5-PFBA_s50` (Fig. 5) | `tumour_runs.tar.zst` | `python scripts/make_movies.py runs/TUM5-PFBA_s50 OUTDIR 8` |
| S3 | Uptake of O2, glucose, glutamine, glycine and lactate secretion (Fig. S9), 0–72 h, every hour, 8 fps | `TUM5-PFBA_s50` | `tumour_runs.tar.zst` | same command as S2 |
| S4 | Consortium, growth rate by species and counters, 0–45 h, every 30 min, 8 fps | `CF3-NC50_s50` (Fig. 6) | `crossfeeding_runs.tar.zst` | `python scripts/make_movie.py runs/CF3-NC50_s50 Video_S4.mp4 2700 8` |

Extract only the run you need, without storing the archive:

    curl -L https://zenodo.org/api/records/24093301/files/<archive>/content | zstd -dc | tar -x usecases/<case>/runs/<run>

The colour scales are fixed over each video. The dashed zone lines of S2 and S3 are recomputed at every snapshot. All four need ffmpeg (H.264, yuv420p).
