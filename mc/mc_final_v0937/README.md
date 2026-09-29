# Final Monte Carlo — xtdpthresh 0.9.37 (frozen)

All folders use the command files of `stata/` (0.9.37).

| Folder | Content |
|---|---|
| `version_check/` | 0.9.35 / 0.9.36 vs 0.9.37 on common samples |
| `final/` | the final registries (point estimation; inference) on the Gong–Seo design and geometry; see `final/FINAL.md` |

`run_tonight.ps1` runs the version check, then each registry (smoke, formal,
verified merge) in sequence; see its log `run_tonight_<date>.log`.
