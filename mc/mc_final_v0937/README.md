# Final Monte Carlo — xtdpthresh 0.9.37 (frozen)

All folders use the same command files (`stata/` of the repository, 0.9.37).

| Folder | Content | Launcher |
|---|---|---|
| `version_check/` | 0.9.35 / 0.9.36 vs 0.9.37 on common samples | `version_check.do` |
| `study_a/` | Study A core (116 cells) and coefficient supplement (32 cells) | `run_study_a.ps1`, `run_supplement.ps1` |
| `study_b_supp/` | kink supplement (8 cells), power FD/FOD (48 cells), power Gong–Seo geometry (12 cells) | `run_supplement.ps1`, `run_supp2.ps1` |

The Study B core (coverage and linearity) and the new blocks of
`../KE_HOACH_MC_CUOI.md` are added in a later step.

Overnight: `run_tonight.ps1` runs the version check and, for each study, a
smoke followed by the formal run and its verified merge; see its log
`run_tonight_<date>.log`. The harnesses are those of the checking runs
(Study A 0.9.35, Study B supplements 0.9.36) with the version contract moved
to 0.9.37; the version check shows that the command's results are unchanged.
