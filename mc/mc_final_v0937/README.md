# Final Monte Carlo — xtdpthresh 0.9.37 (frozen)

All folders use the command files of `stata/` (0.9.37).

| Folder | Content |
|---|---|
| `final/` | the final registries (point estimation; inference) on the Gong–Seo design and geometry; see `final/FINAL.md` |

`run_tonight.ps1` checks that xthenreg/moremata are installed, then runs each registry (smoke, formal,
verified merge) in sequence; see its log `run_tonight_<date>.log`.
