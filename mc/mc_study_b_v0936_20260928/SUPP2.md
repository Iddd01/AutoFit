# Supplement 2 — continuity-test comparison and threshold-test power

Separate block; it does not change the core registry or the kink supplement,
and uses its own bootstrap-seed namespace. Launch with `run_supp2.ps1`.

## Block C — continuity test, 0.9.35 vs 0.9.36 (16 cells, 8,000 fits)

Gong–Seo benchmark DGP (predetermined q), N=400, T=6, balanced, jump model,
kappa in {0, 1, 2, 3}, FD and FOD. Core B1 geometry: grid(199) refine(4)
trim(.15) maxlag(1 3), uniform grid, effective support, B=499; gridci(10)
because the CI is not the target here (linearity and continuity run).

The same outer samples and bootstrap seeds are used for both versions, both
methods and all kappa, so the version comparison is paired:
- 0.9.35: continuity under the first-step weight with kink residuals;
- 0.9.36: Gong–Seo statistic (efficient W2) with jump (unrestricted) residuals.
The worker loads `xtdpthresh_v0935.ado` (the exact 0.9.35 release, SHA-256
7B97571397D0FAAB...) or `xtdpthresh.ado` (0.9.36) per cell; all 0.9.35 cells
come first in the registry, so each shard switches version at most once.
0.9.36 cells also request `citest(.25)`.

Outputs: `supp2_summary.csv` (size/power, MCSE) and `supp2_paired.csv`
(per method and kappa: both rejection rates, their difference, discordant
pairs and an exact McNemar p-value).

## Block P — power of the threshold test (12 cells, 6,000 fits)

B6 geometry (the published Gong–Seo design): grid(46) refine(0) trim(.10)
maxlag(1 5), quantile grid, observed support, B=500, FD, N in {400, 800},
kappa in {0, 1}; `citest(.25 + c)` with c in {.10, .25, .50}, as in Gong–Seo
Table 2. `notest gridci(10)`: `citest()` does not depend on the CI grid.
`citest_reject5` is the rejection rate at the false threshold .25 + c.

## Run

```powershell
.\run_supp2.ps1 -Action Fresh -RunId supp2_smoke -NShard 4 -RepCap 1 -B 19 -Grid 10 -GridCI 10
.\run_supp2.ps1 -Action Status -RunId supp2_smoke
.\run_supp2.ps1 -Action Merge  -RunId supp2_smoke
# formal, after the smoke passes and the kink supplement has finished:
.\run_supp2.ps1 -Action Fresh -RunId supp2_formal -NShard 28
.\run_supp2.ps1 -Action Merge -RunId supp2_formal
```

With `-Grid 10` the smoke uses a 10-point grid in both blocks; it checks
execution and contracts only. Formal results require no overrides.
