# Supplement 2 — power of the threshold test (FD vs FOD)

Separate block; it does not change the core registry or the kink supplement.
Launch with `run_supp2.ps1`. Command: the packaged `xtdpthresh.ado` 0.9.36,
the version of the Study B core run.

## Target

Rejection rate of the grid-bootstrap test of H0: gamma = gamma0 + c at the
5% level (`citest(.25 + c)`), c in {.10, .25, .50}: the power layout of
Gong–Seo (2026, Table 2). Together with c = 0 (the `citest` coverage of the
Study B core) this gives the power curve of the confidence set, which shows
how informative the set is, beyond its coverage.

## Cells (48 cells, 24,000 fits)

| Block | Design | Cells |
|---|---|---:|
| P400 | N=400; balanced / MCAR .30 gaps / monotone attrition .15; kappa 0, 1; c .10, .25, .50; FD, FOD | 36 |
| P800 | N=800; balanced; kappa 0, 1; c .10, .25, .50; FD, FOD | 12 |

Every cell uses 500 replications and the Study B core geometry: jump model,
q predetermined, grid(199) refine(4) trim(.15) maxlag(1 3), uniform grid,
effective support, B=499, `notest gridci(10)` (`citest()` does not depend on
the CI grid).

## Same samples as Study B

The DGP code, the DGP seed and the missingness seed are those of
`study_b_worker.do` (gs26_base_v2). Rep r of a cell is rep r of the Study B
cell in `link_cell` (b1 balanced/gap30, b5 attrition, b2 N=800). Within a
panel type and N, FD/FOD and all kappa/c cells share the outer sample and
the bootstrap draws (own namespace 314159), so FD–FOD and c comparisons are
paired.

## Outputs

- `supp2_summary.csv`: per cell, `citest_reject5` (rejection among evaluable
  tests) with MCSE, the effective rate per requested replication, delivery,
  the two-step rate and the share of mechanical accepts (D = 0).
- `supp2_paired.csv`: FD vs FOD on the same samples and draws: both rates,
  the difference, discordant pairs, exact McNemar p-value.
- `supp2_power_table.csv` (from `supp2_link_check.do`): rejection at
  c = 0, .10, .25, .50 by panel type, kappa and method.

## Run

```powershell
.\run_supp2.ps1 -Action Fresh -RunId supp2_smoke -NShard 4 -RepCap 1 -B 19 -Grid 10 -GridCI 10
.\run_supp2.ps1 -Action Status -RunId supp2_smoke
.\run_supp2.ps1 -Action Merge  -RunId supp2_smoke
# formal:
.\run_supp2.ps1 -Action Fresh -RunId supp2_formal -NShard 28
.\run_supp2.ps1 -Action Status -RunId supp2_formal
.\run_supp2.ps1 -Action Merge -RunId supp2_formal
```

Then, inside `supp2_runs\supp2_formal`, in Stata:

```stata
do supp2_link_check.do "<formal Study B folder>\_merge_stage_<nonce>\study_b_all.dta"
```

It asserts that every replication has the same realized panel and the same
gamma-hat as its Study B counterpart (`SUPP2_LINK_PASS`) and writes the power
table. The smoke run uses a 10-point grid, so its gamma-hat differs from Study
B; run the link check on the formal run only.

Rough cost: 3–5 s per fit at N=400 and 6–8 s at N=800 (from the Study B
`gridci(10)` cells), about 1–1.5 hours on 28 shards.
