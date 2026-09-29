# Final Monte Carlo — xtdpthresh 0.9.37

Main axis: the Gong–Seo (2026) design and geometry — all available lags
(no `maxlag()`), 46-point p5–p95 quantile grid (`grid(46) gridtype(quantile)
gridsample(observed) trim(.10)`), no refinement, wild grid bootstrap with
B = 500. The FD, balanced cells are the benchmark (they should reproduce
Gong–Seo, Tables 1–2); every other cell changes one factor: FOD, MCAR gaps
(.30) or monotone attrition (.15), N, kappa, the DGP, the instrument set, T,
or the model (kink). 500 replications per cell; FD and FOD are paired.

## Registries

`final_point_cells.csv` (mode POINT, `noboot`; 198 cells, 99,000 fits)

| Block | Design |
|---|---|
| P1 | N 400/800/1600; balanced/gap30/attr15; kappa 0/.1/.2/.5/1 |
| P2 | N 200; balanced/gap30; kappa 0/1 |
| P3 | endogenous q, rho = .9, t(5) shocks; N 400; balanced/gap30; kappa 0/1 |
| P4 | lags 1–3 (`maxlag(1 3)`) and `collapse`; N 400; T 6 |
| P5 | T = 10; N 200/400; all lags, lags 1–3, `collapse` |
| P6 | Windmeijer (jump, kappa 1); kink model, robust and Windmeijer |

`final_inf_cells.csv` (202 cells, 101,000 fits)

| Block | Mode | Design |
|---|---|---|
| I1 | FULL: set, `citest(.25)`, linearity | N 400/800; balanced/gap30/attr15; kappa 0/.1/.2/.5/1 |
| I2 | FULL | N 200; endogenous q (N 400); balanced/gap30; kappa 0/1 |
| I3 | LIN: linearity size | linear DGP; N 400/800; balanced/gap30/attr15 |
| I4 | CI | kink model; N 400/800; balanced/gap30 |
| I5 | CI: `citest(.25 + c)`, c .10/.25/.50 | kappa 0/.5/1; N 400 balanced/gap30/attr15; N 800 balanced |
| I6 | FULL and CI | lags 1–3; N 400; balanced/gap30; kappa 0/1; c 0/.25/.50 |
| X1 | XTH: Seo–Shin via `xthenreg` | FD; balanced; N 400/800; kappa 0/.1/.2/.5/1; same samples as I1; `xthenreg y q, endogenous(q) grid_num(46) trim_rate(.1)` (xthenreg has no predetermined option; its default would instrument q by itself, invalid in this DGP); their asymptotic 95% interval for gamma is stored in `ci_lo`/`ci_hi` |

`make_cells.py` generates both files.

## Outputs (per run folder `final_runs/<RunId>`)

- `final_all.dta/csv`: every replication.
- `final_summary.csv`: per cell — delivery, two-step rate, gamma bias/RMSE/
  median absolute error/P(|error| <= .1), Hansen/AR(1)/AR(2) rejection,
  confidence-set delivery, `citest` rejection (at c = 0: one minus the
  pointwise coverage, `citest_coverage`), hull coverage, mean/median hull
  length, boundary rate, linearity delivery and rejection, MCSEs.
- `final_coefficients.csv`: per cell and coefficient — bias, RMSE, SD, mean
  SE, SE/SD, coverage with the reported and the conditional covariance.
- `final_paired.csv`: FD vs FOD on the same samples — MSE of gamma-hat and of
  the lagged-dependent coefficient (difference and MCSE), McNemar tests for
  the `citest` and linearity rejections.

## Seeds

DGP and missingness seeds follow the Study B scheme (method, kappa, c, model,
covariance, instrument set and block are absent): FD/FOD and all cells of a
DGP, N, T and missingness pattern share the outer sample. Base-DGP cells
therefore reuse the samples of the checking runs. The bootstrap seed has its
own namespace (424242) and excludes method, kappa, c, model and covariance.

## Requirements

The X1 cells call `xthenreg` (Seo, Kim and Kim 2019, SSC) and its dependency
`moremata`: `ssc install xthenreg` and `ssc install moremata` before the run.
`run_tonight.ps1` checks both first (`preflight.do`).

## Run

```powershell
.\run_final.ps1 -Action Fresh -Registry final_point_cells.csv -RunId point_smoke -NShard 4 -RepCap 1 -B 19 -Grid 10 -GridCI 10
.\run_final.ps1 -Action Merge -RunId point_smoke
.\run_final.ps1 -Action Fresh -Registry final_point_cells.csv -RunId point_final -NShard 28
.\run_final.ps1 -Action Merge -RunId point_final
```

and the same with `final_inf_cells.csv`. `..\run_tonight.ps1` runs both registries (smoke, formal, merge) in sequence.
Rough cost on 28 cores: POINT about 1–2 hours, INF about 5 hours.
