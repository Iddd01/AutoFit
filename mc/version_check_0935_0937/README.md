# Version check: do the published Monte Carlo numbers hold for 0.9.37?

`version_check.do` fits the same generated samples with the archived releases
and with 0.9.37, then compares the stored results.

- **Part A (Study A):** 0.9.35 vs 0.9.37, `noboot`, Study A options; 24
  configurations (FD/FOD × balanced/MCAR .30 × jump κ=1, jump κ=0, imposed
  kink × `vce(robust)`/`vce(windmeijer)`), 10 samples each. Compares `e(b)`,
  `e(V)`, `e(V_cond)`, γ̂, the objective, Hansen/AR, sample and instrument counts.
- **Part B (Study B and supplements):** 0.9.36 vs 0.9.37 with the wild
  bootstrap (B=99, `gridci(20)`, `citest(.25)`); 8 configurations, 3 samples
  each. Compares `e(b)`, `e(ci_segments)`, the hull, the linearity p-value,
  `citest()` and the continuity p-value (0.9.37 with `conttest`, which now
  requests the continuity test).

Run from this folder in Stata 17 (about 15–25 minutes):

```stata
do version_check.do
```

`VERSION_CHECK_PASS` means every compared quantity agrees within a relative
difference of 1e-10, so Study A (0.9.35) and Study B with its supplements
(0.9.36) describe 0.9.37.

Files: `xtdpthresh_v0935.ado` (Study A release), `xtdpthresh_v0936.ado`
(Study B release), `xtdpthresh_v0937.ado` (current `stata/xtdpthresh.ado`).
