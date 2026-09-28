# Study B — xtdpthresh 0.9.36 (28 September 2026)

New immutable-run package. Earlier folders and results are not modified or
reused. Both FD and FOD must be rerun. The package pins the audited 0.9.36
command and records the pointwise threshold test at the true gamma.

## Scope

Core: 52 cells, 26,000 outer fits, retaining the original DGPs,
registry and seed formulas. B1–B5 use B=499; B6 uses B=500.
Supplement: 8 actual restricted-kink CI cells, 4,000 fits at B=499.
Total planned outer fits: 30,000, in separate core and supplement runs.

Every identified-threshold core cell and every supplement cell calls
`citest(.25)`. The merger treats `citest_accept` as the primary pointwise
coverage indicator and checks its status, draw count and independent seed.
The former union-of-`e(ci_segments)` coverage remains in the outputs as a
descriptive set-geometry measure, not as the primary Gong–Seo coverage score.

The core design is documented in [_DESIGN_B.md](_DESIGN_B.md).
The targeted blocks and their launcher are in [SUPPLEMENT.md](SUPPLEMENT.md).
Do not merge supplement rows into the core summaries.

## Core runbook

Use PowerShell 5.1 and Stata 17/MP on a local SSD, not a live cloud-sync folder.
Each worker uses one Stata processor. Choose shard count after the pilot;
28 is intended for the original 32-core VPS, not a 16-core workstation.

```powershell
# If scripts are blocked, run this once in this PowerShell session only:
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
# One row per registered cell: execution/merge smoke only.
.\run_study_b.ps1 -Fresh -RunId smoke_0936 -NShard 6 -RepCap 1 -B 19 -Grid 10 -GridCI 10
.\monitor_study_b.ps1 -RunId smoke_0936
& ".\runs\smoke_0936\verify_and_merge_study_b.ps1" -RunId smoke_0936 -AllowNonFormal

# Production-geometry pilot (same numerical settings as formal run).
.\run_study_b.ps1 -Fresh -RunId pilot_0936 -NShard 6 -RepCap 1

# Formal run: only after smoke and production-geometry pilot.
.\run_study_b.ps1 -Fresh -RunId formal_b0936 -NShard 28
.\monitor_study_b.ps1 -RunId formal_b0936 -Watch 60
& ".\runs\formal_b0936\verify_and_merge_study_b.ps1" -RunId formal_b0936
```

Resume an interrupted run with `-Resume -RunId <same id>` and the original
settings. Never replace the code inside a staged run. The manifest, version,
schema, raw keys, completion markers and hashes are verified before merge.
Estimation failures are outcomes, not reasons to retry different random draws.

## Updated output contract

Raw rows now retain joint_vce, ar_joint, vce_applied, bwscale, gamma_bw,
q_nvals_bw, dropped/near-dependent instrument diagnostics, conditional AR
statistics/p-values, V_cond coefficient SEs and SE delivery.
Point success is separate from SE delivery. Missing SE does not remove a
valid point estimate from bias/RMSE merely because inference was unavailable.

Cell summaries add separate SE panels:
- j: reported joint covariance;
- f: reported conditional fallback;
- c: V_cond diagnostic, on the current fit (not a reproduction of old code).

For each panel and coefficient, n_*, meanse_* (mean standard error),
esd_* (matched-sample empirical SD), se_sd_*, cov_*, cov_mcse_* and
cov_eff_* are reported. The historical wcov_* fields describe the reported
SE mixture; use the new panels when making a claim specifically about joint SEs.
Joint/AR/SE delivery and near-dependent-IV rates are reported separately.
At kappa=0 the core still fits unrestricted jump: do not relabel it as
restricted-kink inference or infer regular identification from joint_vce=1.

CI inversion now follows the shipped fixed-W2 criterion when two-step
estimation succeeds (fixed-W1 only after one-step fallback), common unit wild
weights across candidates, the current finite-B critical-value rule and
all-draw validity gates. The linearity test retains its one-step criterion;
the continuity test uses the fixed W2 of the jump fit after a two-step estimate
(W1 after fallback) and bootstraps kink fitted values plus jump residuals.
ci_criterion_code is 2 or 1 respectively.
The worker and merger check complete-draw delivery and gamma-hat membership
with floating-point tolerance. The primary threshold coverage is pointwise
`citest_accept`; segment-union coverage is retained as set geometry.
A design-dependent failure is retained and reported; a contract violation
stops the harness for inspection. This remains an evaluation of the shipped
approximation, not exact Gong–Seo Algorithm 1 replication.

The study is not already complete. Formal performance results must come from
new verified runs. Keep old results only as separately labelled historical
evidence. See the project audit folder for migration test logs.
