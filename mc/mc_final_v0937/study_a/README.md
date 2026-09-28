# Study A — xtdpthresh 0.9.35 (27 September 2026)

New immutable-run package. The old 0.9.26 and FOD-only 0.9.27 folders and
results are not modified or reused. Both FD and FOD must be rerun.
The command, predictor and help are the same files as the SSC 0.9.35 release.

## Scope

Core: 116 cells, 48,800 fits, retaining the original DGPs, design cells,
seeds and paired FD/FOD design. The r400 revision increases only the
20 N=1600 cells from R=250 to R=400 before the new formal run.
No bootstrap is run in the core.
Supplement: 32 targeted cells and 12,800 fits for robust/Windmeijer
and actual restricted-kink coefficient inference.
Total planned fits: 61,600, in separate core and supplement runs.

The core design is documented in [_DESIGN.md](_DESIGN.md).
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
.\run_study_a.ps1 -Fresh -RunId smoke_0935_r400 -NShard 6 -RepCap 1
.\monitor_study_a.ps1 -RunId smoke_0935_r400
& ".\runs\smoke_0935_r400\verify_and_merge_study_a.ps1" -RunId smoke_0935_r400

# Production-geometry pilot (same numerical settings as formal run).
.\run_study_a.ps1 -Fresh -RunId pilot_0935_r400 -NShard 6 -RepCap 1

# Formal run: only after smoke and production-geometry pilot.
.\run_study_a.ps1 -Fresh -RunId formal_0935_r400 -NShard 28
.\monitor_study_a.ps1 -RunId formal_0935_r400 -Watch 60
& ".\runs\formal_0935_r400\verify_and_merge_study_a.ps1" -RunId formal_0935_r400
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

The allocation is R=400 in all 92 regular cells, including N=1600, and
R=500 in the 24 persistence/endogeneity/heavy-tail stress cells.
Inference comparisons in the supplement are prespecified, not selected
after looking at formal results.

Use a new RunId for this r400 revision. Its harness identifier is
`xtdpthresh_study_a_v0935_r400`; do not append to or merge with the earlier
R=250 package. Historical archives and staged runs remain unchanged.

The study is not already complete. Formal performance results must come from
new verified runs. Keep old results only as separately labelled historical
evidence. See the project audit folder for migration test logs.
