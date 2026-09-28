# Validation record — Study A 0.9.35

Date: 2026-09-27. Platform: Windows, Stata 17/MP.
This validates the software harness and summary calculations.
It is NOT a completed formal Monte Carlo study or coverage certificate.

## R400 revision checks

The r400 revision changes only the replication allocation for the 20 N=1600
cells from 250 to 400, before any new formal MC. Core: 48,800 fits;
core plus unchanged supplement: 61,600 fits. The harness identifier is
`xtdpthresh_study_a_v0935_r400` and requires new runs.

- Stata registry comparison: exactly 20 R fields changed; every other cell
  field matches the prior release. Total: 92 cells at R=400, 24 at R=500.
- The updated registry probe passed, including FD/FOD pairing and 28-shard
  smoke allocation checks.
- 8,118 per-shard allocation comparisons passed across the actual PowerShell
  launcher, monitor and verifier functions, 12 cap settings (including 0,
  250, 251, 400 and 500), and applicable shard counts up to 256.
- Fresh execution/merge smoke: all 116 cells, one replication each, six
  shards, grid(199), refine(4). Frozen-source verification and merge passed;
  completed-run Resume launched zero jobs.
- 2,911 source/result/archive checks passed. Seeds, point estimates, SEs and
  diagnostics in all 116 rows match the earlier validation at the same seed
  (numeric tolerance 1e-12). The estimator and supplement sources are unchanged;
  worker/shard source changes are limited to the harness identifier.
- Original A and B ZIP hashes remain unchanged. B is not revised.

These checks cover allocation and integration, not execution of all 48,800
formal core fits. Evidence: `_audit_mc_r400_20260927/` and
`mc_study_a_v0935/runs/validation_r400_01/`.

## Prior 0.9.35 migration checks (before the R400 revision)

- Core: all 116 cells, two outer replications each (232 fits), grid(199) refine(4), no bootstrap.
  The frozen-source verifier and merger completed successfully.
- Supplement: all 32 cells, two replications each (64 fits), grid=10; robust and Windmeijer, actual jump and kink.
  Merge and hash attestation completed successfully.
- Production geometry: 32 supplementary cells at two replications each (64 fits) with grid=199; the core check also used its registered numerical geometry.
- Core DGP, missingness and seed executable blocks match the old source;
  no old FD results are reused. The migration retained the old registry;
  the subsequent R400 allocation change is tested separately above.
- Supplement DGP exactly matches core Study A's generated id/t/y/q values
  in eight N/T, missingness and kappa configurations.
- Independent Python/stdlib recomputation of both studies' SE-panel counts,
  mean SEs, matched empirical SDs, coverage and effective coverage:
  13,395 assertions passed across A/B and their supplements.
- A deliberately missing-SE fixture retained its valid point estimate.
- Completed supplement workers replay without changing raw CSV bytes.
  The frozen launcher's Resume path starts zero jobs on a completed run.
- A deliberately duplicated-row fixture is rejected and its stale
  attestation is removed.

## Provenance

Command SHA256:
`7B97571397D0FAABC8DBD4BCF081779705DC771F8D91D7BED76B05F4A0C9FAB3`.

The command/predictor/help match the 0.9.35 SSC snapshot.
The estimator and supplement sources match the preceding release. Core
allocation and harness identifiers differ as documented above; the revised
core sources match the fresh R400 smoke snapshot.
Supplement worker/merger match the final validation snapshot; launcher
hardening was checked through fresh launch, merge and frozen-run resume.
Documentation was finalized after execution without further code changes.

Logs and raw test rows remain in the project workspace, excluded from this
source-only VPS archive:

- `mc_study_a_v0935/runs/validation_v0935_02`
- `mc_study_a_v0935/supplement_runs/supp_final04`
- `mc_study_a_v0935/supplement_runs/supp_prod03`
- `_audit_mc_v0935_20260927/`
- `_audit_mc_r400_20260927/`
- `mc_study_a_v0935/runs/validation_r400_01/`

Earlier failed development-smoke attempts are retained, not presented as
successful tests. The final tests above supersede those attempts.

## Remaining execution

Run the full production-geometry pilot on the target VPS to choose a safe
shard count. Local testing did not certify 28 simultaneous processes on the
VPS. Then create new formal RunIds with no smoke overrides.
All requested replications, including failures, must remain in the relevant
denominators. Do not transfer smoke/pilot estimates into manuscript tables.

