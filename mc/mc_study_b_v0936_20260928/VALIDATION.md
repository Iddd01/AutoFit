# Validation record — Study B 0.9.36

Validated locally on 28 September 2026 with Stata 17/MP.

## Estimator audit inherited by this package

The targeted 0.9.36 audit is stored in `_audit_0936_20260928/REPORT.md` at the
project root. It verified continuity-test batch/scalar agreement, independent
recalculation of `citest()`, seed isolation, FD/FOD/kink edge cases, and exact
0.9.35 equality for outputs that were not intended to change.

## Core harness smoke

- Run ID: `smoke_local_0936b`
- Settings: 52 cells, one replication per cell, six shards, B=19,
  grid=10, gridci=10.
- Result: 52/52 rows completed; the fail-closed verifier and merger passed.
- Identified-threshold cells: `citest()` requested 44/44, returned 44/44,
  evaluable 44/44, and component seed correct 44/44.
- Linear-null cells: 8/8 correctly omitted `citest()`.
- All eight published output files and their attestation were generated.

## Supplement smoke

- Run ID: `supp_smoke_local_0936`
- Settings: eight cells, one replication per cell, four shards, B=19,
  grid=10, gridci=10.
- Result: 8/8 rows completed and verified merge passed.
- `citest()` was evaluable and its component seed was correct in 8/8 rows.

Smoke results test execution, schema, accounting, seed and merge contracts.
They are not Monte Carlo performance evidence. Formal results require new runs
with no smoke overrides.
