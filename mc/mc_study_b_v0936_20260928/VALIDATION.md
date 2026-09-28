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

## Post-smoke change (pending re-smoke)

After the smokes above, the core and supplement workers/mergers gained the
reported-interval columns (`gamma_lo`, `gamma_hi`, `hull_covered`,
`hull_length`; supplement `ci_*` equivalents) and their summaries, and the
supplement no longer passes `citest()` when `B=0`. The estimator is unchanged.
Both smokes must be rerun on this package before any formal run.
Both smokes were rerun and passed (`smoke_b0936`: 52/52 rows, verified
merge; `supp_smoke_0936c`: 8/8 rows, verified merge).

## Post-formal-launch hardening (for later runs)

A code review found that several "all fields empty" checks used Stata's
`missing(a,b,...)`, which is true when ANY argument is missing, so they
accepted rows with some fields filled. They now test each field
(`missing(a) & missing(b) ...`, or `!missing(a) | !missing(b) ...` for the
flag form). Data written by the worker are unaffected; only the strength of
the validation changes. The staged copies inside an already launched run are
not modified.

Smoke results test execution, schema, accounting, seed and merge contracts.
They are not Monte Carlo performance evidence. Formal results require new runs
with no smoke overrides.

## Supplement 2 (redesigned 28 September 2026)

Supplement 2 now measures the power of the threshold test (SUPP2.md): 48 cells,
24,000 fits, `citest(.25 + c)` on the Study B core samples. The earlier
continuity block (0.9.35 vs 0.9.36) and `xtdpthresh_v0935.ado` were removed
before any run, because the continuity test is outside the article's scope.
`supp2_link_check.do` verifies against the formal Study B output that every
replication reproduces the realized panel and gamma-hat of its Study B
counterpart. A smoke (`-RepCap 1 -B 19 -Grid 10 -GridCI 10`) must pass before
the formal run. Core and kink-supplement files are unchanged.

## Launcher date parsing (28 September 2026, after the supp2 formal launch)

`Live-Jobs` in `run_supp2.ps1` and `run_supplement.ps1` parsed the ISO start
time with culture-dependent `[DateTime]::Parse`, which failed on the VPS while
shards were live (`Status` error; running shards unaffected; `Merge` after all
shards finish never reaches the parse). It now parses with the invariant
culture and round-trip kind. Runs launched with the earlier launcher keep it
(frozen by hash); check their progress through the `*_done_SH*.txt` markers.
