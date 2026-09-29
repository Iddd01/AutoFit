# Frozen design: Study B inference validation (`study_b_type4_v0936`)

Version 0.9.36 preserves the 52-cell registry and adds `citest(.25)` to every
identified-threshold cell. Its acceptance indicator is the primary pointwise
threshold-coverage outcome. Coverage computed from the displayed union of
`e(ci_segments)` is retained only as a descriptive set-geometry measure.

## Decision

Study B is the inferential complement to Study A. It evaluates the shipped
`xtdpthresh` 0.9.36 wild-bootstrap procedure for FD and FOD; it does not
re-estimate the broad point-estimator surface already covered by Study A and it
does not claim to reproduce Gong-Seo Algorithm 1 exactly.

The formal registry contains **52 cells and 26,000 outer Monte Carlo
replications**. B1--B5 use **B=499**, fixed `grid(199)`, `refine(4)`,
`coefboot(none)`, `trim(.15)`, `maxlag(1 3)`, `gridtype(uniform)`, and
`gridsample(effective)`. B6 is a four-cell FD-only calibration bridge using
**B=500**, `grid(46)`, `gridci(46)`, `refine(0)`, `trim(.10)`,
`maxlag(1 5)`, `gridtype(quantile)`, and `gridsample(observed)`. Every cell
uses `history(panel)` and one Stata processor per shard. There are 24 paired
FD/FOD designs with identical outer samples and synchronized bootstrap seeds,
plus four FD-only calibration cells. Realized inner bootstrap mappings need not
be identical when FD and FOD retain different effective unit sets.

This design is frozen because it answers the referee's inferential questions
while keeping the bootstrap workload defensible on the 32-core, 24-GB VPS.
The v8r1 execution contract preserves the two targeted additions introduced in
v6—B5 monotone attrition and B6 published-design calibration—while migrating
all estimation to 0.9.36 fixed-grid search. Adaptive-search controls and
System-GMM paths are outside Study B. The finite grid is part of the reported
estimator, as in the comparison commands; `refine(4)` only adds local
resolution after that global search. A refinement-incomplete draw remains a
valid recorded Monte Carlo outcome and is reported diagnostically rather than
silently retried or deleted.
Calibration calls are heavier than their row share because they use all 24
lag columns plus four block constants, so row shares are not runtime shares.

## What answers the referee

1. **Extensive FOD evidence.** All 24 paired substantive designs are estimated by
   both FD and FOD with identical outer DGP/missingness draws and synchronized
   bootstrap seeds.
2. **Bootstrap behavior is evaluated, not asserted.** The outputs report CI
   delivery, complete-set coverage and length, test delivery, valid bootstrap
   draws, size, power, and unconditional effective performance.
3. **Nulls are genuine.** Linearity size uses a true no-threshold DGP;
   continuity size uses a continuous kink (`kappa=0`).
4. **Endogeneity is not confounded with gaps.** The endogenous-q block has
   both balanced and 30% interior-gap panels.
5. **Sample-size behavior is visible without a convergence over-claim.**
   N=400 and N=800 are used for CI behavior and linearity size; continuity
   size and linearity power at `kappa=0` are also repeated at N=800. Because
   the N-specific cells are independent and R=500 gives about 0.97 percentage
   point MCSE near 95% coverage, B2 describes whether conservative coverage
   persists at N=800; it does not claim to identify a small trend toward .95.
6. **Failures cannot disappear.** The merger reports all requested,
   contract-valid, delivered, and valid-test denominators. A one-step fallback
   is retained as a delivered command result but is not labeled two-step
   success.
7. **Unbalanced-panel inference is not represented by interior MCAR alone.**
   B5 evaluates CI and test behavior under monotone attrition at the two
   registered endpoints, `kappa=0` and `kappa=1`. It does not establish
   near-continuity performance under attrition.
8. **An external numerical bridge is explicit.** B6 aligns the FD design more
   closely with Gong-Seo's published Monte Carlo geometry while preserving a
   strict no-exact-replication label.

## Exact implemented bootstrap

The formal design uses `boottype(wild)`, the fast procedure actually shipped
in 0.9.36:

- each bootstrap draw assigns **one Mammen two-point weight to each panel
  unit**, shared by all retained observations/equations of that unit;
- confidence-set inversion holds the reported fit's second-step weight fixed
  when two-step estimation succeeds, or its first-step weight after a fallback;
  each bootstrap profile is solved under that fixed weight, not by rebuilding it;
- the linearity test keeps its fixed first-step criterion;
- the continuity test compares kink and jump under the fixed second-step
  weight of the jump fit after two-step estimation (and W1 after fallback),
  and its null bootstrap outcome combines the kink fit with jump residuals;
- all threshold candidates share the same unit-level wild weights; both
  sample and draw unrestricted minima use the documented common initial grid;
- the finite-B critical value uses ceil(p*(B+1)), and every requested draw must
  be valid before the corresponding test or confidence set is delivered;
- the threshold confidence set inverts the profile GMM-distance statistic at
  each candidate gamma, imposing that candidate in the restricted fit;
- the linearity statistic is the restricted-linear one-step objective minus
  the minimum unrestricted one-step profile objective, bootstrapped under the
  linear null;
- the continuity statistic is the nested kink-versus-jump profile
  GMM-distance on their common feasible grid, bootstrapped under the kink null;
- p-values use the command's add-one correction and are delivered only when
  its valid-draw gate is met.

This statistic is **not** the Seo-Shin sup-Wald statistic. The wild procedure
is related to the fast `xthenreg` convention but is an empirical approximation
to Gong-Seo's test-inversion framework. The package correctly stores
`e(ci_bootstrap_certified)=0`; Study B supplies the finite-sample evidence
needed for that approximation. The experimental `boottype(unit)` path is
outside the formal registry because it is not certified as exact Algorithm 1,
is FD-only, and is more than an order of magnitude slower, so it cannot provide
the required FD/FOD comparison.

## Published benchmark DGP

The balanced benchmark follows the local Gong-Seo (2026) paper and supplement:

- T=6 observed periods, N in {400,800};
- `rho_y=.6`, `beta_q=1`, `q_d=2`, `gamma0=.25`;
- `sigma_e=.5`, `rho_q=.7`, and correlation parameter `rho_eu=.5`;
- no additional x regressor and no individual effect in the standard mode.

The implementation simulates 20 burn-in periods before retaining the six
analysis rows. This is an initialization/stabilization choice; burn-in rows
are dropped and cannot enter `history(panel)` as extra instruments.
The registered inference DGP uses Gaussian homoskedastic innovations. Study B
therefore validates the shipped bootstrap under this benchmark class; it does
not by itself establish heavy-tail or conditional-heteroskedastic robustness.

For predetermined q, the q innovation loads on the preceding standardized
outcome shock. For endogenous q, it loads on the contemporaneous shock. The
endogenous block changes that timing and the matching `endogenous(q)`
declaration only: **it keeps `rho_y=.6` and `rho_q=.7`**, as in the
published appendix. Correspondingly, the transformed-equation q instruments
begin at lag 1 in the predetermined block and lag 2 in the endogenous block.

`T=6` means six observed rows per unit, not six estimating equations. With
the published t0=3 geometry, the first valid differenced equation is t=3.
The main B1--B5 design is a close benchmark, not an exact replication. It
deliberately uses `maxlag(1 3)` rather than Gong-Seo's full 24-column lag
ladder. This is a standard instrument-proliferation control, not a design
weakness: Seo-Shin (2016, footnote 17) likewise cap the maximum lag order at
four and cite Roodman (2009). The worker records the realized command-specific
instrument count for every replication. Study B also adds FOD, interior gaps,
monotone attrition, and the shipped wild approximation.

B6 is a separate FD-only **published-design calibration bridge**. It keeps the
full available lag range `maxlag(1 5)`, which is essential because q is
predetermined and its valid ladder begins at lag 1. `maxlag(2 5)` is therefore
not an admissible way to target a numerical instrument count: it would delete
q's lag-1 moments and silently impose the endogenous-q timing instead.

The command's block constructor starts at `t_min+1` (=2), whereas the
published convention sets `t0=3`. Under the dynamic-FD complete-case rule the
potential t=2 block has no realized transformed row, so its all-zero columns
are pruned. Across the surviving t=3,...,6 blocks, the GS lag ladder contributes
3+5+7+9=24 data-driven columns, but `xtdpthresh` also retains one block
constant at each equation time. Consequently, B6's frozen command-specific
contract is **`e(N_iv)=28`**, and the worker and merger fail closed on any other
count. This is precisely why B6 calibrates against a published design but does
not claim to reproduce the published moment matrix.

Its 46 quantile-grid points span p5--p95, matching the quantile sequence in the
official Monte Carlo code more closely than the 81-point p10--p90 grid used in
the paper's empirical application. B6 still computes the displayed
confidence-set inversion, but 0.9.36 also evaluates `citest(.25)` directly;
the latter is the registered pointwise coverage score.

The authors' public R code also uses a q-innovation loading of
`1/sqrt(2)`, whereas the published DGP states `rho_eu=.5`. B6 keeps `.5` so it
remains comparable with Study A and the main Study B registry. Together with
the shipped wild bootstrap and the command's different first-stage weighting,
this prevents any claim that B6 exactly reproduces Table 1. The defensible
claim is quantitative consistency under a closely aligned published design.

The numerical targets are frozen before the VPS run. From Gong-Seo Table 1,
the Grid-B coverage references used by B6 are:

| N | kappa=0 | kappa=1 |
|---:|---:|---:|
| 400 | 0.992 | 0.966 |
| 800 | 0.986 | 0.955 |

The merger reports pointwise `citest()` coverage and effective pointwise
coverage differences from these references, together with each study's own
Monte Carlo standard error. The segment-union set score remains descriptive.
Because B6 and Gong-Seo do not use identical
DGP/moment/bootstrap/inversion procedures, these differences are descriptive
calibration gaps—not an equivalence test—and no combined-MCSE acceptance rule
is used.
Gong-Seo's main simulation is also broader and more intensive
(N in {400,800,1600}, jump sizes in {0,.1,.2,.5,1}, 2,000 outer repetitions,
and B=500). Study B retains the two core sample sizes, kappa=.1 for the most
near-continuous stress, kappa=.5 as the intermediate response-surface anchor,
kappa=1 as the clear jump, and R=500. The main B1--B5 cells use B=499; B6 uses
B=500. The remaining differences make Study B a validation of the shipped
command under a published benchmark DGP, not a numerical reproduction of
Gong-Seo's tables.
The second published near-continuity value, kappa=.2, is not duplicated here:
kappa=.1 is the more demanding continuity-neighborhood stress, while kappa=.5
adds the missing intermediate power/coverage anchor at lower incremental cost.

The formal standard-DGP label is `gs26_base_v2`, shared with Study A. The
legacy Study-B-only label `gs26_official_v1` denotes the same recursion and is
accepted by the worker for manual compatibility, but it is not used in the v8r1
registry. This prevents a cross-study comparison from mistaking two labels for
two benchmark DGPs.

Study B has no internally generated alternative threshold-location CI such as
Gong-Seo's NP-B, NP-B(S), or the Seo-Shin asymptotic interval. Its FD/FOD
comparison answers the referee's transformation question, while published
comparisons with those threshold-CI methods must be cited from Gong-Seo's
tables and must not be described as reproduced by this registry.

## Missing-panel intervention

The `gap30` condition is generated after burn-in. Conditional on the six
observed calendar periods, each interior row is independently selected for
deletion with probability .30; the first and last rows are protected, and a
unit is removed if fewer than four observations remain. The deletion draw is
independent of outcomes, regressors, innovations, and potential bootstrap
weights, so this is an MCAR stress design with deliberately created interior
gaps. Because endpoints are protected and short units are removed, `.30`
describes the interior-row deletion probability, not the final whole-panel
missing fraction. The worker records both the realized missing rate and the
FD/FOD transformation opportunities for every replication.

The B5 `attr15` condition is monotone. Half of units are selected for
attrition independently of all DGP and bootstrap draws. Conditional on being
selected, a unit loses one or two terminal periods with probabilities chosen
so that the expected loss is 1.8 of six periods; hence the unconditional
expected calendar-period loss is exactly 15%. Every unit retains at least four
observations. The registry records `miss=attrition, missp=.15`, and the worker
records realized missingness and FD/FOD transformation opportunities. B5's
registered evidence is limited to `kappa` in {0,1}; claims about attrition must
not be extended to the near-continuity case.

## Frozen registry

| Block | Scientific target | Cells | R per cell | Rows | Mode |
|---|---|---:|---:|---:|---|
| B1 | N=400 CI coverage/length; linearity power; continuity size/power; balanced and gap30; kappa in {0,.1,.5,1} | 16 | 500 | 8,000 | FULL |
| B2 | Descriptive N=800 CI behavior for kappa in {0,.1,1}; at kappa=0 also continuity size and linearity power | 12 | 500 | 6,000 | FULL at kappa=0; otherwise COVERAGE |
| B3 | Linearity size under a genuine no-threshold DGP; N in {400,800}; balanced and gap30 | 8 | 500 | 4,000 | FULL |
| B4 | Endogenous-q robustness separated from missingness; balanced and gap30; kappa in {0,1} | 8 | 500 | 4,000 | FULL |
| B5 | Monotone-attrition inference at N=400; kappa in {0,1}; paired FD/FOD | 4 | 500 | 2,000 | FULL |
| B6 | FD-only published-design calibration; balanced N in {400,800}; kappa in {0,1}; B=500, 46-point p5--p95 quantile grids; `maxlag(1 5)`; 28 command instruments (24 lag columns + 4 block constants) | 4 | 500 | 2,000 | COVERAGE |

Totals: **52 cells, 26,000 rows, 20,000 FULL rows, 6,000
coverage-only rows; 24 FD/FOD pairs and four FD-only cells**. B3 uses
`grid(199)` for the linearity statistic but
`gridci(10)`: gamma has no population meaning under the no-threshold null, so
B3's automatically generated gamma-CI, continuity, and nuisance point/regime
outputs are not target estimands and are excluded from scientific summaries by
the merger. B1, B2, B4, and B5 use `gridci(100)`; B6 uses `gridci(46)`.

### Estimand map

- Linearity size: B3 only (4,000 rows).
- Linearity power: FULL threshold cells in B1, B2 at kappa=0, B4, and B5
  (16,000 rows).
- Continuity size: FULL threshold cells with kappa=0 (7,000 rows).
- Continuity power: FULL threshold cells with kappa>0 (9,000 rows).
- Gamma-CI coverage/length: all threshold cells; never B3 (22,000 rows).
- Published-design calibration: B6 only (2,000 FD rows); summarize separately
  and do not pool it with the main B1--B5 response surface.

The merger materializes these target flags instead of treating every quantity
returned by a FULL command as a scientific estimand.

## Why B=499/500 and these outer R values

Gong-Seo use B=500, while Davidson-MacKinnon identify 399 as roughly the
minimum useful count for a 5% bootstrap test. B1--B5 use B=499, giving the
conventional `1/(B+1)` resolution and substantially improving on the former
B=299 design. B6 uses B=500 to align that calibration dimension exactly.
With the command's add-one p-value and a strict `p<.05` decision, B=499 has the
finite-B randomization reference 24/500=0.048, not exactly 0.050. Tables report
that finite-B reference alongside the requested nominal level.
Although the help file recommends B>=999 for a single final empirical
application, this is a Monte Carlo performance study: bootstrap noise is
averaged over hundreds of independent outer replications, and matching the
published B approximately 500 convention is the relevant cost/precision
compromise.

At p=.05 or .95, MCSE is `sqrt(p(1-p)/R)`:

| R | MCSE |
|---:|---:|
| 500 | 0.97 percentage points |

Every reported rate must be accompanied by its realized denominator and MCSE.
The xthenreg paper's 500-iteration linearity experiment uses a warp-speed
one-bootstrap-draw-per-outer-iteration construction and is therefore not a
precision or cost precedent for this nested B approximately 500 design. Gong-Seo use R=2,000
with B=500, and their endogenous appendix uses R=1,000 with B=500. Here R=500
is an explicit computational constraint: it estimates each cell-level 5% or
95% rate to about 0.97 percentage points MCSE, which is adequate for coverage,
size, and broad power patterns but not for detecting sub-percentage-point
changes across independent N cells.

The evidence base used to freeze these choices is:

- `ref/1-s2.0-S0304407625002064-main.pdf` and
  `ref/1-s2.0-S0304407625002064-mmc1.pdf` (Gong-Seo main paper and supplement);
- `ref/Estimation of dynamic panel threshold model using Stata.pdf`
  (the published `xthenreg` implementation and its simulation convention);
- `ref/seo2016.pdf` (the published lag-cap precedent in footnote 17);
- `ref/Bootstrap tests How many bootstraps.pdf` (bootstrap replication-count
  guidance);
- `ref/mammen1992 Bootstrap and wild bootstrap for high dimensional linear
  random design models.pdf` (two-point wild weights).
- Gong and Seo's official replication repository at commit
  `7c65e2199d7f05ed489fd7f0533d60deffb384eb`
  (`https://github.com/wsggong/DPTR_bootstrap`), pinned on 2026-08-21 for the
  main-MC grid, resampling algorithm, and code-versus-prose DGP cross-check.

## Success, delivery, and denominators

For each outer replication:

- `contract_ok=1` means the pinned 0.9.36 returned the required command, coefficient,
  sample, and diagnostic contract;
- `success=1` additionally requires `e(estimator_twostep)=1`;
- `one_step_fallback=1` means a contract-valid estimate was delivered after
  the second step failed.

Inference delivery/coverage is reported for all contract-valid delivered
results and again for the strict two-step subset. Effective rates retain all
requested outer replications in the denominator. Incomplete CI inversions are
never converted into empty sets or noncoverage; their formal confidence set is
withdrawn and non-delivery is explicit.

CI length is summarized by its mean, median, 95th percentile, and a Monte
Carlo standard error for the mean, alongside incomplete, empty, disconnected,
and boundary-pinned set rates with explicit denominators. Paired FD/FOD output
reports both joint-delivery contrasts and unconditional effective contrasts,
so a method cannot look better merely by failing to deliver difficult draws.

## Reproducibility and merge contract

- A formal run is immutable by run ID, harness/schema, 0.9.36 code version,
  RNG family, master seed, shard count, registry, and SHA-256 inventory.
- FD/FOD omit method from all three worker seed formulas. The bootstrap hash
  includes numerical tuning fields, preventing B1/B6 stream reuse. The merger
  recomputes the formulas and verifies pair equality and cross-pair uniqueness.
- Failed fits are completed Monte Carlo outcomes and are never selectively
  retried.
- Every replication is a separate open/append/close CSV transaction. Resume
  may repair one demonstrably torn final physical line, but fails closed on
  non-final corruption and revalidates schema, provenance, seeds, outcomes,
  and accounting before skipping an existing key.
- The VPS launcher staggers shard initialization by three seconds and uses a
  600-second hard gate requiring at least 1.5 GB of both available physical RAM
  and Windows commit headroom (`Commit Limit - Committed Bytes`). If the gate
  expires, existing shards continue and no additional shard is opened; a later
  selective resume starts only dead/pending shards. This changes only launch
  timing, never the registry, allocation, or random streams.
- Direct execution of `merge_study_b.do` is denied. The supported
  `verify_and_merge_study_b.ps1` entrypoint executes the registry probe,
  verifies staged/raw hashes and the Stata executable hash, grants a one-use
  nonce, and publishes one complete eight-file output generation through an
  atomic JSON pointer and nonce-specific attestation.
- The verifier recognizes a formal run only when `Master=20260814`, `RepCap=0`,
  and all three B/grid overrides are dots. Non-formal smoke/sensitivity output requires the
  explicit `-AllowNonFormal` switch and is stamped as such in its attestation.
- Formal results must use manifest overrides `B='.'`, `Grid='.'`, and
  `GridCI='.'`. Overrides are for smoke/sensitivity runs and are recorded.

Study A and Study B are independent Monte Carlo runs. Common random numbers are
guaranteed within Study B FD/FOD pairs, not across the two studies.


## Migration and additional restricted-kink block

Both FD and FOD are rerun on 0.9.36. Earlier confidence sets used a different
criterion, critical-value rule and candidate-weight coupling, so the old
rows cannot serve as 0.9.36 inference results. Point/SE/AR telemetry now
distinguishes joint inference from conditional fallback. The merger records
separate coefficient-SE panels with their own valid-SE denominators.
B3 coefficient coverage remains a non-target, as do its threshold sets.

The core's continuous-DGP cells still estimate unrestricted jump.
SUPPLEMENT.md adds eight separately labelled actual-kink CI cells
(N=400/800, balanced/interior gaps, FD/FOD), 4,000 outer fits in total.
These supplement rows do not change the 52-cell core or its target flags.
