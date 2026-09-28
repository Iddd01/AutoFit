# Frozen design: `xtdpthresh_study_a_v0935_r400`

## Purpose and scope

Study A evaluates the finite-sample point-estimation behavior of the FOD
estimator supplied by `xtdpthresh` 0.9.35, using FD as the published-theory
benchmark. It targets the editor's request for substantially broader Monte
Carlo evidence on the FOD extension. System GMM is outside the revised scope.

The design spans sample size, panel length, jump strength, interior gaps,
monotone attrition, persistence, endogenous threshold status, and
non-Gaussian shocks. It does not run threshold confidence sets, bootstrap
tests, or coefficient-bootstrap intervals; those tasks belong to Study B.

The evidence is described as finite-sample convergence patterns, not as a
proof of consistency or an estimate of an asymptotic convergence rate.

## Migration to 0.9.35

The core retains the 116-cell registry, DGP parameters, replication allocation
and seed formulas. Both transformations are rerun; no old FD results are reused.
The new covariance and AR fields distinguish joint inference, conditional
fallback and the V_cond diagnostic. Point success does not require SE delivery.
Additional targeted robust/Windmeijer and restricted-kink runs are documented
separately in SUPPLEMENT.md and never enter the core estimand denominators.

## Data-generating process

The baseline recursion follows the Gong--Seo (2026) Monte Carlo structure:

\[
y_{it}=\rho_y y_{i,t-1}+q_{it}
 +(\delta_0+2q_{it})1\{q_{it}>0.25\}+\varepsilon_{it},
\]

\[
q_{it}=0.7q_{i,t-1}+u_{it},\qquad
\varepsilon_{it}=0.5e^*_{it},\qquad
u_{it}=0.5e^*_{i,t-1}+\sqrt{0.75}\,v_{it},
\]

where `e*` and `v` are independent, mean-zero, variance-one innovations. Thus
`var(u)=1`, `sigma_e=.5`, and the effective lagged shock correlation is `.5`.
The baseline has `rho_y=.6`, 20 burn-in periods, and a predetermined `q` in
the estimation command.

There are `T=6` or `T=10` retained observations, with **no extra
presample \(y_{i0}\) row**: this is the published Gong--Seo sample geometry, in
which the observed panel is \(y_{i1},\dots,y_{iT}\) and the instrument set
\(z_{it}=(y_{i,t-2},\dots,y_{i1},q_{i,t-1},\dots,q_{i1})'\) runs over
\(t=t_0,\dots,T\) with \(t_0=3\), giving \(3+5+7+9=24\) instruments at \(T=6\).
That count describes the uncapped published set. Study A instead fixes
`maxlag(1 3)` in every cell to control instrument growth, so the reported
`N_iv` need not equal 24 or match between FD and FOD. `T` counts observations,
not estimating equations. Study B uses the same retained-panel rule. The
actual jump at the threshold is

\[
\kappa=\delta_0+2(0.25),\qquad \delta_0=\kappa-0.5.
\]

The five jump-model coefficients are `Lag_y_b`, `q_b`, `cons_d`, `q_d`, and
`Lag_y_d`, with truths `rho_y`, `1`, `kappa-.5`, `2`, and `0`. The DGP has no
`x`, so `x_b` and `x_d` are deliberately outside Study A.

`kappa=0` is a continuous threshold DGP.  The estimator remains the
unrestricted jump estimator; the `kink` option is not imposed.  Positive
`kappa` values span near continuity through a strong discontinuity.

### Stress DGPs

- `gs26_base_v2`: baseline Gaussian, weakly exogenous/predetermined `q`.
- `gs26_persist_v2`: changes only `rho_y` to `.9` and uses a 50-period
  burn-in.
- `gs26_endog_v2`: follows Gong--Seo Supplement C.3 by setting
  `u_it=.5 e*_it+sqrt(.75)v_it` and declaring `endogenous(q)`; all other
  baseline parameters, including `rho_y=.6`, remain unchanged.
- `gs26_t5_v2`: replaces `e*` and `v` with independent standardized `t(5)`
  innovations, preserving their declared variances and covariance structure.

## Missing-data mechanisms

Missingness is imposed only after the complete panel has been generated and
is independent of its values.

- `balanced`: no deletion.
- `mcar .15` or `mcar .30`: independently delete eligible interior retained
  observations while protecting each unit's first and last retained rows;
  remove a unit if fewer than four rows remain.
- `attrit .15`: an independent `.15` dropout hazard begins at retained row
  five; once dropout occurs, that row and all subsequent rows are removed.
  Every unit therefore retains at least four rows.

In the `T=10` A4 comparison, the protected-endpoint/unit-retention MCAR rule
with `p=.30` has an expected realized missing fraction of about 24.0%, while
the monotone rule with hazard `.15` has an expected fraction of about 24.7%.
Thus A4 compares mechanisms at approximately matched information loss rather
than confounding pattern with the much larger loss induced by hazard `.30`.
The `.15` and `.30` MCAR cells use the same missingness-uniform stream, so the
mild deletion set is nested inside the severe deletion set; attrition uses a
separate stream.

The probability is a deletion or hazard parameter, not the realized fraction
missing. Every replication records observed rows, units retained/dropped,
realized missing rate, gap events, and gap periods.

## Frozen 116-cell registry

| Block | Cells | Design |
|---|---:|---|
| A1 `GSCORE` | 60 | `N=400,800,1600`; `T=6`; `kappa=0,.1,.2,.5,1`; balanced/MCAR `.30`; FD/FOD |
| A2 `GSSMALL` | 16 | `N=100,200`; `T=10`; `kappa=0,1`; balanced/MCAR `.30`; FD/FOD |
| A3 `PERSIST` | 8 | `N=400`; `T=6`; `rho_y=.9`; `kappa=0,1`; balanced/MCAR `.30`; FD/FOD |
| A4 `MISSPAT` | 16 | `N=400`; `T=10`; `kappa=0,1`; balanced, MCAR `.15/.30`, attrition `.15`; FD/FOD |
| A5 `ENDOG` | 8 | `N=400`; `T=6`; endogenous `q`; balanced/MCAR `.30`; `kappa=0,1`; FD/FOD |
| A6 `HEAVYTAIL` | 8 | `N=400`; `T=6`; standardized `t(5)` shocks; balanced/MCAR `.30`; `kappa=0,1`; FD/FOD |

All 92 regular cells, including the 20 `N=1600` A1 cells, use `R=400`;
the 24 high-variance stress cells in A3, A5, and A6 use `R=500`. This yields
exactly 48,800 method-estimation rows, 24,400 paired DGP replications, and 58
FD/FOD pairs. All cells estimate the unrestricted jump specification and use
`grid(199)`, `refine(4)`, `trim(.15)`, `history(panel)`,
`vce(robust)`, `maxlag(1 3)`, and unconditional
`noboot coefboot(none)`.

The global grid is fixed and disclosed. Four local refinement iterations can
add at most 120 observed support values around its optimum. This is an ex-ante
compute/accuracy compromise for the 48,800-call study, not a claim of
exhaustive support search. A historical 0.9.26 preflight on 96 paired draws found that
the nested 397-point arm cost about 2.92 times as much as 199 points without
improving threshold accuracy in that diagnostic sample. The 397-point arm is
therefore retained only as historical numerical-sensitivity evidence, not as
a new 0.9.35 grid-sensitivity experiment. The grid choice is held fixed across
the migration rather than selected again using the new formal results.

The raw output records requested/effective/admitted grid counts,
stage-specific spans, remaining refinement candidates, neighbourhood
completeness, and whether the reported threshold lies on an admitted-grid
boundary. `refine_complete` never enters the primary-sample selection rule:
a finite-grid estimate remains defined if four local rounds exhaust their
cap. Completion and boundary rates must accompany the loss summaries.

The registry's `design=kink` value is a DGP label for `kappa=0`, not an
estimation option. The worker does not pass `kink`, requires
`e(flag_kink)=0`, and therefore retains the same five estimable coefficients
in every cell. This distinction also makes support-point refinement valid in
all 116 cells.

The merger reports Monte Carlo uncertainty explicitly. Bias MCSE is
`SD(error)/sqrt(R)`. RMSE MCSE is the delta-method quantity
`SD(error^2)/(2*RMSE*sqrt(R))`; proportion MCSEs use the corresponding
binomial expression. Conclusions must be calibrated to these MCSEs.
At a `.05` rejection probability, the
binomial MCSE is approximately 1.09 percentage points for `R=400`, .97 for
`R=500`. Bias MCSE equals 5% and 4.47% of the
corresponding empirical error SD. Under a Gaussian benchmark, the relative
delta-method MCSE of RMSE is approximately 3.54% and 3.16%, respectively.

The allocation is fixed ex ante. Regular designs receive `R=400`; persistence,
threshold endogeneity, and heavy tails receive `R=500` because they are more
variable and more failure-prone. The r400 revision raises `N=1600` from
the historical `R=250` to `R=400`, adding 3,000 fits before the new formal
run; all other design cells, seeds and numerical settings are unchanged.
Common-random-number paired contrasts are
typically more precise than unpaired comparisons. This balances simulation
precision, reviewer credibility, and the need to reserve substantial compute
for bootstrap Study B.

### FD/FOD comparison estimand

Method is excluded from both DGP and missingness seeds. Within each `pair_id`
and replication, FD and FOD receive the identical generated panel and
missingness draw. Paired differences use this common-random-number design.

FD and FOD nevertheless invoke their own transformation, first-stage weight,
usable transformed rows, and instrument stack. With `history(panel)` and the
fixed lag cap, their `N_trans` and `N_iv` can differ even on a balanced panel.
The paired result therefore compares two complete implemented estimator
configurations; it does not hold every moment and weight fixed while changing
only an algebraic transformation.

## Outcomes

For the threshold and each of the five coefficients, cell outputs report
effective counts, bias, RMSE, empirical SD, estimate quantiles, and Monte
Carlo standard errors. Threshold summaries include the 5th, 25th, 50th, 75th,
and 95th percentiles and the concentration probability
`Pr(|gamma_hat-gamma_0| <= .05)`. Coefficient summaries additionally report mean
analytic robust SE, mean-SE/empirical-SD ratio, and analytic 95% Wald
coverage. Hansen and AR(2) rejection rates are descriptive implementation
diagnostics. The paired report evaluates within-replication FOD-minus-FD
differences in signed error, absolute error, and squared error.

The primary patterns are:

1. bias and RMSE across sample size and panel length;
2. paired FD--FOD behavior under balanced panels, interior gaps, and attrition;
3. behavior from continuity (`kappa=0`) to a strong jump (`kappa=1`);
4. sensitivity to persistence, threshold endogeneity, and heavy tails.

The analytic coefficient coverage and specification diagnostics are not
bootstrap validation. No threshold confidence set, threshold coverage,
bootstrap rejection rate, or coefficient-bootstrap claim may be made from
Study A. Those quantities belong to Study B.

### Estimation-failure reporting

Failure and contract-delivery rates are first-class outcomes. A delivered
estimate is an `rc=0`, command/version/option/shape-valid result; delivery does
not silently relabel a one-step fallback as the intended two-step estimator.
The frozen primary success indicator additionally requires
`e(estimator_twostep)=1`. Coefficient bias and RMSE are conditional on this
primary success sample and must always be reported beside their effective
counts. `refine(4)` itself defines the finite search algorithm;
`e(refine_complete)` is reported separately and is not used to select the
loss sample, because such post-estimation conditioning would mechanically
select against denser large-`N` support. Direct FD--FOD efficiency comparisons
use the common subset of paired replications in which both methods meet the
primary contract; marginal method-specific RMSEs may be shown only when
labeled as such. To diagnose selection caused by failure to complete two-step
GMM, cell summaries also provide (i) a delivered-sample panel that mixes the
intended two-step estimate with any contract-valid one-step fallback and (ii)
a fallback-only panel. Both panels report separate effective counts and are
secondary sensitivity analyses, not substitutes for the frozen primary
two-step estimand. The merger separately reports delivery, two-step,
refinement-completion, and admitted-grid-boundary rates. No failed, fallback,
or refinement-capped replication is selectively regenerated.

## Reproducibility contract

The formal harness, worker schema, code version, RNG, master seed, staged file
hashes, shard count and run ID are immutable within a run.  Resume skips
verified `cell_id+rep` keys.  Failed fits remain Monte Carlo outcomes and are
not selectively retried.  The merger fails if any expected cell, replication,
shard, completion marker, provenance field, or unique key is missing or mixed.
