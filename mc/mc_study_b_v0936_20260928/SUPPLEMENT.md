# Study B supplement — xtdpthresh 0.9.36

This is a separate, prespecified block, not a change to the core registry.
It uses the exact command snapshot included in this package.

## Design

8 cells, R=500 each (4,000 fits):
N=400/800, T=6, balanced/30% interior deletion, FD/FOD.
Every cell has kappa=0 and is estimated using the actual `kink` restriction.
The primary target is pointwise acceptance of `citest(.25)` at the true
threshold. The displayed threshold set remains a descriptive geometry check.
Use B=499, grid(199), gridci(100), robust, and notest; the main
Study B already evaluates linearity/continuity tests under the unrestricted model.
No continuity test is requested for a model that already imposes continuity.

## Common DGP and numerical settings

The baseline recursion is the same as the core predetermined-q benchmark:
q_t=.7 q_(t-1)+.5 e*_(t-1)+sqrt(.75) w_t;
y_t=.6 y_(t-1)+q_t+(kappa-.5+2q_t)1(q_t>.25)+.5 e*_t.
Innovations are independent standard normal, with 20 discarded burn-in rows.
There is no extra presample observation or individual effect.
Interior deletion protects first/last periods and removes units with fewer
than four retained observations. It is independent of generated values.

The DGP and missingness seed formulas equal core Study A's baseline formulas.
Method and covariance option never enter those seeds. Supplement B uses its
own bootstrap namespace, synchronized across FD/FOD, and is not paired with
the old core B results. Jump uses grid(199) refine(4); kink omits refine()
and uses the command's automatic continuous refinement. Both use
maxlag(1 3), trim(.15), history(panel), and bwscale(1.5).

## Run

From this package directory, on a local disk:
```powershell
# Small numerical smoke; never use these rows as formal MC evidence.
.\run_supplement.ps1 -Action Fresh -RunId supp_smoke -NShard 4 -RepCap 1 -Grid 10 -B 19 -GridCI 10
.\run_supplement.ps1 -Action Status -RunId supp_smoke
.\run_supplement.ps1 -Action Merge -RunId supp_smoke

# Production-geometry pilot; measures runtime, not coverage.
.\run_supplement.ps1 -Action Fresh -RunId supp_pilot -NShard 4 -RepCap 1
# Formal run only after inspecting the pilot and resource use.
.\run_supplement.ps1 -Action Fresh -RunId supp_formal -NShard 16
.\run_supplement.ps1 -Action Status -RunId supp_formal
.\run_supplement.ps1 -Action Merge -RunId supp_formal
```

Use `-Action Resume` with the same RunId after all workers have stopped.
Resume reads the frozen manifest settings, not new command-line overrides.
It skips completed outcomes, including failures; it does not redraw a failed fit.
If the release-folder launcher has changed, invoke the frozen launcher inside
the corresponding supplement_runs/RunId directory instead.
A malformed/torn CSV fails closed and is not silently repaired.
Do not edit a staged run or combine old and new versions.

## Outputs and interpretation

The verified merger produces row-level CSV/DTA, cell summaries, coefficient
summaries and within-draw FD/FOD and robust/Windmeijer contrasts.
An attestation hashes source, raw input and outputs.
Point success requires the intended two-step fit; covariance availability is
separate. Joint coefficient coverage requires joint_vce=1; conditional SE
coverage is a separate diagnostic. Mean SE and empirical SD use the same
available-SE subset. Coverage denominators and MCSEs accompany each panel.
Supplement paired CSVs contain row-level contrasts, not independent draws
across methods or covariance options.

`citest_coverage` is the pointwise coverage outcome. The older `coverage`
field uses the union of e(ci_segments) and is retained for set geometry;
length sums segment lengths. Report non-delivery, unresolved tests,
incomplete sets, disconnection and boundary contact. Smoke/pilot runs are
labelled formal=0. Neither a smoke pass nor MC performance proves a theorem.

