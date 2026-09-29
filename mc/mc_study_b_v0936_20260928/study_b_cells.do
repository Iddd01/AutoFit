*! study_b_cells.do -- frozen Study B inference registry (study_b_type4_v0936)
*!
*! Caller must define program _cell with this positional interface:
*!   PAIR_ID STUDY DESIGN METHOD N T MISS MISSP QSTATUS R KAPPA DGP_SPEC
*!   BOOTTYPE COEFBOOT GRID GRIDCI B REFINE INFMODE
*!   MAXLAG_LO MAXLAG_HI TRIM GRIDTYPE GRIDSAMPLE PAIRMODE
*!
*! Study B evaluates the INFERENTIAL validity of the wild-bootstrap grid-
*! inversion procedure: threshold-location CI coverage and length, and the
*! empirical size/power of the linearity and continuity tests.  It is the
*! inference counterpart to Study A (point estimation).  Point estimates are
*! distributionally comparable to Study A on the shared point-estimation path.
*!
*! Common to the paired B1--B5 validation cells:
*!   B=499, fixed grid=199 with refine(4), coefboot(none), wild bootstrap,
*!   maxlag(1 3), trim(.15), gridtype(uniform), gridsample(effective),
*!   and FD/FOD common outer DGP/missingness draws with synchronized bootstrap
*!   seeds (realized inner mappings may differ with effective unit sets).
*! B6 is a deliberately separate FD-only published-design calibration bridge:
*!   B=500, grid/gridci=46, refine=0, maxlag(1 5), trim(.10),
*!   gridtype(quantile), and gridsample(observed).  It is not an exact
*!   reproduction of Gong-Seo Algorithm 1. Version 0.9.36 additionally scores
*!   the pointwise test at the true gamma via citest(.25).
*!
*! INFMODE:
*!   FULL     -> threshold CI + linearity + continuity bootstrap tests.
*!   COVERAGE -> threshold CI only (-notest-); used for descriptive N=800 CI
*!               behavior, so the (expensive, grid-profile) tests are not run.
*!   POINT    -> -noboot- (no CI, no tests); reserved, not used in Study B.
*!
*! kappa in {0,.1,.5,1}: continuity null -> near-continuity -> intermediate
*! jump -> clear jump.  The kappa=.5 B1 anchor closes the wide inference gap
*! between .1 and 1 in the published Gong-Seo response surface.  The
*! linearity-size block uses DESIGN=none, a GENUINE no-threshold DGP in which
*! the entire regime term is zero (not the kappa=0 kink), so the reported size
*! is the size of the linearity test under a true linear null.

version 15.0

* ============================ B1 =====================================
* Core inference at N=400: threshold-CI coverage + CI length, linearity and
* continuity size/power, from continuity (k=0) through near-continuity (k=.1)
* and the intermediate anchor (k=.5) to a clear jump (k=1), balanced and with
* 30% interior gaps.  R=500, FULL.
foreach k in 0 .1 .5 1 {
    local d = cond(`k'==0, "kink", "jump")
    local kk : display %03.0f 100*`k'
    foreach m in fd fod {
        _cell b1_bal_n400_k`kk' GSBAL `d' `m' 400 6 balanced 0 pred 500 `k' ///
            gs26_base_v2 wild none 199 100 499 4 FULL ///
            1 3 .15 uniform effective paired
    }
    foreach m in fd fod {
        _cell b1_gap30_n400_k`kk' GSGAP `d' `m' 400 6 mcar .30 pred 500 `k' ///
            gs26_base_v2 wild none 199 100 499 4 FULL ///
            1 3 .15 uniform effective paired
    }
}

* ============================ B2 =====================================
* Descriptive N=800 coverage check.  At kappa=0 use FULL so continuity size
* and linearity power can also be reported at N=800.  At kappa=.1,1 use
* COVERAGE.  R=500 is not precise enough to claim a small convergence trend
* toward .95 across independent N cells; the target is whether conservative
* coverage and delivery remain acceptable at N=800.
foreach k in 0 .1 1 {
    local d = cond(`k'==0, "kink", "jump")
    local kk : display %03.0f 100*`k'
    local mode = cond(`k'==0, "FULL", "COVERAGE")
    foreach m in fd fod {
        _cell b2_bal_n800_k`kk' GSBAL `d' `m' 800 6 balanced 0 pred 500 `k' ///
            gs26_base_v2 wild none 199 100 499 4 `mode' ///
            1 3 .15 uniform effective paired
    }
    foreach m in fd fod {
        _cell b2_gap30_n800_k`kk' GSGAP `d' `m' 800 6 mcar .30 pred 500 `k' ///
            gs26_base_v2 wild none 199 100 499 4 `mode' ///
            1 3 .15 uniform effective paired
    }
}

* ============================ B3 =====================================
* Empirical SIZE of the linearity test under a genuine no-threshold DGP
* (DESIGN=none: regime term identically zero).  N=400 and 800, balanced and
* gap.  R=500 in each cell (MCSE ~0.97pp at nominal 5%).  FULL.  GRIDCI=10 is
* deliberate: gamma-CI quantities are not target estimands under no threshold;
* the linearity test still uses the registered point grid, GRID=199.
foreach nn in 400 800 {
    foreach ms in balanced mcar {
        local p   = cond("`ms'"=="mcar", .30, 0)
        local tag = cond("`ms'"=="mcar", "gap30", "bal")
        foreach m in fd fod {
            _cell b3_`tag'_n`nn' LIN none `m' `nn' 6 `ms' `p' pred 500 0 ///
                gs26_base_v2 wild none 199 10 499 4 FULL ///
                1 3 .15 uniform effective paired
        }
    }
}

* ============================ B4 =====================================
* Endogenous-q robustness: does inference still work when q is
* contemporaneously endogenous (estimated with endogenous(q))?  Balanced and
* interior-gap panels separate endogeneity from missingness.  Apart from the
* contemporaneous (e_t,u_t) correlation, gs26_endog_v2 keeps the published
* benchmark constants.  Continuity null and clear jump; R=500; FULL.
foreach k in 0 1 {
    local d = cond(`k'==0, "kink", "jump")
    local kk : display %03.0f 100*`k'
    foreach ms in balanced mcar {
        local p   = cond("`ms'"=="mcar", .30, 0)
        local tag = cond("`ms'"=="mcar", "gap30", "bal")
        foreach m in fd fod {
            _cell b4_endog_`tag'_n400_k`kk' ENDOG `d' `m' 400 6 `ms' `p' endog 500 `k' ///
                gs26_endog_v2 wild none 199 100 499 4 FULL ///
                1 3 .15 uniform effective paired
        }
    }
}

* ============================ B5 =====================================
* Monotone-attrition inference stress for the estimator's intended
* unbalanced-panel use case.  Attrition targets 15% expected calendar-period
* loss, retains at least four observations per unit, and is independent of the
* DGP innovations.  Continuity null and clear jump; N=400; R=500; FULL.
foreach k in 0 1 {
    local d = cond(`k'==0, "kink", "jump")
    local kk : display %03.0f 100*`k'
    foreach m in fd fod {
        _cell b5_attr15_n400_k`kk' ATTR `d' `m' 400 6 attrition .15 pred 500 `k' ///
            gs26_base_v2 wild none 199 100 499 4 FULL ///
            1 3 .15 uniform effective paired
    }
}

* ============================ B6 =====================================
* FD-only published-design calibration bridge.  These cells move the shipped
* fast wild-bootstrap implementation toward Gong-Seo's main Monte Carlo
* geometry: balanced T=6 panels, N={400,800}, kappa={0,1}, B=500, the full
* available lag range maxlag(1 5) preserving q's predetermined lag-1 moments,
* and 46 empirical-quantile grid points spanning p5--p95.  The published
* t=3..6 lag ladder has 24 data-driven columns; this command additionally
* retains one constant in each of the four realized equation-time blocks, so
* its frozen e(N_iv) contract is 28.  A potential t=2 block is indexed from
* t_min+1 but is all-zero under dynamic FD and is pruned.  COVERAGE avoids
* unrelated tests.  Because the realized moment matrix, bootstrap algorithm,
* first-stage weight, and confidence-set inversion differ from the published
* implementation, this block supports only a quantitative
* calibration/consistency claim, never an exact-replication claim.
foreach nn in 400 800 {
    foreach k in 0 1 {
        local d = cond(`k'==0, "kink", "jump")
        local kk : display %03.0f 100*`k'
        _cell b6_cal_bal_n`nn'_k`kk' GSCAL `d' fd `nn' 6 balanced 0 pred 500 `k' ///
            gs26_base_v2 wild none 46 46 500 0 COVERAGE ///
            1 5 .10 quantile observed fdonly
    }
}
