*! study_a_cells.do -- Study A v6 point-estimation registry
*! Registry version: xtdpthresh_study_a_v0937_r400 (2026-08-17)
*!
*! Caller contract (positional):
*!   _cell PAIR_ID STUDY DESIGN METHOD N T MISS MISSP QSTATUS R KAPPA ///
*!         DGP_SPEC GRID REFINE
*!
*! Six targeted blocks, 116 cells and 48,800 replication rows.
*! FD and FOD in a pair receive identical DGP and missingness streams.
*! This registry is point-only: study_a_worker.do always supplies -noboot-.
*! v6 uses xtdpthresh 0.9.37 fixed grid(199) plus refine(4) in every cell,
*! removes adaptive-search stopping rules, and retains the matched
*! balanced controls for the endogenous-q and heavy-tail blocks. It also
*! matches the monotone-attrition block's expected missing fraction to severe
*! interior MCAR (hazard .15 vs p .30).  DESIGN=kink labels the continuous
*! DGP at kappa=0; it does not request the command's -kink- restriction.

version 15.0

* A1. Gong--Seo benchmark family: continuity through a clear jump, the
* published N/T geometry, balanced panels and 30% interior MCAR gaps.
* All regular cells, including N=1600, use R=400; hard stresses use R=500.
foreach n in 400 800 1600 {
    local RN = 400
    foreach k in 0 .1 .2 .5 1 {
        if `k' == 0   local kk 000
        if `k' == .1  local kk 010
        if `k' == .2  local kk 020
        if `k' == .5  local kk 050
        if `k' == 1   local kk 100
        local d = cond(`k' == 0, "kink", "jump")

        foreach m in fd fod {
            _cell a1_bal_n`n'_k`kk' GSCORE `d' `m' `n' 6 balanced 0 pred ///
                `RN' `k' gs26_base_v2 199 4
            _cell a1_gap30_n`n'_k`kk' GSCORE `d' `m' `n' 6 mcar .30 pred ///
                `RN' `k' gs26_base_v2 199 4
        }
    }
}

* A2. Small-N / longer-panel block.  This prevents the evidence from relying
* only on N>=400 and checks the same estimator at T=10.
foreach n in 100 200 {
    foreach k in 0 1 {
        local kk = cond(`k' == 0, "000", "100")
        local d  = cond(`k' == 0, "kink", "jump")
        foreach m in fd fod {
            _cell a2_bal_n`n'_t10_k`kk' GSSMALL `d' `m' `n' 10 balanced 0 pred ///
                400 `k' gs26_base_v2 199 4
            _cell a2_gap30_n`n'_t10_k`kk' GSSMALL `d' `m' `n' 10 mcar .30 pred ///
                400 `k' gs26_base_v2 199 4
        }
    }
}

* A3. Persistent-y stress: rho_y=.9 and a longer burn-in.  Other parameters
* stay at the Gong--Seo values so the stress isolates persistence.
foreach k in 0 1 {
    local kk = cond(`k' == 0, "000", "100")
    local d  = cond(`k' == 0, "kink", "jump")
    foreach m in fd fod {
        _cell a3_bal_n400_k`kk' PERSIST `d' `m' 400 6 balanced 0 pred ///
            500 `k' gs26_persist_v2 199 4
        _cell a3_gap30_n400_k`kk' PERSIST `d' `m' 400 6 mcar .30 pred ///
            500 `k' gs26_persist_v2 199 4
    }
}

* A4. Missingness-pattern block at T=10: balanced, mild/severe interior MCAR,
* and independent monotone attrition.  Attrition begins only after four rows.
foreach k in 0 1 {
    local kk = cond(`k' == 0, "000", "100")
    local d  = cond(`k' == 0, "kink", "jump")
    foreach m in fd fod {
        _cell a4_bal_n400_t10_k`kk' MISSPAT `d' `m' 400 10 balanced 0 pred ///
            400 `k' gs26_base_v2 199 4
        _cell a4_gap15_n400_t10_k`kk' MISSPAT `d' `m' 400 10 mcar .15 pred ///
            400 `k' gs26_base_v2 199 4
        _cell a4_gap30_n400_t10_k`kk' MISSPAT `d' `m' 400 10 mcar .30 pred ///
            400 `k' gs26_base_v2 199 4
        _cell a4_attr15_n400_t10_k`kk' MISSPAT `d' `m' 400 10 attrit .15 pred ///
            400 `k' gs26_base_v2 199 4
    }
}

* A5. Contemporaneously endogenous q, following Gong--Seo Appendix C.3:
* corr(e_t,u_t)=.5 and endogenous(q), with other parameters unchanged.
foreach k in 0 1 {
    local kk = cond(`k' == 0, "000", "100")
    local d  = cond(`k' == 0, "kink", "jump")
    foreach m in fd fod {
        _cell a5_endog_bal_n400_k`kk' ENDOG `d' `m' 400 6 balanced 0 endog ///
            500 `k' gs26_endog_v2 199 4
        _cell a5_endog_gap30_n400_k`kk' ENDOG `d' `m' 400 6 mcar .30 endog ///
            500 `k' gs26_endog_v2 199 4
    }
}

* A6. Heavy-tail robustness.  Standardised t(5) innovations preserve the
* declared variances and rho_eu=.5 while relaxing joint Gaussianity.
foreach k in 0 1 {
    local kk = cond(`k' == 0, "000", "100")
    local d  = cond(`k' == 0, "kink", "jump")
    foreach m in fd fod {
        _cell a6_t5_bal_n400_k`kk' HEAVYTAIL `d' `m' 400 6 balanced 0 pred ///
            500 `k' gs26_t5_v2 199 4
        _cell a6_t5_gap30_n400_k`kk' HEAVYTAIL `d' `m' 400 6 mcar .30 pred ///
            500 `k' gs26_t5_v2 199 4
    }
}

* Registry total: 60 + 16 + 8 + 16 + 8 + 8 = 116 cells.
* 92 cells use R=400 (including 20 N=1600 cells); 24 stresses use R=500:
* 48,800 method-estimation rows (24,400 pair-replications).
