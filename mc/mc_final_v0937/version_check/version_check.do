*! version_check.do -- are the published Monte Carlo numbers still those of
*! xtdpthresh 0.9.37?
*!   Part A (Study A): 0.9.35 vs 0.9.37, -noboot-, Study A options: e(b), e(V),
*!     e(V_cond), gamma, objective, Hansen/AR, instrument and sample counts.
*!   Part B (Study B and supplements): 0.9.36 vs 0.9.37 with the bootstrap:
*!     threshold set, linearity p-value, citest(), continuity p-value
*!     (0.9.37 with conttest, the option that now requests it).
*! Samples: Gong-Seo benchmark DGP (q predetermined), T=6, N=400; balanced and
*! MCAR .30 gaps; FD and FOD; jump kappa=1, jump kappa=0, imposed kink; vce
*! robust and windmeijer. Run from this folder in Stata 17:
*!     do version_check.do
*! Result: VERSION_CHECK_PASS, or the first difference found.
version 17.0
clear all
set more off
set varabbrev off
set processors 1
set rng mt64
local NREP_A 10
local NREP_B 3

capture program drop _gendata
program define _gendata
    args N missp kappa seed
    clear
    quietly set seed `seed'
    quietly set obs `=`N'*26'
    quietly gen long id = ceil(_n/26)
    quietly bysort id: gen int t = _n
    quietly gen double e  = rnormal()*.5
    quietly gen double es = e/.5
    quietly gen double w  = rnormal()
    quietly gen double u = .
    quietly bysort id (t): replace u = .5*es[_n-1] + sqrt(.75)*w if _n > 1
    quietly bysort id (t): replace u = w if _n == 1
    quietly gen double q = .
    quietly bysort id (t): replace q = rnormal() if _n == 1
    quietly bysort id (t): replace q = .7*q[_n-1] + u if _n > 1
    quietly gen double g = (`kappa' - .5 + 2*q)*(q > .25)
    quietly gen double y = .
    quietly bysort id (t): replace y = rnormal() if _n == 1
    quietly bysort id (t): replace y = .6*y[_n-1] + q + g + e if _n > 1
    quietly drop if t <= 20
    if `missp' > 0 {
        quietly gen double du = runiform()
        quietly drop if du < `missp' & t != 21 & t != 26
        quietly bysort id: gen int nr = _N
        quietly drop if nr < 4
    }
    quietly xtset id t
end

capture program drop _loadver
program define _loadver
    args file
    foreach pr in xtdpthresh _xdpt_tsdouble _xdpt_tsterm {
        capture program drop `pr'
    }
    mata: mata clear
    quietly do `file'
end

* ---------------- configurations ----------------
* Part A: method x panel x spec x vce = 2 x 2 x 3 x 2 = 24 configurations
local ka 0
foreach m in fd fod {
    foreach p in 0 .3 {
        foreach s in jump1 jump0 kink {
            foreach v in robust windmeijer {
                local ++ka
                local A`ka'_m `m'
                local A`ka'_p `p'
                local A`ka'_s `s'
                local A`ka'_v `v'
            }
        }
    }
}
* Part B: method x panel x kappa = 2 x 2 x 2 = 8 configurations
local kb 0
foreach m in fd fod {
    foreach p in 0 .3 {
        foreach k in 0 1 {
            local ++kb
            local B`kb'_m `m'
            local B`kb'_p `p'
            local B`kb'_k `k'
        }
    }
}

* ---------------- fits ----------------
foreach ver in 35 36 37 {
    _loadver xtdpthresh_v09`ver'.ado
    if inlist(`ver', 35, 37) {
        forvalues j = 1/`ka' {
            local kap = cond("`A`j'_s'"=="jump1", 1, 0)
            local mopt = cond("`A`j'_s'"=="kink", "kink", "refine(4)")
            forvalues r = 1/`NREP_A' {
                _gendata 400 `A`j'_p' `kap' `=1000*`j'+`r''
                quietly xtdpthresh y, qx(q) predetermined(q) method(`A`j'_m') ///
                    maxlag(1 3) grid(199) gridtype(uniform) gridsample(effective) ///
                    trim(.15) `mopt' noboot coefboot(none) history(panel) ///
                    vce(`A`j'_v') bwscale(1.5) nowarn
                assert "`e(cmdversion)'" == "0.9.`ver'"
                matrix Ab`ver'_`j'_`r' = e(b)
                matrix AV`ver'_`j'_`r' = e(V)
                capture matrix AC`ver'_`j'_`r' = e(V_cond)
                if _rc matrix AC`ver'_`j'_`r' = J(1,1,0)
                matrix AS`ver'_`j'_`r' = (e(gamma), e(obj), e(hansen), ///
                    e(hansen_p), e(ar1_p), e(ar2_p), e(N), e(N_trans), ///
                    e(N_iv), e(N_units), e(estimator_twostep), e(joint_vce), ///
                    e(vce_applied))
            }
        }
    }
    if inlist(`ver', 36, 37) {
        * 0.9.37 runs the continuity test only on request (conttest)
        local copt = cond(`ver' == 37, "conttest", "")
        forvalues j = 1/`kb' {
            forvalues r = 1/`NREP_B' {
                _gendata 400 `B`j'_p' `B`j'_k' `=5000*`j'+`r''
                quietly xtdpthresh y, qx(q) predetermined(q) method(`B`j'_m') ///
                    maxlag(1 3) grid(199) refine(4) gridtype(uniform) ///
                    gridsample(effective) trim(.15) gridci(20) boot(99) ///
                    boottype(wild) rseed(`=777*`j'+`r'') citest(.25) ///
                    coefboot(none) history(panel) vce(robust) nowarn `copt'
                assert "`e(cmdversion)'" == "0.9.`ver'"
                matrix Bb`ver'_`j'_`r' = e(b)
                capture matrix BG`ver'_`j'_`r' = e(ci_segments)
                if _rc matrix BG`ver'_`j'_`r' = J(1,2,.)
                matrix BS`ver'_`j'_`r' = (e(gamma), e(gamma_lo), e(gamma_hi), ///
                    e(ci_nseg), e(ci_incomplete), e(pval_lin), e(citest_D), ///
                    e(citest_crit), e(citest_accept), e(citest_p), e(citest_status))
                matrix BC`ver'_`j'_`r' = (e(pval_cont))
            }
        }
    }
}

* ---------------- comparisons ----------------
* mreldif of two matrices, treating equal missing entries as equal
capture program drop _cmp
program define _cmp, rclass
    args A B
    confirm matrix `A'
    confirm matrix `B'
    if rowsof(`A') != rowsof(`B') | colsof(`A') != colsof(`B') {
        return scalar d = .
        exit
    }
    local d 0
    forvalues i = 1/`=rowsof(`A')' {
        forvalues k = 1/`=colsof(`A')' {
            local a = `A'[`i',`k']
            local b = `B'[`i',`k']
            if missing(`a') & missing(`b') continue
            if missing(`a') | missing(`b') {
                return scalar d = .
                exit
            }
            local d = max(`d', reldif(`a', `b'))
        }
    }
    return scalar d = `d'
end

local tol 1e-10
local fail 0
local worstA 0
forvalues j = 1/`ka' {
    forvalues r = 1/`NREP_A' {
        foreach mm in Ab AV AC AS {
            _cmp `mm'35_`j'_`r' `mm'37_`j'_`r'
            if missing(r(d)) | r(d) > `tol' {
                di as err "Study A diff: `mm' config `j' (`A`j'_m' missp=`A`j'_p' `A`j'_s' `A`j'_v') rep `r': " r(d)
                local fail 1
            }
            else local worstA = max(`worstA', r(d))
        }
    }
}
local worstB 0
local worstC 0
forvalues j = 1/`kb' {
    forvalues r = 1/`NREP_B' {
        foreach mm in Bb BG BS {
            _cmp `mm'36_`j'_`r' `mm'37_`j'_`r'
            if missing(r(d)) | r(d) > `tol' {
                di as err "Study B diff: `mm' config `j' (`B`j'_m' missp=`B`j'_p' kappa=`B`j'_k') rep `r': " r(d)
                local fail 1
            }
            else local worstB = max(`worstB', r(d))
        }
        _cmp BC36_`j'_`r' BC37_`j'_`r'
        if missing(r(d)) | r(d) > `tol' {
            di as err "continuity p-value differs: config `j' rep `r'"
            local fail 1
        }
        else local worstC = max(`worstC', r(d))
    }
}
di as txt "Part A (0.9.35 vs 0.9.37, `=`ka'*`NREP_A'' fits): max reldif = " %9.2e `worstA'
di as txt "Part B (0.9.36 vs 0.9.37, `=`kb'*`NREP_B'' fits): max reldif = " %9.2e `worstB'
di as txt "continuity p-value (0.9.36 vs 0.9.37): max reldif = " %9.2e `worstC'
if `fail' {
    di as err "VERSION_CHECK_FAIL"
    exit 9
}
di as res "VERSION_CHECK_PASS"
