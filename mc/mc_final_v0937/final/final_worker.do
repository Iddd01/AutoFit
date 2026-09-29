*! Final Monte Carlo worker (xtdpthresh 0.9.37). One shard of a registry
*! (final_point_cells.csv or final_inf_cells.csv, frozen as final_cells.csv).
*! DGP: Gong-Seo (2026) benchmark, as in the checking studies:
*!   y_it = rho y_i,t-1 + q_it + (kappa - .5 + 2 q_it) 1(q_it > .25) + e_it,
*!   q_it = .7 q_i,t-1 + u_it, e = .5 e*, u = .5 e*_{t-1} + sqrt(.75) v;
*! dgp: base; endog (u = .5 e*_t + sqrt(.75) v, q declared endogenous);
*! persist (rho = .9, 50 burn-in periods); heavy (standardized t(5));
*! linear (no regime term). Twenty burn-in periods otherwise; T periods kept.
*! Seeds: DGP and missingness seeds follow the Study B scheme (method, kappa,
*! c, model, vce, instrument set and block are absent), so FD/FOD and all
*! cells of a DGP, N, T and missingness pattern share the outer sample.
*! Bootstrap seed: own namespace (424242); method, kappa, c, model and vce
*! are absent, so FD/FOD share the multiplier draws.
version 17.0
clear all
set more off
set varabbrev off
set processors 1
args SHARD
do final_config.do
set rng mt64
adopath ++ "."
quietly do xtdpthresh.ado
local OUT "final_SH`SHARD'.csv"
local CELLVARS block cell_id pair_id dgp spec method vce iv N T miss missp kappa c mode R B grid gridci gridtype gridsample trim refine
local RESVARS master shard rep dgp_seed missing_seed boot_seed units_realized analysis_observed rc version_ok twostep gamma_hat b_rho b_q b_cons b_qd b_rhod se_rho se_q se_cons se_qd se_rhod sc_rho sc_q sc_cons sc_qd sc_rhod joint_vce vce_applied ar_joint hansen_p ar1_p ar2_p N_obs N_iv N_units ci_delivered ci_lo ci_hi ci_nseg ci_boundary ci_incomplete citest_gamma citest_accept citest_p citest_D citest_crit citest_status citest_draws seed_citest p_lin lin_valid elapsed_s
local HEADER "run_id"
foreach nm in `CELLVARS' `RESVARS' {
    local HEADER "`HEADER',`nm'"
}
local done
capture confirm file "`OUT'"
if !_rc {
    tempname READ
    file open `READ' using "`OUT'", read text
    file read `READ' first
    file close `READ'
    if `"`first'"' != `"`HEADER'"' {
        di as err "final_worker: `OUT' has a different header; use a new RunId"
        exit 459
    }
    import delimited using "`OUT'", clear case(preserve) varnames(1) stringcols(_all)
    if _N {
        destring cell_id rep shard, replace
        isid cell_id rep
        assert shard==`SHARD' & run_id=="${fin_run}"
        forvalues j=1/`=_N' {
            local done "`done' `=cell_id[`j']':`=rep[`j']'"
        }
    }
}
else {
    tempname HH
    file open `HH' using "`OUT'", write text
    file write `HH' "`HEADER'" _n
    file close `HH'
}
import delimited using final_cells.csv, clear case(preserve) asdouble ///
    stringcols(1 3 4 5 6 7 8 11 15 20 21)
local ncells = _N
forvalues j=1/`ncells' {
    foreach nm of local CELLVARS {
        local c`j'_`nm' = `nm'[`j']
    }
}
local expected = 0
local ordinal = 0
forvalues j=1/`ncells' {
    foreach nm of local CELLVARS {
        local `nm' "`c`j'_`nm''"
    }
    if ${fin_repcap}>0 local R=min(`R',${fin_repcap})
    if ${fin_B}>0 & `B'>0 local B=${fin_B}
    if ${fin_grid}>0 local grid=${fin_grid}
    if ${fin_gridci}>0 & `gridci'>0 local gridci=${fin_gridci}
    * DGP constants
    local RHOY  = cond("`dgp'"=="persist", .9, .6)
    local TBURN = cond("`dgp'"=="persist", 50, 20)
    local TMAX  = `T' + `TBURN'
    local T_qd    = cond("`dgp'"=="linear", 0, 2)
    local T_consd = cond("`dgp'"=="linear", 0, `kappa' - .5)
    local GAMMA0  = cond("`dgp'"=="linear", ., .25)
    local SPCODE = cond(inlist("`dgp'","base","linear"),1,cond("`dgp'"=="endog",2,cond("`dgp'"=="persist",3,4)))
    local QCODE  = cond("`dgp'"=="endog",2,1)
    local MCODE  = cond("`miss'"=="balanced",1,cond("`miss'"=="mcar",2,3))
    local PCODE  = round(`missp'*1000000)
    local IVCODE = cond("`iv'"=="all",1,cond("`iv'"=="L3",2,3))
    * estimation options
    local qopt = cond("`dgp'"=="endog", "endogenous(q)", "predetermined(q)")
    local ivopt ""
    if "`iv'"=="L3" local ivopt "maxlag(1 3)"
    if "`iv'"=="collapse" local ivopt "collapse"
    local specopt = cond("`spec'"=="kink", "kink", "")
    local refopt = cond(`refine'>0, "refine(`refine')", "")
    local cit = .
    if inlist("`mode'","CI","FULL") local cit = .25 + `c'
    forvalues rep=1/`R' {
        local ordinal=`ordinal'+1
        if mod(`ordinal'-1,${fin_nshard})+1 != `SHARD' continue
        local expected=`expected'+1
        if strpos(" `done' "," `cell_id':`rep' ") continue
        local _h0 = mod(${fin_master}*1000003 + `rep'*10007 + `N'*101 + ///
            `T'*1009 + `SPCODE'*1000033 + `QCODE'*1000037, 699999999)
        local dgp_seed = 1 + `_h0'
        local missing_seed = 700000001 + ///
            mod(`_h0'*1009 + `MCODE'*10007 + `PCODE'*101 + 7919, 699999999)
        local boot_seed = 1400000001 + mod(`_h0'*1013 + `MCODE'*10007 + ///
            `PCODE'*101 + `B'*131 + `grid'*137 + `gridci'*139 + `IVCODE'*151 + ///
            424242, 699999999)
        timer clear 9
        timer on 9
        *----- generate panel -----
        clear
        quietly set seed `dgp_seed'
        quietly set obs `=`N'*`TMAX''
        quietly gen long id = ceil(_n/`TMAX')
        quietly bysort id: gen int t = _n
        quietly xtset id t
        if "`dgp'"=="heavy" {
            quietly gen double es = rt(5)/sqrt(5/3)
            quietly gen double e  = es*.5
            quietly gen double w  = rt(5)/sqrt(5/3)
        }
        else {
            quietly gen double e  = rnormal()*.5
            quietly gen double es = e/.5
            quietly gen double w  = rnormal()
        }
        quietly gen double u = .
        if "`dgp'"=="endog" {
            quietly replace u = .5*es + sqrt(.75)*w
        }
        else {
            quietly bysort id (t): replace u = .5*es[_n-1] + sqrt(.75)*w if _n > 1
            quietly bysort id (t): replace u = w if _n == 1
        }
        quietly gen double q = .
        quietly bysort id (t): replace q = rnormal() if _n == 1
        quietly bysort id (t): replace q = .7*q[_n-1] + u if _n > 1
        tempvar g
        quietly gen double `g' = (`T_consd' + `T_qd'*q)*(q > .25)
        quietly gen double y = .
        quietly bysort id (t): replace y = rnormal() if _n == 1
        quietly bysort id (t): replace y = `RHOY'*y[_n-1] + q + `g' + e if _n > 1
        quietly gen byte analysis = (t > `TBURN')
        quietly drop if t <= `TBURN'
        quietly xtset id t
        *----- missingness (as study_b_worker.do) -----
        if "`miss'" == "mcar" & `missp' > 0 {
            quietly set seed `missing_seed'
            quietly bysort id (t): egen double _tmin = min(cond(analysis, t, .))
            quietly bysort id (t): egen double _tmax = max(cond(analysis, t, .))
            quietly gen double _du = runiform() if analysis
            quietly drop if analysis & _du < `missp' & t != _tmin & t != _tmax
            quietly drop _du _tmin _tmax
            quietly bysort id: egen int _na = total(analysis)
            quietly drop if _na < 4
            quietly drop _na
            quietly xtset id t
        }
        else if "`miss'" == "attrition" {
            quietly set seed `missing_seed'
            quietly gen double _au = runiform() if t == `=`TBURN'+1'
            quietly bysort id (t): egen double _a = max(_au)
            quietly gen double _du = runiform() if t == `=`TBURN'+1'
            quietly bysort id (t): egen double _d0 = max(_du)
            local DTARGET = 2*`missp'*`T'
            local DLO = floor(`DTARGET')
            local DHI = ceil(`DTARGET')
            local DMAX = `T' - 4
            if `DLO' < 1 local DLO = 1
            if `DHI' < 1 local DHI = 1
            if `DLO' > `DMAX' local DLO = `DMAX'
            if `DHI' > `DMAX' local DHI = `DMAX'
            local PHI = cond(`DHI'==`DLO', 0, (`DTARGET'-`DLO')/(`DHI'-`DLO'))
            quietly gen int _d = `DLO' + (`DHI'-`DLO')*(_d0 < `PHI') if _a < 0.5
            quietly drop if analysis & _a < 0.5 & t > `=`TBURN'+`T'' - _d
            quietly drop _au _a _du _d0 _d
            quietly bysort id: egen int _na = total(analysis)
            quietly drop if _na < 4
            quietly drop _na
            quietly xtset id t
        }
        quietly count if analysis
        local analysis_observed = r(N)
        tempvar _utag
        quietly egen byte `_utag' = tag(id) if analysis
        quietly count if `_utag' == 1
        local units_realized = r(N)
        *----- estimate -----
        if "`mode'"=="POINT" local infopt "noboot"
        else {
            local infopt "boot(`B') gridci(`gridci') boottype(wild) rseed(`boot_seed')"
            if "`mode'"=="CI" local infopt "`infopt' notest"
            if !missing(`cit') local infopt "`infopt' citest(`cit')"
        }
        capture quietly xtdpthresh y if analysis, qx(q) `qopt' method(`method') ///
            `ivopt' grid(`grid') gridtype(`gridtype') gridsample(`gridsample') ///
            trim(`trim') `refopt' `specopt' `infopt' coefboot(none) ///
            history(panel) level(95) vce(`vce') bwscale(1.5) nowarn
        local rc=_rc
        foreach nm in version_ok twostep gamma_hat b_rho b_q b_cons b_qd b_rhod ///
            se_rho se_q se_cons se_qd se_rhod sc_rho sc_q sc_cons sc_qd sc_rhod ///
            joint_vce vce_applied ar_joint hansen_p ar1_p ar2_p N_obs N_iv N_units ///
            ci_delivered ci_lo ci_hi ci_nseg ci_boundary ci_incomplete ///
            citest_gamma citest_accept citest_p citest_D citest_crit ///
            citest_status citest_draws seed_citest p_lin lin_valid {
            local `nm'=.
        }
        if `rc'==0 {
            local version_ok=("`e(cmd)'"=="xtdpthresh" & "`e(cmdversion)'"=="0.9.37")
            if !`version_ok' {
                di as err "final_worker: expected xtdpthresh 0.9.37; got `e(cmd)' `e(cmdversion)'"
                exit 459
            }
            assert "`e(method)'"=="`method'" & e(flag_kink)==("`spec'"=="kink")
            local twostep=e(estimator_twostep)
            local gamma_hat=e(gamma)
            foreach nm in joint_vce vce_applied ar_joint hansen_p ar1_p ar2_p N_iv N_units {
                capture local `nm'=e(`nm')
            }
            local N_obs=e(N)
            * coefficients: jump (5) or kink (3); SE of the reported type and
            * the covariance conditional on gamma-hat
            capture matrix VC=e(V_cond)
            local hasvc=(_rc==0)
            if "`spec'"=="kink" {
                local names "rho q qd"
                local labels "Lag_y_b q_b kink_slope"
            }
            else {
                local names "rho q cons qd rhod"
                local labels "Lag_y_b q_b cons_d q_d Lag_y_d"
            }
            local ii=0
            foreach nm of local names {
                local ++ii
                local label: word `ii' of `labels'
                capture local b_`nm'=_b[`label']
                capture local se_`nm'=_se[`label']
                if `hasvc' {
                    local col=colnumb(VC,"`label'")
                    if !missing(`col') {
                        if VC[`col',`col']>=0 & !missing(VC[`col',`col']) ///
                            local sc_`nm'=sqrt(VC[`col',`col'])
                    }
                }
            }
            if "`mode'"!="POINT" {
                local ci_incomplete=e(ci_incomplete)
                local ci_nseg=e(ci_nseg)
                local ci_boundary=e(boundary_warn)
                capture matrix CI=e(ci_segments)
                if !_rc & `ci_incomplete'==0 {
                    local ci_delivered=0
                    forvalues rr=1/`=rowsof(CI)' {
                        if !missing(CI[`rr',1],CI[`rr',2]) local ci_delivered=1
                    }
                    if `ci_delivered' {
                        local ci_lo=e(gamma_lo)
                        local ci_hi=e(gamma_hi)
                        assert !missing(`ci_lo',`ci_hi') & `ci_lo'<=`ci_hi' & ///
                            `ci_lo'<=`gamma_hat'+1e-12*max(1,abs(`gamma_hat')) & ///
                            `gamma_hat'-1e-12*max(1,abs(`gamma_hat'))<=`ci_hi'
                    }
                }
                else local ci_delivered=0
                foreach nm in citest_gamma citest_accept citest_p citest_D ///
                    citest_crit citest_status citest_draws seed_citest {
                    capture local `nm'=e(`nm')
                }
                if !missing(`cit') {
                    assert abs(`citest_gamma'-`cit')<=1e-12 & ///
                        inrange(`citest_status',1,6) & ///
                        `seed_citest'==mod(`boot_seed'+477377,2147483648)
                    if `citest_status'==1 {
                        assert `citest_draws'==`B' & inlist(`citest_accept',0,1) & ///
                            `citest_accept'==(`citest_D'<=`citest_crit') & ///
                            `citest_accept'==(`citest_p'>.05)
                    }
                    if `citest_status'==2 {
                        assert `citest_D'==0 & `citest_accept'==1 & `citest_p'==1
                    }
                }
                else assert missing(`citest_status')
                if inlist("`mode'","FULL","LIN") {
                    local p_lin=e(pval_lin)
                    capture local lin_valid=e(boot_linearity_valid)
                    if !missing(`p_lin') assert `lin_valid'==`B'
                }
                else assert missing(e(pval_lin))
            }
        }
        timer off 9
        quietly timer list 9
        local elapsed_s=r(t9)
        local master ${fin_master}
        local shard `SHARD'
        local row "${fin_run}"
        foreach nm in `CELLVARS' `RESVARS' {
            local row "`row',``nm''"
        }
        tempname FH
        file open `FH' using "`OUT'", write append text
        file write `FH' "`row'" _n
        file close `FH'
        di "FINAL_ROW cell=`cell_id' rep=`rep' rc=`rc'"
    }
}
tempname DONE
file open `DONE' using "final_done_SH`SHARD'.txt", write replace text
file write `DONE' "${fin_run},`SHARD',`expected'" _n
file close `DONE'
di "FINAL_COMPLETE shard=`SHARD' rows=`expected'"
