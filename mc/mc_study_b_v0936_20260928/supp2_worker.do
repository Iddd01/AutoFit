*! Supplement 2 (xtdpthresh 0.9.36 package): power of the grid-bootstrap
*! threshold test at a false threshold gamma0 + c (Gong-Seo 2026, Table 2
*! layout), FD vs FOD, balanced / interior gaps / monotone attrition.
*! The outer samples are those of the Study B core cells named in link_cell:
*! the DGP code, the DGP seed and the missingness seed below are copied from
*! study_b_worker.do (gs26_base_v2, q predetermined), so rep r of a cell here
*! is rep r of its link_cell there. supp2_link_check.do verifies this through
*! gamma_hat and the realized sample, and adds the c = 0 column (coverage)
*! from Study B. Independent of the core registry and of the kink supplement.
version 17.0
clear all
set more off
set varabbrev off
set processors 1
args SHARD
do supp2_config.do
set rng mt64
adopath ++ "."
local OUT "supp2_SH`SHARD'.csv"
local HEADER "run_id,block,cell_id,pair_id,link_cell,method,N,T,miss,missp,kappa,c,citest_value,R,B,grid,gridci,refine,trim,maxlag_lo,maxlag_hi,gridtype,gridsample,master,shard,rep,dgp_seed,missing_seed,boot_seed,units_realized,analysis_observed,rc,version_ok,twostep,gamma_hat,N_used,citest_gamma,citest_accept,citest_p,citest_D,citest_crit,citest_status,citest_draws,seed_citest,elapsed_s"
local done
capture confirm file "`OUT'"
if !_rc {
    tempname READ
    file open `READ' using "`OUT'", read text
    file read `READ' first
    file close `READ'
    if `"`first'"' != `"`HEADER'"' {
        di as err "supp2_worker: `OUT' has a different header; use a new RunId"
        exit 459
    }
    import delimited using "`OUT'", clear case(preserve) varnames(1) stringcols(_all)
    if _N {
        destring cell_id rep shard, replace
        isid cell_id rep
        assert shard==`SHARD' & run_id=="${supp_run}"
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
local CELLVARS block cell_id pair_id link_cell method N T miss missp kappa c citest_value R B grid gridci refine trim maxlag_lo maxlag_hi gridtype gridsample
import delimited using supp2_cells.csv, clear case(preserve) asdouble ///
    stringcols(1 3 4 5 8 21 22)
local ncells = _N
forvalues j=1/`ncells' {
    foreach nm of local CELLVARS {
        local c`j'_`nm' = `nm'[`j']
    }
}
quietly do xtdpthresh.ado
* Gong-Seo benchmark constants (study_b_worker.do, gs26_base_v2)
local RHOY  = 0.6
local BQ    = 1.0
local D2    = 2.0
local GAMMA = 0.25
local RHOQ  = 0.7
local SIGE  = 0.5
local TBURN = 20
local expected = 0
local ordinal = 0
forvalues j=1/`ncells' {
    foreach nm of local CELLVARS {
        local `nm' "`c`j'_`nm''"
    }
    if ${supp_repcap}>0 local R=min(`R',${supp_repcap})
    if ${supp_B}>0 local B=${supp_B}
    if ${supp_grid}>0 local grid=${supp_grid}
    if ${supp_gridci}>0 local gridci=${supp_gridci}
    local T_qd = `D2'
    local T_consd = `kappa' - `T_qd'*`GAMMA'
    local TMAX = `T' + `TBURN'
    local MCODE = cond("`miss'"=="balanced",1,cond("`miss'"=="mcar",2,3))
    local PCODE = round(`missp'*1000000)
    forvalues rep=1/`R' {
        local ordinal=`ordinal'+1
        if mod(`ordinal'-1,${supp_nshard})+1 != `SHARD' continue
        local expected=`expected'+1
        if strpos(" `done' "," `cell_id':`rep' ") continue
        * Study B seed contract (SPCODE = QCODE = 1): same outer sample and
        * missingness as rep `rep' of `link_cell'.
        local _h0 = mod(${supp_master}*1000003 + `rep'*10007 + `N'*101 + ///
            `T'*1009 + 1000033 + 1000037, 699999999)
        local dgp_seed = 1 + `_h0'
        local missing_seed = 700000001 + ///
            mod(`_h0'*1009 + `MCODE'*10007 + `PCODE'*101 + 7919, 699999999)
        * Own bootstrap namespace (314159). Method, kappa and c are absent:
        * FD/FOD and all kappa/c cells of a panel type share the draws.
        local boot_seed = 1400000001 + mod(`_h0'*1013 + `MCODE'*10007 + ///
            `PCODE'*101 + `B'*131 + `grid'*137 + `gridci'*139 + 314159, 699999999)
        timer clear 9
        timer on 9
        *----- generate panel (verbatim from study_b_worker.do) -----
        clear
        quietly set seed `dgp_seed'
        quietly set obs `=`N'*`TMAX''
        quietly gen long id = ceil(_n/`TMAX')
        bysort id: gen int t = _n
        quietly xtset id t
        quietly gen double e  = rnormal()*`SIGE'
        quietly gen double es = e/`SIGE'
        quietly gen double w  = rnormal()
        quietly gen double u = .
        quietly bysort id (t): replace u = .5*es[_n-1] + sqrt(.75)*w if _n > 1
        quietly bysort id (t): replace u = w if _n == 1
        quietly gen double q = .
        quietly bysort id (t): replace q = rnormal() if _n == 1
        quietly bysort id (t): replace q = `RHOQ'*q[_n-1] + u if _n > 1
        tempvar g
        quietly gen double `g' = (`T_consd' + `T_qd'*q)*(q > `GAMMA')
        quietly gen double y = .
        quietly bysort id (t): replace y = rnormal() if _n == 1
        quietly bysort id (t): replace y = `RHOY'*y[_n-1] + `BQ'*q + `g' + e if _n > 1
        quietly gen byte analysis = (t > `TBURN')
        quietly drop if t <= `TBURN'
        quietly xtset id t
        *----- missingness (verbatim from study_b_worker.do) -----
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
        *----- estimate: Study B core options, -notest-, citest(gamma0 + c) -----
        capture quietly xtdpthresh y if analysis, qx(q) predetermined(q) ///
            method(`method') maxlag(`maxlag_lo' `maxlag_hi') grid(`grid') ///
            gridci(`gridci') gridtype(`gridtype') gridsample(`gridsample') ///
            refine(`refine') trim(`trim') boot(`B') boottype(wild) ///
            rseed(`boot_seed') notest citest(`citest_value') coefboot(none) ///
            history(panel) level(95) vce(robust) nowarn
        local rc=_rc
        foreach nm in version_ok twostep gamma_hat N_used citest_gamma ///
            citest_accept citest_p citest_D citest_crit citest_status ///
            citest_draws seed_citest {
            local `nm'=.
        }
        if `rc'==0 {
            local version_ok=("`e(cmd)'"=="xtdpthresh" & "`e(cmdversion)'"=="0.9.36")
            if !`version_ok' {
                di as err "supp2_worker: expected xtdpthresh 0.9.36; got `e(cmd)' `e(cmdversion)'"
                exit 459
            }
            local twostep=e(estimator_twostep)
            local gamma_hat=e(gamma)
            local N_used=e(N)
            foreach nm in citest_gamma citest_accept citest_p citest_D ///
                citest_crit citest_status citest_draws seed_citest {
                capture local `nm'=e(`nm')
            }
            assert missing(e(pval_lin)) & missing(e(pval_cont))
            assert abs(`citest_gamma'-`citest_value')<=1e-12 & ///
                inrange(`citest_status',1,6) & ///
                `seed_citest'==mod(`boot_seed'+477377,2147483648)
        }
        timer off 9
        quietly timer list 9
        local elapsed_s=r(t9)
        local row "${supp_run}"
        foreach nm of local CELLVARS {
            local row "`row',``nm''"
        }
        local row "`row',${supp_master},`SHARD',`rep',`dgp_seed',`missing_seed',`boot_seed',`units_realized',`analysis_observed',`rc'"
        foreach nm in version_ok twostep gamma_hat N_used citest_gamma ///
            citest_accept citest_p citest_D citest_crit citest_status ///
            citest_draws seed_citest elapsed_s {
            local row "`row',``nm''"
        }
        tempname FH
        file open `FH' using "`OUT'", write append text
        file write `FH' "`row'" _n
        file close `FH'
        di "SUPP2_ROW cell=`cell_id' rep=`rep' rc=`rc'"
    }
}
tempname DONE
file open `DONE' using "supp2_done_SH`SHARD'.txt", write replace text
file write `DONE' "${supp_run},`SHARD',`expected'" _n
file close `DONE'
di "SUPP2_COMPLETE shard=`SHARD' rows=`expected'"
