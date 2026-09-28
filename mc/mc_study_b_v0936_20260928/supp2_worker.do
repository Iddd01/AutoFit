*! Supplement 2 (xtdpthresh 0.9.36 package): block C compares the continuity
*! test of 0.9.35 and 0.9.36 on the same samples and bootstrap seeds
*! (kappa = 0,1,2,3); block P measures the power of the threshold test at
*! gamma0 + c (Gong-Seo 2026, Table 2 layout). Independent of the core
*! registry and of the kink supplement.
version 17.0
clear all
set more off
set varabbrev off
set processors 1
args SHARD
do supp2_config.do
set rng mt64
local OUT "supp2_SH`SHARD'.csv"
local HEADER "run_id,block,cell_id,pair_id,version,method,N,T,missp,kappa,R,B,grid,gridci,refine,trim,maxlag_lo,maxlag_hi,gridtype,gridsample,mode,citest_value,master,shard,rep,dgp_seed,missing_seed,boot_seed,rc,version_ok,twostep,gamma,p_lin,lin_valid,p_cont,cont_valid,cont_common,seed_continuity,citest_gamma,citest_accept,citest_p,citest_D,citest_crit,citest_status,citest_draws,seed_citest,elapsed_s"
local done
capture confirm file "`OUT'"
if !_rc {
    tempname READ
    file open `READ' using "`OUT'", read text
    file read `READ' first
    file close `READ'
    assert `"`first'"' == `"`HEADER'"'
    import delimited using "`OUT'", clear case(preserve) stringcols(_all)
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
local CELLVARS block cell_id pair_id version method N T missp kappa R B grid gridci refine trim maxlag_lo maxlag_hi gridtype gridsample mode citest_value
import delimited using supp2_cells.csv, clear case(preserve) asdouble stringcols(1 3 4 5 18 19 20)
local ncells = _N
forvalues j=1/`ncells' {
    foreach nm of local CELLVARS {
        local c`j'_`nm' = `nm'[`j']
    }
}
* The two command versions share program and Mata names; each shard loads the
* version of the current cell, and the registry lists all 0.9.35 cells first,
* so a shard switches at most once.
local loaded ""
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
    forvalues rep=1/`R' {
        local ordinal=`ordinal'+1
        if mod(`ordinal'-1,${supp_nshard})+1 != `SHARD' continue
        local expected=`expected'+1
        if strpos(" `done' "," `cell_id':`rep' ") continue
        if "`version'" != "`loaded'" {
            capture program drop _all
            mata: mata clear
            if "`version'"=="0.9.35"      quietly do xtdpthresh_v0935.ado
            else if "`version'"=="0.9.36" quietly do xtdpthresh.ado
            else {
                di as err "supp2_worker: unknown version `version'"
                exit 198
            }
            local loaded "`version'"
        }
        local mcode=cond(`missp'==0,1,2)
        local h=mod(${supp_master}*1000003+`rep'*10007+`N'*101+`T'*1009+1000033+1000037,699999999)
        local dgp_seed=1+`h'
        local missing_seed=700000001+mod(`h'*1009+`mcode'*10007+7919,699999999)
        * Own bootstrap namespace (314159); version, method and kappa are absent,
        * so 0.9.35/0.9.36 and FD/FOD share draws and kappa cells share samples.
        local boot_seed=1400000001+mod(`h'*1013+`mcode'*10007+`B'*131+`grid'*137+`gridci'*139+314159,699999999)
        timer clear 9
        timer on 9
        clear
        quietly set seed `dgp_seed'
        quietly set obs `=`N'*(`T'+20)'
        quietly gen long id=ceil(_n/(`T'+20))
        quietly bysort id: gen int t=_n
        quietly gen double e=rnormal()*.5
        quietly gen double es=e/.5
        quietly gen double w=rnormal()
        quietly gen double u=.
        quietly bysort id (t): replace u=.5*es[_n-1]+sqrt(.75)*w if _n>1
        quietly bysort id (t): replace u=w if _n==1
        quietly gen double q=.
        quietly bysort id (t): replace q=rnormal() if _n==1
        quietly bysort id (t): replace q=.7*q[_n-1]+u if _n>1
        * jump at gamma0 = .25 equals kappa: (kappa - .5 + 2q) 1(q > .25)
        quietly gen double ge=(`kappa'-.5+2*q)*(q>.25)
        quietly gen double y=.
        quietly bysort id (t): replace y=rnormal() if _n==1
        quietly bysort id (t): replace y=.6*y[_n-1]+q+ge+e if _n>1
        quietly drop if t<=20
        if `missp'>0 {
            quietly set seed `missing_seed'
            quietly gen double du=runiform()
            quietly drop if du<`missp' & t!=21 & t!=20+`T'
            quietly bysort id: gen int nr=_N
            quietly drop if nr<4
        }
        quietly xtset id t
        local testopt = cond("`mode'"=="FULL", "", "notest")
        local citestopt ""
        if "`citest_value'" != "." local citestopt "citest(`citest_value')"
        capture quietly xtdpthresh y, qx(q) predetermined(q) method(`method') ///
            maxlag(`maxlag_lo' `maxlag_hi') grid(`grid') gridci(`gridci') ///
            gridtype(`gridtype') gridsample(`gridsample') refine(`refine') ///
            trim(`trim') boot(`B') boottype(wild) rseed(`boot_seed') ///
            `testopt' `citestopt' coefboot(none) history(panel) ///
            vce(robust) bwscale(1.5) nowarn
        local rc=_rc
        foreach nm in version_ok twostep gamma p_lin lin_valid p_cont cont_valid ///
            cont_common seed_continuity citest_gamma citest_accept citest_p ///
            citest_D citest_crit citest_status citest_draws seed_citest {
            local `nm'=.
        }
        if `rc'==0 {
            local version_ok=("`e(cmd)'"=="xtdpthresh" & "`e(cmdversion)'"=="`version'")
            if !`version_ok' {
                di as err "supp2_worker: expected xtdpthresh `version'; got `e(cmd)' `e(cmdversion)'"
                exit 459
            }
            local twostep=e(estimator_twostep)
            local gamma=e(gamma)
            local p_lin=e(pval_lin)
            local p_cont=e(pval_cont)
            foreach nm in lin_valid cont_valid cont_common seed_continuity {
                local src=cond("`nm'"=="lin_valid","boot_linearity_valid", ///
                    cond("`nm'"=="cont_valid","boot_continuity_valid", ///
                    cond("`nm'"=="cont_common","continuity_common_grid","seed_continuity")))
                capture local `nm'=e(`src')
            }
            foreach nm in citest_gamma citest_accept citest_p citest_D ///
                citest_crit citest_status citest_draws seed_citest {
                capture local `nm'=e(`nm')
            }
            if "`mode'"=="FULL" {
                assert `seed_continuity'==mod(`boot_seed'+224737,2147483648)
            }
            else {
                assert missing(`p_lin') & missing(`p_cont')
            }
            if "`citest_value'" != "." {
                assert abs(`citest_gamma'-`citest_value')<=1e-12 & ///
                    inrange(`citest_status',1,6) & ///
                    `seed_citest'==mod(`boot_seed'+477377,2147483648)
            }
            else {
                assert missing(`citest_gamma') & missing(`citest_status') & ///
                    missing(`seed_citest')
            }
        }
        timer off 9
        quietly timer list 9
        local elapsed_s=r(t9)
        local row "${supp_run}"
        foreach nm of local CELLVARS {
            local row "`row',``nm''"
        }
        local row "`row',${supp_master},`SHARD',`rep',`dgp_seed',`missing_seed',`boot_seed',`rc'"
        foreach nm in version_ok twostep gamma p_lin lin_valid p_cont cont_valid ///
            cont_common seed_continuity citest_gamma citest_accept citest_p ///
            citest_D citest_crit citest_status citest_draws seed_citest elapsed_s {
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
