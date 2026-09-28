*! Targeted 0.9.36 supplement; independent of the frozen core registry.
version 17.0
clear all
set more off
set varabbrev off
set processors 1
args SHARD
do supplement_config.do
adopath ++ "."
quietly do xtdpthresh.ado
set rng mt64
local OUT "supplement_SH`SHARD'.csv"
local HEADER "run_id,study,cell_id,pair_id,method,spec,vce,N,T,missp,R,B,grid,gridci,master,shard,rep,dgp_seed,missing_seed,boot_seed,rc,success,twostep,joint_vce,vce_applied,ar_joint,gamma,bwscale,gamma_bw,b_rho,b_q,b_delta,b_cons,b_lagd,se_rho,se_q,se_delta,se_cons,se_lagd,sc_rho,sc_q,sc_delta,sc_cons,sc_lagd,ar2_p,ar2_p_cond,N_iv,N_iv_dep_near,ci_delivered,ci_covered,ci_length,ci_nseg,ci_boundary,ci_incomplete,ci_valid,ci_criterion,gammahat_in_set,citest_gamma,citest_accept,citest_p,citest_D,citest_crit,citest_status,citest_draws,seed_citest,ci_gamma_lo,ci_gamma_hi,ci_hull_covered,ci_hull_length,elapsed_s"
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
import delimited using supplement_cells.csv, clear case(preserve) asdouble
local ncells = _N
forvalues j=1/`ncells' {
    foreach nm in study cell_id pair_id method spec vce N T missp R B grid gridci {
        local c`j'_`nm' = `nm'[`j']
    }
}
local expected = 0
local ordinal = 0
forvalues j=1/`ncells' {
    foreach nm in study cell_id pair_id method spec vce N T missp R B grid gridci {
        local `nm' "`c`j'_`nm''"
    }
    if ${supp_repcap}>0 local R=min(`R',${supp_repcap})
    if ${supp_B}>0 & `B'>0 local B=${supp_B}
    if ${supp_grid}>0 local grid=${supp_grid}
    if ${supp_gridci}>0 & `B'>0 local gridci=${supp_gridci}
    forvalues rep=1/`R' {
        local ordinal=`ordinal'+1
        if mod(`ordinal'-1,${supp_nshard})+1 != `SHARD' continue
        local expected=`expected'+1
        if strpos(" `done' "," `cell_id':`rep' ") continue
        local mcode=cond(`missp'==0,1,2)
        local h=mod(${supp_master}*1000003+`rep'*10007+`N'*101+`T'*1009+1000033+1000037,699999999)
        local dgp_seed=1+`h'
        local missing_seed=700000001+mod(`h'*1009+`mcode'*10007+7919,699999999)
        * Separate supplement bootstrap namespace; method and vce are absent.
        local boot_seed=1400000001+mod(`h'*1013+`mcode'*10007+`B'*131+`grid'*137+`gridci'*139+271828,699999999)
        local kappa=cond("`spec'"=="jump",1,0)
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
        local modelopt "refine(4)"
        if "`spec'"=="kink" local modelopt "kink"
        local infopt "noboot"
        if `B'>0 local infopt "boot(`B') gridci(`gridci') boottype(wild) rseed(`boot_seed') notest"
        * citest() requires the grid bootstrap: never with noboot
        local citestopt ""
        if `B'>0 local citestopt "citest(.25)"
        capture quietly xtdpthresh y, qx(q) predetermined(q) method(`method') ///
            maxlag(1 3) grid(`grid') trim(.15) history(panel) ///
            `modelopt' `infopt' `citestopt' coefboot(none) ///
            vce(`vce') bwscale(1.5) nowarn
        local rc=_rc
        local success=0
        local ci_delivered=0
        foreach nm in twostep joint_vce vce_applied ar_joint gamma bwscale gamma_bw ///
            b_rho b_q b_delta b_cons b_lagd se_rho se_q se_delta se_cons se_lagd ///
            sc_rho sc_q sc_delta sc_cons sc_lagd ar2_p ar2_p_cond N_iv N_iv_dep_near ///
            ci_covered ci_length ci_nseg ci_boundary ci_incomplete ci_valid ///
            ci_criterion gammahat_in_set citest_gamma citest_accept citest_p ///
            citest_D citest_crit citest_status citest_draws seed_citest ///
            ci_gamma_lo ci_gamma_hi ci_hull_covered ci_hull_length {
            local `nm'=.
        }
        if `rc'==0 {
            assert "`e(cmdversion)'"=="0.9.36" & "`e(method)'"=="`method'"
            assert "`e(vce_requested)'"=="`vce'" & e(flag_kink)==("`spec'"=="kink")
            foreach nm in joint_vce vce_applied ar_joint gamma bwscale gamma_bw ar2_p N_iv N_iv_dep_near {
                capture local `nm'=e(`nm')
            }
            local twostep=e(estimator_twostep)
            local success=(`twostep'==1)
            local ar2_p_cond=2*normal(-abs(e(ar2_cond)))
            foreach nm in citest_gamma citest_accept citest_p citest_D ///
                citest_crit citest_status citest_draws seed_citest {
                capture local `nm'=e(`nm')
            }
            if `B'==0 {
                assert missing(`citest_gamma',`citest_status',`seed_citest')
            }
            else {
            assert abs(`citest_gamma'-.25)<=1e-12 & ///
                inrange(`citest_status',1,6) & ///
                `seed_citest'==mod(`boot_seed'+477377,2147483648)
            }
            if `B'>0 & `citest_status'==1 {
                assert `citest_draws'==`B' & !missing(`citest_crit') & ///
                    inlist(`citest_accept',0,1) & inrange(`citest_p',0,1) & ///
                    `citest_accept'==(`citest_D'<=`citest_crit') & ///
                    `citest_accept'==(`citest_p'>.05)
            }
            else if `B'>0 & `citest_status'==2 {
                assert `citest_D'==0 & `citest_accept'==1 & `citest_p'==1 & ///
                    missing(`citest_crit',`citest_draws')
            }
            capture matrix VC=e(V_cond)
            local hasvc=(_rc==0)
            local names "rho q delta cons lagd"
            local labels "Lag_y_b q_b q_d cons_d Lag_y_d"
            if "`spec'"=="kink" {
                local names "rho q delta"
                local labels "Lag_y_b q_b kink_slope"
            }
            local ii=0
            foreach nm of local names {
                local ++ii
                local label: word `ii' of `labels'
                capture local b_`nm'=_b[`label']
                if missing(`b_`nm'') local success=0
                capture local se_`nm'=_se[`label']
                if `hasvc' {
                    local col=colnumb(VC,"`label'")
                    if !missing(`col') & `col'>0 {
                        if VC[`col',`col']>=0 & !missing(VC[`col',`col']) ///
                            local sc_`nm'=sqrt(VC[`col',`col'])
                    }
                }
            }
            if `B'>0 {
                local ci_incomplete=e(ci_incomplete)
                local ci_nseg=e(ci_nseg)
                local ci_boundary=e(boundary_warn)
                local ci_valid=e(gridboot_min_draws)
                local ci_criterion=cond("`e(ci_criterion)'"=="twostep",2,1)
                capture matrix CI=e(ci_segments)
                if !_rc & `ci_incomplete'==0 {
                    local ci_covered=0
                    local gammahat_in_set=0
                    local ci_length=0
                    forvalues rr=1/`=rowsof(CI)' {
                        if !missing(CI[`rr',1],CI[`rr',2]) {
                            local ci_delivered=1
                            local ci_length=`ci_length'+CI[`rr',2]-CI[`rr',1]
                            if CI[`rr',1]<=.25 & .25<=CI[`rr',2] local ci_covered=1
                            if CI[`rr',1]<=e(gamma)+1e-12*max(1,abs(e(gamma))) & ///
                                e(gamma)-1e-12*max(1,abs(e(gamma)))<=CI[`rr',2] local gammahat_in_set=1
                        }
                    }
                }
                if `ci_delivered' {
                    assert `ci_valid'==`B' & `gammahat_in_set'==1
                    assert `ci_criterion'==cond(`twostep'==1,2,1)
                    * reported interval: hull of the accepted grid points
                    local ci_gamma_lo=e(gamma_lo)
                    local ci_gamma_hi=e(gamma_hi)
                    local nrci=rowsof(CI)
                    assert !missing(`ci_gamma_lo',`ci_gamma_hi') & ///
                        `ci_gamma_lo'<=`ci_gamma_hi' & ///
                        reldif(`ci_gamma_lo',CI[1,1])<=1e-12 & ///
                        reldif(`ci_gamma_hi',CI[`nrci',2])<=1e-12
                    local ci_hull_length=`ci_gamma_hi'-`ci_gamma_lo'
                    local ci_hull_covered=(`ci_gamma_lo'<=.25 & .25<=`ci_gamma_hi')
                    assert `ci_covered'!=1 | `ci_hull_covered'==1
                }
            }
        }
        timer off 9
        quietly timer list 9
        local elapsed_s=r(t9)
        local row "${supp_run},`study',`cell_id',`pair_id',`method',`spec',`vce',`N',`T',`missp',`R',`B',`grid',`gridci',${supp_master},`SHARD',`rep',`dgp_seed',`missing_seed',`boot_seed',`rc',`success'"
        foreach nm in twostep joint_vce vce_applied ar_joint gamma bwscale gamma_bw ///
            b_rho b_q b_delta b_cons b_lagd se_rho se_q se_delta se_cons se_lagd ///
            sc_rho sc_q sc_delta sc_cons sc_lagd ar2_p ar2_p_cond N_iv N_iv_dep_near ///
            ci_delivered ci_covered ci_length ci_nseg ci_boundary ci_incomplete ///
            ci_valid ci_criterion gammahat_in_set citest_gamma citest_accept ///
            citest_p citest_D citest_crit citest_status citest_draws ///
            seed_citest ci_gamma_lo ci_gamma_hi ci_hull_covered ci_hull_length ///
            elapsed_s {
            local row "`row',``nm''"
        }
        tempname FH
        file open `FH' using "`OUT'", write append text
        file write `FH' "`row'" _n
        file close `FH'
        di "SUPP_ROW cell=`cell_id' rep=`rep' rc=`rc' success=`success'"
    }
}
tempname DONE
file open `DONE' using "supplement_done_SH`SHARD'.txt", write replace text
file write `DONE' "${supp_run},`SHARD',`expected'" _n
file close `DONE'
di "SUPPLEMENT_COMPLETE shard=`SHARD' rows=`expected'"
