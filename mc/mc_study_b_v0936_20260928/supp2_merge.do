*! Supplement 2 merge; launch through run_supp2.ps1 -Action Merge.
version 17.0
clear all
set more off
set varabbrev off
do supp2_config.do
tempfile registry pooled
local CELLVARS block pair_id version method N T missp kappa R B grid gridci refine trim maxlag_lo maxlag_hi gridtype gridsample mode citest_value
import delimited using supp2_cells.csv, clear case(preserve) asdouble stringcols(1 3 4 5 18 19 20)
isid cell_id
if ${supp_repcap}>0 replace R=min(R,${supp_repcap})
if ${supp_B}>0 replace B=${supp_B}
if ${supp_grid}>0 replace grid=${supp_grid}
if ${supp_gridci}>0 replace gridci=${supp_gridci}
egen double totalrows=total(R)
assert totalrows==${supp_expected}
drop totalrows
sort cell_id
gen double offset=sum(R)-R
foreach nm of local CELLVARS {
    rename `nm' reg_`nm'
}
save `registry'
local first=1
local nshards=${supp_nshard}
forvalues sh=1/`nshards' {
    import delimited using "supp2_SH`sh'.csv", clear case(preserve) stringcols(_all)
    local expected run_id block cell_id pair_id version method N T missp kappa R B grid gridci refine trim maxlag_lo maxlag_hi gridtype gridsample mode citest_value master shard rep dgp_seed missing_seed boot_seed rc version_ok twostep gamma p_lin lin_valid p_cont cont_valid cont_common seed_continuity citest_gamma citest_accept citest_p citest_D citest_crit citest_status citest_draws seed_citest elapsed_s
    unab actual: _all
    assert "`actual'"=="`expected'"
    local strings run_id block pair_id version method gridtype gridsample mode
    ds
    local allvars `r(varlist)'
    local numeric: list allvars - strings
    destring `numeric', replace
    assert run_id=="${supp_run}" & shard==`sh'
    if `first' {
        save `pooled'
        local first=0
    }
    else {
        append using `pooled'
        save `pooled', replace
    }
}
use `pooled', clear
isid cell_id rep
assert _N==${supp_expected}
merge m:1 cell_id using `registry', assert(match) nogen
foreach nm of local CELLVARS {
    capture confirm string variable `nm'
    if !_rc assert `nm'==reg_`nm'
    else assert `nm'==reg_`nm' | (missing(`nm') & missing(reg_`nm'))
    drop reg_`nm'
}
bysort cell_id (rep): assert _N==R & rep==_n
assert shard==mod(offset+rep-1,${supp_nshard})+1
assert master==${supp_master}
gen double hh=mod(master*1000003+rep*10007+N*101+T*1009+1000033+1000037,699999999)
gen byte mc=cond(missp==0,1,2)
assert dgp_seed==1+hh
assert missing_seed==700000001+mod(hh*1009+mc*10007+7919,699999999)
assert boot_seed==1400000001+mod(hh*1013+mc*10007+B*131+grid*137+gridci*139+314159,699999999)
drop hh mc offset
* contracts (each "all missing" condition is tested field by field)
assert inlist(rc,0) | missing(version_ok)
assert version_ok==1 & inlist(twostep,0,1) if rc==0
assert seed_continuity==mod(boot_seed+224737,2147483648) if rc==0 & mode=="FULL"
assert missing(p_lin) & missing(p_cont) if mode=="CI"
assert inrange(p_lin,0,1) & lin_valid==B if !missing(p_lin)
assert inrange(p_cont,0,1) & cont_valid==B if !missing(p_cont)
gen byte cit_req=!missing(citest_value)
assert cit_req==0 if version=="0.9.35"
assert missing(citest_gamma) & missing(citest_accept) & missing(citest_p) & ///
    missing(citest_status) & missing(seed_citest) if !cit_req
assert abs(citest_gamma-citest_value)<=1e-12 & inrange(citest_status,1,6) & ///
    seed_citest==mod(boot_seed+477377,2147483648) if cit_req & rc==0
assert citest_draws==B & !missing(citest_crit) & inlist(citest_accept,0,1) & ///
    inrange(citest_p,0,1) & citest_accept==(citest_D<=citest_crit) & ///
    citest_accept==(citest_p>.05) if cit_req & rc==0 & citest_status==1
assert citest_D==0 & citest_accept==1 & citest_p==1 & missing(citest_crit) & ///
    missing(citest_draws) if cit_req & rc==0 & citest_status==2
* paired samples: same dgp/boot seeds across versions and kappa (block C)
bysort block N rep: assert dgp_seed==dgp_seed[1] & boot_seed==boot_seed[1] if block=="C"
gen byte formal=${supp_formal}
sort cell_id rep
save supp2_all.dta, replace
export delimited using supp2_all.csv, replace

* ---- per-cell summary ----
preserve
    gen byte ok=(rc==0)
    gen byte lin_del=!missing(p_lin) if mode=="FULL"
    gen byte cont_del=!missing(p_cont) if mode=="FULL"
    gen byte lin_rej=(p_lin<.05) if !missing(p_lin)
    gen byte cont_rej=(p_cont<.05) if !missing(p_cont)
    gen byte cit_eval=inlist(citest_status,1,2) if cit_req & ok
    gen byte cit_rej=(citest_accept==0) if cit_eval==1
    gen byte cit_acc=(citest_accept==1) if cit_eval==1
    collapse (count) n_rep=rep n_lin=lin_rej n_cont=cont_rej n_cit=cit_rej ///
        (sum) n_ok=ok ///
        (mean) twostep_rate=twostep lin_delivery=lin_del cont_delivery=cont_del ///
            lin_reject5=lin_rej cont_reject5=cont_rej ///
            citest_reject5=cit_rej citest_accept_rate=cit_acc ///
            mean_p_cont=p_cont mean_elapsed_s=elapsed_s ///
        (first) N missp kappa B grid gridci citest_value, ///
        by(run_id block cell_id pair_id version method mode formal)
    gen double lin_reject5_mcse=sqrt(lin_reject5*(1-lin_reject5)/n_lin)
    gen double cont_reject5_mcse=sqrt(cont_reject5*(1-cont_reject5)/n_cont)
    gen double citest_reject5_mcse=sqrt(citest_reject5*(1-citest_reject5)/n_cit)
    gen double citest_c=citest_value-.25
    order run_id block cell_id pair_id version method mode N kappa citest_value ///
        citest_c n_rep n_ok twostep_rate cont_delivery n_cont cont_reject5 ///
        cont_reject5_mcse lin_delivery n_lin lin_reject5 lin_reject5_mcse ///
        n_cit citest_reject5 citest_reject5_mcse citest_accept_rate
    sort block version method kappa N citest_value
    export delimited using supp2_summary.csv, replace
restore

* ---- block C: 0.9.35 vs 0.9.36 on the same samples and draws ----
preserve
    keep if block=="C"
    gen byte cont_rej=(p_cont<.05) if !missing(p_cont)
    gen vtag=cond(version=="0.9.35","v35","v36")
    keep method N kappa rep vtag p_cont cont_rej
    reshape wide p_cont cont_rej, i(method N kappa rep) j(vtag) string
    gen byte both=!missing(cont_rejv35) & !missing(cont_rejv36)
    gen byte only35=(cont_rejv35==1 & cont_rejv36==0) if both
    gen byte only36=(cont_rejv35==0 & cont_rejv36==1) if both
    collapse (sum) n_both=both n_only35=only35 n_only36=only36 ///
        (mean) reject_v35=cont_rejv35 reject_v36=cont_rejv36, by(method N kappa)
    gen double diff_v36_v35=reject_v36-reject_v35
    * exact McNemar-type count of discordant pairs; p from the binomial(1/2)
    gen double n_disc=n_only35+n_only36
    gen double mcnemar_p=min(1,2*binomial(n_disc,min(n_only35,n_only36),.5)) if n_disc>0
    sort method kappa
    export delimited using supp2_paired.csv, replace
restore
tempname DONE
file open `DONE' using supp2_merge.ok, write replace text
file write `DONE' "${supp_run},${supp_expected}" _n
file close `DONE'
di "SUPP2_MERGE_PASS rows=${supp_expected}"
