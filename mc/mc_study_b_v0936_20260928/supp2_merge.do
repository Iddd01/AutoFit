*! Supplement 2 merge (threshold-test power); launch through run_supp2.ps1 -Action Merge.
version 17.0
clear all
set more off
set varabbrev off
do supp2_config.do
tempfile registry pooled
local CELLVARS block pair_id link_cell method N T miss missp kappa c citest_value R B grid gridci refine trim maxlag_lo maxlag_hi gridtype gridsample
import delimited using supp2_cells.csv, clear case(preserve) asdouble ///
    stringcols(1 3 4 5 8 21 22)
isid cell_id
isid pair_id method
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
    import delimited using "supp2_SH`sh'.csv", clear case(preserve) varnames(1) stringcols(_all)
    local expected run_id block cell_id pair_id link_cell method N T miss missp kappa c citest_value R B grid gridci refine trim maxlag_lo maxlag_hi gridtype gridsample master shard rep dgp_seed missing_seed boot_seed units_realized analysis_observed rc version_ok twostep gamma_hat N_used citest_gamma citest_accept citest_p citest_D citest_crit citest_status citest_draws seed_citest elapsed_s
    unab actual: _all
    assert "`actual'"=="`expected'"
    local strings run_id block pair_id link_cell method miss gridtype gridsample
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
    else assert abs(`nm'-reg_`nm')<=1e-12
    drop reg_`nm'
}
bysort cell_id (rep): assert _N==R & rep==_n
assert shard==mod(offset+rep-1,${supp_nshard})+1
assert master==${supp_master}
* Study B seed contract (outer sample) and the supplement's bootstrap namespace
gen double hh=mod(master*1000003+rep*10007+N*101+T*1009+1000033+1000037,699999999)
gen byte mc=cond(miss=="balanced",1,cond(miss=="mcar",2,3))
gen double pc=round(missp*1000000)
assert dgp_seed==1+hh
assert missing_seed==700000001+mod(hh*1009+mc*10007+pc*101+7919,699999999)
assert boot_seed==1400000001+mod(hh*1013+mc*10007+pc*101+B*131+grid*137+gridci*139+314159,699999999)
drop hh mc pc offset
* contracts (each "all missing" condition is tested field by field)
assert abs(citest_value-(.25+c))<=1e-12
assert rc==0 | (missing(version_ok) & missing(twostep) & missing(gamma_hat) & ///
    missing(citest_status) & missing(seed_citest))
assert version_ok==1 & inlist(twostep,0,1) & !missing(gamma_hat) if rc==0
assert abs(citest_gamma-citest_value)<=1e-12 & inrange(citest_status,1,6) & ///
    seed_citest==mod(boot_seed+477377,2147483648) if rc==0
assert citest_draws==B & !missing(citest_crit) & inlist(citest_accept,0,1) & ///
    inrange(citest_p,0,1) & citest_accept==(citest_D<=citest_crit) & ///
    citest_accept==(citest_p>.05) if rc==0 & citest_status==1
assert citest_D==0 & citest_accept==1 & citest_p==1 & missing(citest_crit) & ///
    missing(citest_draws) if rc==0 & citest_status==2
* paired samples: FD/FOD and every kappa/c cell of a panel type and N share
* the outer sample; the realized panel must therefore coincide
bysort block miss N rep (cell_id): assert dgp_seed==dgp_seed[1] & ///
    missing_seed==missing_seed[1] & boot_seed==boot_seed[1] & ///
    units_realized==units_realized[1] & analysis_observed==analysis_observed[1]
gen byte formal=${supp_formal}
sort cell_id rep
save supp2_all.dta, replace
export delimited using supp2_all.csv, replace

* ---- per-cell summary: rejection of the false threshold .25 + c ----
preserve
    gen byte ok=(rc==0)
    gen byte cit_eval=inlist(citest_status,1,2) if ok
    gen byte cit_rej=(citest_accept==0) if cit_eval==1
    gen byte cit_rej_eff=(cit_eval==1 & citest_accept==0)
    gen byte cit_mech=(citest_status==2) if ok
    collapse (count) n_rep=rep n_cit=cit_rej ///
        (sum) n_ok=ok ///
        (mean) twostep_rate=twostep citest_delivery=cit_eval ///
            citest_reject5=cit_rej citest_reject5_eff=cit_rej_eff ///
            citest_mech_rate=cit_mech mean_citest_p=citest_p ///
            mean_elapsed_s=elapsed_s ///
        (first) N missp kappa c citest_value B grid gridci link_cell, ///
        by(run_id block cell_id pair_id method miss formal)
    gen double citest_reject5_mcse=sqrt(citest_reject5*(1-citest_reject5)/n_cit)
    gen double citest_reject5_eff_mcse=sqrt(citest_reject5_eff*(1-citest_reject5_eff)/n_rep)
    order run_id block cell_id pair_id link_cell method miss N kappa c ///
        citest_value n_rep n_ok twostep_rate citest_delivery n_cit ///
        citest_reject5 citest_reject5_mcse citest_reject5_eff ///
        citest_reject5_eff_mcse citest_mech_rate mean_citest_p mean_elapsed_s
    sort block miss kappa c method
    export delimited using supp2_summary.csv, replace
restore

* ---- FD vs FOD on the same samples and draws (registries with both) ----
quietly levelsof method
local nmeth : word count `r(levels)'
preserve
if `nmeth' == 2 {
    gen byte cit_rej=(citest_accept==0) if inlist(citest_status,1,2) & rc==0
    keep block miss N kappa c rep method cit_rej
    reshape wide cit_rej, i(block miss N kappa c rep) j(method) string
    gen byte both=!missing(cit_rejfd) & !missing(cit_rejfod)
    gen byte only_fd=(cit_rejfd==1 & cit_rejfod==0) if both
    gen byte only_fod=(cit_rejfd==0 & cit_rejfod==1) if both
    * rates on the pairs where both tests are evaluable, as McNemar
    gen byte rfd=cit_rejfd if both
    gen byte rfod=cit_rejfod if both
    collapse (sum) n_both=both n_only_fd=only_fd n_only_fod=only_fod ///
        (mean) reject_fd=rfd reject_fod=rfod, by(block miss N kappa c)
    gen double diff_fod_fd=reject_fod-reject_fd
    * exact McNemar test on the discordant pairs
    gen double n_disc=n_only_fd+n_only_fod
    gen double mcnemar_p=min(1,2*binomial(n_disc,min(n_only_fd,n_only_fod),.5)) if n_disc>0
    sort block miss kappa c
}
else {
    clear
    set obs 1
    gen str40 note = "single-method registry: no FD-FOD pairs"
}
export delimited using supp2_paired.csv, replace
restore
tempname DONE
file open `DONE' using supp2_merge.ok, write replace text
file write `DONE' "${supp_run},${supp_expected}" _n
file close `DONE'
di "SUPP2_MERGE_PASS rows=${supp_expected}"
