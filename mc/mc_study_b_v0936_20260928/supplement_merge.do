*! Supplement merge; launch through run_supplement.ps1 -Action Merge.
version 17.0
clear all
set more off
set varabbrev off
do supplement_config.do
tempfile registry pooled
import delimited using supplement_cells.csv, clear case(preserve) asdouble
isid cell_id
if ${supp_repcap}>0 replace R=min(R,${supp_repcap})
if ${supp_B}>0 replace B=${supp_B} if B>0
if ${supp_grid}>0 replace grid=${supp_grid}
if ${supp_gridci}>0 replace gridci=${supp_gridci} if B>0
egen double totalrows=total(R)
assert totalrows==${supp_expected}
drop totalrows
sort cell_id
gen double offset=sum(R)-R
foreach nm in study pair_id method spec vce N T missp R B grid gridci {
    rename `nm' reg_`nm'
}
save `registry'
local first=1
local nshards=${supp_nshard}
forvalues sh=1/`nshards' {
    import delimited using "supplement_SH`sh'.csv", clear case(preserve) stringcols(_all)
    local expected run_id study cell_id pair_id method spec vce N T missp R B grid gridci master shard rep dgp_seed missing_seed boot_seed rc success twostep joint_vce vce_applied ar_joint gamma bwscale gamma_bw b_rho b_q b_delta b_cons b_lagd se_rho se_q se_delta se_cons se_lagd sc_rho sc_q sc_delta sc_cons sc_lagd ar2_p ar2_p_cond N_iv N_iv_dep_near ci_delivered ci_covered ci_length ci_nseg ci_boundary ci_incomplete ci_valid ci_criterion gammahat_in_set citest_gamma citest_accept citest_p citest_D citest_crit citest_status citest_draws seed_citest elapsed_s
    unab actual: _all
    assert "`actual'"=="`expected'"
    local strings run_id study pair_id method spec vce
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
foreach nm in study pair_id method spec vce N T missp R B grid gridci {
    assert `nm'==reg_`nm'
    drop reg_`nm'
}
bysort cell_id (rep): assert _N==R & rep==_n
assert shard==mod(offset+rep-1,${supp_nshard})+1
assert master==${supp_master}
gen double hh=mod(master*1000003+rep*10007+N*101+T*1009+1000033+1000037,699999999)
gen byte mc=cond(missp==0,1,2)
assert dgp_seed==1+hh
assert missing_seed==700000001+mod(hh*1009+mc*10007+7919,699999999)
assert boot_seed==1400000001+mod(hh*1013+mc*10007+B*131+grid*137+gridci*139+271828,699999999)
assert seed_citest==mod(boot_seed+477377,2147483648) if rc==0
drop hh mc offset
assert inlist(success,0,1) & inlist(ci_delivered,0,1)
assert success==0 if rc!=0
assert inlist(joint_vce,0,1) & inlist(ar_joint,0,1) & inlist(vce_applied,0,1) if rc==0
assert abs(bwscale-1.5)<1e-12 if rc==0
assert success==1 if rc==0 & twostep==1 & !missing(gamma,b_rho,b_q,b_delta) & spec=="kink"
assert success==1 if rc==0 & twostep==1 & !missing(gamma,b_rho,b_q,b_delta,b_cons,b_lagd) & spec=="jump"
assert ci_valid==B & gammahat_in_set==1 & ci_incomplete==0 if ci_delivered
assert ci_criterion==cond(twostep==1,2,1) if rc==0 & B>0
assert ci_length>=0 & !missing(ci_length,ci_covered) if ci_delivered
assert ci_delivered==0 if B==0
assert abs(citest_gamma-.25)<=1e-12 & inrange(citest_status,1,6) if rc==0
assert citest_draws==B & !missing(citest_crit,citest_accept,citest_p,citest_D) & ///
    inlist(citest_accept,0,1) & inrange(citest_p,0,1) & ///
    citest_accept==(citest_D<=citest_crit) & citest_accept==(citest_p>.05) ///
    if rc==0 & citest_status==1
assert citest_D==0 & citest_accept==1 & citest_p==1 & ///
    missing(citest_crit,citest_draws) if rc==0 & citest_status==2
bysort pair_id rep: assert dgp_seed==dgp_seed[1] & missing_seed==missing_seed[1] & boot_seed==boot_seed[1]
gen byte formal=${supp_formal}
sort cell_id rep
save supplement_all.dta, replace
export delimited using supplement_all.csv, replace
preserve
    gen double eg=gamma-.25 if success
    gen double sg=eg^2
    gen byte delivered=rc==0
    gen byte joint_ok=joint_vce if delivered
    gen byte wind_ok=vce_applied if delivered & vce=="windmeijer"
    gen byte ar_ok=ar_joint if delivered
    gen byte cov=ci_covered if ci_delivered
    gen byte ceff=ci_delivered & ci_covered==1 if B>0
    gen double clen=ci_length if ci_delivered
    gen byte cdisc=ci_nseg>1 if ci_delivered
    gen byte citest_evaluable=inlist(citest_status,1,2) if rc==0
    gen byte citest_cov=citest_accept if citest_evaluable
    gen byte citest_eff=(citest_evaluable==1 & citest_accept==1)
    collapse (count) n_rep=rep n_gamma=eg n_ci=cov n_citest=citest_cov ///
        n_joint=joint_ok n_wind=wind_ok ///
        (sum) n_success=success n_delivered=delivered ///
        (mean) bias_gamma=eg mse_gamma=sg joint_rate=joint_ok wind_applied_rate=wind_ok ///
        ar_joint_rate=ar_ok ci_delivery_rate=ci_delivered coverage=cov effective_coverage=ceff ///
        citest_delivery=citest_evaluable citest_coverage=citest_cov ///
        citest_effective_coverage=citest_eff mean_citest_p=citest_p ///
        mean_length=clen disconnected_rate=cdisc boundary_rate=ci_boundary ///
        incomplete_rate=ci_incomplete elapsed_s ///
        (sd) sd_gamma=eg sd_sq_gamma=sg sd_length=clen, ///
        by(run_id study cell_id pair_id method spec vce N T missp R B grid gridci master formal)
    gen double rmse_gamma=sqrt(mse_gamma)
    gen double bias_mcse=sd_gamma/sqrt(n_gamma)
    gen double rmse_mcse=sd_sq_gamma/(2*rmse_gamma*sqrt(n_gamma))
    gen double coverage_mcse=sqrt(coverage*(1-coverage)/n_ci)
    gen double citest_coverage_mcse=sqrt(citest_coverage*(1-citest_coverage)/n_citest)
    gen double citest_effective_coverage_mcse=sqrt(citest_effective_coverage* ///
        (1-citest_effective_coverage)/n_rep)
    gen double effective_coverage_mcse=sqrt(effective_coverage*(1-effective_coverage)/n_rep)
    gen double length_mcse=sd_length/sqrt(n_ci)
    export delimited using supplement_summary.csv, replace
restore
keep run_id study cell_id pair_id method spec vce N T missp R B grid gridci master formal rep ///
    success joint_vce vce_applied b_* se_* sc_*
reshape long b_ se_ sc_, i(cell_id rep) j(parameter) string
drop if spec=="kink" & inlist(parameter,"cons","lagd")
gen double truth=cond(parameter=="rho",.6,cond(parameter=="q",1, ///
    cond(parameter=="delta",2,cond(parameter=="cons",.5,0))))
gen double err=b_-truth if success
gen double sq=err^2
gen byte covered=abs(err)<=invnormal(.975)*se_ if success & joint_vce==1 & !missing(se_) & se_>=0
gen byte covered_cond=abs(err)<=invnormal(.975)*sc_ if success & !missing(sc_) & sc_>=0
gen double err_se=err if !missing(covered)
gen double err_sc=err if !missing(covered_cond)
gen double joint_se=se_ if !missing(covered)
gen double cond_se=sc_ if !missing(covered_cond)
preserve
    keep pair_id rep parameter method vce err sq joint_se covered
    gen str30 arm=method+"_"+vce
    drop method vce
    reshape wide err sq joint_se covered, i(pair_id rep parameter) j(arm) string
    foreach kind in robust windmeijer {
        capture confirm variable errfod_`kind'
        if !_rc {
            gen double d_error_`kind'=errfod_`kind'-errfd_`kind'
            gen double d_mse_`kind'=sqfod_`kind'-sqfd_`kind'
            gen double d_cov_`kind'=coveredfod_`kind'-coveredfd_`kind'
        }
    }
    foreach method in fd fod {
        capture confirm variable covered`method'_windmeijer
        if !_rc gen double d_cov_wind_`method'=covered`method'_windmeijer-covered`method'_robust
    }
    export delimited using supplement_paired.csv, replace
restore
collapse (count) n_rep=rep n_est=err n_se=covered n_cond=covered_cond ///
    (mean) bias=err mse=sq mean_se=joint_se mean_cond_se=cond_se coverage=covered coverage_cond=covered_cond ///
    (sd) empirical_sd=err matched_sd=err_se matched_cond_sd=err_sc sd_sq=sq, ///
    by(run_id study cell_id pair_id method spec vce N T missp R B grid gridci master formal parameter truth)
gen double rmse=sqrt(mse)
gen double bias_mcse=empirical_sd/sqrt(n_est)
gen double rmse_mcse=sd_sq/(2*rmse*sqrt(n_est))
gen double se_sd_ratio=mean_se/matched_sd
gen double cond_se_sd_ratio=mean_cond_se/matched_cond_sd
gen double coverage_mcse=sqrt(coverage*(1-coverage)/n_se)
gen double coverage_cond_mcse=sqrt(coverage_cond*(1-coverage_cond)/n_cond)
gen double effective_coverage=cond(n_se>0,coverage*n_se,0)/n_rep
export delimited using supplement_coefficients.csv, replace
tempname DONE
file open `DONE' using supplement_merge.ok, write replace text
file write `DONE' "${supp_run},${supp_expected}" _n
file close `DONE'
di "SUPPLEMENT_MERGE_PASS rows=${supp_expected}"
