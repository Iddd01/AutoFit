*! Final Monte Carlo merge (xtdpthresh 0.9.37); launch through run_final.ps1 -Action Merge.
version 17.0
clear all
set more off
set varabbrev off
do final_config.do
tempfile registry pooled
local CELLVARS block pair_id dgp spec method vce iv N T miss missp kappa c mode R B grid gridci gridtype gridsample trim refine
local STRS block pair_id dgp spec method vce iv miss mode gridtype gridsample
import delimited using final_cells.csv, clear case(preserve) asdouble ///
    stringcols(1 3 4 5 6 7 8 11 15 20 21)
isid cell_id
isid pair_id method
if ${fin_repcap}>0 replace R=min(R,${fin_repcap})
if ${fin_B}>0 replace B=${fin_B} if B>0
if ${fin_grid}>0 replace grid=${fin_grid}
if ${fin_gridci}>0 replace gridci=${fin_gridci} if gridci>0
egen double totalrows=total(R)
assert totalrows==${fin_expected}
drop totalrows
sort cell_id
gen double offset=sum(R)-R
foreach nm of local CELLVARS {
    rename `nm' reg_`nm'
}
save `registry'
local first=1
forvalues sh=1/${fin_nshard} {
    import delimited using "final_SH`sh'.csv", clear case(preserve) varnames(1) stringcols(_all)
    local expected run_id `CELLVARS' master shard rep dgp_seed missing_seed boot_seed units_realized analysis_observed rc version_ok twostep gamma_hat b_rho b_q b_cons b_qd b_rhod se_rho se_q se_cons se_qd se_rhod sc_rho sc_q sc_cons sc_qd sc_rhod joint_vce vce_applied ar_joint hansen_p ar1_p ar2_p N_obs N_iv N_units ci_delivered ci_lo ci_hi ci_nseg ci_boundary ci_incomplete citest_gamma citest_accept citest_p citest_D citest_crit citest_status citest_draws seed_citest p_lin lin_valid elapsed_s
    local expected : subinstr local expected "block pair_id" "block cell_id pair_id"
    unab actual: _all
    assert "`actual'"=="`expected'"
    ds
    local allvars `r(varlist)'
    local numeric: list allvars - STRS
    local numeric: list numeric - run_id
    destring `numeric', replace
    assert run_id=="${fin_run}" & shard==`sh'
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
assert _N==${fin_expected}
merge m:1 cell_id using `registry', assert(match) nogen
foreach nm of local CELLVARS {
    capture confirm string variable `nm'
    if !_rc assert `nm'==reg_`nm'
    else assert abs(`nm'-reg_`nm')<=1e-12
    drop reg_`nm'
}
bysort cell_id (rep): assert _N==R & rep==_n
assert shard==mod(offset+rep-1,${fin_nshard})+1
assert master==${fin_master}
* seed contract
gen byte spc=cond(inlist(dgp,"base","linear"),1,cond(dgp=="endog",2,cond(dgp=="persist",3,4)))
gen byte qc=cond(dgp=="endog",2,1)
gen byte mc=cond(miss=="balanced",1,cond(miss=="mcar",2,3))
gen double pc=round(missp*1000000)
gen byte ivc=cond(iv=="all",1,cond(iv=="L3",2,3))
gen double hh=mod(master*1000003+rep*10007+N*101+T*1009+spc*1000033+qc*1000037,699999999)
assert dgp_seed==1+hh
assert missing_seed==700000001+mod(hh*1009+mc*10007+pc*101+7919,699999999)
assert boot_seed==1400000001+mod(hh*1013+mc*10007+pc*101+B*131+grid*137+gridci*139+ivc*151+424242,699999999)
drop spc qc mc pc ivc hh offset
* result contracts
assert rc==0 | (missing(version_ok) & missing(twostep) & missing(gamma_hat))
assert version_ok==1 & inlist(twostep,0,1) & !missing(gamma_hat) if rc==0
gen double cit=cond(inlist(mode,"CI","FULL") & dgp!="linear",.25+c,.)
assert missing(citest_status) if missing(cit) | rc!=0
assert abs(citest_gamma-cit)<=1e-12 & inrange(citest_status,1,6) & ///
    seed_citest==mod(boot_seed+477377,2147483648) if !missing(cit) & rc==0
assert citest_draws==B & inlist(citest_accept,0,1) & citest_accept==(citest_p>.05) ///
    if rc==0 & citest_status==1
assert citest_D==0 & citest_accept==1 & citest_p==1 if rc==0 & citest_status==2
assert missing(p_lin) if !inlist(mode,"FULL","LIN")
assert lin_valid==B if !missing(p_lin)
assert inlist(ci_delivered,0,1) & !missing(ci_incomplete) if rc==0 & mode!="POINT"
assert missing(ci_delivered) if mode=="POINT"
assert ci_lo<=ci_hi if ci_delivered==1
* shared outer samples: FD/FOD and every cell of a DGP, N, T and missingness pattern
bysort dgp N T miss missp rep (cell_id): assert dgp_seed==dgp_seed[1] & ///
    missing_seed==missing_seed[1] & units_realized==units_realized[1] & ///
    analysis_observed==analysis_observed[1]
* one FD and one FOD row per pair and replication
bysort pair_id rep: assert _N==2 & method[1]!=method[2]
gen byte formal=${fin_formal}
* truths
gen double t_rho=cond(dgp=="persist",.9,.6)
gen double t_q=1
gen double t_cons=cond(dgp=="linear" | spec=="kink",cond(spec=="kink",.,0),kappa-.5)
gen double t_qd=cond(dgp=="linear",0,2)
gen double t_rhod=cond(spec=="kink",.,0)
gen double gamma0=cond(dgp=="linear",.,.25)
sort cell_id rep
save final_all.dta, replace
export delimited using final_all.csv, replace

* ---- cell summary: delivery, threshold estimate, diagnostics ----
preserve
    gen byte ok=(rc==0)
    gen double eg=gamma_hat-gamma0 if ok
    gen double sg=eg^2
    gen double ag=abs(eg)
    gen byte gclose=(ag<=.1) if !missing(ag)
    gen byte hj=(hansen_p<.05) if !missing(hansen_p)
    gen byte a1=(ar1_p<.05) if !missing(ar1_p)
    gen byte a2=(ar2_p<.05) if !missing(ar2_p)
    gen byte ciok=ci_delivered if mode!="POINT" & ok
    gen byte cit_eval=inlist(citest_status,1,2) if !missing(cit) & ok
    gen byte cit_rej=(citest_accept==0) if cit_eval==1
    gen byte hcov=(ci_lo<=gamma0 & gamma0<=ci_hi) if ci_delivered==1 & !missing(gamma0)
    gen double hlen=ci_hi-ci_lo if ci_delivered==1
    gen byte bnd=ci_boundary if ci_delivered==1
    gen byte lin_rej=(p_lin<.05) if !missing(p_lin)
    gen byte lin_del=!missing(p_lin) if inlist(mode,"FULL","LIN") & ok
    collapse (count) n_rep=rep n_gamma=eg n_hansen=hj n_ar2=a2 n_cit=cit_rej ///
        n_ci=hcov n_lin=lin_rej ///
        (sum) n_ok=ok ///
        (mean) twostep_rate=twostep joint_vce_rate=joint_vce ///
            vce_applied_rate=vce_applied bias_gamma=eg mse_gamma=sg ///
            p_gamma_close=gclose hansen_reject5=hj ar1_reject5=a1 ///
            ar2_reject5=a2 ci_delivery=ciok citest_delivery=cit_eval ///
            citest_reject5=cit_rej hull_coverage=hcov mean_hull_length=hlen ///
            boundary_rate=bnd lin_delivery=lin_del lin_reject5=lin_rej ///
            mean_N_iv=N_iv mean_units=units_realized mean_elapsed_s=elapsed_s ///
        (median) mad_gamma=ag median_hull_length=hlen ///
        (sd) sd_gamma=eg, ///
        by(block cell_id pair_id dgp spec method vce iv N T miss missp kappa c mode B grid gridci formal)
    gen double rmse_gamma=sqrt(mse_gamma)
    gen double bias_gamma_mcse=sd_gamma/sqrt(n_gamma)
    gen double citest_reject5_mcse=sqrt(citest_reject5*(1-citest_reject5)/n_cit)
    gen double hull_coverage_mcse=sqrt(hull_coverage*(1-hull_coverage)/n_ci)
    gen double lin_reject5_mcse=sqrt(lin_reject5*(1-lin_reject5)/n_lin)
    gen double ar2_reject5_mcse=sqrt(ar2_reject5*(1-ar2_reject5)/n_ar2)
    gen double hansen_reject5_mcse=sqrt(hansen_reject5*(1-hansen_reject5)/n_hansen)
    * citest_reject5 at c = 0 is one minus the pointwise coverage of the set
    gen double citest_coverage=1-citest_reject5 if c==0
    sort cell_id
    export delimited using final_summary.csv, replace
restore

* ---- coefficient summary (reported and conditional SE) ----
preserve
    keep if rc==0
    keep cell_id pair_id block dgp spec method vce iv N T miss missp kappa c mode rep ///
        b_* se_* sc_* t_*
    rename (t_rho t_q t_cons t_qd t_rhod) (truth_rho truth_q truth_cons truth_qd truth_rhod)
    reshape long b_ se_ sc_ truth_, i(cell_id rep) j(param) string
    drop if missing(truth_)
    gen double err=b_-truth_
    gen double sq=err^2
    gen byte cov=abs(err)<=invnormal(.975)*se_ if !missing(se_,err) & se_>=0
    gen byte covc=abs(err)<=invnormal(.975)*sc_ if !missing(sc_,err) & sc_>=0
    collapse (count) n_est=err n_se=cov n_sc=covc ///
        (mean) truth=truth_ bias=err mse=sq mean_se=se_ mean_sc=sc_ ///
            coverage=cov coverage_cond=covc ///
        (sd) sd=err, by(block cell_id pair_id dgp spec method vce iv N T miss missp kappa c mode param)
    gen double rmse=sqrt(mse)
    gen double bias_mcse=sd/sqrt(n_est)
    gen double se_sd_ratio=mean_se/sd
    gen double coverage_mcse=sqrt(coverage*(1-coverage)/n_se)
    gen double coverage_cond_mcse=sqrt(coverage_cond*(1-coverage_cond)/n_sc)
    sort cell_id param
    export delimited using final_coefficients.csv, replace
restore

* ---- FD vs FOD on the same samples (and bootstrap draws) ----
preserve
    keep if rc==0
    gen double sq_gamma=(gamma_hat-gamma0)^2
    gen double sq_rho=(b_rho-t_rho)^2
    gen byte cit_rej=(citest_accept==0) if inlist(citest_status,1,2)
    gen byte lin_rej=(p_lin<.05) if !missing(p_lin)
    keep block pair_id dgp spec vce iv N T miss missp kappa c mode rep method ///
        sq_gamma sq_rho cit_rej lin_rej
    reshape wide sq_gamma sq_rho cit_rej lin_rej, i(pair_id rep) j(method) string
    gen double d_sq_gamma=sq_gammafod-sq_gammafd
    gen double d_sq_rho=sq_rhofod-sq_rhofd
    foreach v in cit lin {
        gen byte both_`v'=!missing(`v'_rejfd,`v'_rejfod)
        gen byte only_fd_`v'=(`v'_rejfd==1 & `v'_rejfod==0) if both_`v'
        gen byte only_fod_`v'=(`v'_rejfd==0 & `v'_rejfod==1) if both_`v'
    }
    collapse (count) n_pairs_gamma=d_sq_gamma n_pairs_rho=d_sq_rho ///
        (mean) mse_gamma_fd=sq_gammafd mse_gamma_fod=sq_gammafod ///
            mse_rho_fd=sq_rhofd mse_rho_fod=sq_rhofod ///
            d_mse_gamma=d_sq_gamma d_mse_rho=d_sq_rho ///
        (sd) sd_d_gamma=d_sq_gamma sd_d_rho=d_sq_rho ///
        (sum) n_both_cit=both_cit n_only_fd_cit=only_fd_cit n_only_fod_cit=only_fod_cit ///
            n_both_lin=both_lin n_only_fd_lin=only_fd_lin n_only_fod_lin=only_fod_lin, ///
        by(block pair_id dgp spec vce iv N T miss missp kappa c mode)
    gen double d_mse_gamma_mcse=sd_d_gamma/sqrt(n_pairs_gamma)
    gen double d_mse_rho_mcse=sd_d_rho/sqrt(n_pairs_rho)
    foreach v in cit lin {
        gen double n_disc_`v'=n_only_fd_`v'+n_only_fod_`v'
        gen double mcnemar_p_`v'=min(1,2*binomial(n_disc_`v',min(n_only_fd_`v',n_only_fod_`v'),.5)) if n_disc_`v'>0
    }
    sort block pair_id
    export delimited using final_paired.csv, replace
restore
tempname DONE
file open `DONE' using final_merge.ok, write replace text
file write `DONE' "${fin_run},${fin_expected}" _n
file close `DONE'
di "FINAL_MERGE_PASS rows=${fin_expected}"
