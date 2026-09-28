*! Supplement 2 link check and power table. Run in the merged supplement-2
*! run folder, after the merge:
*!     do supp2_link_check.do "<path to the formal Study B study_b_all.dta>"
*! 1. Verifies that rep r of every supplement cell is rep r of its Study B
*!    link_cell: same realized panel and the same threshold estimate (gamma
*!    does not depend on the bootstrap seed).
*! 2. Writes supp2_power_table.csv: rejection rate of H0: gamma = .25 + c,
*!    with c = 0 (one minus the Study B citest coverage at the true
*!    threshold, on the same samples) and c = .10, .25, .50 from this run.
version 17.0
args BALL
if `"`BALL'"' == "" {
    di as err "usage: do supp2_link_check.do <path to study_b_all.dta>"
    exit 198
}
confirm file `"`BALL'"'
use supp2_all.dta, clear
keep block cell_id link_cell method miss N kappa c rep rc gamma_hat ///
    units_realized analysis_observed citest_accept citest_status formal
tempfile supp sb
save `supp'
use cell_id rep rc gamma_hat units_realized analysis_observed citest_accept ///
    citest_status citest_evaluable using `"`BALL'"', clear
rename cell_id link_cell
foreach v in rc gamma_hat units_realized analysis_observed citest_accept citest_status citest_evaluable {
    rename `v' b_`v'
}
isid link_cell rep
save `sb'
use `supp', clear
merge m:1 link_cell rep using `sb', keep(master match)
assert _merge==3
drop _merge
* same outer sample
assert units_realized==b_units_realized & analysis_observed==b_analysis_observed
* same point estimate when both fits succeeded (Study B's gridci differs,
* which changes only the bootstrap seed, never gamma-hat)
assert reldif(gamma_hat,b_gamma_hat)<=1e-9 if rc==0 & b_rc==0
count if rc==0 & b_rc==0
local nmatch=r(N)
count
di as res "SUPP2_LINK_PASS rows=`r(N)' gamma-matched=`nmatch'"

* ---- power table ----
gen byte rej=(citest_accept==0) if rc==0 & inlist(citest_status,1,2)
gen byte rej0=(b_citest_accept==0) if b_rc==0 & b_citest_evaluable==1
preserve
    keep if c==.1 | abs(c-.1)<1e-12
    collapse (mean) r=rej0 (count) n=rej0, by(block miss N kappa method)
    gen double c=0
    tempfile zero
    save `zero'
restore
collapse (mean) r=rej (count) n=rej, by(block miss N kappa method c)
append using `zero'
gen double mcse=sqrt(r*(1-r)/n)
gen int ck=round(100*c)
drop c
reshape wide r n mcse, i(block miss N kappa method) j(ck)
rename (r0 r10 r25 r50) (reject_c000 reject_c010 reject_c025 reject_c050)
rename (mcse0 mcse10 mcse25 mcse50) (mcse_c000 mcse_c010 mcse_c025 mcse_c050)
rename (n0 n10 n25 n50) (n_c000 n_c010 n_c025 n_c050)
order block miss N kappa method reject_c* mcse_c* n_c*
sort block miss kappa method
export delimited using supp2_power_table.csv, replace
list block miss N kappa method reject_c*, noobs sepby(block miss kappa) abbreviate(12)
