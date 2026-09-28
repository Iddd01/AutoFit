version 15.0
clear all
set more off
set varabbrev off

capture erase "_registry_probe.ok"

tempfile registry
tempname PP
postfile `PP' str80 cell_id str72 pair_id str12 study str12 design ///
    str4 method int N int T str12 miss double missp str8 qstatus long R ///
    double kappa str24 dgp_spec str8 boottype str10 coefboot int grid ///
    int gridci long B byte refine str8 infmode byte maxlag_lo byte maxlag_hi ///
    double trim str10 gridtype str10 gridsample str8 pairmode ///
    using `registry', replace
global REGPOST "`PP'"

capture program drop _cell
program define _cell
    version 15.0
    args PAIR_ID STUDY DESIGN METHOD N T MISS MISSP QSTATUS R KAPPA ///
         DGP_SPEC BOOTTYPE COEFBOOT GRID GRIDCI B REFINE INFMODE ///
         MAXLAG_LO MAXLAG_HI TRIM GRIDTYPE GRIDSAMPLE PAIRMODE
    local CELL_ID "`PAIR_ID'_`METHOD'"
    post ${REGPOST} ("`CELL_ID'") ("`PAIR_ID'") ("`STUDY'") ///
        ("`DESIGN'") ("`METHOD'") (`N') (`T') ("`MISS'") (`MISSP') ///
        ("`QSTATUS'") (`R') (`KAPPA') ("`DGP_SPEC'") ("`BOOTTYPE'") ///
        ("`COEFBOOT'") (`GRID') (`GRIDCI') (`B') (`REFINE') ("`INFMODE'") ///
        (`MAXLAG_LO') (`MAXLAG_HI') (`TRIM') ("`GRIDTYPE'") ///
        ("`GRIDSAMPLE'") ("`PAIRMODE'")
end

do study_b_cells.do
postclose `PP'
macro drop REGPOST
use `registry', clear

assert _N == 52
isid cell_id
isid pair_id method
assert inlist(pairmode, "paired", "fdonly")
sort pair_id method
by pair_id: assert pairmode == pairmode[1]
by pair_id: assert (_N == 2 & method[1] == "fd" & method[2] == "fod") ///
    if pairmode == "paired"
by pair_id: assert (_N == 1 & method == "fd") if pairmode == "fdonly"
egen byte pair_tag = tag(pair_id)
count if pair_tag & pairmode == "paired"
assert r(N) == 24
count if pair_tag & pairmode == "fdonly"
assert r(N) == 4

assert T == 6 & inlist(N, 400, 800) & R == 500
assert boottype == "wild" & coefboot == "none"
assert inlist(infmode, "FULL", "COVERAGE")
assert (miss == "balanced" & missp == 0) | ///
    (miss == "mcar" & missp == .30) | ///
    (miss == "attrition" & missp == .15)
assert inlist(kappa, 0, .1, .5, 1)
assert (study == "LIN" & design == "none" & kappa == 0) | ///
    (study != "LIN" & kappa == 0 & design == "kink") | ///
    (study != "LIN" & kappa > 0 & design == "jump")

gen byte block = real(substr(pair_id, 2, 1))
assert inrange(block, 1, 6)

* Paired main validation design (B1--B5).
assert block > 5 | (pairmode == "paired" & B == 499 & grid == 199 & ///
    refine == 4 & maxlag_lo == 1 & maxlag_hi == 3 & abs(trim-.15)<1e-12 & ///
    gridtype == "uniform" & gridsample == "effective")
assert block != 1 | (N == 400 & inlist(study, "GSBAL", "GSGAP") & ///
    qstatus == "pred" & inlist(kappa, 0, .1, .5, 1) & ///
    dgp_spec == "gs26_base_v2" & gridci == 100 & infmode == "FULL")
assert block != 2 | (N == 800 & inlist(study, "GSBAL", "GSGAP") & ///
    qstatus == "pred" & inlist(kappa, 0, .1, 1) & ///
    dgp_spec == "gs26_base_v2" & gridci == 100 & ///
    ((kappa == 0 & infmode == "FULL") | ///
     (kappa > 0 & infmode == "COVERAGE")))
assert block != 3 | (study == "LIN" & inlist(N, 400, 800) & ///
    qstatus == "pred" & dgp_spec == "gs26_base_v2" & ///
    gridci == 10 & infmode == "FULL")
assert block != 4 | (study == "ENDOG" & N == 400 & ///
    inlist(kappa, 0, 1) & qstatus == "endog" & ///
    dgp_spec == "gs26_endog_v2" & gridci == 100 & infmode == "FULL")
assert block != 5 | (study == "ATTR" & N == 400 & miss == "attrition" & ///
    missp == .15 & qstatus == "pred" & inlist(kappa, 0, 1) & ///
    dgp_spec == "gs26_base_v2" & gridci == 100 & infmode == "FULL")

* FD-only published-design calibration bridge (not exact table replication).
assert block != 6 | (study == "GSCAL" & method == "fd" & ///
    pairmode == "fdonly" & inlist(N, 400, 800) & miss == "balanced" & ///
    missp == 0 & qstatus == "pred" & inlist(kappa, 0, 1) & ///
    dgp_spec == "gs26_base_v2" & B == 500 & grid == 46 & gridci == 46 & ///
    refine == 0 & infmode == "COVERAGE" & maxlag_lo == 1 & ///
    maxlag_hi == 5 & abs(trim-.10)<1e-12 & gridtype == "quantile" & ///
    gridsample == "observed")

forvalues j = 1/6 {
    count if block == `j'
    local NB`j' = r(N)
}
assert `NB1' == 16 & `NB2' == 12 & `NB3' == 8 & `NB4' == 8 & ///
    `NB5' == 4 & `NB6' == 4
foreach k in 0 .1 .5 1 {
    count if block == 1 & kappa == `k'
    assert r(N) == 4
}

count if R == 500
assert r(N) == 52
count if infmode == "FULL"
assert r(N) == 40
count if infmode == "COVERAGE"
assert r(N) == 12
count if B == 499
assert r(N) == 48
count if B == 500
assert r(N) == 4
count if grid == 199
assert r(N) == 48
count if grid == 46
assert r(N) == 4
count if refine == 4
assert r(N) == 48
count if refine == 0
assert r(N) == 4
count if gridci == 100
assert r(N) == 40
count if gridci == 46
assert r(N) == 4
count if gridci == 10
assert r(N) == 8

egen long total_reps = total(R)
assert total_reps[1] == 26000
egen long full_reps = total(cond(infmode=="FULL", R, 0))
assert full_reps[1] == 20000
egen long paired_reps = total(cond(pairmode=="paired", R, 0))
assert paired_reps[1] == 24000

* Preflight the canonical formal bootstrap namespace.  Methods in a registered
* FD/FOD pair must share a seed, while no two substantive pair/rep keys may
* reuse one.  Numerical tuning fields distinguish B6 from coincident main-DGP
* cells without breaking method pairing.
expand R
bysort cell_id: gen long rep = _n
gen byte dcode = cond(design=="kink",1,cond(design=="jump",2, ///
    cond(design=="weakkink",3,4)))
gen byte mcode = cond(miss=="balanced",1,cond(miss=="mcar",2,3))
gen byte qcode = cond(qstatus=="pred",1,2)
gen byte spcode = cond(inlist(dgp_spec,"gs26_base_v2","gs26_official_v1"),1, ///
    cond(dgp_spec=="gs26_endog_v2",2,3))
gen byte btcode = cond(boottype=="wild",1,2)
gen byte cbcode = cond(coefboot=="none",1,cond(coefboot=="onestep",2,3))
gen double kcode = round((kappa+5)*1000)
gen double pcode = round(missp*1000000)
gen double h0 = mod(20260814*1000003 + rep*10007 + N*101 + T*1009 + ///
    spcode*1000033 + qcode*1000037, 699999999)
gen double hb = mod(h0*1013 + dcode*100003 + mcode*10007 + ///
    pcode*101 + qcode*1009 + btcode*17 + cbcode*29 + kcode*37 + ///
    B*131 + grid*137 + gridci*139 + refine*149 + ///
    maxlag_lo*151 + maxlag_hi*157 + 104729, 699999999)
gen double boot_seed = 1400000001 + hb
sort pair_id rep method
by pair_id rep: assert boot_seed == boot_seed[1]
by pair_id rep: keep if _n == 1
isid boot_seed
assert inrange(boot_seed,1400000001,2100000000)

file open ok using "_registry_probe.ok", write replace text
file write ok "cells=52;pairs=24;fdonly=4;reps=26000;R500=52;FULL=40;COVERAGE=12;B499=48;B500=4;G199=48;G46=4;GC100=40;GC46=4;GC10=8;RF4=48;RF0=4" _n
file close ok
exit 0, clear
