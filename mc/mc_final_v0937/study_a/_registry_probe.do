version 15.0
clear all
set more off
set varabbrev off

tempfile registry
tempname PP
postfile `PP' str72 pair_id str12 study str12 design str4 method ///
    int N int T str12 miss double missp str8 qstatus int R double kappa ///
    str24 dgp_spec int grid byte refine using `registry', replace
global REGPOST "`PP'"

capture program drop _cell
program define _cell
    version 15.0
    args PAIR_ID STUDY DESIGN METHOD N T MISS MISSP QSTATUS R KAPPA ///
         DGP_SPEC GRID REFINE
    post ${REGPOST} ("`PAIR_ID'") ("`STUDY'") ("`DESIGN'") ///
        ("`METHOD'") (`N') (`T') ("`MISS'") (`MISSP') ///
        ("`QSTATUS'") (`R') (`KAPPA') ("`DGP_SPEC'") (`GRID') (`REFINE')
end

do study_a_cells.do
postclose `PP'
macro drop REGPOST
use `registry', clear

assert _N == 116
gen long registry_order = _n
gen byte smoke_shard28 = mod(registry_order-1, 28)+1
bysort smoke_shard28: gen byte smoke_cells28 = _N
egen byte smoke_shard_tag = tag(smoke_shard28)
assert inlist(smoke_cells28, 4, 5)
count if smoke_shard_tag & smoke_cells28 == 5
assert r(N) == 4
count if smoke_shard_tag & smoke_cells28 == 4
assert r(N) == 24
isid pair_id method
sort pair_id method
by pair_id: assert _N == 2 & method[1] == "fd" & method[2] == "fod"
egen byte pair_tag = tag(pair_id)
count if pair_tag
assert r(N) == 58

assert inlist(study, "GSCORE", "GSSMALL", "PERSIST", "MISSPAT", ///
    "ENDOG", "HEAVYTAIL")
assert inlist(method, "fd", "fod")
assert inlist(N, 100, 200, 400, 800, 1600) & inlist(T, 6, 10)
assert (inlist(study, "PERSIST", "ENDOG", "HEAVYTAIL") & R == 500) | ///
    (!inlist(study, "PERSIST", "ENDOG", "HEAVYTAIL") & R == 400)
assert grid == 199 & refine == 4
assert inlist(kappa, 0, .1, .2, .5, 1)
assert (kappa == 0 & design == "kink") | (kappa > 0 & design == "jump")
assert (miss == "balanced" & missp == 0) | ///
    (miss == "mcar" & inlist(missp, .15, .30)) | ///
    (miss == "attrit" & missp == .15)
assert (dgp_spec == "gs26_endog_v2" & qstatus == "endog") | ///
    (dgp_spec != "gs26_endog_v2" & qstatus == "pred")

count if study == "GSCORE"
assert r(N) == 60
count if study == "GSSMALL"
assert r(N) == 16
count if study == "PERSIST"
assert r(N) == 8
count if study == "MISSPAT"
assert r(N) == 16
count if study == "ENDOG"
assert r(N) == 8
count if study == "HEAVYTAIL"
assert r(N) == 8
count if study == "ENDOG" & miss == "balanced"
assert r(N) == 4
count if study == "ENDOG" & miss == "mcar" & missp == .30
assert r(N) == 4
count if study == "HEAVYTAIL" & miss == "balanced"
assert r(N) == 4
count if study == "HEAVYTAIL" & miss == "mcar" & missp == .30
assert r(N) == 4
count if R == 400
assert r(N) == 92
count if R == 500
assert r(N) == 24
count if N == 1600 & R == 400
assert r(N) == 20
egen long total_reps = total(R)
assert total_reps[1] == 48800

file open ok using "_registry_probe.ok", write replace text
file write ok "cells=116;pairs=58;reps=48800;R400=92;R500=24;N1600_R400=20;grid=199;refine=4;smoke28=4-5" _n
file close ok
exit 0, clear
