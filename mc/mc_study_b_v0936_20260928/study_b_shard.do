*!==========================================================================
*! study_b_shard.do -- resumable shard of the frozen Type-4 FOD Monte Carlo
*! Harness: study_b_type4_v0936 (2026-09-28)
*!--------------------------------------------------------------------------
*! Args:
*!   SHARD NSHARD MASTER [B_OVERRIDE GRID_OVERRIDE GRIDCI_OVERRIDE REP_CAP]
*!   [RUN_ID HARNESS_VERSION CODE_VERSION_EXPECTED RNG_KIND]
*!
*! Use a dot (.) for an unset numeric override.  The formal registry values
*! are cell-specific: B1--B5 use B=499/grid=199/refine=4, whereas B6 uses
*! B=500/grid=46/refine=0. Overrides and REP_CAP exist only for controlled
*! sensitivity/smoke runs and are recorded in every output row.
*!==========================================================================
version 15.0
clear all
set more off
set varabbrev off

args SHARD NSHARD MASTER Bov Gov GCov RCAP RUN_ID HARNESS_VERSION ///
     CODE_VERSION_EXPECTED RNG_KIND

if "`SHARD'"                == "" local SHARD 1
if "`NSHARD'"               == "" local NSHARD 1
if "`MASTER'"               == "" local MASTER 20260814
if "`Bov'"                  == "" local Bov .
if "`Gov'"                  == "" local Gov .
if "`GCov'"                 == "" local GCov .
if "`RCAP'"                 == "" local RCAP 0
if "`RUN_ID'"               == "" local RUN_ID legacy_manual
if "`HARNESS_VERSION'"      == "" local HARNESS_VERSION study_b_type4_v0936
if "`CODE_VERSION_EXPECTED'"== "" local CODE_VERSION_EXPECTED 0.9.36
if "`RNG_KIND'"             == "" local RNG_KIND mt64

foreach z in SHARD NSHARD MASTER RCAP {
    capture confirm integer number ``z''
    if _rc {
        di as err "study_b_shard: `z' must be an integer"
        exit 198
    }
}
foreach z in Bov Gov GCov {
    if "``z''" != "." {
        capture confirm integer number ``z''
        if _rc {
            di as err "study_b_shard: `z' must be an integer or ."
            exit 198
        }
    }
}
if `NSHARD' < 1 | `SHARD' < 1 | `SHARD' > `NSHARD' | ///
        `MASTER' < 0 | `MASTER' > 2147483000 | `RCAP' < 0 | `RCAP' > 500 {
    di as err "study_b_shard: invalid shard/seed/rep-cap contract"
    exit 198
}
if ("`Bov'" != "." & (`Bov' < 10 | `Bov' > 100000)) | ///
   ("`Gov'" != "." & (`Gov' < 10 | `Gov' > 10000)) | ///
   ("`GCov'" != "." & (`GCov' < 10 | `GCov' > 10000)) {
    di as err "study_b_shard: B/grid override outside supported range"
    exit 198
}
foreach tok in RUN_ID HARNESS_VERSION CODE_VERSION_EXPECTED RNG_KIND {
    if !regexm("``tok''", "^[A-Za-z0-9][A-Za-z0-9_.-]*$") {
        di as err "study_b_shard: `tok' is empty or contains an unsafe character"
        exit 198
    }
}
if !inlist("`RNG_KIND'", "mt64", "kiss32") {
    di as err "study_b_shard: RNG_KIND must be mt64 or kiss32"
    exit 198
}
if "`HARNESS_VERSION'" != "study_b_type4_v0936" | ///
   "`CODE_VERSION_EXPECTED'" != "0.9.36" | "`RNG_KIND'" != "mt64" {
    di as err "study_b_shard: unsupported frozen Study-B contract"
    exit 198
}
foreach f in study_b_worker.do study_b_cells.do xtdpthresh.ado xtdpthresh_p.ado {
    capture confirm file "`f'"
    if _rc {
        di as err "study_b_shard: required staged file is missing: `f'"
        exit 601
    }
}

global mc_shard        `SHARD'
global mc_nshard       `NSHARD'
global mc_master       `MASTER'
global mc_B_override   `Bov'
global mc_G_override   `Gov'
global mc_GC_override  `GCov'
global mc_rep_cap      `RCAP'
global mc_run_id       "`RUN_ID'"
global mc_harness      "`HARNESS_VERSION'"
global mc_code_version "`CODE_VERSION_EXPECTED'"
global mc_rng          "`RNG_KIND'"
global mc_outstem      "study_b_results"
global mc_cells_seen   0
global mc_reps_seen    0
global mc_cell_index   0
global mc_last_cell    ""
global mc_last_rep     0

capture program drop _cell
program define _cell
    version 15.0
    args PAIR_ID STUDY DESIGN METHOD N T MISS MISSP QSTATUS R KAPPA DGP_SPEC ///
         BOOTTYPE COEFBOOT GRID GRIDCI B REFINE INFMODE ///
         MAXLAG_LO MAXLAG_HI TRIM GRIDTYPE GRIDSAMPLE PAIRMODE

    if "`INFMODE'" == "" local INFMODE FULL

    * Registry identity is immutable; METHOD completes the primary cell key.
    local CELL_ID "`PAIR_ID'_`METHOD'"
    global mc_cell_index = ${mc_cell_index} + 1
    local RR = `R'
    if ${mc_rep_cap} > 0 & `RR' > ${mc_rep_cap} local RR = ${mc_rep_cap}
    local BB = `B'
    local GG = `GRID'
    local GC = `GRIDCI'
    if "${mc_B_override}"  != "." local BB = ${mc_B_override}
    if "${mc_G_override}"  != "." local GG = ${mc_G_override}
    if "${mc_GC_override}" != "." local GC = ${mc_GC_override}

    * The one-replication smoke assigns whole cells round-robin so 28 shards
    * exercise the complete 52-cell registry with only 52 expensive fits.
    if ${mc_rep_cap} == 1 & ${mc_nshard} > 1 {
        local owner = mod(${mc_cell_index}-1, ${mc_nshard}) + 1
        if `owner' != ${mc_shard} exit
        local rs 1
        local re 1
    }
    else {
        local rs = floor((${mc_shard}-1)*`RR'/${mc_nshard}) + 1
        local re = floor(${mc_shard}*`RR'/${mc_nshard})
    }
    if `re' < `rs' exit

    global mc_cells_seen = ${mc_cells_seen} + 1
    global mc_reps_seen  = ${mc_reps_seen} + (`re' - `rs' + 1)
    di as txt ">> `CELL_ID': reps `rs'-`re' / `RR'; B=`BB', grid=`GG'/`GC', maxlag=`MAXLAG_LO':`MAXLAG_HI'"

    do study_b_worker.do `STUDY' `DESIGN' `METHOD' `N' `T' `MISS' `MISSP' ///
        `QSTATUS' `RR' `BB' `GG' `GC' `REFINE' ${mc_master} `rs' `re' ///
        ${mc_shard} ${mc_outstem} `BOOTTYPE' `COEFBOOT' ${mc_run_id} ///
        ${mc_harness} ${mc_nshard} ${mc_code_version} ${mc_rng} ///
        `DGP_SPEC' `KAPPA' `CELL_ID' `PAIR_ID' `INFMODE' ///
        `MAXLAG_LO' `MAXLAG_HI' `TRIM' `GRIDTYPE' `GRIDSAMPLE' `PAIRMODE'

    global mc_last_cell "`CELL_ID'"
    global mc_last_rep  `re'
    tempname PROG
    local PROGTMP "study_b_progress_SH${mc_shard}.csv.tmp"
    local PROGOUT "study_b_progress_SH${mc_shard}.csv"
    file open `PROG' using "`PROGTMP'", write replace text
    file write `PROG' "run_id,harness_version,code_version_expected,master,shard,nshard,cells_completed,reps_completed,last_cell_id,last_rep" _n
    file write `PROG' "${mc_run_id},${mc_harness},${mc_code_version},${mc_master},${mc_shard},${mc_nshard},${mc_cells_seen},${mc_reps_seen},${mc_last_cell},${mc_last_rep}" _n
    file close `PROG'
    quietly copy "`PROGTMP'" "`PROGOUT'", replace
    quietly erase "`PROGTMP'"
end

capture noisily do study_b_cells.do
local registry_rc = _rc
if `registry_rc' {
    di as err "study_b_shard: registry/worker aborted with rc=`registry_rc'"
    exit `registry_rc'
}

* An atomic-enough completion marker: it is created only after every cell has
* returned successfully.  Runner/monitor/merge validate its exact contents.
tempname DONE
file open `DONE' using "study_b_complete_SH`SHARD'.csv", write replace text
file write `DONE' "run_id,harness_version,code_version_expected,master,shard,nshard,cells,reps,exit_code" _n
file write `DONE' "`RUN_ID',`HARNESS_VERSION',`CODE_VERSION_EXPECTED',`MASTER',`SHARD',`NSHARD',${mc_cells_seen},${mc_reps_seen},0" _n
file close `DONE'

di as res "SHARD `SHARD'/`NSHARD' COMPLETE: cells=${mc_cells_seen}; reps=${mc_reps_seen}"
