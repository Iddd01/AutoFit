*!==========================================================================
*! study_a_shard.do -- resumable shard for point-estimation Study A
*! Harness: xtdpthresh_study_a_v0937_r400 (2026-08-17)
*!--------------------------------------------------------------------------
*! Args: SHARD NSHARD MASTER REP_CAP RUN_ID HARNESS_VERSION
*!       CODE_VERSION_EXPECTED RNG_KIND
*!==========================================================================
version 15.0
clear all
set more off
set varabbrev off

args SHARD NSHARD MASTER RCAP RUN_ID HARNESS_VERSION ///
     CODE_VERSION_EXPECTED RNG_KIND

if "`SHARD'"                 == "" local SHARD 1
if "`NSHARD'"                == "" local NSHARD 1
if "`MASTER'"                == "" local MASTER 20260814
if "`RCAP'"                  == "" local RCAP 0
if "`RUN_ID'"                == "" local RUN_ID legacy_manual
if "`HARNESS_VERSION'"       == "" local HARNESS_VERSION xtdpthresh_study_a_v0937_r400
if "`CODE_VERSION_EXPECTED'" == "" local CODE_VERSION_EXPECTED 0.9.37
if "`RNG_KIND'"              == "" local RNG_KIND mt64

foreach z in SHARD NSHARD MASTER RCAP {
    capture confirm integer number ``z''
    if _rc {
        di as err "study_a_shard: `z' must be an integer"
        exit 198
    }
}
if `NSHARD' < 1 | `SHARD' < 1 | `SHARD' > `NSHARD' | ///
        `MASTER' < 0 | `MASTER' > 2147483000 | `RCAP' < 0 | `RCAP' > 500 {
    di as err "study_a_shard: invalid shard/seed/rep-cap contract"
    exit 198
}
foreach tok in RUN_ID HARNESS_VERSION CODE_VERSION_EXPECTED RNG_KIND {
    if !regexm("``tok''", "^[A-Za-z0-9][A-Za-z0-9_.-]*$") {
        di as err "study_a_shard: `tok' is empty or unsafe"
        exit 198
    }
}
if !inlist("`RNG_KIND'", "mt64", "kiss32") {
    di as err "study_a_shard: RNG_KIND must be mt64 or kiss32"
    exit 198
}
if "`HARNESS_VERSION'" != "xtdpthresh_study_a_v0937_r400" | ///
   "`CODE_VERSION_EXPECTED'" != "0.9.37" | "`RNG_KIND'" != "mt64" {
    di as err "study_a_shard: unsupported frozen Study-A contract"
    exit 198
}
foreach f in study_a_worker.do study_a_cells.do xtdpthresh.ado xtdpthresh_p.ado {
    capture confirm file "`f'"
    if _rc {
        di as err "study_a_shard: required staged file is missing: `f'"
        exit 601
    }
}

global sta_shard        `SHARD'
global sta_nshard       `NSHARD'
global sta_master       `MASTER'
global sta_rep_cap      `RCAP'
global sta_run_id       "`RUN_ID'"
global sta_harness      "`HARNESS_VERSION'"
global sta_code_version "`CODE_VERSION_EXPECTED'"
global sta_rng          "`RNG_KIND'"
global sta_outstem      "study_a_results"
global sta_cells_seen   0
global sta_reps_seen    0
global sta_cell_index   0
global sta_last_cell    ""
global sta_last_rep     0

capture program drop _cell
program define _cell
    version 15.0
    args PAIR_ID STUDY DESIGN METHOD N T MISS MISSP QSTATUS R KAPPA ///
         DGP_SPEC GRID REFINE

    local CELL_ID "`PAIR_ID'_`METHOD'"
    global sta_cell_index = ${sta_cell_index} + 1
    local RR = `R'
    if ${sta_rep_cap} > 0 & `RR' > ${sta_rep_cap} local RR = ${sta_rep_cap}

    * A one-replication smoke distributes registry cells round-robin so all
    * shards do useful work.  Formal/capped runs otherwise retain the frozen
    * within-cell replication-range partition.
    if ${sta_rep_cap} == 1 & ${sta_nshard} > 1 {
        local owner = mod(${sta_cell_index}-1, ${sta_nshard}) + 1
        if `owner' != ${sta_shard} exit
        local rs 1
        local re 1
    }
    else {
        local rs = floor((${sta_shard}-1)*`RR'/${sta_nshard}) + 1
        local re = floor(${sta_shard}*`RR'/${sta_nshard})
    }
    if `re' < `rs' exit

    global sta_cells_seen = ${sta_cells_seen} + 1
    global sta_reps_seen  = ${sta_reps_seen} + (`re' - `rs' + 1)
    di as txt ">> `CELL_ID': reps `rs'-`re' / `RR'; point-only grid=`GRID' refine=`REFINE'"

    do study_a_worker.do `STUDY' `DESIGN' `METHOD' `N' `T' `MISS' `MISSP' ///
        `QSTATUS' `RR' `GRID' `REFINE' ${sta_master} `rs' `re' ///
        ${sta_shard} ${sta_outstem} ${sta_run_id} ${sta_harness} ///
        ${sta_nshard} ${sta_code_version} ${sta_rng} `DGP_SPEC' `KAPPA' ///
        `CELL_ID' `PAIR_ID'

    global sta_last_cell "`CELL_ID'"
    global sta_last_rep  `re'
    tempname PROG
    local PROGTMP "study_a_progress_SH${sta_shard}.csv.tmp"
    local PROGOUT "study_a_progress_SH${sta_shard}.csv"
    file open `PROG' using "`PROGTMP'", write replace text
    file write `PROG' "run_id,harness_version,code_version_expected,master,shard,nshard,cells_completed,reps_completed,last_cell_id,last_rep" _n
    file write `PROG' "${sta_run_id},${sta_harness},${sta_code_version},${sta_master},${sta_shard},${sta_nshard},${sta_cells_seen},${sta_reps_seen},${sta_last_cell},${sta_last_rep}" _n
    file close `PROG'
    quietly copy "`PROGTMP'" "`PROGOUT'", replace
    quietly erase "`PROGTMP'"
end

capture noisily do study_a_cells.do
local registry_rc = _rc
if `registry_rc' {
    di as err "study_a_shard: registry/worker aborted with rc=`registry_rc'"
    exit `registry_rc'
}

tempname DONE
file open `DONE' using "study_a_complete_SH`SHARD'.csv", write replace text
file write `DONE' "run_id,harness_version,code_version_expected,master,shard,nshard,cells,reps,exit_code" _n
file write `DONE' "`RUN_ID',`HARNESS_VERSION',`CODE_VERSION_EXPECTED',`MASTER',`SHARD',`NSHARD',${sta_cells_seen},${sta_reps_seen},0" _n
file close `DONE'

di as res "STUDY A SHARD `SHARD'/`NSHARD' COMPLETE: cells=${sta_cells_seen}; reps=${sta_reps_seen}"
