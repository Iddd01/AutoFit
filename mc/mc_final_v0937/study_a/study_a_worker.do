*!==========================================================================
*! study_a_worker.do -- Study A: point-estimator properties (NO BOOTSTRAP)
*! Schema: study_a_worker_v0937 (2026-08-17)
*!--------------------------------------------------------------------------
*! Gong--Seo point-estimator family with targeted small-N/T, persistence,
*! missing-pattern, endogenous-q, and heavy-tail stresses.  Every DGP keeps
*! the same five estimable coefficients (Lag_y_b, q_b, cons_d, q_d, Lag_y_d).
*!
*! Args:
*!   STUDY DESIGN METHOD N T MISS MISSP QSTATUS R GRID REFINE MASTER ///
*!   REP_START REP_END SHARD OUTSTEM RUN_ID HARNESS_VERSION NSHARD ///
*!   CODE_VERSION_EXPECTED RNG_KIND DGP_SPEC KAPPA CELL_ID PAIR_ID
*!
*! Invariants:
*!   - -noboot- is literal and unconditional in the only estimation call.
*!   - METHOD is excluded from both seeds, giving exact FD/FOD CRN pairing.
*!   - completed failures are resumable outcomes and are never retried.
*!   - DESIGN=kink labels the continuous kappa=0 DGP; every cell estimates
*!     the unrestricted jump model (e(flag_kink)=0) with five coefficients.
*!   - success is the frozen primary estimator: contract-valid two-step GMM
*!     over a fixed 199-point global grid followed by refine(4) under one W2.
*!     Retired adaptive-search fields are retained only as missing/zero schema
*!     compatibility fields; they never select or retry replications.
*!     Contract-valid one-step fallback remains delivered but not successful.
*!==========================================================================
version 15.0
clear all
set more off
set varabbrev off

args STUDY DESIGN METHOD N T MISS MISSP QSTATUS R GRID REFINE MASTER ///
     REP_START REP_END SHARD OUTSTEM RUN_ID HARNESS_VERSION NSHARD ///
     CODE_VERSION_EXPECTED RNG_KIND DGP_SPEC KAPPA CELL_ID PAIR_ID

if "`RUN_ID'"                == "" local RUN_ID legacy_manual
if "`HARNESS_VERSION'"       == "" local HARNESS_VERSION xtdpthresh_study_a_v0937_r400
if "`NSHARD'"                == "" local NSHARD 1
if "`CODE_VERSION_EXPECTED'" == "" local CODE_VERSION_EXPECTED 0.9.37
if "`RNG_KIND'"              == "" local RNG_KIND mt64

local SCHEMA_VERSION study_a_worker_v0937
local SEED_SCHEME    study_a_crn_v3
local INFERENCE_MODE point
local FORMAL_HARNESS  = ("`HARNESS_VERSION'" == "xtdpthresh_study_a_v0937_r400")

*----- fail-fast argument contract -----------------------------------------
if !`FORMAL_HARNESS' {
    di as err "study_a_worker: unsupported harness version"
    exit 198
}
if !inlist("`STUDY'", "GSCORE", "GSSMALL", "PERSIST", "MISSPAT", ///
        "ENDOG", "HEAVYTAIL") {
    di as err "study_a_worker: STUDY is outside the v6 registry"
    exit 198
}
if !inlist("`DESIGN'", "kink", "jump") {
    di as err "study_a_worker: DESIGN must be kink or jump"
    exit 198
}
if !inlist("`METHOD'", "fd", "fod") {
    di as err "study_a_worker: METHOD must be fd or fod"
    exit 198
}
if !inlist("`MISS'", "balanced", "mcar", "attrit") {
    di as err "study_a_worker: MISS must be balanced, mcar, or attrit"
    exit 198
}
if !inlist("`QSTATUS'", "pred", "endog") {
    di as err "study_a_worker: QSTATUS must be pred or endog"
    exit 198
}
if !inlist("`DGP_SPEC'", "gs26_base_v2", "gs26_persist_v2", ///
        "gs26_endog_v2", "gs26_t5_v2") {
    di as err "study_a_worker: DGP_SPEC is outside the v2 registry"
    exit 198
}
if !inlist("`RNG_KIND'", "mt64", "kiss32") {
    di as err "study_a_worker: RNG_KIND must be mt64 or kiss32"
    exit 198
}
foreach tok in STUDY RUN_ID HARNESS_VERSION CODE_VERSION_EXPECTED DGP_SPEC ///
                   CELL_ID PAIR_ID {
    if !regexm("``tok''", "^[A-Za-z0-9][A-Za-z0-9_.-]*$") {
        di as err "study_a_worker: `tok' is empty or contains an unsafe character"
        exit 198
    }
}
if "`OUTSTEM'" == "" | strpos("`OUTSTEM'", ",") | ///
        strpos("`OUTSTEM'", char(34)) | strpos("`OUTSTEM'", char(10)) | ///
        strpos("`OUTSTEM'", char(13)) {
    di as err "study_a_worker: OUTSTEM is empty or unsafe for CSV output"
    exit 198
}
foreach z in N T R GRID REFINE MASTER REP_START REP_END SHARD NSHARD {
    capture confirm integer number ``z''
    if _rc {
        di as err "study_a_worker: `z' must be an integer; received [``z'']"
        exit 198
    }
}
foreach z in MISSP KAPPA {
    capture confirm number ``z''
    if _rc {
        di as err "study_a_worker: `z' must be numeric; received [``z'']"
        exit 198
    }
}
if !inlist(`N', 100, 200, 400, 800, 1600) | !inlist(`T', 6, 10) | ///
        `R' < 1 | `R' > 500 {
    di as err "study_a_worker: N/T/R lies outside the supported Study A range"
    exit 198
}
if `GRID' != 199 | `REFINE' != 4 {
    di as err "study_a_worker: formal Study A v6 is frozen at grid(199), refine(4)"
    exit 198
}
if `MASTER' < 0 | `MASTER' > 2147483000 | ///
        `REP_START' < 1 | `REP_END' < `REP_START' | `REP_END' > `R' | ///
        `NSHARD' < 1 | `SHARD' < 1 | `SHARD' > `NSHARD' {
    di as err "study_a_worker: invalid seed, replication interval, or shard contract"
    exit 198
}
if ("`MISS'" == "balanced" & `MISSP' != 0) | ///
        ("`MISS'" == "mcar" & !inlist(`MISSP', .15, .30)) | ///
        ("`MISS'" == "attrit" & `MISSP' != .15) {
    di as err "study_a_worker: invalid MISS/MISSP combination"
    exit 198
}
if !inlist(`KAPPA', 0, .1, .2, .5, 1) | ///
        (`KAPPA' == 0 & "`DESIGN'" != "kink") | ///
        (`KAPPA' > 0 & "`DESIGN'" != "jump") {
    di as err "study_a_worker: DESIGN/KAPPA is outside the frozen registry"
    exit 198
}
if ("`DGP_SPEC'" == "gs26_endog_v2") != ("`QSTATUS'" == "endog") {
    di as err "study_a_worker: gs26_endog_v2 and QSTATUS=endog must occur together"
    exit 198
}
if "`DGP_SPEC'" != "gs26_endog_v2" & "`QSTATUS'" != "pred" {
    di as err "study_a_worker: non-endogenous v2 DGPs require QSTATUS=pred"
    exit 198
}
if "`DGP_SPEC'" == "gs26_persist_v2" & "`STUDY'" != "PERSIST" {
    di as err "study_a_worker: persistent DGP is reserved for PERSIST cells"
    exit 198
}
if "`DGP_SPEC'" == "gs26_t5_v2" & "`STUDY'" != "HEAVYTAIL" {
    di as err "study_a_worker: t(5) DGP is reserved for HEAVYTAIL cells"
    exit 198
}

capture set rng `RNG_KIND'
if _rc {
    di as err "study_a_worker: RNG_KIND=`RNG_KIND' is unavailable"
    exit 198
}
local RNG_ACTUAL "`c(rng_current)'"
local STATA_VERSION "`c(stata_version)'"
capture set processors 1
adopath ++ "."

*----- DGP constants -------------------------------------------------------
local RHOY   = cond("`DGP_SPEC'" == "gs26_persist_v2", .9, .6)
local BQ     1
local D2     2
local GAMMA0 .25
local RHOQ   .7
local RHOEU  = .5
local SIGE   .5
local TBURN  = cond("`DGP_SPEC'" == "gs26_persist_v2", 50, 20)
local TMAX   = `T' + `TBURN'

local T_Lagyb `RHOY'
local T_qb    `BQ'
local T_consd = `KAPPA' - `D2'*`GAMMA0'
local T_qd    `D2'
local T_Lagyd 0

* Numeric namespaces used by both fresh-run and resume seed checks.
local SPCODE = cond("`DGP_SPEC'" == "gs26_base_v2", 1, ///
               cond("`DGP_SPEC'" == "gs26_persist_v2", 2, ///
               cond("`DGP_SPEC'" == "gs26_endog_v2", 3, 4)))
local QCODE  = cond("`QSTATUS'" == "pred", 1, 2)
local MCODE  = cond("`MISS'" == "balanced", 1, ///
               cond("`MISS'" == "mcar", 2, 3))

*----- output and fail-closed resume contract ------------------------------
local CONFIG_KEY "`STUDY'|`DESIGN'|`METHOD'|N=`N'|T=`T'|`MISS'=`MISSP'|`QSTATUS'|R=`R'|G=`GRID'|RF=`REFINE'|K=`KAPPA'|`DGP_SPEC'|point"
local OUT "`OUTSTEM'_SH`SHARD'.csv"
capture confirm file "`OUT'"
local newfile = _rc

local EXPECT_HEADER "run_id,harness_version,schema_version,code_version_expected,code_version_actual,stata_version,rng_kind,rng_actual,seed_scheme,dgp_spec,inference_mode,cell_id,pair_id,config_key,full_key,study,design,method,N,T,miss,missp,qstatus,kappa,R,grid,refine,master,shard,nshard,rep,dgp_seed,missing_seed,rc,fit_rc0,cmd_ok,version_ok,contract_ok,success,gamma0,gamma_hat,g_bias,g_sqerr,obj,Lagyb,qb,consd,qd,Lagyd,b_Lagyb,b_qb,b_consd,b_qd,b_Lagyd,err_Lagyb,sqerr_Lagyb,err_qb,sqerr_qb,err_consd,sqerr_consd,err_qd,sqerr_qd,err_Lagyd,sqerr_Lagyd,se_Lagyb,se_qb,se_consd,se_qd,se_Lagyd,hansen_p,ar1_p,ar2_p,N_used,N_stack,N_trans,N_iv,N_units,estimator_twostep,flag_kink,refine_complete,refine_iterations,refine_added,refine_remaining,refine_neigh_unevaluated,searchmode,searchtol,search_converged,search_incomplete,search_hit_max,search_W2_builds,search_stage1_level,search_stage2_level,search_stage1_points,search_stage2_points,search_stage1_same_split,search_stage2_same_split,search_stage1_rel_gain,search_stage2_rel_gain,grid_max_requested,point_boundary,q_lo,q_hi,gamma_grid1_lo,gamma_grid1_hi,gamma_grid2_lo,gamma_grid2_hi,grid_requested,grid_effective,grid_admitted,grid_structural,grid_twostep_admitted,minregime_applied,units_target,units_realized,units_dropped,analysis_potential,analysis_observed,missing_n,missing_rate,gap_events,gap_periods,fd_pair_rows_potential,fod_rows_potential,rhoy_true,bq_true,consd_true,qd_true,lagyd_true,rhoq,rhoeu_effective,sige,tburn,maxlag_lo,maxlag_hi,trim,history,gridtype,gridsample,vce,elapsed_s,joint_vce,ar_joint,vce_applied,bwscale,gamma_bw,q_nvals_bw,N_iv_dep,N_iv_dep_near,iv_dep_res,ar1_cond,ar2_cond,ar1_p_cond,ar2_p_cond,se_cond_Lagyb,se_cond_qb,se_cond_consd,se_cond_qd,se_cond_Lagyd,se_delivered"
local EXPECT_COMMAS = length("`EXPECT_HEADER'") - ///
    length(subinstr("`EXPECT_HEADER'", ",", "", .))

local done_reps
if !`newfile' {
    tempname RH
    file open `RH' using "`OUT'", read text
    file read `RH' firstline
    if `"`firstline'"' != `"`EXPECT_HEADER'"' {
        file close `RH'
        di as err "study_a_worker: existing `OUT' has an incompatible schema"
        di as err "Archive it or select a new RUN_ID/run directory."
        exit 459
    }

    * A crash can interrupt only the final append.  Auto-repair exactly that
    * case: every prior physical line must have the frozen separator count and
    * only the final physical line may differ.  Any non-final corruption stays
    * fail-closed.  The imported semantic checks below still decide what skips.
    tempname WH
    local REPAIRTMP "`OUT'.repair.tmp"
    file open `WH' using "`REPAIRTMP'", write replace text
    file write `WH' `"`firstline'"' _n
    local physical_line 1
    local malformed_line 0
    local complete_rows 0
    file read `RH' thisline
    while r(eof) == 0 {
        local ++physical_line
        local ncomma = length(`"`thisline'"') - ///
            length(subinstr(`"`thisline'"', ",", "", .))
        if `ncomma' != `EXPECT_COMMAS' {
            if `malformed_line' {
                file close `RH'
                file close `WH'
                capture erase "`REPAIRTMP'"
                di as err "study_a_worker: multiple malformed rows in resume CSV"
                exit 459
            }
            local malformed_line `physical_line'
        }
        else {
            if `malformed_line' {
                file close `RH'
                file close `WH'
                capture erase "`REPAIRTMP'"
                di as err "study_a_worker: malformed non-final row in resume CSV"
                exit 459
            }
            file write `WH' `"`thisline'"' _n
            local ++complete_rows
        }
        file read `RH' thisline
    }
    file close `RH'
    file close `WH'
    if `malformed_line' {
        quietly copy "`REPAIRTMP'" "`OUT'", replace
        di as txt "study_a_worker: removed demonstrably torn final CSV line `malformed_line'; replication will be rerun"
    }
    capture erase "`REPAIRTMP'"

    if `complete_rows' {
        capture quietly import delimited using "`OUT'", clear varnames(1) ///
            case(preserve) stringcols(_all) bindquote(strict)
        if _rc {
            di as err "study_a_worker: cannot parse existing resume file `OUT'"
            exit 459
        }
        * import all columns as string (so numeric-looking string columns such
        * as stata_version stay string), then restore the genuinely-numeric
        * columns; this makes the string/numeric confirms below type-consistent
        * regardless of the values in any given resume file.
        quietly destring N T missp kappa R grid refine master shard nshard ///
            rep dgp_seed missing_seed rc fit_rc0 cmd_ok version_ok ///
            contract_ok success gamma0 gamma_hat g_bias g_sqerr obj ///
            Lagyb qb consd qd Lagyd b_Lagyb b_qb b_consd b_qd b_Lagyd ///
            err_Lagyb sqerr_Lagyb err_qb sqerr_qb err_consd sqerr_consd ///
            err_qd sqerr_qd err_Lagyd sqerr_Lagyd ///
            se_Lagyb se_qb se_consd se_qd se_Lagyd ///
            hansen_p ar1_p ar2_p N_used N_stack N_trans N_iv N_units ///
            estimator_twostep flag_kink refine_complete refine_iterations ///
            refine_added refine_remaining refine_neigh_unevaluated ///
            searchtol search_converged search_incomplete search_hit_max ///
            search_W2_builds search_stage1_level search_stage2_level ///
            search_stage1_points search_stage2_points ///
            search_stage1_same_split search_stage2_same_split ///
            search_stage1_rel_gain search_stage2_rel_gain ///
            grid_max_requested point_boundary q_lo q_hi gamma_grid1_lo ///
            gamma_grid1_hi gamma_grid2_lo gamma_grid2_hi grid_requested ///
            grid_effective grid_admitted grid_structural ///
            grid_twostep_admitted minregime_applied units_target ///
            units_realized units_dropped analysis_potential ///
            analysis_observed missing_n missing_rate gap_events gap_periods ///
            fd_pair_rows_potential fod_rows_potential rhoy_true bq_true ///
            consd_true qd_true lagyd_true rhoq rhoeu_effective sige tburn ///
            maxlag_lo maxlag_hi trim elapsed_s joint_vce ar_joint vce_applied bwscale gamma_bw q_nvals_bw N_iv_dep N_iv_dep_near iv_dep_res ar1_cond ar2_cond ar1_p_cond ar2_p_cond se_cond_Lagyb se_cond_qb se_cond_consd se_cond_qd se_cond_Lagyd se_delivered, replace
    foreach vv in run_id harness_version schema_version code_version_expected ///
                  code_version_actual stata_version rng_kind rng_actual ///
                  seed_scheme dgp_spec inference_mode cell_id pair_id ///
                  config_key full_key study design method N T miss missp ///
                  qstatus kappa R grid refine master shard nshard rep ///
                  dgp_seed missing_seed rc fit_rc0 cmd_ok version_ok ///
                   contract_ok success N_used N_stack N_trans N_iv N_units ///
                   estimator_twostep flag_kink refine_complete ///
                   refine_iterations ///
                   refine_added refine_remaining refine_neigh_unevaluated ///
                   searchmode obj searchtol search_converged ///
                   search_incomplete search_hit_max search_W2_builds ///
                   search_stage1_level search_stage2_level ///
                   search_stage1_points search_stage2_points ///
                   search_stage1_same_split search_stage2_same_split ///
                   search_stage1_rel_gain search_stage2_rel_gain ///
                   grid_max_requested ///
                   point_boundary q_lo q_hi gamma_grid1_lo gamma_grid1_hi ///
                   gamma_grid2_lo gamma_grid2_hi grid_requested ///
                   grid_effective grid_admitted grid_structural ///
                   grid_twostep_admitted minregime_applied ///
                   units_target units_realized ///
                  units_dropped analysis_potential analysis_observed ///
                  missing_n missing_rate gap_events gap_periods ///
                  maxlag_lo maxlag_hi trim history gridtype gridsample ///
                  vce elapsed_s {
        capture confirm variable `vv'
        if _rc {
            di as err "study_a_worker: resume file lacks mandatory column `vv'"
            exit 459
        }
    }
    foreach vv in run_id harness_version schema_version code_version_expected ///
                  code_version_actual stata_version rng_kind rng_actual ///
                  seed_scheme dgp_spec inference_mode cell_id pair_id ///
                  config_key full_key study design method miss qstatus ///
                  searchmode history gridtype gridsample vce {
        capture confirm string variable `vv'
        if _rc {
            di as err "study_a_worker: incomplete/truncated resume row or incompatible string column `vv'"
            exit 459
        }
    }
    foreach vv in N T missp kappa R grid refine master shard nshard rep ///
                  dgp_seed missing_seed rc fit_rc0 cmd_ok version_ok ///
                  contract_ok success gamma0 gamma_hat g_bias g_sqerr obj ///
                  Lagyb qb consd qd Lagyd b_Lagyb b_qb b_consd b_qd ///
                  b_Lagyd err_Lagyb sqerr_Lagyb err_qb sqerr_qb ///
                  err_consd sqerr_consd err_qd sqerr_qd err_Lagyd ///
                  sqerr_Lagyd se_Lagyb se_qb se_consd se_qd se_Lagyd ///
                  hansen_p ar1_p ar2_p N_used N_stack N_trans N_iv N_units ///
                  estimator_twostep flag_kink refine_complete ///
                  refine_iterations refine_added refine_remaining ///
                  refine_neigh_unevaluated searchtol search_converged ///
                  search_incomplete search_hit_max search_W2_builds ///
                  search_stage1_level search_stage2_level ///
                  search_stage1_points search_stage2_points ///
                  search_stage1_same_split search_stage2_same_split ///
                  search_stage1_rel_gain search_stage2_rel_gain ///
                  grid_max_requested point_boundary q_lo q_hi ///
                  gamma_grid1_lo gamma_grid1_hi gamma_grid2_lo ///
                  gamma_grid2_hi grid_requested grid_effective ///
                  grid_admitted grid_structural grid_twostep_admitted ///
                  minregime_applied units_target units_realized ///
                  units_dropped analysis_potential analysis_observed ///
                  missing_n missing_rate gap_events gap_periods ///
                  fd_pair_rows_potential fod_rows_potential rhoy_true ///
                  bq_true consd_true qd_true lagyd_true rhoq ///
                  rhoeu_effective sige tburn maxlag_lo maxlag_hi trim ///
                  elapsed_s {
        capture confirm numeric variable `vv'
        if _rc {
            di as err "study_a_worker: incomplete/truncated resume row or incompatible numeric column `vv'"
            exit 459
        }
    }
    quietly count if full_key == "" | cell_id == "" | pair_id == "" | ///
        config_key == "" | code_version_actual == "" | stata_version == "" | ///
        rng_actual == "" | method == "" | history == "" | ///
        gridtype == "" | gridsample == "" | vce == "" | ///
        missing(N) | missing(T) | missing(kappa) | missing(R) | ///
        missing(grid) | missing(refine) | missing(master) | ///
        missing(shard) | missing(nshard) | missing(rep) | ///
        missing(dgp_seed) | missing(missing_seed) | missing(rc) | ///
        missing(fit_rc0) | missing(cmd_ok) | missing(version_ok) | ///
        missing(contract_ok) | missing(success) | missing(units_target) | ///
        missing(units_realized) | missing(units_dropped) | ///
        missing(analysis_potential) | missing(analysis_observed) | ///
        missing(missing_n) | missing(missing_rate) | ///
        missing(gap_events) | missing(gap_periods) | ///
        missing(fd_pair_rows_potential) | missing(fod_rows_potential) | ///
        missing(gamma0) | missing(Lagyb) | missing(qb) | ///
        missing(consd) | missing(qd) | missing(Lagyd) | ///
        missing(rhoy_true) | missing(bq_true) | missing(consd_true) | ///
        missing(qd_true) | missing(lagyd_true) | missing(rhoq) | ///
        missing(rhoeu_effective) | missing(sige) | missing(tburn) | ///
        missing(maxlag_lo) | missing(maxlag_hi) | missing(trim) | ///
        missing(elapsed_s)
    if r(N) {
        di as err "study_a_worker: incomplete/truncated row in existing `OUT'"
        exit 459
    }
    capture isid full_key
    if _rc {
        di as err "study_a_worker: duplicate full_key in existing `OUT'"
        exit 459
    }
    tempvar _expected_key
    quietly gen strL `_expected_key' = run_id + "|" + cell_id + ///
        "|rep=" + trim(string(rep, "%12.0f"))
    quietly count if full_key != `_expected_key'
    if r(N) {
        di as err "study_a_worker: inconsistent full_key in existing `OUT'"
        exit 459
    }
    quietly count if run_id != "`RUN_ID'" | ///
        harness_version != "`HARNESS_VERSION'" | ///
        schema_version != "`SCHEMA_VERSION'" | ///
        code_version_expected != "`CODE_VERSION_EXPECTED'" | ///
        rng_kind != "`RNG_KIND'" | seed_scheme != "`SEED_SCHEME'" | ///
        inference_mode != "`INFERENCE_MODE'" | ///
        master != `MASTER' | shard != `SHARD' | nshard != `NSHARD'
    if r(N) {
        di as err "study_a_worker: run-level provenance mismatch in `OUT'"
        exit 459
    }
    * DGP_SPEC varies by registered cell (A3/A5/A6 differ from A1/A2), so it
    * belongs to the cell contract and must never be imposed on all prior rows
    * in the shard-level CSV.
    quietly count if cell_id == "`CELL_ID'" & ///
        (study != "`STUDY'" | design != "`DESIGN'" | ///
         dgp_spec != "`DGP_SPEC'" | ///
         method != "`METHOD'" | N != `N' | T != `T' | ///
         miss != "`MISS'" | abs(missp-`MISSP') > 1e-12 | ///
         qstatus != "`QSTATUS'" | abs(kappa-`KAPPA') > 1e-12 | ///
         R != `R' | grid != `GRID' | refine != `REFINE' | ///
         history != "panel" | gridtype != "uniform" | ///
         gridsample != "effective" | vce != "robust" | ///
         units_target != `N' | analysis_potential != `N'*`T' | ///
         abs(gamma0-`GAMMA0') > 1e-12 | ///
         abs(Lagyb-`T_Lagyb') > 1e-12 | abs(qb-`T_qb') > 1e-12 | ///
         abs(consd-`T_consd') > 1e-12 | abs(qd-`T_qd') > 1e-12 | ///
         abs(Lagyd-`T_Lagyd') > 1e-12 | ///
         abs(rhoy_true-`RHOY') > 1e-12 | abs(bq_true-`BQ') > 1e-12 | ///
         abs(consd_true-`T_consd') > 1e-12 | ///
         abs(qd_true-`T_qd') > 1e-12 | ///
         abs(lagyd_true-`T_Lagyd') > 1e-12 | ///
         abs(rhoq-`RHOQ') > 1e-12 | ///
         abs(rhoeu_effective-`RHOEU') > 1e-12 | ///
         abs(sige-`SIGE') > 1e-12 | tburn != `TBURN' | ///
         maxlag_lo != 1 | maxlag_hi != 3 | abs(trim-.15) > 1e-12)
    if r(N) {
        di as err "study_a_worker: completed row has an invalid cell-level contract"
        exit 459
    }
    quietly count if cell_id == "`CELL_ID'" & ///
        (config_key != "`CONFIG_KEY'" | pair_id != "`PAIR_ID'")
    if r(N) {
        di as err "study_a_worker: CELL_ID reused with different configuration"
        exit 459
    }
    tempvar _h_expected _dseed_expected _mseed_expected
    quietly gen double `_h_expected' = mod(`MASTER'*1000003 + ///
        rep*10007 + `N'*101 + `T'*1009 + `SPCODE'*1000033 + ///
        `QCODE'*1000037, 699999999) ///
        if cell_id == "`CELL_ID'"
    quietly gen double `_dseed_expected' = 1 + `_h_expected' ///
        if cell_id == "`CELL_ID'"
    quietly gen double `_mseed_expected' = 700000001 + ///
        mod(`_h_expected'*1009 + `MCODE'*10007 + 7919, 699999999) ///
        if cell_id == "`CELL_ID'"
    quietly count if cell_id == "`CELL_ID'" & ///
        (dgp_seed != `_dseed_expected' | missing_seed != `_mseed_expected' | ///
         rep < `REP_START' | rep > `REP_END' | rep > `R' | ///
         !inlist(fit_rc0,0,1) | !inlist(cmd_ok,0,1) | ///
         !inlist(version_ok,0,1) | !inlist(contract_ok,0,1) | ///
         !inlist(success,0,1) | fit_rc0 != (rc == 0) | ///
         success != (contract_ok == 1 & estimator_twostep == 1) | ///
         contract_ok > fit_rc0 | ///
         contract_ok > cmd_ok | contract_ok > version_ok | ///
         units_realized + units_dropped != units_target | ///
         analysis_observed + missing_n != analysis_potential | ///
         abs(missing_rate-missing_n/analysis_potential) > 1e-12)
    if r(N) {
        di as err "study_a_worker: completed row has invalid seed/outcome/accounting semantics"
        exit 459
    }
    quietly count if cell_id == "`CELL_ID'" & contract_ok == 1 & ///
        (missing(N_used) | missing(N_stack) | missing(N_trans) | ///
         missing(N_iv) | missing(N_units) | missing(estimator_twostep) | ///
         missing(flag_kink) | ///
         missing(refine_complete) | missing(refine_iterations) | ///
         missing(refine_added) | missing(refine_remaining) | ///
         missing(refine_neigh_unevaluated) | missing(obj) | ///
         missing(gamma_hat) | missing(g_bias) | missing(g_sqerr) | ///
         missing(b_Lagyb) | missing(b_qb) | missing(b_consd) | ///
         missing(b_qd) | missing(b_Lagyd) | ///
         missing(err_Lagyb) | missing(sqerr_Lagyb) | ///
         missing(err_qb) | missing(sqerr_qb) | ///
         missing(err_consd) | missing(sqerr_consd) | ///
         missing(err_qd) | missing(sqerr_qd) | ///
         missing(err_Lagyd) | missing(sqerr_Lagyd) | ///
         se_Lagyb < 0 | se_qb < 0 | se_consd < 0 | ///
         se_qd < 0 | se_Lagyd < 0 | ///
         abs(g_bias-(gamma_hat-gamma0)) > 1e-8 | ///
         abs(g_sqerr-g_bias^2) > 1e-8 | ///
         abs(err_Lagyb-(b_Lagyb-Lagyb)) > 1e-8 | ///
         abs(sqerr_Lagyb-err_Lagyb^2) > 1e-8 | ///
         abs(err_qb-(b_qb-qb)) > 1e-8 | ///
         abs(sqerr_qb-err_qb^2) > 1e-8 | ///
         abs(err_consd-(b_consd-consd)) > 1e-8 | ///
         abs(sqerr_consd-err_consd^2) > 1e-8 | ///
         abs(err_qd-(b_qd-qd)) > 1e-8 | ///
         abs(sqerr_qd-err_qd^2) > 1e-8 | ///
         abs(err_Lagyd-(b_Lagyd-Lagyd)) > 1e-8 | ///
         abs(sqerr_Lagyd-err_Lagyd^2) > 1e-8 | ///
         searchmode == "" | missing(search_hit_max) | ///
         missing(search_W2_builds) | missing(search_stage1_level) | ///
         missing(search_stage2_level) | missing(search_stage1_points) | ///
         missing(search_stage2_points) | ///
         missing(grid_max_requested) | missing(point_boundary) | ///
         missing(q_lo) | missing(q_hi) | missing(gamma_grid1_lo) | ///
         missing(gamma_grid1_hi) | missing(grid_requested) | ///
         missing(grid_effective) | missing(grid_admitted) | ///
         missing(grid_structural) | missing(minregime_applied) | ///
         !inlist(estimator_twostep,0,1) | flag_kink != 0 | ///
         !inlist(point_boundary,0,1,2,3) | ///
         N_stack != N_trans | N_units < 1 | N_iv < 1 | ///
         q_lo >= q_hi | gamma_grid1_lo > gamma_grid1_hi | ///
         gamma_grid1_lo < q_lo | gamma_grid1_hi > q_hi | ///
         grid_requested != `GRID' | grid_max_requested != `GRID' | ///
         grid_effective != `GRID' + refine_added | ///
         grid_admitted < 2 | grid_structural < grid_admitted | ///
         grid_structural > grid_effective | minregime_applied < 1 | ///
         refine != `REFINE' | !inrange(refine_iterations,0,`REFINE') | ///
         refine_added < 0 | refine_remaining < 0 | ///
         refine_neigh_unevaluated < 0 | ///
         !inlist(refine_complete,0,1) | ///
         (refine_complete == 1 & ///
              (refine_remaining != 0 | refine_neigh_unevaluated != 0)) | ///
         searchmode != "fixed" | ///
         !missing(searchtol) | !missing(search_converged) | ///
         !missing(search_incomplete) | search_hit_max != 0 | ///
         search_W2_builds != 1 | search_stage1_level != 1 | ///
         !inlist(search_stage2_level,0,1) | ///
         search_stage1_points != `GRID' | ///
         (search_stage2_level==0 & search_stage2_points != 0) | ///
         (search_stage2_level==1 & search_stage2_points != `GRID') | ///
         !missing(search_stage1_same_split) | ///
         !missing(search_stage2_same_split) | ///
         !missing(search_stage1_rel_gain) | ///
         !missing(search_stage2_rel_gain) | ///
         (estimator_twostep == 1 & ///
             (missing(gamma_grid2_lo) | missing(gamma_grid2_hi) | ///
              missing(grid_twostep_admitted) | ///
              gamma_grid2_lo > gamma_grid2_hi | ///
              gamma_grid2_lo < q_lo | gamma_grid2_hi > q_hi | ///
              grid_twostep_admitted < 2 | ///
              grid_twostep_admitted > grid_effective)) | ///
         (estimator_twostep == 0 & ///
             ((missing(grid_twostep_admitted) & ///
                  (!missing(gamma_grid2_lo) | !missing(gamma_grid2_hi))) | ///
              (!missing(grid_twostep_admitted) & ///
                  (grid_twostep_admitted < 0 | ///
                   grid_twostep_admitted > grid_effective | ///
                   (grid_twostep_admitted == 0 & ///
                       (!missing(gamma_grid2_lo) | ///
                        !missing(gamma_grid2_hi))) | ///
                   (grid_twostep_admitted > 0 & ///
                       (missing(gamma_grid2_lo) | ///
                        missing(gamma_grid2_hi) | ///
                        gamma_grid2_lo > gamma_grid2_hi)))))))
    if r(N) {
        di as err "study_a_worker: completed contract-valid row has invalid search diagnostics"
        exit 459
    }
    capture quietly levelsof rep if cell_id == "`CELL_ID'", local(done_reps)
        clear
    }
}

* Open/close the CSV for every durable write.  This makes each completed row
* visible to the monitor, releases the Windows file lock between estimates,
* and limits an interrupted write to the final line (which resume rejects).
if `newfile' {
    tempname HH
    file open `HH' using "`OUT'", write replace text
    file write `HH' "`EXPECT_HEADER'" _n
    file close `HH'
}

forvalues rep = `REP_START'/`REP_END' {
    * A recorded failure is a completed MC result; never condition on retry.
    if strpos(" `done_reps' ", " `rep' ") {
        di as txt "  resume: skip `CELL_ID' rep `rep'"
        continue
    }

    * METHOD, DESIGN and KAPPA are deliberately absent from the base-shock
    * namespace.  Thus FD/FOD and kappa designs share the same innovations
    * within a DGP family; distinct stress DGPs have disjoint streams.
    local _h0 = mod(`MASTER'*1000003 + `rep'*10007 + `N'*101 + ///
        `T'*1009 + `SPCODE'*1000033 + `QCODE'*1000037, 699999999)
    local dgp_seed = 1 + `_h0'
    * v3: probability is deliberately absent.  MCAR .15 and .30 therefore
    * use the same uniforms, making the mild deletion set a strict subset of
    * the severe deletion set.  MCODE keeps attrition on a disjoint stream.
    local _hm = mod(`_h0'*1009 + `MCODE'*10007 + 7919, ///
        699999999)
    local missing_seed = 700000001 + `_hm'
    local FULL_KEY "`RUN_ID'|`CELL_ID'|rep=`rep'"

    timer clear 8
    timer on 8

    *----- Gong--Seo weak-exogeneity recursion ----------------------------
    clear
    quietly set seed `dgp_seed'
    quietly set obs `=`N'*`TMAX''
    quietly gen long id = ceil(_n/`TMAX')
    bysort id: gen int t = _n
    quietly xtset id t

    if "`DGP_SPEC'" == "gs26_t5_v2" {
        * Standardise t(5) by sqrt(5/3), so es and w each have variance one.
        quietly gen double es = rt(5)/sqrt(5/3)
        quietly gen double e  = es*`SIGE'
        quietly gen double w  = rt(5)/sqrt(5/3)
    }
    else {
        quietly gen double e  = rnormal()*`SIGE'
        quietly gen double es = e/`SIGE'
        quietly gen double w  = rnormal()
    }
    quietly gen double u  = .
    if "`QSTATUS'" == "endog" {
        * Gong--Seo Appendix C.3: contemporaneous corr(u_t,e*_t)=.5.
        quietly replace u = .5*es + sqrt(.75)*w
    }
    else {
        * Section-5 weak-exogeneity timing: corr(u_t,e*_{t-1})=.5.
        bysort id (t): replace u = .5*es[_n-1] + sqrt(.75)*w if _n > 1
        bysort id (t): replace u = w if _n == 1
    }

    quietly gen double q = .
    bysort id (t): replace q = rnormal() if _n == 1
    bysort id (t): replace q = `RHOQ'*q[_n-1] + u if _n > 1

    tempvar regime_effect
    quietly gen double `regime_effect' = ///
        (`T_consd' + `T_qd'*q)*(q > `GAMMA0')
    quietly gen double y = .
    bysort id (t): replace y = rnormal() if _n == 1
    bysort id (t): replace y = `RHOY'*y[_n-1] + `BQ'*q + ///
        `regime_effect' + e if _n > 1

    * Expose exactly T observations per unit, with NO extra presample row.  This
    * matches the Gong--Seo sample geometry: the observed panel is y_i1..y_iT,
    * so T counts observations rather than estimating equations.  Their uncapped
    * instrument set would contain 3+5+7+9=24 columns at T=6.  Study A instead
    * applies the command's fixed maxlag(1 3) cap to every cell; N_iv is therefore
    * an implemented-configuration diagnostic and need not equal 24 (or match
    * between FD and FOD).
    * Retaining t=TBURN as a y_i0 row would hand history(panel) an extra lag and
    * buy one more equation per unit than the published design has -- a gratuitous
    * departure from the benchmark we are being compared against.  Kept identical
    * to study_b_worker.do so the two studies share one sample geometry.
    quietly gen byte analysis = (t > `TBURN')
    quietly drop if t <= `TBURN'
    quietly xtset id t

    *----- balanced, interior MCAR gaps, or monotone attrition -------------
    if "`MISS'" == "mcar" {
        quietly set seed `missing_seed'
        bysort id (t): egen double _tmin = min(cond(analysis, t, .))
        bysort id (t): egen double _tmax = max(cond(analysis, t, .))
        quietly gen double _du = runiform() if analysis
        quietly drop if analysis & _du < `MISSP' & ///
            t != _tmin & t != _tmax
        quietly drop _du _tmin _tmax
        bysort id: egen int _na = total(analysis)
        quietly drop if _na < 4
        quietly drop _na
        quietly xtset id t
    }
    else if "`MISS'" == "attrit" {
        * Independent monotone dropout.  The first eligible hazard is the
        * fifth retained observation; dropping that row onward leaves four
        * observations and therefore never selects units by fit success.
        quietly set seed `missing_seed'
        bysort id (t): gen int _arank = sum(analysis)
        quietly gen double _du = runiform() if analysis & _arank >= 5
        bysort id: egen double _drop_t = min(cond(_du < `MISSP', t, .))
        quietly drop if analysis & !missing(_drop_t) & t >= _drop_t
        quietly drop _arank _du _drop_t
        quietly xtset id t
    }

    *----- realised panel diagnostics before estimation -------------------
    local analysis_potential = `N'*`T'
    quietly count if analysis
    local analysis_observed = r(N)
    local missing_n = `analysis_potential' - `analysis_observed'
    local missing_rate = `missing_n'/`analysis_potential'

    tempvar _utag _gap _fdpot _complete _lastcomplete _fodpot
    quietly egen byte `_utag' = tag(id) if analysis
    quietly count if `_utag' == 1
    local units_realized = r(N)
    local units_dropped = `N' - `units_realized'

    bysort id (t): gen int `_gap' = max(t-t[_n-1]-1, 0) ///
        if analysis & analysis[_n-1]
    quietly count if `_gap' > 0 & !missing(`_gap')
    local gap_events = r(N)
    quietly summarize `_gap' if analysis, meanonly
    local gap_periods = cond(r(N), r(sum), 0)

    bysort id (t): gen byte `_fdpot' = analysis & _n > 2 & ///
        t[_n-1] == t-1 & t[_n-2] == t-2
    quietly count if `_fdpot' == 1
    local fd_pair_rows_potential = r(N)

    bysort id (t): gen byte `_complete' = analysis & _n > 1 & ///
        t[_n-1] == t-1
    bysort id: egen int `_lastcomplete' = max(cond(`_complete', t, .))
    quietly gen byte `_fodpot' = `_complete' & t < `_lastcomplete'
    quietly count if `_fodpot' == 1
    local fod_rows_potential = r(N)

    *----- point estimation only: noboot is unconditional -----------------
    local qopt = cond("`QSTATUS'" == "endog", "endogenous(q)", "predetermined(q)")
    capture quietly xtdpthresh y if analysis, qx(q) `qopt' ///
        method(`METHOD') maxlag(1 3) grid(`GRID') gridtype(uniform) ///
        gridsample(effective) refine(`REFINE') trim(0.15) ///
        noboot coefboot(none) history(panel) vce(robust) nowarn
    local rc = _rc
    timer off 8
    quietly timer list 8
    local el = r(t8)

    *----- collect point estimates and loss inputs -------------------------
    foreach z in gh gbias gsq hp a1 a2 nu nstk ntr niv nun ///
                 est2 fkink rcomp rit radd rrem rneigh pedge ///
                 obj stol sconv sinc shit sw2 s1lev s2lev s1n s2n ///
                 s1same s2same s1gain s2gain gmax ///
                 qlo qhi g1lo g1hi g2lo g2hi ///
                 greq geff gadm gstruct g2adm minreg ///
                 b_Lagyb b_qb b_consd b_qd b_Lagyd ///
                 err_Lagyb err_qb err_consd err_qd err_Lagyd ///
                 sqerr_Lagyb sqerr_qb sqerr_consd sqerr_qd sqerr_Lagyd ///
                 se_Lagyb se_qb se_consd se_qd se_Lagyd {
        local `z' = .
    }
    local code_actual NA
    local smode NA
    local fit_rc0 = (`rc' == 0)
    local cmdok = 0
    local vok = 0
    local contract = 0
    local succ = 0

    if `rc' == 0 {
        local code_actual "`e(cmdversion)'"
        local cmdok = ("`e(cmd)'" == "xtdpthresh")
        local vok = ("`e(cmdversion)'" == "`CODE_VERSION_EXPECTED'")
        local contract = (`cmdok' & `vok')

        capture local hp      = e(hansen_p)
        capture local a1      = e(ar1_p)
        capture local a2      = e(ar2_p)
        capture local nu      = e(N)
        capture local nstk    = e(N_stack)
        capture local ntr     = e(N_trans)
        capture local niv     = e(N_iv)
        capture local nun     = e(N_units)
        capture local est2    = e(estimator_twostep)
        capture local fkink   = e(flag_kink)
        capture local rcomp   = e(refine_complete)
        capture local rit     = e(refine_iterations)
        capture local radd    = e(refine_added)
        capture local rrem    = e(refine_remaining)
        capture local rneigh  = e(refine_neigh_unevaluated)
        capture local obj     = e(obj)
        capture local smode   "`e(searchmode)'"
        capture local stol    = e(searchtol)
        capture local sconv   = e(search_converged)
        capture local sinc    = e(search_incomplete)
        capture local shit    = e(search_hit_max)
        capture local sw2     = e(search_W2_builds)
        capture local s1lev   = e(search_stage1_level)
        capture local s2lev   = e(search_stage2_level)
        capture local s1n     = e(search_stage1_points)
        capture local s2n     = e(search_stage2_points)
        capture local s1same  = e(search_stage1_same_split)
        capture local s2same  = e(search_stage2_same_split)
        capture local s1gain  = e(search_stage1_rel_gain)
        capture local s2gain  = e(search_stage2_rel_gain)
        capture local gmax    = e(grid_max_requested)
        capture local qlo     = e(q_lo)
        capture local qhi     = e(q_hi)
        capture local g1lo    = e(gamma_grid1_lo)
        capture local g1hi    = e(gamma_grid1_hi)
        capture local g2lo    = e(gamma_grid2_lo)
        capture local g2hi    = e(gamma_grid2_hi)
        capture local greq    = e(grid_requested)
        capture local geff    = e(grid_effective)
        capture local gadm    = e(grid_admitted)
        capture local gstruct = e(grid_structural)
        capture local g2adm   = e(grid_twostep_admitted)
        capture local minreg  = e(minregime_applied)

        capture matrix _b_contract = e(b)
        if _rc local contract = 0
        capture assert !missing(e(N), e(N_stack), e(N_trans), e(N_iv), ///
            e(N_units), e(gamma), e(estimator_twostep), e(flag_kink), ///
            e(q_lo), e(q_hi), e(gamma_grid1_lo), e(gamma_grid1_hi), ///
            e(grid_requested), e(grid_effective), e(grid_admitted), ///
            e(grid_structural), e(refine_requested), e(refine_complete), ///
            e(refine_iterations), e(refine_added), e(refine_remaining), ///
            e(refine_neigh_unevaluated), e(minregime_applied), e(obj), ///
            e(search_hit_max), e(search_W2_builds), ///
            e(search_stage1_level), e(search_stage2_level), ///
            e(search_stage1_points), e(search_stage2_points), ///
            e(grid_max_requested))
        if _rc local contract = 0
        capture assert e(N_stack) == e(N_trans) & e(N_units) >= 1 & ///
            e(N_iv) >= 1 & e(q_lo) < e(q_hi) & ///
            e(gamma_grid1_lo) <= e(gamma_grid1_hi) & ///
            e(grid_requested) == `GRID' & ///
            e(grid_max_requested) == `GRID' & ///
            e(grid_effective) == `GRID' + e(refine_added) & ///
            e(grid_admitted) >= 2 & ///
            e(grid_admitted) <= e(grid_structural) & ///
            e(grid_structural) <= e(grid_effective) & ///
            e(minregime_applied) >= 1 & e(flag_kink) == 0 & ///
            inlist(e(estimator_twostep),0,1)
        if _rc local contract = 0
        if "`e(method)'" != "`METHOD'" | "`e(history)'" != "panel" | ///
                "`e(searchmode)'" != "fixed" | ///
                "`e(gridtype)'" != "uniform" | ///
                "`e(gridsample)'" != "effective" | ///
                "`e(vce_requested)'" != "robust" {
            local contract = 0
        }
        capture assert e(flag_static) == 0 & abs(e(trim)-.15) <= 1e-12
        if _rc local contract = 0
        capture assert e(refine_requested) == `REFINE' & ///
            inrange(e(refine_iterations),0,`REFINE') & ///
            e(refine_added) >= 0 & e(refine_remaining) >= 0 & ///
            e(refine_neigh_unevaluated) >= 0 & ///
            inlist(e(refine_complete),0,1) & ///
            (e(refine_complete) == 0 | ///
                (e(refine_remaining) == 0 & ///
                 e(refine_neigh_unevaluated) == 0)) & ///
            missing(e(searchtol)) & missing(e(search_converged)) & ///
            missing(e(search_incomplete)) & e(search_hit_max) == 0 & ///
            e(search_W2_builds) == 1 & ///
            e(search_stage1_level) == 1 & ///
            inlist(e(search_stage2_level),0,1) & ///
            e(search_stage1_points) == `GRID' & ///
            e(search_stage2_points) == ///
                cond(e(search_stage2_level)==1,`GRID',0) & ///
            missing(e(search_stage1_same_split)) & ///
            missing(e(search_stage2_same_split)) & ///
            missing(e(search_stage1_rel_gain)) & ///
            missing(e(search_stage2_rel_gain)) & ///
            missing(e(search_level2_points)) & ///
            missing(e(search_level3_points))
        if _rc local contract = 0
        if `est2' == 1 {
            capture assert !missing(e(gamma_grid2_lo), e(gamma_grid2_hi), ///
                e(grid_twostep_admitted)) & ///
                e(gamma_grid2_lo) <= e(gamma_grid2_hi) & ///
                e(grid_twostep_admitted) >= 2 & ///
                e(grid_twostep_admitted) <= e(grid_effective)
            if _rc local contract = 0
        }
        else if `est2' == 0 {
            capture assert ///
                (missing(e(grid_twostep_admitted)) & ///
                    missing(e(gamma_grid2_lo), e(gamma_grid2_hi))) | ///
                (!missing(e(grid_twostep_admitted)) & ///
                    inrange(e(grid_twostep_admitted),0,e(grid_effective)) & ///
                    ((e(grid_twostep_admitted) == 0 & ///
                        missing(e(gamma_grid2_lo),e(gamma_grid2_hi))) | ///
                     (e(grid_twostep_admitted) > 0 & ///
                        !missing(e(gamma_grid2_lo),e(gamma_grid2_hi)) & ///
                        e(gamma_grid2_lo) <= e(gamma_grid2_hi))))
            if _rc local contract = 0
        }
        capture assert e(boot_threshold_requested) == 0 & ///
            e(boot_linearity_requested) == 0 & ///
            e(boot_continuity_requested) == 0 & ///
            e(boot_coef_requested) == 0
        if _rc local contract = 0
        foreach lab in Lag_y_b q_b cons_d q_d Lag_y_d {
            capture scalar __contract = _b[`lab']
            if _rc local contract = 0
            else if missing(__contract) local contract = 0
            * Point delivery is independent of covariance availability.
        }

        local gh = e(gamma)
        if !missing(`gh') {
            local gbias = `gh' - `GAMMA0'
            local gsq = (`gh' - `GAMMA0')^2
        }
        foreach nm in Lagyb qb consd qd Lagyd {
            local lab = cond("`nm'" == "Lagyb", "Lag_y_b", ///
                        cond("`nm'" == "qb", "q_b", ///
                        cond("`nm'" == "consd", "cons_d", ///
                        cond("`nm'" == "qd", "q_d", "Lag_y_d"))))
            capture local b_`nm' = _b[`lab']
            capture local se_`nm' = _se[`lab']
            if !missing(`b_`nm'') {
                local err_`nm' = `b_`nm'' - `T_`nm''
                local sqerr_`nm' = (`b_`nm'' - `T_`nm'')^2
            }
        }
        if `contract' {
            local _blo = cond(`est2' == 1, `g2lo', `g1lo')
            local _bhi = cond(`est2' == 1, `g2hi', `g1hi')
            local _btol = 1e-10*max(1,abs(`_blo'),abs(`_bhi'))
            capture assert `gh' >= `_blo'-`_btol' & ///
                `gh' <= `_bhi'+`_btol'
            if _rc local contract = 0
            else local pedge = (abs(`gh'-`_blo') <= `_btol') + ///
                2*(abs(`gh'-`_bhi') <= `_btol')
        }
        local succ = (`contract' & `est2' == 1)
    }


    * 0.9.37 telemetry: missing SE is not a point-estimation failure.
    foreach zz in joint_vce ar_joint vce_applied bwscale gamma_bw q_nvals_bw N_iv_dep N_iv_dep_near iv_dep_res ar1_cond ar2_cond ar1_p_cond ar2_p_cond se_cond_Lagyb se_cond_qb se_cond_consd se_cond_qd se_cond_Lagyd se_delivered {
        local `zz' = .
    }
    local se_delivered = 0
    if `rc' == 0 {
        foreach zz in joint_vce ar_joint vce_applied bwscale gamma_bw q_nvals_bw N_iv_dep N_iv_dep_near iv_dep_res ar1_cond ar2_cond {
            capture local `zz' = e(`zz')
        }
        local ar1_p_cond = 2*normal(-abs(`ar1_cond'))
        local ar2_p_cond = 2*normal(-abs(`ar2_cond'))
        local se_delivered = 1
        capture matrix _mc_vcond = e(V_cond)
        local _mc_has_vcond = (_rc == 0)
        foreach nm in Lagyb qb consd qd Lagyd {
            local lab = cond("`nm'"=="Lagyb","Lag_y_b", ///
                cond("`nm'"=="qb","q_b",cond("`nm'"=="consd","cons_d", ///
                cond("`nm'"=="qd","q_d","Lag_y_d"))))
            if missing(`se_`nm'') | `se_`nm'' < 0 local se_delivered = 0
            if `_mc_has_vcond' {
                local cc = colnumb(_mc_vcond, "`lab'")
                if !missing(`cc') {
                    if !missing(_mc_vcond[`cc',`cc']) & _mc_vcond[`cc',`cc'] >= 0 ///
                        local se_cond_`nm' = sqrt(_mc_vcond[`cc',`cc'])
                }
            }
        }
    }

    tempname FH
    file open `FH' using "`OUT'", write append text
    file write `FH' ///
        "`RUN_ID',`HARNESS_VERSION',`SCHEMA_VERSION',`CODE_VERSION_EXPECTED',`code_actual',`STATA_VERSION',`RNG_KIND',`RNG_ACTUAL',`SEED_SCHEME',`DGP_SPEC',`INFERENCE_MODE'," ///
        "`CELL_ID',`PAIR_ID',`CONFIG_KEY',`FULL_KEY',`STUDY',`DESIGN',`METHOD',`N',`T',`MISS',`MISSP',`QSTATUS',`KAPPA',`R',`GRID',`REFINE',`MASTER',`SHARD',`NSHARD',`rep'," ///
        "`dgp_seed',`missing_seed',`rc',`fit_rc0',`cmdok',`vok',`contract',`succ'," ///
        "`GAMMA0',`gh',`gbias',`gsq',`obj'," ///
        "`T_Lagyb',`T_qb',`T_consd',`T_qd',`T_Lagyd'," ///
        "`b_Lagyb',`b_qb',`b_consd',`b_qd',`b_Lagyd'," ///
        "`err_Lagyb',`sqerr_Lagyb',`err_qb',`sqerr_qb',`err_consd',`sqerr_consd',`err_qd',`sqerr_qd',`err_Lagyd',`sqerr_Lagyd'," ///
        "`se_Lagyb',`se_qb',`se_consd',`se_qd',`se_Lagyd'," ///
        "`hp',`a1',`a2',`nu',`nstk',`ntr',`niv',`nun',`est2',`fkink',`rcomp',`rit',`radd',`rrem',`rneigh'," ///
        "`smode',`stol',`sconv',`sinc',`shit',`sw2',`s1lev',`s2lev',`s1n',`s2n',`s1same',`s2same',`s1gain',`s2gain',`gmax',`pedge'," ///
        "`qlo',`qhi',`g1lo',`g1hi',`g2lo',`g2hi',`greq',`geff',`gadm',`gstruct',`g2adm',`minreg'," ///
        "`N',`units_realized',`units_dropped',`analysis_potential',`analysis_observed',`missing_n',`missing_rate'," ///
        "`gap_events',`gap_periods',`fd_pair_rows_potential',`fod_rows_potential'," ///
        "`RHOY',`BQ',`T_consd',`T_qd',`T_Lagyd',`RHOQ',`RHOEU',`SIGE',`TBURN'," ///
        "1,3,0.15,panel,uniform,effective,robust,`el',`joint_vce',`ar_joint',`vce_applied',`bwscale',`gamma_bw',`q_nvals_bw',`N_iv_dep',`N_iv_dep_near',`iv_dep_res',`ar1_cond',`ar2_cond',`ar1_p_cond',`ar2_p_cond',`se_cond_Lagyb',`se_cond_qb',`se_cond_consd',`se_cond_qd',`se_cond_Lagyd',`se_delivered'" _n
    file close `FH'

    if mod(`rep', 25) == 0 {
        di as txt "  [`CELL_ID'] rep `rep': rc=`rc', gamma=`gh', 2step=`est2', fixed-grid=`s1n'/`s2n', refine=`rit', success=`succ'"
    }
}

di as res "WORKER DONE: `CELL_ID' reps `REP_START'-`REP_END' shard `SHARD'"
