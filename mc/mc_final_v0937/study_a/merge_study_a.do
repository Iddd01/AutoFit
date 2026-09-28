*!==========================================================================
*! merge_study_a.do -- fail-closed merger for Study A (point estimates only)
*!--------------------------------------------------------------------------
*! Harness: xtdpthresh_study_a_v0937_r400
*! Worker schema: study_a_worker_v0937
*! Registry: 116 targeted cells; 48,800 formal replication rows.
*!--------------------------------------------------------------------------
*! No bootstrap estimand is constructed here.  Analytic-SE calibration and
*! diagnostic rejection rates are descriptive point-estimator diagnostics.
*! Outputs are published only after the manifest, all shards, completion
*! markers, the exact registry, row keys, CRN pairs, and result schema pass.
*!==========================================================================
version 15.0
clear all
set more off
set varabbrev off

args RUN_ID MERGE_NONCE
local HARNESS_EXPECTED "xtdpthresh_study_a_v0937_r400"
local SCHEMA_EXPECTED  "study_a_worker_v0937"
local CODE_EXPECTED    "0.9.37"
local RNG_EXPECTED     "mt64"
local SEEDS_EXPECTED   "study_a_crn_v3"

* The PowerShell verifier recomputes every frozen SHA-256 before granting a
* one-use merge authorization.  Direct invocation of this .do file is denied.
if !regexm("`RUN_ID'", "^[A-Za-z0-9][A-Za-z0-9_.-]*$") | ///
   !regexm("`RUN_ID'", "[A-Za-z]") | ///
   !regexm("`MERGE_NONCE'", "^[0-9A-Fa-f]+$") | ///
   strlen("`MERGE_NONCE'") != 32 {
    di as err "merge_study_a: valid RUN_ID and verifier nonce are required"
    exit 198
}
capture confirm file "_merge_authorized.txt"
if _rc {
    di as err "merge_study_a: verifier authorization is required"
    exit 601
}
tempname AUTH
file open `AUTH' using "_merge_authorized.txt", read text
file read `AUTH' AUTH_RUN
file read `AUTH' AUTH_HASH
file read `AUTH' AUTH_NONCE
file read `AUTH' AUTH_EXTRA
file close `AUTH'
local AUTH_RUN = subinstr(`"`AUTH_RUN'"', char(13), "", .)
local AUTH_HASH = subinstr(`"`AUTH_HASH'"', char(13), "", .)
local AUTH_NONCE = subinstr(`"`AUTH_NONCE'"', char(13), "", .)
local AUTH_EXTRA = subinstr(`"`AUTH_EXTRA'"', char(13), "", .)
if `"`AUTH_RUN'"' != "run_id=`RUN_ID'" | ///
   !regexm(`"`AUTH_HASH'"', "^manifest_sha256=[0-9A-Fa-f]+$") | ///
   strlen(`"`AUTH_HASH'"') != 80 | ///
   `"`AUTH_NONCE'"' != "nonce=`MERGE_NONCE'" | `"`AUTH_EXTRA'"' != "" {
    di as err "merge_study_a: malformed or mismatched verifier authorization"
    exit 459
}
local MANIFEST_SHA = substr(`"`AUTH_HASH'"', 17, 64)
capture erase "_merge_authorized.txt"
if _rc {
    di as err "merge_study_a: cannot consume verifier authorization"
    exit 603
}

* An optional RUN_ID lets this file be launched from the source directory.
if "`RUN_ID'" != "" {
    local here : dir "." files "study_a_results_SH*.csv"
    if `: word count `here'' == 0 {
        capture cd "runs/`RUN_ID'"
        if _rc {
            di as err "merge_study_a: cannot find runs/`RUN_ID'"
            exit 601
        }
    }
}

* A prior attestation must never survive a new merge attempt.  The marker is
* recreated only after every validated output has been copied into place.
local MERGE_MARKER "study_a_merge_complete.csv"
capture confirm file "`MERGE_MARKER'"
if !_rc {
    capture erase "`MERGE_MARKER'"
    if _rc {
        di as err "merge_study_a: cannot remove stale merge attestation"
        exit 603
    }
}

foreach f in study_a_cells.do study_a_worker.do {
    capture confirm file "`f'"
    if _rc {
        di as err "merge_study_a: required frozen file is missing: `f'"
        exit 601
    }
}

* -------------------------------------------------------------------------
* Immutable launcher manifest.
* REP_CAP is read only from this contract; it is never inferred from output.
* -------------------------------------------------------------------------
capture confirm file "manifest.json"
if _rc {
    di as err "merge_study_a: manifest.json is required"
    exit 601
}
local MJ = fileread("manifest.json")

if !regexm(`"`MJ'"', `""schemaVersion"[ ]*:[ ]*([0-9]+)"') {
    di as err "merge_study_a: manifest lacks schemaVersion"
    exit 459
}
local M_Schema = real(regexs(1))
foreach q in NShard Master RepCap ExpectedCells ExpectedReplications {
    if !regexm(`"`MJ'"', `""`q'"[ ]*:[ ]*([0-9]+)"') {
        di as err "merge_study_a: manifest lacks numeric config `q'"
        exit 459
    }
    local M_`q' = real(regexs(1))
}
foreach q in RunId HarnessVersion CodeVersionExpected RngKind {
    if !regexm(`"`MJ'"', `""`q'"[ ]*:[ ]*"([^"]+)""') {
        di as err "merge_study_a: manifest lacks string config `q'"
        exit 459
    }
    local M_`q' "`=regexs(1)'"
}
if !regexm("`M_RunId'", "[A-Za-z]") | ///
   !regexm("`M_RunId'", "^[A-Za-z0-9][A-Za-z0-9_.-]*$") {
    di as err "merge_study_a: manifest RunId is unsafe"
    exit 459
}
if `M_Schema' != 1 | `M_NShard' < 1 | `M_Master' < 1 | ///
   `M_RepCap' < 0 | `M_RepCap' > 500 | ///
   (`M_RepCap' > 0 & `M_RepCap' < `M_NShard' & `M_RepCap' != 1) | ///
   (`M_RepCap' == 1 & `M_NShard' > 116) | ///
   `M_ExpectedCells' != 116 | ///
   "`M_HarnessVersion'" != "`HARNESS_EXPECTED'" | ///
   "`M_CodeVersionExpected'" != "`CODE_EXPECTED'" | ///
   "`M_RngKind'" != "`RNG_EXPECTED'" {
    di as err "merge_study_a: unsupported/inconsistent manifest contract"
    exit 459
}
if "`RUN_ID'" != "" & "`RUN_ID'" != "`M_RunId'" {
    di as err "merge_study_a: requested RUN_ID differs from manifest"
    exit 459
}

* Presence and syntax of the core frozen-source SHA-256 declarations.
* run_study_a.ps1 is responsible for recomputing these hashes before launch.
foreach f in study_a_cells.do _registry_probe.do ///
             study_a_worker.do study_a_shard.do ///
             merge_study_a.do xtdpthresh.ado xtdpthresh_p.ado ///
             run_study_a.ps1 monitor_study_a.ps1 ///
             verify_and_merge_study_a.ps1 README.md _DESIGN.md {
    capture confirm file "`f'"
    if _rc {
        di as err "merge_study_a: frozen manifest member is missing: `f'"
        exit 601
    }
    local fp = subinstr("`f'", ".", "[.]", .)
    if !regexm(`"`MJ'"', `""`fp'"[ ]*:[ ]*"([0-9A-Fa-f]+)""') {
        di as err "merge_study_a: manifest lacks SHA-256 for `f'"
        exit 459
    }
    local hf "`=regexs(1)'"
    if strlen("`hf'") != 64 {
        di as err "merge_study_a: malformed SHA-256 for `f'"
        exit 459
    }
}

* Two replication groups: 24 hard stresses @ R=500 and 92 regular cells
* @ R=400, including the 20 N=1600 cells.
local R500 = cond(`M_RepCap' > 0, min(500, `M_RepCap'), 500)
local R400 = cond(`M_RepCap' > 0, min(400, `M_RepCap'), 400)
local G_EFF = 199
local TOTAL_EFF = 24*`R500' + 92*`R400'
if `M_ExpectedReplications' != `TOTAL_EFF' {
    di as err "merge_study_a: manifest ExpectedReplications is not the v6 registry total"
    exit 459
}

* -------------------------------------------------------------------------
* Import every shard with an exact, ordered schema.  No extra columns are
* accepted, so an inference-mode or legacy file cannot enter this study.
* -------------------------------------------------------------------------
local EXPECT_VARS run_id harness_version schema_version ///
    code_version_expected code_version_actual stata_version rng_kind ///
    rng_actual seed_scheme dgp_spec inference_mode cell_id pair_id ///
    config_key full_key study design method N T miss missp qstatus kappa ///
    R grid refine master shard nshard rep dgp_seed missing_seed rc fit_rc0 ///
    cmd_ok version_ok contract_ok success gamma0 gamma_hat g_bias g_sqerr obj ///
    Lagyb qb consd qd Lagyd b_Lagyb b_qb b_consd b_qd b_Lagyd ///
    err_Lagyb sqerr_Lagyb err_qb sqerr_qb err_consd sqerr_consd ///
    err_qd sqerr_qd err_Lagyd sqerr_Lagyd se_Lagyb se_qb se_consd se_qd ///
    se_Lagyd hansen_p ar1_p ar2_p N_used N_stack N_trans N_iv N_units ///
    estimator_twostep flag_kink refine_complete refine_iterations ///
    refine_added ///
    refine_remaining refine_neigh_unevaluated searchmode searchtol ///
    search_converged search_incomplete search_hit_max search_W2_builds ///
    search_stage1_level search_stage2_level search_stage1_points ///
    search_stage2_points search_stage1_same_split ///
    search_stage2_same_split search_stage1_rel_gain ///
    search_stage2_rel_gain grid_max_requested point_boundary q_lo q_hi ///
    gamma_grid1_lo gamma_grid1_hi gamma_grid2_lo gamma_grid2_hi ///
    grid_requested grid_effective grid_admitted grid_structural ///
    grid_twostep_admitted minregime_applied units_target ///
    units_realized units_dropped analysis_potential analysis_observed ///
    missing_n missing_rate gap_events gap_periods fd_pair_rows_potential ///
    fod_rows_potential rhoy_true bq_true consd_true qd_true lagyd_true ///
    rhoq rhoeu_effective sige tburn maxlag_lo maxlag_hi trim history ///
    gridtype gridsample vce elapsed_s ///
    joint_vce ar_joint vce_applied bwscale gamma_bw q_nvals_bw N_iv_dep N_iv_dep_near iv_dep_res ar1_cond ar2_cond ar1_p_cond ar2_p_cond se_cond_Lagyb se_cond_qb se_cond_consd se_cond_qd se_cond_Lagyd se_delivered
local EXPECT_VARS : list retokenize EXPECT_VARS

local files : dir "." files "study_a_results_SH*.csv"
local nfiles : word count `files'
if `nfiles' == 0 {
    di as err "merge_study_a: no result shards found"
    exit 601
}
tempfile pooled one
local first 1
foreach f of local files {
    local fl = lower("`f'")
    if !regexm("`fl'", "^study_a_results_sh([0-9]+)[.]csv$") {
        di as err "merge_study_a: unexpected result filename `f'"
        exit 459
    }
    local fsh = real(regexs(1))
    * The worker writes decimal truths and estimator metadata to CSV as
    * doubles.  Stata's default delimited import may demote short decimals
    * such as .7 and .15 to float, which can falsely fail the strict frozen-
    * design assertions below.  Preserve the posted numeric precision.
    capture noisily import delimited using "`f'", clear varnames(1) ///
        case(preserve) bindquote(strict) asdouble
    if _rc | _N == 0 {
        di as err "merge_study_a: unreadable or empty result shard `f'"
        exit 459
    }
    unab GOT_VARS : _all
    if `"`GOT_VARS'"' != `"`EXPECT_VARS'"' {
        di as err "merge_study_a: exact schema/order mismatch in `f'"
        exit 459
    }
    foreach v in run_id harness_version schema_version ///
        code_version_expected code_version_actual rng_kind rng_actual ///
        seed_scheme dgp_spec inference_mode cell_id pair_id config_key ///
        full_key study design method miss qstatus history gridtype ///
        gridsample vce searchmode {
        capture confirm string variable `v'
        if _rc {
            di as err "merge_study_a: wrong string type for `v' in `f'"
            exit 459
        }
    }
    foreach v in stata_version N T missp kappa R grid refine master ///
        shard nshard rep ///
        dgp_seed missing_seed rc fit_rc0 cmd_ok version_ok contract_ok ///
        success gamma0 gamma_hat g_bias g_sqerr obj Lagyb qb consd qd Lagyd ///
        b_Lagyb b_qb b_consd b_qd b_Lagyd err_Lagyb sqerr_Lagyb ///
        err_qb sqerr_qb err_consd sqerr_consd err_qd sqerr_qd ///
        err_Lagyd sqerr_Lagyd se_Lagyb se_qb se_consd se_qd se_Lagyd ///
        hansen_p ar1_p ar2_p N_used N_stack N_trans N_iv N_units ///
        estimator_twostep flag_kink refine_complete refine_iterations ///
        refine_added refine_remaining refine_neigh_unevaluated searchtol ///
        search_converged search_incomplete search_hit_max search_W2_builds ///
        search_stage1_level search_stage2_level search_stage1_points ///
        search_stage2_points search_stage1_same_split ///
        search_stage2_same_split search_stage1_rel_gain ///
        search_stage2_rel_gain grid_max_requested point_boundary q_lo q_hi ///
        gamma_grid1_lo gamma_grid1_hi gamma_grid2_lo gamma_grid2_hi ///
        grid_requested grid_effective grid_admitted grid_structural ///
        grid_twostep_admitted minregime_applied units_target ///
        units_realized units_dropped analysis_potential analysis_observed ///
        missing_n missing_rate gap_events gap_periods ///
        fd_pair_rows_potential fod_rows_potential rhoy_true bq_true ///
        consd_true qd_true lagyd_true rhoq rhoeu_effective sige tburn ///
        maxlag_lo maxlag_hi trim elapsed_s joint_vce ar_joint vce_applied bwscale gamma_bw q_nvals_bw N_iv_dep N_iv_dep_near iv_dep_res ar1_cond ar2_cond ar1_p_cond ar2_p_cond se_cond_Lagyb se_cond_qb se_cond_consd se_cond_qd se_cond_Lagyd se_delivered {
        capture confirm numeric variable `v'
        if _rc {
            di as err "merge_study_a: wrong numeric type for `v' in `f'"
            exit 459
        }
    }
    quietly count if shard != `fsh'
    if r(N) {
        di as err "merge_study_a: filename/data shard mismatch in `f'"
        exit 459
    }
    quietly save `one', replace
    if `first' {
        quietly save `pooled', replace
        local first 0
    }
    else {
        quietly append using `pooled'
        quietly save `pooled', replace
    }
}

use `pooled', clear
foreach v in run_id harness_version schema_version code_version_expected ///
                 stata_version rng_kind rng_actual seed_scheme master ///
                 nshard grid refine inference_mode {
    capture assert `v' == `v'[1]
    if _rc {
        di as err "merge_study_a: mixed run provenance/configuration in `v'"
        exit 459
    }
}
if run_id[1] == "" | harness_version[1] != "`HARNESS_EXPECTED'" | ///
   schema_version[1] != "`SCHEMA_EXPECTED'" | ///
   code_version_expected[1] != "`CODE_EXPECTED'" | ///
   rng_kind[1] != "`RNG_EXPECTED'" | rng_actual[1] != "`RNG_EXPECTED'" | ///
   seed_scheme[1] != "`SEEDS_EXPECTED'" | inference_mode[1] != "point" {
    di as err "merge_study_a: wrong/blank harness, schema, RNG, seed scheme, or mode"
    exit 459
}
if "`RUN_ID'" != "" & run_id[1] != "`RUN_ID'" {
    di as err "merge_study_a: requested RUN_ID differs from data"
    exit 459
}
local ACT_RUN "`=run_id[1]'"
local ACT_HARNESS "`=harness_version[1]'"
local ACT_CODE "`=code_version_expected[1]'"
local ACT_MASTER = master[1]
local ACT_NSHARD = nshard[1]
local ACT_GRID = grid[1]
local ACT_ROWS = _N
if "`ACT_RUN'" != "`M_RunId'" | `ACT_NSHARD' != `M_NShard' | ///
   `ACT_MASTER' != `M_Master' | `ACT_GRID' != `G_EFF' | ///
   "`ACT_HARNESS'" != "`M_HarnessVersion'" | ///
   "`ACT_CODE'" != "`M_CodeVersionExpected'" {
    di as err "merge_study_a: results differ from immutable manifest config"
    exit 459
}
capture assert `ACT_NSHARD' >= 1 & shard >= 1 & shard <= `ACT_NSHARD'
if _rc | `nfiles' != `ACT_NSHARD' {
    di as err "merge_study_a: result-file count/range differs from NShard"
    exit 459
}
forvalues k = 1/`ACT_NSHARD' {
    local exact : dir "." files "study_a_results_SH`k'.csv"
    if `: word count `exact'' != 1 {
        di as err "merge_study_a: missing exact result shard `k'"
        exit 601
    }
}

* Primary keys, exact shard allocation, and point-mode row contract.
capture isid full_key
if _rc {
    di as err "merge_study_a: duplicate full_key"
    exit 459
}
capture isid cell_id rep
if _rc {
    di as err "merge_study_a: duplicate cell_id/rep"
    exit 459
}
gen str244 __key = run_id + "|" + cell_id + "|rep=" + ///
    trim(string(rep, "%12.0f"))
quietly count if full_key != __key | cell_id != pair_id + "_" + method
if r(N) {
    di as err "merge_study_a: malformed full_key or cell_id"
    exit 459
}
drop __key
* Use the worker's allocation inequalities directly.  The algebraically
* equivalent ceil(rep*nshard/R) inverse is unsafe at exact shard boundaries
* because Stata can round an integer quotient infinitesimally upward.
capture assert rep >= 1 & rep <= R & R == floor(R) & ///
    ((`M_RepCap' == 1 & nshard > 1) | ///
     (rep*nshard > (shard-1)*R & rep*nshard <= shard*R)) & ///
    ((inlist(study, "PERSIST", "ENDOG", "HEAVYTAIL") & R == `R500') | ///
     (!inlist(study, "PERSIST", "ENDOG", "HEAVYTAIL") & R == `R400')) & ///
    grid == `G_EFF' & refine == 4 & ///
    inlist(T, 6, 10) & inference_mode == "point"
if _rc {
    di as err "merge_study_a: row outside formal/capped allocation"
    exit 459
}
sort cell_id rep
capture by cell_id: assert rep == _n & _N == R
if _rc {
    di as err "merge_study_a: incomplete replication sequence"
    exit 459
}
foreach v in fit_rc0 cmd_ok contract_ok success {
    capture assert inlist(`v', 0, 1)
    if _rc {
        di as err "merge_study_a: nonbinary flag `v'"
        exit 459
    }
}
capture assert inlist(version_ok, 0, 1, .) & ///
    (!missing(version_ok) | fit_rc0 == 0) & ///
    success == (contract_ok == 1 & estimator_twostep == 1) & ///
    contract_ok <= fit_rc0 & ///
    contract_ok <= cmd_ok & contract_ok <= version_ok & ///
    fit_rc0 == (rc == 0)
if _rc {
    di as err "merge_study_a: inconsistent fit/contract/success accounting"
    exit 459
}
quietly count if rc == 0 & code_version_actual != code_version_expected
if r(N) {
    di as err "merge_study_a: a fitted row used the wrong ado version"
    exit 459
}
capture assert dgp_seed == floor(dgp_seed) & dgp_seed > 0 & ///
    missing_seed == floor(missing_seed) & missing_seed > 0
if _rc {
    di as err "merge_study_a: invalid DGP/missingness seed"
    exit 459
}
tempvar __spcode __qcode __mcode __h_expected __dseed_expected ///
    __mseed_expected
gen byte `__spcode' = cond(dgp_spec == "gs26_base_v2", 1, ///
    cond(dgp_spec == "gs26_persist_v2", 2, ///
    cond(dgp_spec == "gs26_endog_v2", 3, 4)))
gen byte `__qcode' = cond(qstatus == "pred", 1, 2)
gen byte `__mcode' = cond(miss == "balanced", 1, ///
    cond(miss == "mcar", 2, 3))
gen double `__h_expected' = mod(master*1000003 + rep*10007 + ///
    N*101 + T*1009 + `__spcode'*1000033 + `__qcode'*1000037, ///
    699999999)
gen double `__dseed_expected' = 1 + `__h_expected'
gen double `__mseed_expected' = 700000001 + ///
    mod(`__h_expected'*1009 + `__mcode'*10007 + 7919, 699999999)
capture assert dgp_seed == `__dseed_expected' & ///
    missing_seed == `__mseed_expected'
if _rc {
    di as err "merge_study_a: row violates the frozen v6 registry CRN seed map"
    exit 459
}
drop `__spcode' `__qcode' `__mcode' `__h_expected' ///
    `__dseed_expected' `__mseed_expected'

* Audit worker-posted truth and errors before using them in any summary.
capture assert inlist(dgp_spec, "gs26_base_v2", "gs26_persist_v2", ///
        "gs26_endog_v2", "gs26_t5_v2") & ///
    inlist(study, "GSCORE", "GSSMALL", "PERSIST", "MISSPAT", ///
        "ENDOG", "HEAVYTAIL") & ///
    inlist(qstatus, "pred", "endog") & inlist(method, "fd", "fod") & ///
    inlist(N, 100, 200, 400, 800, 1600) & inlist(T, 6, 10) & ///
    (abs(kappa) < 1e-6 | abs(kappa-.1) < 1e-6 | ///
     abs(kappa-.2) < 1e-6 | abs(kappa-.5) < 1e-6 | ///
     abs(kappa-1) < 1e-6) & ///
    ((miss == "balanced" & missp == 0) | ///
     (miss == "mcar" & (abs(missp-.15) < 1e-6 | abs(missp-.30) < 1e-6)) | ///
     (miss == "attrit" & abs(missp-.15) < 1e-6)) & ///
    ((abs(kappa) < 1e-6 & design == "kink") | ///
     (kappa > 1e-6 & design == "jump")) & ///
    ((dgp_spec == "gs26_endog_v2" & qstatus == "endog") | ///
     (dgp_spec != "gs26_endog_v2" & qstatus == "pred"))
if _rc {
    di as err "merge_study_a: row lies outside the Study-A design"
    exit 459
}
capture assert abs(Lagyb-rhoy_true) < 1e-10 & ///
    abs(qb-bq_true) < 1e-10 & abs(consd-consd_true) < 1e-10 & ///
    abs(qd-qd_true) < 1e-10 & abs(Lagyd-lagyd_true) < 1e-10
if _rc {
    di as err "merge_study_a: duplicated true coefficients disagree"
    exit 459
}
capture assert abs(gamma0-.25) < 1e-10 & ///
    abs(rhoy_true-cond(dgp_spec == "gs26_persist_v2", .9, .6)) < 1e-10 & ///
    abs(bq_true-1) < 1e-10 & abs(consd_true-(kappa-.5)) < 1e-10 & ///
    abs(qd_true-2) < 1e-10 & abs(lagyd_true) < 1e-10 & ///
    abs(rhoq-.7) < 1e-10 & abs(rhoeu_effective-.5) < 1e-10 & ///
    abs(sige-.5) < 1e-10 & ///
    tburn == cond(dgp_spec == "gs26_persist_v2", 50, 20) & ///
    maxlag_lo == 1 & maxlag_hi == 3 & abs(trim-.15) < 1e-10 & ///
    history == "panel" & gridtype == "uniform" & ///
    gridsample == "effective" & vce == "robust"
if _rc {
    di as err "merge_study_a: DGP truth or estimator configuration is not frozen v6"
    exit 459
}
capture assert !missing(units_target, units_realized, units_dropped, ///
    analysis_potential, analysis_observed, missing_n, missing_rate, ///
    gap_events, gap_periods, elapsed_s) & ///
    units_target == N & units_realized >= 0 & ///
    units_dropped == units_target-units_realized & ///
    analysis_potential >= analysis_observed & ///
    missing_n == analysis_potential-analysis_observed & ///
    analysis_potential > 0 & ///
    abs(missing_rate-missing_n/analysis_potential) < 1e-6 & ///
    gap_events >= 0 & gap_periods >= gap_events
if _rc {
    di as err "merge_study_a: invalid missingness/gap accounting"
    exit 459
}
capture assert missing_n == 0 & missing_rate == 0 & ///
    gap_events == 0 & gap_periods == 0 if miss == "balanced"
if _rc {
    di as err "merge_study_a: balanced rows contain missingness/gaps"
    exit 459
}
capture assert !missing(N_used, N_stack, N_trans, N_iv, N_units, ///
        estimator_twostep, flag_kink, refine_complete, ///
        refine_iterations, refine_added, refine_remaining, ///
        refine_neigh_unevaluated, obj, search_hit_max, search_W2_builds, ///
        search_stage1_level, search_stage2_level, ///
        search_stage1_points, search_stage2_points, ///
        grid_max_requested, point_boundary, q_lo, q_hi, ///
        gamma_grid1_lo, gamma_grid1_hi, grid_requested, ///
        grid_effective, grid_admitted, grid_structural, ///
        minregime_applied) & ///
    inlist(estimator_twostep,0,1) & flag_kink == 0 & ///
    inlist(point_boundary,0,1,2,3) & ///
    N_used >= 1 & N_used <= analysis_observed & ///
    N_stack == N_trans & N_trans >= 1 & N_iv >= 1 & ///
    N_units >= 1 & N_units <= units_realized & ///
    q_lo < q_hi & gamma_grid1_lo <= gamma_grid1_hi & ///
    grid_requested == `G_EFF' & grid_max_requested == grid_requested & ///
    grid_effective == grid_requested+refine_added & ///
    grid_admitted >= 2 & grid_admitted <= grid_structural & ///
    grid_structural <= grid_effective & ///
    minregime_applied >= 1 & refine == 4 & ///
    inrange(refine_iterations,0,refine) & refine_added >= 0 & ///
    refine_remaining >= 0 & refine_neigh_unevaluated >= 0 & ///
    inlist(refine_complete,0,1) & ///
    (refine_complete == 0 | ///
        (refine_remaining == 0 & refine_neigh_unevaluated == 0)) ///
    if contract_ok
if _rc {
    di as err "merge_study_a: contract-valid row has invalid search/sample diagnostics"
    exit 459
}
* Covariance delivery is tracked separately from point success.
capture assert inlist(se_delivered,0,1)
if _rc exit 459
capture assert searchmode == "fixed" & missing(searchtol) & ///
    missing(search_converged) & missing(search_incomplete) & ///
    search_hit_max == 0 & search_W2_builds == 1 & ///
    search_stage1_level == 1 & inlist(search_stage2_level,0,1) & ///
    search_stage1_points == grid_requested & ///
    search_stage2_points == ///
        cond(search_stage2_level==1,grid_requested,0) & ///
    missing(search_stage1_same_split) & ///
    missing(search_stage2_same_split) & ///
    missing(search_stage1_rel_gain) & ///
    missing(search_stage2_rel_gain) ///
    if contract_ok
if _rc {
    di as err "merge_study_a: delivered row has invalid fixed-grid diagnostics"
    exit 459
}
capture assert !missing(gamma_grid2_lo, gamma_grid2_hi, ///
        grid_twostep_admitted) & ///
    gamma_grid2_lo <= gamma_grid2_hi & ///
    grid_twostep_admitted >= 2 & ///
    grid_twostep_admitted <= grid_effective ///
    if contract_ok & estimator_twostep == 1
if _rc {
    di as err "merge_study_a: two-step row has invalid stage-2 search space"
    exit 459
}
capture assert ///
    (missing(grid_twostep_admitted) & ///
        missing(gamma_grid2_lo, gamma_grid2_hi)) | ///
    (!missing(grid_twostep_admitted) & ///
        inrange(grid_twostep_admitted,0,grid_effective) & ///
        ((grid_twostep_admitted == 0 & ///
            missing(gamma_grid2_lo,gamma_grid2_hi)) | ///
         (grid_twostep_admitted > 0 & ///
            !missing(gamma_grid2_lo,gamma_grid2_hi) & ///
            gamma_grid2_lo <= gamma_grid2_hi))) ///
    if contract_ok & estimator_twostep == 0
if _rc {
    di as err "merge_study_a: one-step fallback has incoherent stage-2 diagnostics"
    exit 459
}
tempvar __search_lo __search_hi __edge_tol __edge_expected
gen double `__search_lo' = cond(estimator_twostep == 1, ///
    gamma_grid2_lo, gamma_grid1_lo) if contract_ok
gen double `__search_hi' = cond(estimator_twostep == 1, ///
    gamma_grid2_hi, gamma_grid1_hi) if contract_ok
gen double `__edge_tol' = 1e-10*max(1,abs(`__search_lo'), ///
    abs(`__search_hi')) if contract_ok
gen byte `__edge_expected' = ///
    (abs(gamma_hat-`__search_lo') <= `__edge_tol') + ///
    2*(abs(gamma_hat-`__search_hi') <= `__edge_tol') if contract_ok
capture assert gamma_hat >= `__search_lo'-`__edge_tol' & ///
    gamma_hat <= `__search_hi'+`__edge_tol' & ///
    point_boundary == `__edge_expected' if contract_ok
if _rc {
    di as err "merge_study_a: reported gamma/boundary flag violates its search span"
    exit 459
}
drop `__search_lo' `__search_hi' `__edge_tol' `__edge_expected'
foreach v in gamma0 gamma_hat g_bias g_sqerr obj b_Lagyb b_qb b_consd ///
    b_qd b_Lagyd err_Lagyb sqerr_Lagyb err_qb sqerr_qb err_consd ///
    sqerr_consd err_qd sqerr_qd err_Lagyd sqerr_Lagyd N_used ///
    N_trans N_iv {
    quietly count if contract_ok & missing(`v')
    if r(N) {
        di as err "merge_study_a: contract-valid row lacks `v'"
        exit 459
    }
}
capture assert abs(g_bias-(gamma_hat-gamma0)) < 1e-8 & ///
    abs(g_sqerr-g_bias^2) < 1e-8 if contract_ok
if _rc {
    di as err "merge_study_a: posted gamma error is inconsistent"
    exit 459
}
local P Lagyb qb consd qd Lagyd
foreach p of local P {
    capture assert abs(err_`p'-(b_`p'-`p')) < 1e-8 & ///
        abs(sqerr_`p'-err_`p'^2) < 1e-8 if contract_ok
    if _rc {
        di as err "merge_study_a: posted `p' error is inconsistent"
        exit 459
    }
}
quietly save `pooled', replace

* -------------------------------------------------------------------------
* Dry-evaluate the authoritative registry under the manifest overrides.
* -------------------------------------------------------------------------
tempfile expected
tempname PP
postfile `PP' str80 cell_id str72 pair_id str12 study str12 design ///
    str4 method int N int T str12 miss double missp str8 qstatus long R ///
    double kappa str24 dgp_spec int grid byte refine ///
    str8 inference_mode str244 config_key long cell_order ///
    using `expected', replace
global STA_POST "`PP'"
global STA_CAP `M_RepCap'
global STA_CELL_ORDER 0
capture program drop _cell
program define _cell
    version 15.0
    args PAIR_ID STUDY DESIGN METHOD N T MISS MISSP QSTATUS R KAPPA ///
         DGP_SPEC GRID REFINE
    global STA_CELL_ORDER = ${STA_CELL_ORDER} + 1
    local CELL_ID "`PAIR_ID'_`METHOD'"
    local RR = `R'
    if ${STA_CAP} > 0 & `RR' > ${STA_CAP} local RR = ${STA_CAP}
    local GG = `GRID'
    local CK "`STUDY'|`DESIGN'|`METHOD'|N=`N'|T=`T'|`MISS'=`MISSP'|`QSTATUS'|R=`RR'|G=`GG'|RF=`REFINE'|K=`KAPPA'|`DGP_SPEC'|point"
    post ${STA_POST} ("`CELL_ID'") ("`PAIR_ID'") ("`STUDY'") ///
        ("`DESIGN'") ("`METHOD'") (`N') (`T') ("`MISS'") (`MISSP') ///
        ("`QSTATUS'") (`RR') (`KAPPA') ("`DGP_SPEC'") (`GG') ///
        (`REFINE') ("point") ("`CK'") (${STA_CELL_ORDER})
end
capture noisily do study_a_cells.do
local registry_rc = _rc
postclose `PP'
capture program drop _cell
macro drop STA_POST STA_CAP STA_CELL_ORDER
if `registry_rc' {
    di as err "merge_study_a: registry evaluation failed"
    exit `registry_rc'
}
use `expected', clear
capture isid cell_id
if _rc | _N != 116 {
    di as err "merge_study_a: registry is not exactly 116 unique cells"
    exit 459
}
capture assert ((inlist(study, "PERSIST", "ENDOG", "HEAVYTAIL") & R == `R500') | ///
    (!inlist(study, "PERSIST", "ENDOG", "HEAVYTAIL") & R == `R400')) & ///
    grid == `G_EFF' & refine == 4 & inlist(T, 6, 10) & ///
    inference_mode == "point" & inlist(method, "fd", "fod") & ///
    inlist(N, 100, 200, 400, 800, 1600) & ///
    inlist(dgp_spec, "gs26_base_v2", "gs26_persist_v2", ///
        "gs26_endog_v2", "gs26_t5_v2") & ///
    (abs(kappa) < 1e-6 | abs(kappa-.1) < 1e-6 | ///
     abs(kappa-.2) < 1e-6 | abs(kappa-.5) < 1e-6 | ///
     abs(kappa-1) < 1e-6) & ///
    ((miss == "balanced" & missp == 0) | ///
     (miss == "mcar" & (abs(missp-.15) < 1e-6 | abs(missp-.30) < 1e-6)) | ///
     (miss == "attrit" & abs(missp-.15) < 1e-6)) & ///
    ((abs(kappa) < 1e-6 & design == "kink") | ///
     (kappa > 1e-6 & design == "jump")) & ///
    ((dgp_spec == "gs26_endog_v2" & qstatus == "endog") | ///
     (dgp_spec != "gs26_endog_v2" & qstatus == "pred"))
if _rc {
    di as err "merge_study_a: registry does not encode the frozen Study-A design"
    exit 459
}
sort pair_id method
capture by pair_id: assert _N == 2 & method[1] == "fd" & method[2] == "fod"
if _rc {
    di as err "merge_study_a: registry lacks an exact FD/FOD pair"
    exit 459
}
egen byte __pairtag = tag(pair_id)
quietly count if __pairtag
if r(N) != 58 {
    di as err "merge_study_a: registry does not contain exactly 58 FD/FOD pairs"
    exit 459
}
drop __pairtag
egen double __expected_total = total(R)
local EXPECTED_TOTAL = __expected_total[1]
drop __expected_total
if `EXPECTED_TOTAL' != `M_ExpectedReplications' | ///
   `ACT_ROWS' != `EXPECTED_TOTAL' {
    di as err "merge_study_a: expected-replication contract mismatch"
    di as err "manifest=`M_ExpectedReplications'; registry=`EXPECTED_TOTAL'; observed=`ACT_ROWS'"
    exit 459
}
quietly save `expected', replace

* The dedicated RepCap=1 smoke assigns whole registry cells round-robin.
* Verify that assignment exactly; formal runs retain within-cell range splits.
tempfile assignment
preserve
keep cell_id cell_order
quietly save `assignment', replace
restore
use `pooled', clear
merge m:1 cell_id using `assignment'
capture assert _merge == 3
if _rc {
    di as err "merge_study_a: result cell is absent from assignment registry"
    exit 459
}
if `M_RepCap' == 1 & `ACT_NSHARD' > 1 {
    capture assert shard == mod(cell_order-1, `ACT_NSHARD')+1
    if _rc {
        di as err "merge_study_a: parallel-smoke cell/shard assignment mismatch"
        exit 459
    }
}
drop _merge cell_order
quietly save `pooled', replace

* Match every observed cell/configuration/count to the dry-run registry.
use `pooled', clear
sort cell_id rep
foreach v in pair_id config_key study design method N T miss missp qstatus ///
                 R kappa dgp_spec grid refine inference_mode {
    capture by cell_id: assert `v' == `v'[1]
    if _rc {
        di as err "merge_study_a: cell_id reused with mixed `v'"
        exit 459
    }
}
by cell_id: gen long n_actual = _N
by cell_id: keep if _n == 1
keep cell_id pair_id config_key study design method N T miss missp qstatus ///
    R kappa dgp_spec grid refine inference_mode n_actual
tempfile actual_cells
quietly save `actual_cells', replace
use `expected', clear
foreach v in pair_id config_key study design method N T miss missp qstatus ///
                 R kappa dgp_spec grid refine inference_mode {
    rename `v' exp_`v'
}
merge 1:1 cell_id using `actual_cells'
quietly count if _merge != 3
if r(N) {
    di as err "merge_study_a: missing or extra registry cells"
    list cell_id _merge if _merge != 3, noobs abbrev(32)
    exit 459
}
drop _merge
quietly count if pair_id != exp_pair_id | config_key != exp_config_key | ///
    study != exp_study | design != exp_design | method != exp_method | ///
    N != exp_N | T != exp_T | miss != exp_miss | ///
    abs(missp-exp_missp) > 1e-6 | qstatus != exp_qstatus | R != exp_R | ///
    abs(kappa-exp_kappa) > 1e-6 | dgp_spec != exp_dgp_spec | ///
    grid != exp_grid | refine != exp_refine | ///
    inference_mode != exp_inference_mode | n_actual != exp_R
if r(N) {
    di as err "merge_study_a: cell arguments/counts differ from registry"
    exit 459
}

* -------------------------------------------------------------------------
* Exact completion-marker and per-shard count contract.
* -------------------------------------------------------------------------
local done : dir "." files "study_a_complete_SH*.csv"
if `: word count `done'' != `ACT_NSHARD' {
    di as err "merge_study_a: completion-marker count differs from NShard"
    exit 459
}
tempfile markers markerone shard_contract
local first 1
forvalues k = 1/`ACT_NSHARD' {
    local mfile : dir "." files "study_a_complete_SH`k'.csv"
    if `: word count `mfile'' != 1 {
        di as err "merge_study_a: missing exact completion marker `k'"
        exit 601
    }
    local mf : word 1 of `mfile'
    quietly import delimited using "`mf'", clear varnames(1) ///
        case(preserve) bindquote(strict)
    local MARKER_VARS run_id harness_version code_version_expected master ///
        shard nshard cells reps exit_code
    unab GOT_MARKER : _all
    if `"`GOT_MARKER'"' != `"`MARKER_VARS'"' | _N != 1 {
        di as err "merge_study_a: marker `k' has wrong schema or row count"
        exit 459
    }
    if shard[1] != `k' | exit_code[1] != 0 | ///
       run_id[1] != "`ACT_RUN'" | harness_version[1] != "`ACT_HARNESS'" | ///
       code_version_expected[1] != "`ACT_CODE'" | ///
       master[1] != `ACT_MASTER' | nshard[1] != `ACT_NSHARD' {
        di as err "merge_study_a: invalid/mixed completion marker `k'"
        exit 459
    }
    quietly save `markerone', replace
    if `first' {
        quietly save `markers', replace
        local first 0
    }
    else {
        quietly append using `markers'
        quietly save `markers', replace
    }
}
use `expected', clear
if `M_RepCap' == 1 & `ACT_NSHARD' > 1 {
    gen int shard = mod(cell_order-1, `ACT_NSHARD')+1
    gen long expected_cells = 1
    gen long expected_reps = R
    collapse (sum) expected_cells expected_reps, by(shard)
}
else {
    expand `ACT_NSHARD'
    bysort cell_id: gen int shard = _n
    gen long __lo = floor((shard-1)*R/`ACT_NSHARD') + 1
    gen long __hi = floor(shard*R/`ACT_NSHARD')
    gen long __n = max(0, __hi-__lo+1)
    gen byte __cell = (__n > 0)
    collapse (sum) expected_cells=__cell expected_reps=__n, by(shard)
}
quietly save `shard_contract', replace
use `pooled', clear
bysort shard cell_id: gen byte __cell = (_n == 1)
collapse (count) actual_reps=rep (sum) actual_cells=__cell, by(shard)
merge 1:1 shard using `shard_contract'
capture assert _merge == 3
if _rc {
    di as err "merge_study_a: shard range does not match registry"
    exit 459
}
drop _merge
merge 1:1 shard using `markers'
capture assert _merge == 3
if _rc {
    di as err "merge_study_a: marker shard set does not match results"
    exit 459
}
drop _merge
quietly count if actual_reps != expected_reps | ///
    actual_cells != expected_cells | reps != expected_reps | ///
    cells != expected_cells
if r(N) {
    di as err "merge_study_a: marker/result counts violate shard contract"
    exit 459
}

* -------------------------------------------------------------------------
* All validation has passed.  Build outputs in temporary files first.
* -------------------------------------------------------------------------
tempfile OUTALL OUTALLCSV OUTCELL OUTCELLCSV OUTPAIR OUTPAIRCSV ///
         OUTPAIRSUM OUTPAIRSUMCSV
use `pooled', clear
sort cell_id rep
compress
* Added 0.9.37 telemetry contract (does not select point estimates).
capture assert inlist(joint_vce,0,1) & inlist(ar_joint,0,1) & ///
    inlist(vce_applied,0,1) & abs(bwscale-1.5)<1e-12 if contract_ok
if _rc exit 459
foreach vv in se_cond_Lagyb se_cond_qb se_cond_consd se_cond_qd se_cond_Lagyd {
    capture assert `vv' >= 0 | missing(`vv')
    if _rc exit 459
}
save `OUTALL', replace
export delimited using "`OUTALLCSV'", replace

* Replication-level FD/FOD paired contrasts.  Missingness/shocks must match.
preserve
gen double err_gamma = g_bias
gen double sqerr_gamma = g_sqerr
gen double abs_err_gamma = abs(err_gamma)
foreach p of local P {
    gen double abs_err_`p' = abs(err_`p')
}
foreach p in gamma Lagyb qb consd qd Lagyd {
    replace err_`p' = . if !success
    replace sqerr_`p' = . if !success
    replace abs_err_`p' = . if !success
}
keep run_id harness_version schema_version code_version_expected master ///
    pair_id study design method N T miss missp qstatus kappa R grid refine ///
    dgp_spec inference_mode rep cell_id success rc fit_rc0 contract_ok ///
    estimator_twostep refine_complete point_boundary ///
    gamma0 gamma_hat err_gamma sqerr_gamma abs_err_gamma ///
    Lagyb qb consd qd Lagyd b_Lagyb b_qb b_consd b_qd b_Lagyd ///
    err_Lagyb sqerr_Lagyb abs_err_Lagyb ///
    err_qb sqerr_qb abs_err_qb err_consd sqerr_consd abs_err_consd ///
    err_qd sqerr_qd abs_err_qd err_Lagyd sqerr_Lagyd abs_err_Lagyd ///
    N_used N_trans N_iv units_realized units_dropped analysis_potential ///
    analysis_observed missing_n missing_rate gap_events gap_periods ///
    elapsed_s dgp_seed missing_seed
gen byte present = 1
reshape wide cell_id present success rc fit_rc0 contract_ok ///
    estimator_twostep refine_complete point_boundary gamma_hat ///
    b_Lagyb b_qb b_consd b_qd b_Lagyd ///
    err_gamma sqerr_gamma abs_err_gamma ///
    err_Lagyb sqerr_Lagyb abs_err_Lagyb err_qb sqerr_qb abs_err_qb ///
    err_consd sqerr_consd abs_err_consd err_qd sqerr_qd abs_err_qd ///
    err_Lagyd sqerr_Lagyd abs_err_Lagyd N_used N_trans N_iv ///
    units_realized units_dropped analysis_potential analysis_observed ///
    missing_n missing_rate gap_events gap_periods elapsed_s dgp_seed ///
    missing_seed, i(run_id pair_id rep) j(method) string
capture confirm variable presentfd
if _rc {
    di as err "merge_study_a: paired reshape lacks FD"
    exit 459
}
capture confirm variable presentfod
if _rc {
    di as err "merge_study_a: paired reshape lacks FOD"
    exit 459
}
capture assert presentfd == 1 & presentfod == 1 & ///
    dgp_seedfd == dgp_seedfod & missing_seedfd == missing_seedfod
if _rc {
    di as err "merge_study_a: incomplete pair or FD/FOD CRN seeds differ"
    exit 459
}
foreach v in units_realized units_dropped analysis_potential ///
             analysis_observed missing_n gap_events gap_periods {
    capture assert `v'fd == `v'fod
    if _rc {
        di as err "merge_study_a: FD/FOD input accounting differs in `v'"
        exit 459
    }
}
capture assert abs(missing_ratefd-missing_ratefod) < 1e-12
if _rc {
    di as err "merge_study_a: FD/FOD realized missing rates differ"
    exit 459
}
gen byte both_success = successfd == 1 & successfod == 1
gen double d_success_fod_fd = successfod-successfd
gen double d_failure_fod_fd = (1-successfod)-(1-successfd)
gen byte delivered_fd = contract_okfd == 1
gen byte delivered_fod = contract_okfod == 1
gen byte twostep_fd = delivered_fd & estimator_twostepfd == 1
gen byte twostep_fod = delivered_fod & estimator_twostepfod == 1
gen byte refine_applicable_fd = delivered_fd
gen byte refine_applicable_fod = delivered_fod
gen byte refine_complete_fd = refine_applicable_fd & refine_completefd == 1
gen byte refine_complete_fod = refine_applicable_fod & refine_completefod == 1
gen byte boundary_fd = delivered_fd & point_boundaryfd > 0 & ///
    !missing(point_boundaryfd)
gen byte boundary_fod = delivered_fod & point_boundaryfod > 0 & ///
    !missing(point_boundaryfod)
foreach p in gamma Lagyb qb consd qd Lagyd {
    gen double d_err_`p'_fod_fd = err_`p'fod-err_`p'fd
    gen double d_abs_err_`p'_fod_fd = abs_err_`p'fod-abs_err_`p'fd
    gen double d_sqerr_`p'_fod_fd = sqerr_`p'fod-sqerr_`p'fd
    gen double err_`p'_fd_common = err_`p'fd if both_success
    gen double err_`p'_fod_common = err_`p'fod if both_success
    gen double sqerr_`p'_fd_common = sqerr_`p'fd if both_success
    gen double sqerr_`p'_fod_common = sqerr_`p'fod if both_success
}
foreach v in N_used N_trans N_iv elapsed_s {
    gen double d_`v'_fod_fd = `v'fod-`v'fd
}
sort pair_id rep
compress
save `OUTPAIR', replace
export delimited using "`OUTPAIRCSV'", replace

* Pair-level summaries.  Paired MSE differences are the efficiency contrast;
* their MCSE uses the within-replication FD/FOD covariance automatically.
gen double success_fd = successfd
gen double success_fod = successfod
local PCOLL ""
foreach p in gamma Lagyb qb consd qd Lagyd {
    local PCOLL `PCOLL' n_fd_`p'=sqerr_`p'fd ///
        n_fod_`p'=sqerr_`p'fod n_both_`p'=d_sqerr_`p'_fod_fd ///
        n_common_`p'=sqerr_`p'_fd_common
}
foreach v in N_used N_trans N_iv elapsed_s {
    local PCOLL `PCOLL' n_d_`v'=d_`v'_fod_fd
}
local PMEANS ""
local PSDS ""
foreach p in gamma Lagyb qb consd qd Lagyd {
    local PMEANS `PMEANS' bias_`p'_fd_marginal=err_`p'fd ///
        bias_`p'_fod_marginal=err_`p'fod ///
        mse_`p'_fd_marginal=sqerr_`p'fd ///
        mse_`p'_fod_marginal=sqerr_`p'fod ///
        bias_`p'_fd_common=err_`p'_fd_common ///
        bias_`p'_fod_common=err_`p'_fod_common ///
        mse_`p'_fd_common=sqerr_`p'_fd_common ///
        mse_`p'_fod_common=sqerr_`p'_fod_common ///
        d_bias_`p'=d_err_`p'_fod_fd ///
        d_mae_`p'=d_abs_err_`p'_fod_fd ///
        d_mse_`p'=d_sqerr_`p'_fod_fd
    local PSDS `PSDS' sd_err_`p'_fd_marginal=err_`p'fd ///
        sd_err_`p'_fod_marginal=err_`p'fod ///
        sd_sqerr_`p'_fd_marginal=sqerr_`p'fd ///
        sd_sqerr_`p'_fod_marginal=sqerr_`p'fod ///
        sd_err_`p'_fd_common=err_`p'_fd_common ///
        sd_err_`p'_fod_common=err_`p'_fod_common ///
        sd_sqerr_`p'_fd_common=sqerr_`p'_fd_common ///
        sd_sqerr_`p'_fod_common=sqerr_`p'_fod_common ///
        sd_d_bias_`p'=d_err_`p'_fod_fd ///
        sd_d_mae_`p'=d_abs_err_`p'_fod_fd ///
        sd_d_mse_`p'=d_sqerr_`p'_fod_fd
}
collapse (count) n_pair=rep `PCOLL' ///
    (sum) n_both_success=both_success ///
        n_delivered_fd=delivered_fd n_delivered_fod=delivered_fod ///
        n_twostep_fd=twostep_fd n_twostep_fod=twostep_fod ///
        n_refine_applicable_fd=refine_applicable_fd ///
        n_refine_applicable_fod=refine_applicable_fod ///
        n_refine_complete_fd=refine_complete_fd ///
        n_refine_complete_fod=refine_complete_fod ///
        n_boundary_fd=boundary_fd n_boundary_fod=boundary_fod ///
    (mean) success_fd success_fod d_success=d_success_fod_fd ///
        d_failure=d_failure_fod_fd `PMEANS' ///
        d_N_used=d_N_used_fod_fd d_N_trans=d_N_trans_fod_fd ///
        d_N_iv=d_N_iv_fod_fd d_elapsed_s=d_elapsed_s_fod_fd ///
        mean_missing_rate=missing_ratefd mean_gap_events=gap_eventsfd ///
        mean_gap_periods=gap_periodsfd ///
    (sd) sd_d_success=d_success_fod_fd ///
        sd_d_failure=d_failure_fod_fd `PSDS' ///
        sd_d_N_used=d_N_used_fod_fd sd_d_N_trans=d_N_trans_fod_fd ///
        sd_d_N_iv=d_N_iv_fod_fd sd_d_elapsed_s=d_elapsed_s_fod_fd ///
        sd_missing_rate=missing_ratefd sd_gap_events=gap_eventsfd ///
        sd_gap_periods=gap_periodsfd, ///
    by(run_id harness_version schema_version code_version_expected master ///
       pair_id study design N T miss missp qstatus kappa R grid refine ///
       dgp_spec inference_mode)
gen double d_success_mcse = sd_d_success/sqrt(n_pair)
gen double d_failure_mcse = sd_d_failure/sqrt(n_pair)
gen double success_fd_mcse = sqrt(success_fd*(1-success_fd)/n_pair)
gen double success_fod_mcse = sqrt(success_fod*(1-success_fod)/n_pair)
gen double delivery_rate_fd = n_delivered_fd/n_pair
gen double delivery_rate_fod = n_delivered_fod/n_pair
gen double twostep_rate_fd = n_twostep_fd/n_pair
gen double twostep_rate_fod = n_twostep_fod/n_pair
gen double twostep_given_delivery_fd = n_twostep_fd/n_delivered_fd ///
    if n_delivered_fd > 0
gen double twostep_given_delivery_fod = n_twostep_fod/n_delivered_fod ///
    if n_delivered_fod > 0
gen double refine_complete_rate_fd = ///
    n_refine_complete_fd/n_refine_applicable_fd ///
    if n_refine_applicable_fd > 0
gen double refine_complete_rate_fod = ///
    n_refine_complete_fod/n_refine_applicable_fod ///
    if n_refine_applicable_fod > 0
gen double boundary_rate_fd = n_boundary_fd/n_delivered_fd ///
    if n_delivered_fd > 0
gen double boundary_rate_fod = n_boundary_fod/n_delivered_fod ///
    if n_delivered_fod > 0
gen double delivery_rate_fd_mcse = ///
    sqrt(delivery_rate_fd*(1-delivery_rate_fd)/n_pair)
gen double delivery_rate_fod_mcse = ///
    sqrt(delivery_rate_fod*(1-delivery_rate_fod)/n_pair)
gen double twostep_rate_fd_mcse = ///
    sqrt(twostep_rate_fd*(1-twostep_rate_fd)/n_pair)
gen double twostep_rate_fod_mcse = ///
    sqrt(twostep_rate_fod*(1-twostep_rate_fod)/n_pair)
gen double twostep_given_delivery_fd_mcse = ///
    sqrt(twostep_given_delivery_fd*(1-twostep_given_delivery_fd) ///
         /n_delivered_fd) if n_delivered_fd > 0
gen double twostep_given_delivery_fod_mcse = ///
    sqrt(twostep_given_delivery_fod*(1-twostep_given_delivery_fod) ///
         /n_delivered_fod) if n_delivered_fod > 0
gen double refine_complete_rate_fd_mcse = ///
    sqrt(refine_complete_rate_fd*(1-refine_complete_rate_fd) ///
         /n_refine_applicable_fd) if n_refine_applicable_fd > 0
gen double refine_complete_rate_fod_mcse = ///
    sqrt(refine_complete_rate_fod*(1-refine_complete_rate_fod) ///
         /n_refine_applicable_fod) if n_refine_applicable_fod > 0
gen double boundary_rate_fd_mcse = ///
    sqrt(boundary_rate_fd*(1-boundary_rate_fd)/n_delivered_fd) ///
    if n_delivered_fd > 0
gen double boundary_rate_fod_mcse = ///
    sqrt(boundary_rate_fod*(1-boundary_rate_fod)/n_delivered_fod) ///
    if n_delivered_fod > 0
foreach p in gamma Lagyb qb consd qd Lagyd {
    capture assert n_common_`p' == n_both_`p'
    if _rc {
        di as err "merge_study_a: inconsistent common-success count for `p'"
        exit 459
    }
    capture assert abs(d_bias_`p' - ///
        (bias_`p'_fod_common-bias_`p'_fd_common)) < 1e-10 ///
        if n_common_`p' > 0
    if _rc {
        di as err "merge_study_a: paired/common bias identity failed for `p'"
        exit 459
    }
    capture assert abs(d_mse_`p' - ///
        (mse_`p'_fod_common-mse_`p'_fd_common)) < 1e-10 ///
        if n_common_`p' > 0
    if _rc {
        di as err "merge_study_a: paired/common MSE identity failed for `p'"
        exit 459
    }

    gen double rmse_`p'_fd_marginal = sqrt(mse_`p'_fd_marginal)
    gen double rmse_`p'_fod_marginal = sqrt(mse_`p'_fod_marginal)
    gen double rmse_ratio_`p'_fod_fd_marginal = ///
        rmse_`p'_fod_marginal/rmse_`p'_fd_marginal
    gen double bias_`p'_fd_marginal_mcse = ///
        sd_err_`p'_fd_marginal/sqrt(n_fd_`p')
    gen double bias_`p'_fod_marginal_mcse = ///
        sd_err_`p'_fod_marginal/sqrt(n_fod_`p')
    gen double rmse_`p'_fd_marginal_mcse = ///
        sd_sqerr_`p'_fd_marginal/sqrt(n_fd_`p') ///
        /(2*rmse_`p'_fd_marginal) if rmse_`p'_fd_marginal > 0
    gen double rmse_`p'_fod_marginal_mcse = ///
        sd_sqerr_`p'_fod_marginal/sqrt(n_fod_`p') ///
        /(2*rmse_`p'_fod_marginal) if rmse_`p'_fod_marginal > 0

    gen double rmse_`p'_fd_common = sqrt(mse_`p'_fd_common)
    gen double rmse_`p'_fod_common = sqrt(mse_`p'_fod_common)
    gen double rmse_ratio_`p'_fod_fd_common = ///
        rmse_`p'_fod_common/rmse_`p'_fd_common
    gen double bias_`p'_fd_common_mcse = ///
        sd_err_`p'_fd_common/sqrt(n_common_`p')
    gen double bias_`p'_fod_common_mcse = ///
        sd_err_`p'_fod_common/sqrt(n_common_`p')
    gen double rmse_`p'_fd_common_mcse = ///
        sd_sqerr_`p'_fd_common/sqrt(n_common_`p') ///
        /(2*rmse_`p'_fd_common) if rmse_`p'_fd_common > 0
    gen double rmse_`p'_fod_common_mcse = ///
        sd_sqerr_`p'_fod_common/sqrt(n_common_`p') ///
        /(2*rmse_`p'_fod_common) if rmse_`p'_fod_common > 0

    gen double d_bias_`p'_mcse = sd_d_bias_`p'/sqrt(n_common_`p')
    gen double d_mae_`p'_mcse = sd_d_mae_`p'/sqrt(n_common_`p')
    gen double d_mse_`p'_mcse = sd_d_mse_`p'/sqrt(n_common_`p')
    drop mse_`p'_fd_marginal mse_`p'_fod_marginal ///
        mse_`p'_fd_common mse_`p'_fod_common ///
        sd_err_`p'_fd_marginal sd_err_`p'_fod_marginal ///
        sd_sqerr_`p'_fd_marginal sd_sqerr_`p'_fod_marginal ///
        sd_err_`p'_fd_common sd_err_`p'_fod_common ///
        sd_sqerr_`p'_fd_common sd_sqerr_`p'_fod_common ///
        sd_d_bias_`p' sd_d_mae_`p' sd_d_mse_`p'
}
foreach v in N_used N_trans N_iv elapsed_s {
    gen double d_`v'_mcse = sd_d_`v'/sqrt(n_d_`v')
    drop sd_d_`v'
}
* These design diagnostics were asserted identical within every method pair.
gen double mean_missing_rate_mcse = sd_missing_rate/sqrt(n_pair)
gen double mean_gap_events_mcse = sd_gap_events/sqrt(n_pair)
gen double mean_gap_periods_mcse = sd_gap_periods/sqrt(n_pair)
drop sd_missing_rate sd_gap_events sd_gap_periods
drop sd_d_success sd_d_failure
sort design N miss kappa
compress
save `OUTPAIRSUM', replace
export delimited using "`OUTPAIRSUMCSV'", replace
restore

* Cell-level point-estimation summaries.  The frozen primary loss sample is
* contract-valid two-step; delivered and fallback panels are secondary views
* of the same rows and never trigger additional estimation.
gen byte failure = 1-success
gen byte rc_nonzero = (rc != 0)
gen byte contract_failure = (fit_rc0 == 1 & contract_ok == 0)
gen byte delivered = contract_ok == 1
gen byte twostep_delivered = delivered & estimator_twostep == 1
gen byte one_step_fallback = delivered & estimator_twostep == 0
gen byte refine_applicable = delivered
gen byte refine_complete_ok = refine_applicable & refine_complete == 1
gen byte refine_incomplete = refine_applicable & refine_complete == 0
gen byte boundary_hit = delivered & point_boundary > 0 & ///
    !missing(point_boundary)
gen double err_gamma = g_bias if success
gen double sqerr_gamma = g_sqerr if success
gen double err_gamma_delivered = g_bias if delivered
gen double sqerr_gamma_delivered = g_sqerr if delivered
gen double err_gamma_fallback = g_bias if one_step_fallback
gen double sqerr_gamma_fallback = g_sqerr if one_step_fallback
gen double gamma_hat_ok = gamma_hat if success
gen byte gamma_within_005 = (abs(g_bias) <= .05) if success
gen byte hansen_reject5 = (hansen_p < .05) if success & !missing(hansen_p)
gen byte ar2_reject5    = (ar2_p < .05)    if success & !missing(ar2_p)
gen double hansen_p_ok = hansen_p if success
gen double ar1_p_ok = ar1_p if success
gen double ar2_p_ok = ar2_p if success
foreach p of local P {
    gen double err_`p'_delivered = err_`p' if delivered
    gen double sqerr_`p'_delivered = sqerr_`p' if delivered
    gen double err_`p'_fallback = err_`p' if one_step_fallback
    gen double sqerr_`p'_fallback = sqerr_`p' if one_step_fallback
    replace err_`p' = . if !success
    replace sqerr_`p' = . if !success
    replace se_`p' = . if !success
    gen byte wcov_`p' = (abs(err_`p') <= invnormal(.975)*se_`p') ///
        if success & !missing(err_`p', se_`p') & se_`p' >= 0
}
local CCOLL n_gamma=err_gamma ///
    n_gamma_delivered=err_gamma_delivered ///
    n_gamma_fallback=err_gamma_fallback ///
    n_gamma_within_005=gamma_within_005 ///
    n_n_used=N_used n_n_trans=N_trans n_n_iv=N_iv n_elapsed=elapsed_s ///
    n_hansen=hansen_reject5 n_ar2=ar2_reject5
local CMEANS bias_gamma=err_gamma mse_gamma=sqerr_gamma ///
    bias_gamma_delivered=err_gamma_delivered ///
    mse_gamma_delivered=sqerr_gamma_delivered ///
    bias_gamma_fallback=err_gamma_fallback ///
    mse_gamma_fallback=sqerr_gamma_fallback ///
    prob_gamma_within_005=gamma_within_005 ///
    mean_hansen_p=hansen_p_ok hansen_reject5=hansen_reject5 ///
    mean_ar1_p=ar1_p_ok mean_ar2_p=ar2_p_ok ar2_reject5=ar2_reject5
local CSDS sd_err_gamma=err_gamma sd_sqerr_gamma=sqerr_gamma ///
    sd_err_gamma_delivered=err_gamma_delivered ///
    sd_sqerr_gamma_delivered=sqerr_gamma_delivered ///
    sd_err_gamma_fallback=err_gamma_fallback ///
    sd_sqerr_gamma_fallback=sqerr_gamma_fallback
local CP05 p05_gamma=gamma_hat_ok
local CP25 p25_gamma=gamma_hat_ok
local CP50 median_gamma=gamma_hat_ok
local CP75 p75_gamma=gamma_hat_ok
local CP95 p95_gamma=gamma_hat_ok
foreach p of local P {
    local pl = lower("`p'")
    local CCOLL `CCOLL' n_`pl'=err_`p' ///
        n_`pl'_delivered=err_`p'_delivered ///
        n_`pl'_fallback=err_`p'_fallback n_wcov_`pl'=wcov_`p'
    local CMEANS `CMEANS' bias_`pl'=err_`p' mse_`pl'=sqerr_`p' ///
        bias_`pl'_delivered=err_`p'_delivered ///
        mse_`pl'_delivered=sqerr_`p'_delivered ///
        bias_`pl'_fallback=err_`p'_fallback ///
        mse_`pl'_fallback=sqerr_`p'_fallback ///
        mean_se_`pl'=se_`p' wcov_`pl'=wcov_`p'
    local CSDS `CSDS' sd_err_`pl'=err_`p' sd_sqerr_`pl'=sqerr_`p' ///
        sd_err_`pl'_delivered=err_`p'_delivered ///
        sd_sqerr_`pl'_delivered=sqerr_`p'_delivered ///
        sd_err_`pl'_fallback=err_`p'_fallback ///
        sd_sqerr_`pl'_fallback=sqerr_`p'_fallback
    local CP25 `CP25' p25_`pl'=b_`p'
    local CP50 `CP50' median_`pl'=b_`p'
    local CP75 `CP75' p75_`pl'=b_`p'
}

* SE panels use the SAME error subset for mean(SE) and empirical SD.
* j=reported joint, f=reported conditional fallback, c=V_cond diagnostic.
local XCNT
local XMEAN
local XSD
foreach nm in Lagyb qb consd qd Lagyd {
    local nn = lower("`nm'")
    foreach kind in j f c {
        local ss = cond("`kind'"=="c","se_cond_`nm'","se_`nm'")
        local gate = cond("`kind'"=="j","joint_vce==1", ///
            cond("`kind'"=="f","joint_vce==0","1"))
        gen double xe_`kind'_`nn' = b_`nm'-`nm' ///
            if success & (1) & (`gate') & !missing(`ss') & `ss'>=0
        gen double xs_`kind'_`nn' = `ss' if !missing(xe_`kind'_`nn')
        gen byte xc_`kind'_`nn' = abs(xe_`kind'_`nn')<=invnormal(.975)*`ss' ///
            if !missing(xe_`kind'_`nn')
        local XCNT `XCNT' n_`kind'_`nn'=xc_`kind'_`nn'
        local XMEAN `XMEAN' meanse_`kind'_`nn'=xs_`kind'_`nn' ///
            cov_`kind'_`nn'=xc_`kind'_`nn'
        local XSD `XSD' esd_`kind'_`nn'=xe_`kind'_`nn'
    }
}
gen byte x_joint = joint_vce==1 if contract_ok
gen byte x_ar_joint = ar_joint==1 if contract_ok
gen byte x_se = se_delivered if contract_ok
gen byte x_near = N_iv_dep_near>0 if contract_ok & !missing(N_iv_dep_near)
gen byte x_ar2c = ar2_p_cond<.05 if success & !missing(ar2_p_cond)
local XCNT `XCNT' n_joint_eval=x_joint n_ar_joint_eval=x_ar_joint ///
    n_se_eval=x_se n_iv_near_eval=x_near n_ar2_cond=x_ar2c
local XMEAN `XMEAN' joint_rate=x_joint ar_joint_rate=x_ar_joint ///
    se_delivery_rate=x_se iv_near_rate=x_near ar2_cond_reject5=x_ar2c
collapse (count) n_rep=rep `CCOLL' `XCNT' ///
    (sum) n_success=success n_failure=failure n_rc_nonzero=rc_nonzero ///
        n_fit_rc0=fit_rc0 n_contract_failure=contract_failure ///
        n_delivered=delivered n_twostep=twostep_delivered ///
        n_one_step_fallback=one_step_fallback ///
        n_refine_applicable=refine_applicable ///
        n_refine_complete=refine_complete_ok ///
        n_refine_incomplete=refine_incomplete n_boundary=boundary_hit ///
    (mean) `CMEANS' `XMEAN' mean_gamma_hat=gamma_hat_ok ///
        mean_N_used=N_used mean_N_trans=N_trans mean_N_iv=N_iv ///
        mean_units_realized=units_realized ///
        mean_units_dropped=units_dropped ///
        mean_analysis_potential=analysis_potential ///
        mean_analysis_observed=analysis_observed mean_missing_n=missing_n ///
        mean_missing_rate=missing_rate mean_gap_events=gap_events ///
        mean_gap_periods=gap_periods mean_elapsed_s=elapsed_s ///
    (sd) `CSDS' `XSD' sd_N_used=N_used sd_N_trans=N_trans sd_N_iv=N_iv ///
        sd_units_realized=units_realized sd_units_dropped=units_dropped ///
        sd_analysis_potential=analysis_potential ///
        sd_analysis_observed=analysis_observed sd_missing_n=missing_n ///
        sd_missing_rate=missing_rate sd_gap_events=gap_events ///
        sd_gap_periods=gap_periods sd_elapsed_s=elapsed_s ///
    (p5) `CP05' (p25) `CP25' (p50) `CP50' (p75) `CP75' ///
    (p95) `CP95', ///
    by(run_id harness_version schema_version code_version_expected master ///
       cell_id pair_id study design method N T miss missp qstatus kappa R ///
       grid refine dgp_spec inference_mode gamma0 rhoy_true bq_true ///
       consd_true qd_true lagyd_true rhoq rhoeu_effective sige tburn ///
       maxlag_lo maxlag_hi trim history gridtype gridsample vce)
capture assert n_gamma == n_success & ///
    n_gamma_delivered == n_delivered & ///
    n_gamma_fallback == n_one_step_fallback & ///
    n_gamma_within_005 == n_gamma
if _rc {
    di as err "merge_study_a: inconsistent gamma loss denominators"
    exit 459
}
foreach p of local P {
    local pl = lower("`p'")
    capture assert n_`pl' == n_success & ///
        n_`pl'_delivered == n_delivered & ///
        n_`pl'_fallback == n_one_step_fallback
    if _rc {
        di as err "merge_study_a: inconsistent `pl' loss denominators"
        exit 459
    }
}
gen double success_rate = n_success/n_rep
gen double success_mcse = sqrt(success_rate*(1-success_rate)/n_rep)
gen double failure_rate = n_failure/n_rep
gen double failure_mcse = sqrt(failure_rate*(1-failure_rate)/n_rep)
gen double fit_rate = n_fit_rc0/n_rep
gen double fit_rate_mcse = sqrt(fit_rate*(1-fit_rate)/n_rep)
gen double delivery_rate = n_delivered/n_rep
gen double delivery_rate_mcse = ///
    sqrt(delivery_rate*(1-delivery_rate)/n_rep)
gen double twostep_rate = n_twostep/n_rep
gen double twostep_rate_mcse = ///
    sqrt(twostep_rate*(1-twostep_rate)/n_rep)
gen double twostep_given_delivery = n_twostep/n_delivered ///
    if n_delivered > 0
gen double twostep_given_delivery_mcse = ///
    sqrt(twostep_given_delivery*(1-twostep_given_delivery)/n_delivered) ///
    if n_delivered > 0
gen double refine_complete_rate = n_refine_complete/n_refine_applicable ///
    if n_refine_applicable > 0
gen double refine_complete_rate_mcse = ///
    sqrt(refine_complete_rate*(1-refine_complete_rate) ///
         /n_refine_applicable) if n_refine_applicable > 0
gen double boundary_rate = n_boundary/n_delivered if n_delivered > 0
gen double boundary_rate_mcse = ///
    sqrt(boundary_rate*(1-boundary_rate)/n_delivered) ///
    if n_delivered > 0
gen double rmse_gamma = sqrt(mse_gamma) if n_gamma > 0
gen double bias_gamma_mcse = sd_err_gamma/sqrt(n_gamma) if n_gamma > 1
gen double rmse_gamma_mcse = sd_sqerr_gamma/sqrt(n_gamma)/(2*rmse_gamma) ///
    if n_gamma > 1 & rmse_gamma > 0
gen double empirical_sd_gamma = sd_err_gamma
foreach sample in delivered fallback {
    gen double rmse_gamma_`sample' = sqrt(mse_gamma_`sample') ///
        if n_gamma_`sample' > 0
    gen double bias_gamma_`sample'_mcse = ///
        sd_err_gamma_`sample'/sqrt(n_gamma_`sample') ///
        if n_gamma_`sample' > 1
    gen double rmse_gamma_`sample'_mcse = ///
        sd_sqerr_gamma_`sample'/sqrt(n_gamma_`sample') ///
        /(2*rmse_gamma_`sample') ///
        if n_gamma_`sample' > 1 & rmse_gamma_`sample' > 0
}
gen double prob_gamma_within_005_mcse = ///
    sqrt(prob_gamma_within_005*(1-prob_gamma_within_005) ///
         /n_gamma_within_005) if n_gamma_within_005 > 0
drop mse_gamma sd_err_gamma sd_sqerr_gamma
drop mse_gamma_delivered sd_err_gamma_delivered ///
    sd_sqerr_gamma_delivered mse_gamma_fallback ///
    sd_err_gamma_fallback sd_sqerr_gamma_fallback
foreach p of local P {
    local pl = lower("`p'")
    gen double rmse_`pl' = sqrt(mse_`pl') if n_`pl' > 0
    gen double bias_`pl'_mcse = sd_err_`pl'/sqrt(n_`pl') if n_`pl' > 1
    gen double rmse_`pl'_mcse = sd_sqerr_`pl'/sqrt(n_`pl')/(2*rmse_`pl') ///
        if n_`pl' > 1 & rmse_`pl' > 0
    gen double empirical_sd_`pl' = sd_err_`pl'
    foreach sample in delivered fallback {
        gen double rmse_`pl'_`sample' = sqrt(mse_`pl'_`sample') ///
            if n_`pl'_`sample' > 0
        gen double bias_`pl'_`sample'_mcse = ///
            sd_err_`pl'_`sample'/sqrt(n_`pl'_`sample') ///
            if n_`pl'_`sample' > 1
        gen double rmse_`pl'_`sample'_mcse = ///
            sd_sqerr_`pl'_`sample'/sqrt(n_`pl'_`sample') ///
            /(2*rmse_`pl'_`sample') ///
            if n_`pl'_`sample' > 1 & rmse_`pl'_`sample' > 0
    }
    gen double se_sd_ratio_`pl' = mean_se_`pl'/empirical_sd_`pl' ///
        if empirical_sd_`pl' > 0
    gen double wcov_`pl'_mcse = sqrt(wcov_`pl'*(1-wcov_`pl')/n_wcov_`pl') ///
        if n_wcov_`pl' > 0
    drop mse_`pl' sd_err_`pl' sd_sqerr_`pl' ///
        mse_`pl'_delivered sd_err_`pl'_delivered ///
        sd_sqerr_`pl'_delivered mse_`pl'_fallback ///
        sd_err_`pl'_fallback sd_sqerr_`pl'_fallback
}
gen double hansen_reject5_mcse = sqrt(hansen_reject5*(1-hansen_reject5)/n_hansen) ///
    if n_hansen > 0
gen double ar2_reject5_mcse = sqrt(ar2_reject5*(1-ar2_reject5)/n_ar2) ///
    if n_ar2 > 0
foreach v in N_used N_trans N_iv {
    local vl = lower("`v'")
    gen double mean_`v'_mcse = sd_`v'/sqrt(n_`vl')
    drop sd_`v'
}
gen double mean_elapsed_s_mcse = sd_elapsed_s/sqrt(n_elapsed)
drop sd_elapsed_s
foreach v in units_realized units_dropped analysis_potential ///
             analysis_observed missing_n missing_rate gap_events gap_periods {
    gen double mean_`v'_mcse = sd_`v'/sqrt(n_rep)
    drop sd_`v'
}
order run_id harness_version schema_version code_version_expected master ///
    cell_id pair_id study design method N T miss missp qstatus kappa R ///
    grid refine dgp_spec inference_mode gamma0 rhoy_true bq_true ///
    consd_true qd_true lagyd_true rhoq rhoeu_effective sige tburn ///
    maxlag_lo maxlag_hi trim history gridtype gridsample vce ///
    n_rep n_fit_rc0 fit_rate fit_rate_mcse ///
    n_delivered delivery_rate delivery_rate_mcse ///
    n_twostep twostep_rate twostep_rate_mcse ///
    twostep_given_delivery twostep_given_delivery_mcse ///
    n_refine_applicable n_refine_complete n_refine_incomplete ///
    refine_complete_rate refine_complete_rate_mcse ///
    n_boundary boundary_rate boundary_rate_mcse ///
    n_success n_failure success_rate success_mcse failure_rate failure_mcse ///
    n_rc_nonzero n_contract_failure n_one_step_fallback ///
    n_gamma bias_gamma bias_gamma_mcse rmse_gamma ///
    rmse_gamma_mcse empirical_sd_gamma mean_gamma_hat ///
    p05_gamma p25_gamma median_gamma p75_gamma p95_gamma ///
    n_gamma_within_005 prob_gamma_within_005 ///
    prob_gamma_within_005_mcse ///
    n_gamma_delivered bias_gamma_delivered ///
    bias_gamma_delivered_mcse rmse_gamma_delivered ///
    rmse_gamma_delivered_mcse ///
    n_gamma_fallback bias_gamma_fallback bias_gamma_fallback_mcse ///
    rmse_gamma_fallback rmse_gamma_fallback_mcse
sort design N miss kappa method
compress

foreach nm in lagyb qb consd qd lagyd {
    foreach kind in j f c {
        gen double se_sd_`kind'_`nm' = meanse_`kind'_`nm'/esd_`kind'_`nm' ///
            if esd_`kind'_`nm'>0 & !missing(esd_`kind'_`nm')
        gen double cov_mcse_`kind'_`nm' = ///
            sqrt(cov_`kind'_`nm'*(1-cov_`kind'_`nm')/n_`kind'_`nm') ///
            if n_`kind'_`nm'>0
        * Effective joint delivery+coverage, NOT coverage conditional on delivery.
        gen double cov_eff_`kind'_`nm' = ///
            cond(n_`kind'_`nm'>0,cov_`kind'_`nm'*n_`kind'_`nm',0)/n_rep
    }
}
save `OUTCELL', replace
export delimited using "`OUTCELLCSV'", replace

* Publish only after every check, reshape, and summary has completed.
copy "`OUTALL'" "study_a_all.dta", replace
copy "`OUTALLCSV'" "study_a_all.csv", replace
copy "`OUTCELL'" "study_a_summary.dta", replace
copy "`OUTCELLCSV'" "study_a_summary.csv", replace
copy "`OUTPAIR'" "study_a_paired_fd_fod.dta", replace
copy "`OUTPAIRCSV'" "study_a_paired_fd_fod.csv", replace
copy "`OUTPAIRSUM'" "study_a_paired_summary.dta", replace
copy "`OUTPAIRSUMCSV'" "study_a_paired_summary.csv", replace

* This is the sole success attestation.  Its absence means that one or more
* validations/builds/copies did not complete, even if partial outputs exist.
tempfile ATTEST
tempname AH
file open `AH' using "`ATTEST'", write replace text
file write `AH' "run_id,manifest_sha256,rows,cells,nshard,exit_code" _n
file write `AH' "`ACT_RUN',`MANIFEST_SHA',`EXPECTED_TOTAL',116,`ACT_NSHARD',0" _n
file close `AH'
copy "`ATTEST'" "`MERGE_MARKER'", replace

di as res "STUDY A MERGE COMPLETE AND VALIDATED: `ACT_RUN'"
di as txt "116 cells; `EXPECTED_TOTAL' rows; `ACT_NSHARD' shards"
di as txt "Point-only MC: gamma + five structural coefficients; no bootstrap"
di as txt "Wrote study_a_all.*, study_a_summary.*, and paired FD/FOD files"
exit 0
