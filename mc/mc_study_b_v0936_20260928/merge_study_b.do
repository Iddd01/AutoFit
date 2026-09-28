*!==========================================================================
*! merge_study_b.do -- fail-closed merger for the frozen Type-4 Monte Carlo
*!--------------------------------------------------------------------------
*! Harness study_b_type4_v0936; worker schema study_b_worker_v0936 (2026-09-28).
*! Validates the exact registry/shard contract before creating any report.
*!==========================================================================
version 15.0
clear all
set more off
set varabbrev off

* -------------------------------------------------------------------------
* v8r1 fail-closed merger.  The legacy implementation remains below this block
* only as inert history; exit 0 at the end of this block prevents execution.
* -------------------------------------------------------------------------
args RUN_ID MERGE_NONCE
local HARNESS_EXPECTED "study_b_type4_v0936"
local SCHEMA_EXPECTED  "study_b_worker_v0936"
local CODE_EXPECTED    "0.9.36"
local RNG_EXPECTED     "mt64"
local SEEDS_EXPECTED   "hash131_mod3_v2"

* Only the PowerShell verifier may authorize a merge.  It recomputes every
* staged-source and raw-input SHA-256 before creating this one-use nonce.
if !regexm("`RUN_ID'", "^[A-Za-z0-9][A-Za-z0-9_.-]*$") | ///
   !regexm("`RUN_ID'", "[A-Za-z]") | ///
   !regexm("`MERGE_NONCE'", "^[0-9A-Fa-f]+$") | ///
   strlen("`MERGE_NONCE'") != 32 {
    di as err "merge_study_b: valid RUN_ID and verifier nonce are required"
    exit 198
}
capture confirm file "_merge_authorized.txt"
if _rc {
    di as err "merge_study_b: verifier authorization is required"
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
    di as err "merge_study_b: malformed or mismatched verifier authorization"
    exit 459
}
local MANIFEST_SHA = substr(`"`AUTH_HASH'"', 17, 64)
capture erase "_merge_authorized.txt"
if _rc {
    di as err "merge_study_b: cannot consume verifier authorization"
    exit 603
}

local MERGE_MARKER "study_b_merge_complete.csv"
capture confirm file "`MERGE_MARKER'"
if !_rc {
    capture erase "`MERGE_MARKER'"
    if _rc {
        di as err "merge_study_b: cannot remove stale merge marker"
        exit 603
    }
}

if "`RUN_ID'" != "" {
    if !regexm("`RUN_ID'", "^[A-Za-z0-9][A-Za-z0-9_.-]*$") {
        di as err "merge_study_b: unsafe RUN_ID"
        exit 198
    }
    local here : dir "." files "study_b_results_SH*.csv"
    if `: word count `here'' == 0 {
        capture cd "runs/`RUN_ID'"
        if _rc {
            di as err "merge_study_b: cannot find runs/`RUN_ID'"
            exit 601
        }
    }
}
foreach f in study_b_cells.do study_b_worker.do {
    capture confirm file "`f'"
    if _rc {
        di as err "merge_study_b: missing staged file `f'"
        exit 601
    }
}

* The launcher manifest is the immutable run contract.  In particular REP_CAP
* must never be inferred from realized output: doing so could mistake a
* uniformly truncated run for an intentional smoke run.
capture confirm file "manifest.json"
if _rc {
    di as err "merge_study_b: manifest.json is required"
    exit 601
}
local MJ = fileread("manifest.json")
if !regexm(`"`MJ'"', `""schemaVersion"[ ]*:[ ]*([0-9]+)"') {
    di as err "merge_study_b: manifest lacks schemaVersion"
    exit 459
}
local M_Schema = real(regexs(1))
foreach q in NShard Master RepCap ExpectedCells ExpectedReplications {
    if !regexm(`"`MJ'"', `""`q'"[ ]*:[ ]*([0-9]+)"') {
        di as err "merge_study_b: manifest lacks numeric config `q'"
        exit 459
    }
    local M_`q' = real(regexs(1))
}
foreach q in RunId B Grid GridCI HarnessVersion CodeVersionExpected RngKind {
    if !regexm(`"`MJ'"', `""`q'"[ ]*:[ ]*"([^"]+)""') {
        di as err "merge_study_b: manifest lacks string config `q'"
        exit 459
    }
    local M_`q' "`=regexs(1)'"
}
if !regexm("`M_RunId'", "[A-Za-z]") | ///
   !regexm("`M_RunId'", "^[A-Za-z0-9][A-Za-z0-9_.-]*$") {
    di as err "merge_study_b: manifest RunId must contain a letter"
    exit 459
}
* Arbitrary numeric overrides are valid smoke settings; a dot means registry.
foreach q in B Grid GridCI {
    if "`M_`q''" != "." {
        capture confirm integer number `M_`q''
        if _rc {
            di as err "merge_study_b: invalid manifest override `q'"
            exit 459
        }
    }
}
if ("`M_B'" != "." & (real("`M_B'") < 10 | real("`M_B'") > 100000)) | ///
   ("`M_Grid'" != "." & (real("`M_Grid'") < 10 | real("`M_Grid'") > 10000)) | ///
   ("`M_GridCI'" != "." & (real("`M_GridCI'") < 10 | real("`M_GridCI'") > 10000)) {
    di as err "merge_study_b: manifest B/grid override is out of range"
    exit 459
}
if `M_Schema' != 1 | `M_NShard' < 1 | `M_Master' < 1 | `M_RepCap' < 0 | ///
   (`M_RepCap' > 0 & `M_RepCap' < `M_NShard' & `M_RepCap' != 1) | ///
   (`M_RepCap' == 1 & `M_NShard' > 52) | ///
   `M_ExpectedCells' != 52 | ///
   "`M_HarnessVersion'" != "`HARNESS_EXPECTED'" | ///
   "`M_CodeVersionExpected'" != "`CODE_EXPECTED'" | ///
   "`M_RngKind'" != "`RNG_EXPECTED'" {
    di as err "merge_study_b: unsupported/inconsistent manifest contract"
    exit 459
}

* Frozen-source inventory.  The verifier recomputes these values; the merger
* independently requires a complete and syntactically valid declaration.
foreach f in study_b_cells.do _registry_probe.do study_b_worker.do ///
             study_b_shard.do merge_study_b.do xtdpthresh.ado ///
             xtdpthresh_p.ado run_study_b.ps1 ///
             monitor_study_b.ps1 verify_and_merge_study_b.ps1 ///
             README.md _DESIGN_B.md {
    capture confirm file "`f'"
    if _rc {
        di as err "merge_study_b: frozen manifest member is missing: `f'"
        exit 601
    }
    local fp = subinstr("`f'", ".", "[.]", .)
    if !regexm(`"`MJ'"', `""`fp'"[ ]*:[ ]*"([0-9A-Fa-f]+)""') {
        di as err "merge_study_b: manifest lacks SHA-256 for `f'"
        exit 459
    }
    local hf "`=regexs(1)'"
    if strlen("`hf'") != 64 {
        di as err "merge_study_b: malformed SHA-256 for `f'"
        exit 459
    }
}

local R500 = cond(`M_RepCap' > 0, min(500, `M_RepCap'), 500)
local TOTAL_EFF = 52*`R500'
if `M_ExpectedReplications' != `TOTAL_EFF' {
    di as err "merge_study_b: manifest ExpectedReplications is not the v8r1 registry total"
    exit 459
}
if "`RUN_ID'" != "" & "`RUN_ID'" != "`M_RunId'" {
    di as err "merge_study_b: requested RUN_ID differs from manifest"
    exit 459
}

local files : dir "." files "study_b_results_SH*.csv"
local nfiles : word count `files'
if `nfiles' == 0 {
    di as err "merge_study_b: no result shards found"
    exit 601
}
tempfile pooled one
local first 1
foreach f of local files {
    local fl = lower("`f'")
    if !regexm("`fl'", "^study_b_results_sh([0-9]+)[.]csv$") {
        di as err "merge_study_b: unexpected result filename `f'"
        exit 459
    }
    local fsh = real(regexs(1))
    capture noisily import delimited using "`f'", clear varnames(1) ///
        case(preserve) bindquote(strict) stringcols(_all)
    if _rc | _N == 0 {
        di as err "merge_study_b: unreadable or empty result shard `f'"
        exit 459
    }
    local EXPECT_VARS run_id harness_version schema_version ///
        code_version_expected code_version_actual stata_version rng_kind ///
        rng_actual seed_scheme dgp_spec has_x infmode cell_id pair_id pairmode ///
        config_key full_key study design method N T miss missp qstatus ///
        boottype coefboot kappa R B grid gridci refine master shard nshard ///
        rep dgp_seed missing_seed boot_seed seed_threshold seed_linearity ///
        seed_continuity seed_coefficient rc fit_rc0 cmd_ok version_ok ///
        contract_ok success estimator_twostep gamma0 gamma_hat g_bias ///
        g_sqerr obj Lagyb xb qb consd qd Lagyd xd b_Lagyb b_xb b_qb b_consd ///
        b_qd b_Lagyd b_xd se_Lagyb se_xb se_qb se_consd se_qd se_Lagyd ///
        se_xd ci_requested ci_delivered ci_incomplete ci_empty n_seg ///
        boundary_warn covered set_length lin_requested lin_requested_B ///
        lin_valid_draws lin_delivered p_lin lin_reject5_cond ///
        lin_reject5_uncond cont_requested cont_requested_B cont_valid_draws ///
        cont_delivered p_cont cont_reject5_cond cont_reject5_uncond ///
        hansen_p ar1_p ar2_p N_used N_trans N_iv N_units refine_complete ///
        refine_iterations refine_added refine_remaining ///
        refine_neigh_unevaluated searchmode searchtol search_converged ///
        search_incomplete search_hit_max search_W2_builds ///
        search_stage1_level search_stage2_level search_stage1_points ///
        search_stage2_points search_stage1_same_split ///
        search_stage2_same_split search_stage1_rel_gain ///
        search_stage2_rel_gain grid_max_requested ///
        grid_requested grid_effective grid_admitted ///
        grid_structural grid_twostep_admitted gridci_requested ///
        gridci_effective gridci_admitted gridci_evaluated ci_unresolved ///
        gridboot_min_draws threshold_requested_B continuity_common_grid ///
        ci_bootstrap_certified units_target units_realized units_dropped ///
        analysis_potential analysis_observed missing_n missing_rate ///
        gap_events gap_periods fd_pair_rows_potential fod_rows_potential ///
        cb_requested cb_delivered cbcov_Lagyb cbcov_xb cbcov_qb ///
        cbcov_consd cbcov_qd cbcov_Lagyd cbcov_xd rhoy_true bx_true ///
        bq_true rhoq rhoeu_effective sige sigeta tburn maxlag_lo maxlag_hi ///
        trim history gridtype gridsample level vce elapsed_s ///
    joint_vce ar_joint vce_applied bwscale gamma_bw q_nvals_bw N_iv_dep N_iv_dep_near iv_dep_res ar1_cond ar2_cond ar1_p_cond ar2_p_cond se_cond_Lagyb se_cond_qb se_cond_consd se_cond_qd se_cond_Lagyd se_delivered ci_criterion_code citest_requested citest_returned citest_evaluable citest_gamma citest_accept citest_p citest_D citest_crit citest_status citest_draws seed_citest gammahat_in_set
    * A /// continuation keeps the continuation line's leading blanks inside the
    * macro, so the literal EXPECT_VARS text carries runs of spaces while -unab-
    * always returns single-spaced names.  Without retokenize the exact string
    * comparison below can never be true and every merge aborts on shard 1.
    local EXPECT_VARS : list retokenize EXPECT_VARS
    unab GOT_VARS : _all
    if `"`GOT_VARS'"' != `"`EXPECT_VARS'"' {
        di as err "merge_study_b: `f' has an incompatible or reordered schema"
        exit 459
    }
    local STRVARS run_id harness_version schema_version ///
        code_version_expected code_version_actual stata_version rng_kind ///
        rng_actual seed_scheme dgp_spec infmode cell_id pair_id pairmode ///
        config_key full_key study design method miss qstatus boottype ///
        coefboot searchmode history gridtype gridsample vce
    local NUMVARS : list EXPECT_VARS - STRVARS
    capture quietly destring `NUMVARS', replace
    if _rc {
        di as err "merge_study_b: `f' contains a nonnumeric token in a numeric column"
        exit 459
    }
    foreach v in run_id harness_version schema_version code_version_expected ///
        code_version_actual stata_version rng_kind rng_actual seed_scheme ///
        dgp_spec has_x infmode cell_id pair_id pairmode config_key full_key study design method ///
        N T miss missp qstatus boottype coefboot kappa R B grid gridci refine ///
        master shard nshard rep dgp_seed missing_seed boot_seed ///
        seed_threshold seed_linearity seed_continuity seed_coefficient ///
        rc fit_rc0 ///
        cmd_ok version_ok contract_ok success estimator_twostep gamma0 gamma_hat g_bias ///
        g_sqerr obj Lagyb xb qb ///
        consd qd Lagyd xd b_Lagyb b_xb b_qb b_consd b_qd b_Lagyd b_xd ///
        se_Lagyb se_xb se_qb se_consd se_qd se_Lagyd se_xd ci_requested ///
        ci_delivered ci_incomplete ci_empty n_seg boundary_warn covered ///
        set_length lin_requested lin_requested_B lin_valid_draws ///
        lin_delivered p_lin lin_reject5_cond lin_reject5_uncond ///
        cont_requested cont_requested_B cont_valid_draws cont_delivered ///
        p_cont cont_reject5_cond cont_reject5_uncond hansen_p ar1_p ar2_p ///
        N_used N_trans N_iv N_units refine_complete refine_iterations ///
        refine_added refine_remaining refine_neigh_unevaluated ///
        searchmode searchtol search_converged search_incomplete ///
        search_hit_max search_W2_builds search_stage1_level ///
        search_stage2_level search_stage1_points search_stage2_points ///
        search_stage1_same_split search_stage2_same_split ///
        search_stage1_rel_gain search_stage2_rel_gain grid_max_requested ///
        grid_requested grid_effective grid_admitted grid_structural ///
        grid_twostep_admitted gridci_requested gridci_effective ///
        gridci_admitted gridci_evaluated ci_unresolved gridboot_min_draws ///
        threshold_requested_B continuity_common_grid ci_bootstrap_certified ///
        units_target units_realized units_dropped analysis_potential ///
        analysis_observed missing_n missing_rate gap_events gap_periods ///
        fd_pair_rows_potential fod_rows_potential ///
        cb_requested cb_delivered cbcov_Lagyb ///
        cbcov_xb cbcov_qb cbcov_consd cbcov_qd cbcov_Lagyd cbcov_xd ///
        rhoy_true bx_true bq_true rhoq rhoeu_effective sige sigeta tburn ///
        maxlag_lo maxlag_hi trim history gridtype gridsample level vce elapsed_s {
        capture confirm variable `v'
        if _rc {
            di as err "merge_study_b: `f' lacks mandatory column `v'"
            exit 459
        }
    }
    foreach v in run_id harness_version schema_version code_version_expected ///
        code_version_actual stata_version rng_kind rng_actual seed_scheme ///
        dgp_spec infmode cell_id pair_id pairmode config_key full_key study design method miss ///
        qstatus boottype coefboot searchmode history gridtype gridsample vce {
        capture confirm string variable `v'
        if _rc {
            di as err "merge_study_b: `f' has wrong type for string column `v'"
            exit 459
        }
    }
    foreach v of local NUMVARS {
        capture confirm numeric variable `v'
        if _rc {
            di as err "merge_study_b: `f' has wrong type for numeric column `v'"
            exit 459
        }
    }
    foreach v in has_x N T missp kappa R B grid gridci refine master shard ///
        nshard rep dgp_seed missing_seed boot_seed seed_threshold ///
        seed_linearity seed_continuity seed_coefficient rc fit_rc0 cmd_ok ///
        version_ok contract_ok success estimator_twostep gamma0 gamma_hat ///
        g_bias g_sqerr obj N_units refine_complete refine_iterations ///
        refine_added refine_remaining refine_neigh_unevaluated ///
        searchtol search_converged search_incomplete search_hit_max ///
        search_W2_builds search_stage1_level search_stage2_level ///
        search_stage1_points search_stage2_points ///
        search_stage1_same_split search_stage2_same_split ///
        search_stage1_rel_gain search_stage2_rel_gain grid_max_requested ///
        grid_requested grid_effective grid_admitted grid_structural ///
        grid_twostep_admitted gridci_requested gridci_effective ///
        gridci_admitted gridci_evaluated ci_unresolved gridboot_min_draws ///
        threshold_requested_B continuity_common_grid ci_bootstrap_certified ///
        rhoy_true bx_true bq_true rhoq rhoeu_effective sige sigeta tburn ///
        maxlag_lo maxlag_hi trim level elapsed_s {
        capture confirm numeric variable `v'
        if _rc {
            di as err "merge_study_b: `f' has wrong type for numeric column `v'"
            exit 459
        }
    }
    quietly count if shard != `fsh'
    if r(N) {
        di as err "merge_study_b: filename/data shard mismatch in `f'"
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
                 stata_version rng_kind rng_actual seed_scheme master nshard {
    capture assert `v' == `v'[1]
    if _rc {
        di as err "merge_study_b: mixed run provenance in `v'"
        exit 459
    }
}
if run_id[1] == "" | harness_version[1] != "`HARNESS_EXPECTED'" | ///
   schema_version[1] != "`SCHEMA_EXPECTED'" | ///
   code_version_expected[1] != "`CODE_EXPECTED'" | ///
   rng_kind[1] != "`RNG_EXPECTED'" | rng_actual[1] != "`RNG_EXPECTED'" | ///
   seed_scheme[1] != "`SEEDS_EXPECTED'" {
    di as err "merge_study_b: wrong/blank harness, schema, or seed scheme"
    exit 459
}
if "`RUN_ID'" != "" & run_id[1] != "`RUN_ID'" {
    di as err "merge_study_b: requested RUN_ID differs from data"
    exit 459
}
local ACT_RUN "`=run_id[1]'"
local ACT_HARNESS "`=harness_version[1]'"
local ACT_CODE "`=code_version_expected[1]'"
local ACT_MASTER = master[1]
local ACT_NSHARD = nshard[1]
local ACT_ROWS = _N
if "`ACT_RUN'" != "`M_RunId'" | `ACT_NSHARD' != `M_NShard' | ///
   `ACT_MASTER' != `M_Master' | ///
   "`ACT_HARNESS'" != "`M_HarnessVersion'" | ///
   "`ACT_CODE'" != "`M_CodeVersionExpected'" {
    di as err "merge_study_b: results differ from immutable manifest config"
    exit 459
}
foreach q in B Grid GridCI {
    if "`M_`q''" != "." {
        local column = cond("`q'"=="B","B",lower("`q'"))
        quietly count if `column' != real("`M_`q''")
        if r(N) {
            di as err "merge_study_b: `q' differs from immutable manifest override"
            exit 459
        }
    }
}
capture assert `ACT_NSHARD' >= 1 & shard >= 1 & shard <= `ACT_NSHARD'
if _rc | `nfiles' != `ACT_NSHARD' {
    di as err "merge_study_b: result-file count/range differs from nshard"
    exit 459
}
forvalues k = 1/`ACT_NSHARD' {
    local exact : dir "." files "study_b_results_SH`k'.csv"
    if `: word count `exact'' != 1 {
        di as err "merge_study_b: missing result shard `k'"
        exit 601
    }
}

capture isid full_key
if _rc {
    di as err "merge_study_b: duplicate full_key"
    exit 459
}
capture isid cell_id rep
if _rc {
    di as err "merge_study_b: duplicate cell_id/rep"
    exit 459
}
gen str244 __key = run_id + "|" + cell_id + "|rep=" + ///
    trim(string(rep, "%12.0f"))
quietly count if full_key != __key | cell_id != pair_id + "_" + method
if r(N) {
    di as err "merge_study_b: malformed primary/registry key"
    exit 459
}
drop __key
capture assert rep >= 1 & rep <= R & R == floor(R) & ///
    ((`M_RepCap' == 1 & nshard > 1) | shard == ceil(rep*nshard/R))
if _rc {
    di as err "merge_study_b: row outside shard allocation"
    exit 459
}
sort cell_id rep
capture by cell_id: assert rep == _n & _N == R
if _rc {
    di as err "merge_study_b: incomplete replication sequence"
    exit 459
}
foreach v in fit_rc0 cmd_ok contract_ok success ci_requested ///
                 ci_delivered lin_requested lin_delivered cont_requested ///
                 cont_delivered cb_requested cb_delivered citest_requested ///
                 citest_returned citest_evaluable {
    capture assert inlist(`v', 0, 1)
    if _rc {
        di as err "merge_study_b: nonbinary flag `v'"
        exit 459
    }
}
capture assert fit_rc0 == (rc == 0) & cmd_ok <= fit_rc0 & ///
    contract_ok <= cmd_ok & contract_ok <= version_ok
if _rc {
    di as err "merge_study_b: rc/command/contract flags disagree"
    exit 459
}
capture assert inlist(version_ok, 0, 1, .) & ///
    (!missing(version_ok) | fit_rc0 == 0)
if _rc {
    di as err "merge_study_b: invalid version_ok accounting"
    exit 459
}
capture assert inlist(estimator_twostep, 0, 1, .) & ///
    (!missing(estimator_twostep) | fit_rc0 == 0) & ///
    success == (contract_ok == 1 & estimator_twostep == 1) & ///
    ci_delivered <= ci_requested & ci_delivered <= contract_ok & ///
    lin_delivered <= lin_requested & lin_delivered <= contract_ok & ///
    cont_delivered <= cont_requested & cont_delivered <= contract_ok & ///
    cb_delivered <= cb_requested & cb_delivered <= contract_ok
if _rc {
    di as err "merge_study_b: inconsistent success/delivery flags"
    exit 459
}
capture assert ci_requested == (infmode != "POINT") & ///
    lin_requested == (infmode == "FULL") & ///
    cont_requested == (infmode == "FULL") & ///
    cb_requested == (coefboot != "none") & ///
    citest_requested == (study != "LIN" & infmode != "POINT") & ///
    lin_requested_B == cond(lin_requested, B, 0) & ///
    cont_requested_B == cond(cont_requested, B, 0)
if _rc {
    di as err "merge_study_b: requested-inference flags disagree with registry mode"
    exit 459
}
capture assert citest_returned <= citest_requested & ///
    citest_evaluable <= citest_returned & ///
    (citest_requested == 0 | abs(citest_gamma-gamma0) <= ///
        1e-12*max(1,abs(gamma0))) & ///
    (citest_returned == 0 | inrange(citest_status,1,6)) & ///
    (citest_evaluable == inlist(citest_status,1,2))
if _rc {
    di as err "merge_study_b: citest request/delivery accounting is inconsistent"
    exit 459
}
quietly count if citest_requested == 1 & fit_rc0 == 1 & ///
    (citest_returned != 1 | ///
     seed_citest != mod(boot_seed+477377,2147483648) | ///
     (citest_status==1 & ///
        (citest_draws!=B | missing(citest_crit,citest_accept,citest_p,citest_D) | ///
         !inlist(citest_accept,0,1) | !inrange(citest_p,0,1) | ///
         citest_accept!=(citest_D<=citest_crit) | ///
         citest_accept!=(citest_p>.05))) | ///
     (citest_status==2 & ///
        (citest_D!=0 | citest_accept!=1 | citest_p!=1 | ///
         !missing(citest_crit) | !missing(citest_draws))) | ///
     (inrange(citest_status,3,6) & citest_evaluable!=0))
if r(N) {
    di as err "merge_study_b: citest value/status/seed contract failed"
    exit 459
}
quietly count if citest_requested == 0 & ///
    (citest_returned!=0 | citest_evaluable!=0 | ///
     !missing(citest_gamma,citest_accept,citest_p,citest_D,citest_crit, ///
              citest_status,citest_draws,seed_citest))
if r(N) {
    di as err "merge_study_b: unrequested citest returned state"
    exit 459
}
capture assert (missing(lin_valid_draws) | ///
        (inrange(lin_valid_draws, 0, lin_requested_B) & ///
         lin_valid_draws == floor(lin_valid_draws))) & ///
    (missing(cont_valid_draws) | ///
        (inrange(cont_valid_draws, 0, cont_requested_B) & ///
         cont_valid_draws == floor(cont_valid_draws))) & ///
    (missing(p_lin) | inrange(p_lin, 0, 1)) & ///
    (missing(p_cont) | inrange(p_cont, 0, 1)) & ///
    lin_delivered == !missing(p_lin) & cont_delivered == !missing(p_cont) & ///
    (lin_delivered == 0 | inrange(lin_valid_draws, 1, lin_requested_B)) & ///
    (cont_delivered == 0 | inrange(cont_valid_draws, 1, cont_requested_B))
if _rc {
    di as err "merge_study_b: invalid bootstrap-test delivery/accounting"
    exit 459
}
capture assert inlist(ci_incomplete, 0, 1, .) & ///
    inlist(ci_empty, 0, 1, .) & ///
    (missing(boundary_warn) | inlist(boundary_warn, 0, 1, 2, 3)) & ///
    inlist(covered, 0, 1, .) & ///
    ((ci_delivered == 1 & !missing(gamma0) & inlist(covered, 0, 1)) | ///
     ((ci_delivered == 0 | missing(gamma0)) & missing(covered)))
if _rc {
    di as err "merge_study_b: invalid confidence-set outcome accounting"
    exit 459
}
foreach z in lin cont {
    capture assert inlist(`z'_reject5_cond, 0, 1, .) & ///
        inlist(`z'_reject5_uncond, 0, 1, .)
    if _rc {
        di as err "merge_study_b: nonbinary `z' rejection outcome"
        exit 459
    }
    quietly count if ///
        (`z'_requested == 0 & ///
         (`z'_delivered != 0 | !missing(`z'_valid_draws) | ///
          !missing(p_`z') | !missing(`z'_reject5_cond) | ///
          !missing(`z'_reject5_uncond))) | ///
        (`z'_requested == 1 & `z'_delivered == 0 & ///
         (!missing(p_`z') | !missing(`z'_reject5_cond) | ///
          `z'_reject5_uncond != 0)) | ///
        (`z'_delivered == 1 & ///
         (missing(`z'_reject5_cond, `z'_reject5_uncond) | ///
          `z'_reject5_cond != (p_`z' < .05) | ///
          `z'_reject5_uncond != `z'_reject5_cond))
    if r(N) {
        di as err "merge_study_b: `z' p-value/rejection flags disagree"
        exit 459
    }
}
capture assert ci_incomplete != 1 | ci_delivered == 0
if _rc {
    di as err "merge_study_b: an incomplete inversion was reported as delivered"
    exit 459
}
capture assert ci_empty != 1 | ci_delivered == 0
if _rc {
    di as err "merge_study_b: an empty inversion was reported as delivered"
    exit 459
}
capture assert ci_delivered == 0 | ///
    (ci_incomplete == 0 & ci_empty == 0 & n_seg >= 1 & ///
     !missing(set_length) & set_length >= 0)
if _rc {
    di as err "merge_study_b: delivered-CI fields are internally inconsistent"
    exit 459
}
quietly count if contract_ok == 1 & ci_requested == 1 & ///
    (missing(ci_incomplete, ci_unresolved, boundary_warn) | ///
     ci_incomplete != (ci_unresolved > 0) | ///
     (ci_incomplete == 0 & !inlist(ci_empty, 0, 1)) | ///
     (ci_incomplete == 1 & !missing(ci_empty)))
if r(N) {
    di as err "merge_study_b: confidence-inversion completeness is inconsistent"
    exit 459
}
capture assert missing(refine_complete) | inlist(refine_complete, 0, 1)
if _rc {
    di as err "merge_study_b: invalid refinement accounting"
    exit 459
}
quietly count if contract_ok == 1 & ///
    (missing(N_units, refine_complete, refine_iterations, refine_added, ///
             refine_remaining, refine_neigh_unevaluated, grid_requested, ///
             grid_effective, grid_admitted, grid_structural, ///
             search_hit_max, search_W2_builds, search_stage1_level, ///
             search_stage2_level, search_stage1_points, ///
             search_stage2_points, grid_max_requested, ///
             gridci_requested, gridci_effective, gridci_admitted, ///
             gridci_evaluated, threshold_requested_B, ///
             ci_bootstrap_certified) | ///
     N_units < 1 | refine != cond(study=="GSCAL",0,4) | ///
     !inrange(refine_iterations,0,refine) | refine_added < 0 | ///
     refine_remaining < 0 | refine_neigh_unevaluated < 0 | ///
     !inlist(refine_complete,0,1) | ///
     (refine_complete == 1 & ///
         (refine_remaining != 0 | refine_neigh_unevaluated != 0)) | ///
     searchmode != "fixed" | !missing(searchtol) | ///
     !missing(search_converged) | !missing(search_incomplete) | ///
     search_hit_max != 0 | search_W2_builds != 1 | ///
     search_stage1_level != 1 | !inlist(search_stage2_level,0,1) | ///
     search_stage1_points != grid | ///
     (search_stage2_level==0 & search_stage2_points!=0) | ///
     (search_stage2_level==1 & search_stage2_points!=grid) | ///
     !missing(search_stage1_same_split) | ///
     !missing(search_stage2_same_split) | ///
     !missing(search_stage1_rel_gain) | ///
     !missing(search_stage2_rel_gain) | ///
     grid_requested != grid | gridci_requested != gridci | ///
     grid_max_requested != grid | grid_effective != grid+refine_added | ///
     grid_admitted < 2 | grid_structural < grid_admitted | ///
     grid_structural > grid_effective | ///
     (estimator_twostep==1 & ///
         (missing(grid_twostep_admitted) | ///
          !inrange(grid_twostep_admitted,2,grid_effective))) | ///
     (estimator_twostep==0 & !missing(grid_twostep_admitted) & ///
         !inrange(grid_twostep_admitted,0,grid_effective)) | ///
     gridci_effective < 1 | gridci_admitted < 0 | gridci_evaluated < 0 | ///
     gridci_evaluated > gridci_admitted | ///
     threshold_requested_B != B | ci_bootstrap_certified != 0 | ///
     (infmode=="FULL" & (missing(continuity_common_grid) | ///
                         !inrange(continuity_common_grid,0,grid_effective))) | ///
     (!missing(gridboot_min_draws) & !inrange(gridboot_min_draws, 0, B)) | ///
     (ci_delivered == 1 & ///
         (missing(gridboot_min_draws) | gridboot_min_draws < 1 | ///
          gridboot_min_draws != floor(gridboot_min_draws))))
if r(N) {
    di as err "merge_study_b: estimator grid/unit diagnostics violate their contract"
    exit 459
}
quietly count if rc == 0 & code_version_actual != code_version_expected
if r(N) {
    di as err "merge_study_b: successful fit used wrong ado version"
    exit 459
}

* Recompute all worker seed hashes from recorded design fields.  METHOD is
* intentionally absent, yielding paired outer samples and synchronized
* bootstrap seeds; tuning fields prevent cross-block inner-stream reuse.
gen byte __dcode = cond(design=="kink",1,cond(design=="jump",2, ///
    cond(design=="weakkink",3,4)))
gen byte __mcode = cond(miss=="balanced",1,cond(miss=="mcar",2,3))
gen byte __qcode = cond(qstatus=="pred",1,2)
gen byte __spcode = cond(inlist(dgp_spec,"gs26_base_v2","gs26_official_v1"),1, ///
    cond(dgp_spec=="gs26_endog_v2",2,3))
gen byte __btcode = cond(boottype=="wild",1,2)
gen byte __cbcode = cond(coefboot=="none",1,cond(coefboot=="onestep",2,3))
gen double __kcode = round((kappa+5)*1000)
gen double __pcode = round(missp*1000000)
gen double __h0 = mod(master*1000003 + rep*10007 + N*101 + T*1009 + ///
    __spcode*1000033 + __qcode*1000037, 699999999)
gen double __dgp_seed = 1 + __h0
gen double __hm = mod(__h0*1009 + __mcode*10007 + __pcode*101 + 7919, ///
    699999999)
gen double __missing_seed = 700000001 + __hm
gen double __hb = mod(__h0*1013 + __dcode*100003 + __mcode*10007 + ///
    __pcode*101 + __qcode*1009 + __btcode*17 + __cbcode*29 + ///
    __kcode*37 + B*131 + grid*137 + gridci*139 + refine*149 + ///
    maxlag_lo*151 + maxlag_hi*157 + 104729, 699999999)
gen double __boot_seed = 1400000001 + __hb
quietly count if dgp_seed != __dgp_seed | missing_seed != __missing_seed | ///
    boot_seed != __boot_seed
if r(N) {
    di as err "merge_study_b: recorded worker seed differs from frozen formula"
    exit 459
}
quietly count if fit_rc0 == 1 & ///
    (seed_threshold != boot_seed | ///
     (infmode=="FULL" & seed_linearity != mod(boot_seed+104729,2147483648)) | ///
     (infmode=="FULL" & seed_continuity != mod(boot_seed+224737,2147483648)) | ///
     (infmode!="FULL" & (!missing(seed_linearity) | !missing(seed_continuity))) | ///
     !missing(seed_coefficient) | ///
     (citest_requested==1 & seed_citest!=mod(boot_seed+477377,2147483648)) | ///
     (citest_requested==0 & !missing(seed_citest)))
if r(N) {
    di as err "merge_study_b: component-specific bootstrap seeds are inconsistent"
    exit 459
}
drop __dcode __mcode __qcode __spcode __btcode __cbcode __kcode __pcode ///
    __h0 __dgp_seed __hm __missing_seed __hb __boot_seed

* Published-benchmark DGP and estimator contract.  B1--B5 use the main
* lag-capped design; B6 is the registered FD-only calibration bridge.
quietly count if has_x != 0 | abs(rhoy_true-.6) > 1e-12 | ///
    !missing(bx_true) | abs(bq_true-1) > 1e-12 | ///
    abs(rhoq-.7) > 1e-12 | abs(rhoeu_effective-.5) > 1e-12 | ///
    abs(sige-.5) > 1e-12 | sigeta != 0 | tburn != 20 | ///
    history != "panel" | level != 95 | vce != "robust" | ///
    boottype != "wild" | coefboot != "none" | ///
    !inlist(pairmode, "paired", "fdonly") | ///
    (pairmode == "paired" & (maxlag_lo != 1 | maxlag_hi != 3 | ///
        abs(trim-.15) > 1e-12 | gridtype != "uniform" | ///
        gridsample != "effective")) | ///
    (pairmode == "fdonly" & (study != "GSCAL" | method != "fd" | ///
        maxlag_lo != 1 | maxlag_hi != 5 | abs(trim-.10) > 1e-12 | ///
        gridtype != "quantile" | gridsample != "observed"))
if r(N) {
    di as err "merge_study_b: DGP/estimator constants differ from the frozen design"
    exit 459
}
quietly count if ///
    (dgp_spec=="gs26_base_v2" & qstatus!="pred") | ///
    (dgp_spec=="gs26_endog_v2" & qstatus!="endog") | ///
    (study=="ENDOG" & dgp_spec!="gs26_endog_v2") | ///
    (study!="ENDOG" & dgp_spec!="gs26_base_v2") | ///
    (design=="none" & (!missing(gamma0) | qd!=0 | consd!=0)) | ///
    (design!="none" & (abs(gamma0-.25)>1e-12 | qd!=2 | ///
                       abs(consd-(kappa-.5))>1e-12)) | ///
    abs(Lagyb-.6)>1e-12 | abs(qb-1)>1e-12 | Lagyd!=0 | ///
    !missing(xb) | !missing(xd)
if r(N) {
    di as err "merge_study_b: row-level true values do not match its registered DGP"
    exit 459
}
quietly count if study == "GSCAL" & contract_ok == 1 & N_iv != 28
if r(N) {
    di as err "merge_study_b: contract-valid calibration row does not have 28 instruments"
    exit 459
}
local ACT_RCAP = `M_RepCap'
quietly save `pooled', replace

* Dry-evaluate the authoritative registry under recorded global overrides.
tempfile expected
tempname PP
postfile `PP' str80 cell_id str72 pair_id str12 study str12 design ///
    str4 method int N int T str12 miss double missp str8 qstatus long R ///
    double kappa str24 dgp_spec str8 boottype str10 coefboot int grid ///
    int gridci long B byte refine str8 infmode byte maxlag_lo byte maxlag_hi ///
    double trim str10 gridtype str10 gridsample str8 pairmode ///
    str244 config_key using `expected', replace
global MCB_POST "`PP'"
global MCB_BOV "`M_B'"
global MCB_GOV "`M_Grid'"
global MCB_GCOV "`M_GridCI'"
global MCB_CAP `ACT_RCAP'
capture program drop _cell
program define _cell
    version 15.0
    args PAIR_ID STUDY DESIGN METHOD N T MISS MISSP QSTATUS R KAPPA DGP_SPEC ///
         BOOTTYPE COEFBOOT GRID GRIDCI B REFINE INFMODE ///
         MAXLAG_LO MAXLAG_HI TRIM GRIDTYPE GRIDSAMPLE PAIRMODE
    if "`INFMODE'" == "" local INFMODE FULL
    local CELL_ID "`PAIR_ID'_`METHOD'"
    local RR = `R'
    if ${MCB_CAP} > 0 & `RR' > ${MCB_CAP} local RR = ${MCB_CAP}
    local BB = `B'
    local GG = `GRID'
    local GC = `GRIDCI'
    if "${MCB_BOV}" != "." local BB = ${MCB_BOV}
    if "${MCB_GOV}" != "." local GG = ${MCB_GOV}
    if "${MCB_GCOV}" != "." local GC = ${MCB_GCOV}
    local CK "`STUDY'|`DESIGN'|`METHOD'|N=`N'|T=`T'|`MISS'=`MISSP'|`QSTATUS'|R=`RR'|B=`BB'|G=`GG'|GC=`GC'|RF=`REFINE'|`BOOTTYPE'|`COEFBOOT'|K=`KAPPA'|`DGP_SPEC'|INF=`INFMODE'|ML=`MAXLAG_LO':`MAXLAG_HI'|TR=`TRIM'|GT=`GRIDTYPE'|GS=`GRIDSAMPLE'|PM=`PAIRMODE'"
    post ${MCB_POST} ("`CELL_ID'") ("`PAIR_ID'") ("`STUDY'") ///
        ("`DESIGN'") ("`METHOD'") (`N') (`T') ("`MISS'") (`MISSP') ///
        ("`QSTATUS'") (`RR') (`KAPPA') ("`DGP_SPEC'") ("`BOOTTYPE'") ///
        ("`COEFBOOT'") (`GG') (`GC') (`BB') (`REFINE') ("`INFMODE'") ///
        (`MAXLAG_LO') (`MAXLAG_HI') (`TRIM') ("`GRIDTYPE'") ///
        ("`GRIDSAMPLE'") ("`PAIRMODE'") ("`CK'")
end
capture noisily do study_b_cells.do
local registry_rc = _rc
postclose `PP'
capture program drop _cell
macro drop MCB_POST MCB_BOV MCB_GOV MCB_GCOV MCB_CAP
if `registry_rc' {
    di as err "merge_study_b: registry evaluation failed"
    exit `registry_rc'
}
use `expected', clear
gen long cell_order = _n
capture isid cell_id
if _rc | _N != 52 {
    di as err "merge_study_b: registry is not exactly 52 unique cells"
    exit 459
}
sort pair_id method
capture by pair_id: assert pairmode == pairmode[1]
if _rc {
    di as err "merge_study_b: mixed pair modes within a registry ID"
    exit 459
}
capture by pair_id: assert (_N == 2 & method[1] == "fd" & ///
    method[2] == "fod") if pairmode == "paired"
if _rc {
    di as err "merge_study_b: incomplete paired FD/FOD registry entry"
    exit 459
}
capture by pair_id: assert (_N == 1 & method == "fd") if pairmode == "fdonly"
if _rc {
    di as err "merge_study_b: malformed FD-only calibration registry entry"
    exit 459
}
bysort pair_id: gen byte __pair_tag = (_n == 1)
quietly count if __pair_tag & pairmode == "paired"
local EXPECTED_PAIRS = r(N)
quietly count if __pair_tag & pairmode == "fdonly"
local EXPECTED_FDONLY = r(N)
drop __pair_tag
if `EXPECTED_PAIRS' != 24 | `EXPECTED_FDONLY' != 4 {
    di as err "merge_study_b: registry is not 24 FD/FOD pairs plus 4 FD-only cells"
    exit 459
}
* Independently mirror the frozen registry contract instead of treating the
* registry file as its own oracle.  Numeric smoke/sensitivity overrides remain
* explicit in the manifest and are applied uniformly where allowed.
local ER = cond(`M_RepCap'>0, min(500,`M_RepCap'), 500)
local EBMAIN = cond("`M_B'"==".",499,real("`M_B'"))
local EBCAL  = cond("`M_B'"==".",500,real("`M_B'"))
local EGMAIN = cond("`M_Grid'"==".",199,real("`M_Grid'"))
local EGCAL  = cond("`M_Grid'"==".",46,real("`M_Grid'"))
local EGCMAIN = cond("`M_GridCI'"==".",100,real("`M_GridCI'"))
local EGCLIN  = cond("`M_GridCI'"==".",10,real("`M_GridCI'"))
local EGCCAL  = cond("`M_GridCI'"==".",46,real("`M_GridCI'"))
gen byte block = real(substr(cell_id,2,1))
capture assert T==6 & inlist(N,400,800) & R==`ER' & ///
    boottype=="wild" & coefboot=="none" & ///
    inlist(infmode,"FULL","COVERAGE") & ///
    ((miss=="balanced" & missp==0) | (miss=="mcar" & missp==.30) | ///
     (miss=="attrition" & missp==.15)) & ///
    ((study=="LIN" & design=="none" & kappa==0) | ///
     (study!="LIN" & kappa==0 & design=="kink") | ///
     (study!="LIN" & kappa>0 & design=="jump")) & inrange(block,1,6)
if _rc {
    di as err "merge_study_b: registry violates the frozen global design"
    exit 459
}
capture assert block>5 | (pairmode=="paired" & B==`EBMAIN' & ///
    grid==`EGMAIN' & refine==4 & maxlag_lo==1 & maxlag_hi==3 & ///
    abs(trim-.15)<1e-12 & gridtype=="uniform" & gridsample=="effective")
if _rc {
    di as err "merge_study_b: B1--B5 estimator contract differs from the frozen design"
    exit 459
}
capture assert block!=1 | (N==400 & inlist(study,"GSBAL","GSGAP") & ///
    qstatus=="pred" & inlist(kappa,0,.1,.5,1) & ///
    dgp_spec=="gs26_base_v2" & gridci==`EGCMAIN' & infmode=="FULL")
if _rc exit 459
capture assert block!=2 | (N==800 & inlist(study,"GSBAL","GSGAP") & ///
    qstatus=="pred" & inlist(kappa,0,.1,1) & dgp_spec=="gs26_base_v2" & ///
    gridci==`EGCMAIN' & ((kappa==0 & infmode=="FULL") | ///
                         (kappa>0 & infmode=="COVERAGE")))
if _rc exit 459
capture assert block!=3 | (study=="LIN" & inlist(N,400,800) & ///
    qstatus=="pred" & dgp_spec=="gs26_base_v2" & ///
    gridci==`EGCLIN' & infmode=="FULL")
if _rc exit 459
capture assert block!=4 | (study=="ENDOG" & N==400 & ///
    inlist(kappa,0,1) & qstatus=="endog" & dgp_spec=="gs26_endog_v2" & ///
    gridci==`EGCMAIN' & infmode=="FULL")
if _rc exit 459
capture assert block!=5 | (study=="ATTR" & N==400 & miss=="attrition" & ///
    missp==.15 & qstatus=="pred" & inlist(kappa,0,1) & ///
    dgp_spec=="gs26_base_v2" & gridci==`EGCMAIN' & infmode=="FULL")
if _rc exit 459
capture assert block!=6 | (study=="GSCAL" & method=="fd" & ///
    pairmode=="fdonly" & inlist(N,400,800) & miss=="balanced" & ///
    missp==0 & qstatus=="pred" & inlist(kappa,0,1) & ///
    dgp_spec=="gs26_base_v2" & B==`EBCAL' & grid==`EGCAL' & ///
    gridci==`EGCCAL' & refine==0 & infmode=="COVERAGE" & ///
    maxlag_lo==1 & maxlag_hi==5 & abs(trim-.10)<1e-12 & ///
    gridtype=="quantile" & gridsample=="observed")
if _rc {
    di as err "merge_study_b: B6 calibration registry differs from its four frozen targets"
    exit 459
}
quietly count if regexm(cell_id, "^b1_")
local NB1 = r(N)
quietly count if regexm(cell_id, "^b2_")
local NB2 = r(N)
quietly count if regexm(cell_id, "^b3_")
local NB3 = r(N)
quietly count if regexm(cell_id, "^b4_")
local NB4 = r(N)
quietly count if regexm(cell_id, "^b5_")
local NB5 = r(N)
quietly count if regexm(cell_id, "^b6_")
local NB6 = r(N)
if `NB1' != 16 | `NB2' != 12 | `NB3' != 8 | `NB4' != 8 | ///
   `NB5' != 4 | `NB6' != 4 {
    di as err "merge_study_b: registry block counts differ from 16/12/8/8/4/4"
    exit 459
}
quietly count if infmode=="FULL"
if r(N)!=40 {
    di as err "merge_study_b: registry does not contain exactly 40 FULL cells"
    exit 459
}
quietly count if infmode=="COVERAGE"
if r(N)!=12 {
    di as err "merge_study_b: registry does not contain exactly 12 COVERAGE cells"
    exit 459
}
foreach kk in 0 .1 .5 1 {
    quietly count if regexm(cell_id, "^b1_") & abs(kappa-`kk') < 1e-12
    if r(N) != 4 {
        di as err "merge_study_b: B1 lacks a complete kappa=`kk' FD/FOD-by-missingness anchor"
        exit 459
    }
}
egen double __expected_total = total(R)
local EXPECTED_TOTAL = __expected_total[1]
drop __expected_total
egen double __expected_paired = total(cond(pairmode=="paired", R, 0))
local EXPECTED_PAIRED_ROWS = __expected_paired[1]
drop __expected_paired
if `EXPECTED_PAIRED_ROWS' != 48*`R500' {
    di as err "merge_study_b: paired-registry row contract is inconsistent"
    exit 459
}
if `EXPECTED_TOTAL' != `M_ExpectedReplications' | ///
   _N != `M_ExpectedCells' | ///
   `ACT_ROWS' != `EXPECTED_TOTAL' {
    di as err "merge_study_b: expected-replication contract mismatch"
    di as err "manifest=`M_ExpectedReplications'; registry=`EXPECTED_TOTAL'; observed=`ACT_ROWS'"
    exit 459
}
quietly save `expected', replace

* The dedicated RepCap=1 smoke assigns whole registry cells round-robin.
* Verify the exact assignment before any summaries are materialized.
tempfile assignment
preserve
keep cell_id cell_order
quietly save `assignment', replace
restore
use `pooled', clear
merge m:1 cell_id using `assignment'
capture assert _merge == 3
if _rc {
    di as err "merge_study_b: result cell is absent from assignment registry"
    exit 459
}
if `M_RepCap' == 1 & `ACT_NSHARD' > 1 {
    capture assert shard == mod(cell_order-1, `ACT_NSHARD')+1
    if _rc {
        di as err "merge_study_b: parallel-smoke cell/shard assignment mismatch"
        exit 459
    }
}
drop _merge cell_order
quietly save `pooled', replace

use `pooled', clear
sort cell_id rep
foreach v in pair_id config_key study design method N T miss missp qstatus R ///
                 kappa dgp_spec boottype coefboot grid gridci B refine infmode ///
                 maxlag_lo maxlag_hi trim gridtype gridsample pairmode {
    capture by cell_id: assert `v' == `v'[1]
    if _rc {
        di as err "merge_study_b: cell_id reused with mixed `v'"
        exit 459
    }
}
by cell_id: gen long n_actual = _N
by cell_id: keep if _n == 1
keep cell_id pair_id config_key study design method N T miss missp qstatus R ///
    kappa dgp_spec boottype coefboot grid gridci B refine infmode ///
    maxlag_lo maxlag_hi trim gridtype gridsample pairmode n_actual
tempfile actual_cells
quietly save `actual_cells', replace
use `expected', clear
foreach v in pair_id config_key study design method N T miss missp qstatus R ///
                 kappa dgp_spec boottype coefboot grid gridci B refine infmode ///
                 maxlag_lo maxlag_hi trim gridtype gridsample pairmode {
    rename `v' exp_`v'
}
merge 1:1 cell_id using `actual_cells'
quietly count if _merge != 3
if r(N) {
    di as err "merge_study_b: missing or extra registry cells"
    list cell_id _merge if _merge != 3, noobs abbrev(32)
    exit 459
}
drop _merge
quietly count if pair_id != exp_pair_id | config_key != exp_config_key | ///
    study != exp_study | design != exp_design | method != exp_method | ///
    N != exp_N | T != exp_T | miss != exp_miss | ///
    abs(missp-exp_missp)>1e-6 | qstatus != exp_qstatus | R != exp_R | ///
    abs(kappa-exp_kappa)>1e-6 | dgp_spec != exp_dgp_spec | ///
    boottype != exp_boottype | coefboot != exp_coefboot | ///
    grid != exp_grid | gridci != exp_gridci | B != exp_B | ///
    refine != exp_refine | infmode != exp_infmode | ///
    maxlag_lo != exp_maxlag_lo | maxlag_hi != exp_maxlag_hi | ///
    abs(trim-exp_trim)>1e-12 | gridtype != exp_gridtype | ///
    gridsample != exp_gridsample | pairmode != exp_pairmode | ///
    n_actual != exp_R
if r(N) {
    di as err "merge_study_b: cell arguments/counts differ from registry"
    exit 459
}

* Completion-marker contract and exact per-shard counts.
local done : dir "." files "study_b_complete_SH*.csv"
if `: word count `done'' != `ACT_NSHARD' {
    di as err "merge_study_b: completion-marker count differs from nshard"
    exit 459
}
tempfile markers markerone shard_contract
local first 1
forvalues k = 1/`ACT_NSHARD' {
    local mfile : dir "." files "study_b_complete_SH`k'.csv"
    if `: word count `mfile'' != 1 {
        di as err "merge_study_b: missing completion marker `k'"
        exit 601
    }
    local mf : word 1 of `mfile'
    quietly import delimited using "`mf'", clear ///
        varnames(1) case(preserve) bindquote(strict)
    local MARKER_VARS run_id harness_version code_version_expected master ///
        shard nshard cells reps exit_code
    local MARKER_VARS : list retokenize MARKER_VARS
    unab GOT_MARKER : _all
    if `"`GOT_MARKER'"' != `"`MARKER_VARS'"' | _N != 1 {
        di as err "merge_study_b: marker `k' has wrong schema or row count"
        exit 459
    }
    if shard[1] != `k' | exit_code[1] != 0 | ///
       run_id[1] != "`ACT_RUN'" | harness_version[1] != "`ACT_HARNESS'" | ///
       code_version_expected[1] != "`ACT_CODE'" | ///
       master[1] != `ACT_MASTER' | nshard[1] != `ACT_NSHARD' {
        di as err "merge_study_b: invalid/mixed completion marker `k'"
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
if _rc exit 459
drop _merge
merge 1:1 shard using `markers'
capture assert _merge == 3
if _rc exit 459
drop _merge
quietly count if actual_reps != expected_reps | actual_cells != expected_cells | ///
    reps != expected_reps | cells != expected_cells
if r(N) {
    di as err "merge_study_b: marker/result counts violate shard contract"
    exit 459
}

* Only now are pooled data and downstream reports materialized.
tempfile OUTALL OUTALLCSV OUTPAIR OUTPAIRCSV OUTPAIRSUM OUTPAIRSUMCSV ///
    OUTSUMMARY OUTSUMMARYCSV
use `pooled', clear
sort cell_id rep
compress
* Added 0.9.36 telemetry contract (does not select point estimates).
capture assert inlist(joint_vce,0,1) & inlist(ar_joint,0,1) & ///
    inlist(vce_applied,0,1) & abs(bwscale-1.5)<1e-12 if contract_ok
if _rc exit 459
foreach vv in se_cond_Lagyb se_cond_qb se_cond_consd se_cond_qd se_cond_Lagyd {
    capture assert `vv' >= 0 | missing(`vv')
    if _rc exit 459
}
capture assert ci_criterion_code==cond(estimator_twostep==1,2,1) if contract_ok
if _rc exit 459
capture assert gridboot_min_draws==B & gammahat_in_set==1 if ci_delivered
if _rc exit 459
capture assert lin_valid_draws==lin_requested_B if lin_delivered
if _rc exit 459
capture assert cont_valid_draws==cont_requested_B if cont_delivered
if _rc exit 459
save `OUTALL', replace
export delimited using "`OUTALLCSV'", replace

* One bootstrap stream per substantive pair/replication, shared only by its
* FD/FOD methods.  This also catches accidental reuse between B6 and a main
* cell with coincident DGP parameters.
preserve
sort pair_id rep method
capture by pair_id rep: assert boot_seed == boot_seed[1]
if _rc {
    di as err "merge_study_b: methods within a pair do not share a bootstrap seed"
    exit 459
}
by pair_id rep: keep if _n==1
capture isid boot_seed
if _rc {
    di as err "merge_study_b: accidental bootstrap-seed reuse across pair/replication keys"
    exit 459
}
restore

* Replication-level paired FD/FOD differences.  The registered B6 calibration
* cells are FD-only and are intentionally excluded.  Equality of outer-sample
* and synchronized bootstrap seeds is mandatory; realized inner mappings may
* still differ when FD and FOD retain different effective unit sets.
preserve
keep if pairmode == "paired"
keep run_id harness_version schema_version code_version_expected master ///
    pair_id pairmode study design method N T miss missp qstatus boottype coefboot infmode ///
    kappa R B grid gridci refine dgp_spec rep cell_id contract_ok success ///
    estimator_twostep rc gamma_hat ///
    g_bias g_sqerr ci_delivered covered set_length lin_reject5_cond ///
    cont_reject5_cond citest_evaluable citest_accept ///
    N_used N_trans N_iv elapsed_s dgp_seed missing_seed boot_seed
gen byte present = 1
reshape wide cell_id present contract_ok success estimator_twostep rc ///
    gamma_hat g_bias g_sqerr ///
    ci_delivered covered set_length lin_reject5_cond cont_reject5_cond ///
    citest_evaluable citest_accept N_used N_trans N_iv elapsed_s ///
    dgp_seed missing_seed boot_seed, ///
    i(run_id pair_id rep) j(method) string
capture assert presentfd == 1 & presentfod == 1
if _rc | _N != `EXPECTED_PAIRED_ROWS'/2 {
    di as err "merge_study_b: incomplete or extra FD/FOD replication pairs"
    exit 459
}
capture assert dgp_seedfd == dgp_seedfod & ///
    missing_seedfd == missing_seedfod & boot_seedfd == boot_seedfod
if _rc {
    di as err "merge_study_b: paired FD/FOD seeds differ"
    exit 459
}
capture isid boot_seedfd
if _rc {
    di as err "merge_study_b: accidental bootstrap-seed collision across pairs"
    exit 459
}
gen double d_success_fod_fd = successfod-successfd
gen double d_contract_fod_fd = contract_okfod-contract_okfd
gen double d_fallback_fod_fd = ///
    (contract_okfod==1 & estimator_twostepfod==0) - ///
    (contract_okfd==1 & estimator_twostepfd==0)
gen byte target_ci = (study != "LIN")
gen byte target_lin = (infmode == "FULL")
gen byte target_cont = (study != "LIN" & infmode == "FULL")
gen double d_N_trans_fod_fd = N_transfod-N_transfd ///
    if contract_okfd==1 & contract_okfod==1
gen double d_abs_gerr_fod_fd = abs(g_biasfod)-abs(g_biasfd) ///
    if successfd==1 & successfod==1
gen double d_sq_gerr_fod_fd = g_sqerrfod-g_sqerrfd ///
    if successfd==1 & successfod==1
gen double d_ci_delivery_fod_fd = ci_deliveredfod-ci_deliveredfd if target_ci
gen double d_coverage_fod_fd = coveredfod-coveredfd if target_ci
gen double d_citest_eval_fod_fd = citest_evaluablefod-citest_evaluablefd if target_ci
gen double d_citest_accept_fod_fd = citest_acceptfod-citest_acceptfd if target_ci
gen double d_set_length_fod_fd = set_lengthfod-set_lengthfd if target_ci
gen double d_ci_effcov_fod_fd = ///
    (contract_okfod==1 & ci_deliveredfod==1 & coveredfod==1) - ///
    (contract_okfd==1 & ci_deliveredfd==1 & coveredfd==1) if target_ci
gen double d_ci_effcov_2s_fod_fd = ///
    (successfod==1 & ci_deliveredfod==1 & coveredfod==1) - ///
    (successfd==1 & ci_deliveredfd==1 & coveredfd==1) if target_ci
gen double d_lin_reject_fod_fd = ///
    lin_reject5_condfod-lin_reject5_condfd if target_lin
gen double d_lin_effrej_fod_fd = ///
    (contract_okfod==1 & lin_reject5_condfod==1) - ///
    (contract_okfd==1 & lin_reject5_condfd==1) if target_lin
gen double d_lin_effrej_2s_fod_fd = ///
    (successfod==1 & lin_reject5_condfod==1) - ///
    (successfd==1 & lin_reject5_condfd==1) if target_lin
gen double d_cont_reject_fod_fd = ///
    cont_reject5_condfod-cont_reject5_condfd if target_cont
gen double d_cont_effrej_fod_fd = ///
    (contract_okfod==1 & cont_reject5_condfod==1) - ///
    (contract_okfd==1 & cont_reject5_condfd==1) if target_cont
gen double d_cont_effrej_2s_fod_fd = ///
    (successfod==1 & cont_reject5_condfod==1) - ///
    (successfd==1 & cont_reject5_condfd==1) if target_cont
gen double d_elapsed_fod_fd = elapsed_sfod-elapsed_sfd
sort pair_id rep
compress
save `OUTPAIR', replace
export delimited using "`OUTPAIRCSV'", replace

local DPAIRS d_success_fod_fd d_contract_fod_fd d_fallback_fod_fd ///
    d_N_trans_fod_fd d_abs_gerr_fod_fd d_sq_gerr_fod_fd ///
    d_ci_delivery_fod_fd d_coverage_fod_fd d_set_length_fod_fd ///
    d_citest_eval_fod_fd d_citest_accept_fod_fd ///
    d_ci_effcov_fod_fd d_ci_effcov_2s_fod_fd ///
    d_lin_reject_fod_fd d_lin_effrej_fod_fd ///
    d_lin_effrej_2s_fod_fd d_cont_reject_fod_fd ///
    d_cont_effrej_fod_fd d_cont_effrej_2s_fod_fd d_elapsed_fod_fd
local PCNT
local PMEAN
local PSD
foreach v of local DPAIRS {
    local PCNT "`PCNT' n_`v'=`v'"
    local PMEAN "`PMEAN' mean_`v'=`v'"
    local PSD "`PSD' sd_`v'=`v'"
}
local PBY run_id harness_version schema_version code_version_expected master ///
    pair_id pairmode study design N T miss missp qstatus boottype coefboot kappa R B ///
    grid gridci refine dgp_spec infmode target_ci target_lin target_cont
collapse (count) n_pair=rep `PCNT' (mean) `PMEAN' (sd) `PSD', by(`PBY')
foreach v of local DPAIRS {
    gen double mcse_`v' = sd_`v'/sqrt(n_`v')
}
sort study design N T miss missp kappa pair_id
compress
save `OUTPAIRSUM', replace
export delimited using "`OUTPAIRSUMCSV'", replace
restore

* Explicit estimand map.  The no-threshold B3 block targets only linearity
* size: gamma-CI and continuity outputs generated on that shared FULL command
* path are deliberately excluded from scientific summaries.
gen byte target_ci = (study != "LIN" & ci_requested == 1)
gen byte target_citest = (study != "LIN" & citest_requested == 1)
gen byte target_lin = (lin_requested == 1)
gen byte target_cont = (study != "LIN" & cont_requested == 1)
gen byte ok = (contract_ok == 1)
gen byte ok_twostep = (success == 1)
gen byte one_step_fallback = (contract_ok == 1 & estimator_twostep == 0)
gen byte ci_req_i = target_ci
gen byte ci_sreq_i = (target_ci & ok)
gen byte ci_del_i = (target_ci & ok & ci_delivered == 1)
gen byte ci_sreq_2s_i = (target_ci & ok_twostep)
gen byte ci_del_2s_i = (target_ci & ok_twostep & ci_delivered == 1)
gen byte ci_cov_obs = covered if ci_del_i & !missing(gamma0, covered)
gen byte ci_cov_hit = (ci_cov_obs == 1) if !missing(ci_cov_obs)
gen byte ci_cov_2s_obs = covered if target_ci & ok_twostep & ///
    ci_delivered == 1 & !missing(gamma0, covered)
gen byte ci_cov_2s_hit = (ci_cov_2s_obs == 1) if !missing(ci_cov_2s_obs)
gen byte ci_eff_req = target_ci if !missing(gamma0)
gen byte ci_eff_hit = (ok & ci_delivered == 1 & covered == 1) ///
    if ci_eff_req == 1
gen byte ci_eff_2s_hit = (ok_twostep & ci_delivered == 1 & covered == 1) ///
    if ci_eff_req == 1
gen byte ci_inc_obs = ci_incomplete if ci_sreq_i
gen byte ci_empty_obs = ci_empty if ci_sreq_i & ci_incomplete == 0
gen byte ci_disconnected_obs = (n_seg > 1) if ci_del_i
gen byte ci_boundary_obs = (boundary_warn > 0) if ci_sreq_i & ///
    !missing(boundary_warn)
gen double ci_len_obs = set_length if ci_del_i
gen byte citest_req_i = target_citest
gen byte citest_sreq_i = (target_citest & ok)
gen byte citest_ret_i = (target_citest & ok & citest_returned==1)
gen byte citest_eval_i = (target_citest & ok & citest_evaluable==1)
gen byte citest_eval_2s_i = (target_citest & ok_twostep & citest_evaluable==1)
gen byte citest_acc_obs = citest_accept if citest_eval_i
gen byte citest_acc_hit = (citest_acc_obs==1) if !missing(citest_acc_obs)
gen byte citest_acc_2s_obs = citest_accept if citest_eval_2s_i
gen byte citest_acc_2s_hit = (citest_acc_2s_obs==1) if !missing(citest_acc_2s_obs)
gen byte citest_eff_hit = (ok & citest_evaluable==1 & citest_accept==1) ///
    if target_citest
gen byte citest_eff_2s_hit = ///
    (ok_twostep & citest_evaluable==1 & citest_accept==1) if target_citest
gen byte citest_unresolved_obs = inrange(citest_status,3,6) ///
    if target_citest & ok & citest_returned==1
foreach z in lin cont {
    gen byte `z'_req_i = target_`z'
    gen byte `z'_sreq_i = (target_`z' & ok)
    gen byte `z'_del_i = (target_`z' & ok & `z'_delivered == 1)
    gen byte `z'_sreq_2s_i = (target_`z' & ok_twostep)
    gen byte `z'_del_2s_i = (target_`z' & ok_twostep & `z'_delivered == 1)
    gen byte `z'_rej_obs = `z'_reject5_cond if `z'_del_i
    gen byte `z'_rej_2s_obs = `z'_reject5_cond if target_`z' & ///
        ok_twostep & `z'_delivered == 1
    gen byte `z'_eff_hit = (ok & `z'_delivered == 1 & ///
        `z'_reject5_cond == 1) if `z'_req_i
    gen byte `z'_eff_2s_hit = (ok_twostep & `z'_delivered == 1 & ///
        `z'_reject5_cond == 1) if `z'_req_i
}
gen byte hansen_rej = (hansen_p < .05) if ok_twostep & !missing(hansen_p)
gen byte ar1_rej = (ar1_p < .05) if ok_twostep & !missing(ar1_p)
gen byte ar2_rej = (ar2_p < .05) if ok_twostep & !missing(ar2_p)
gen byte cb_req_i = (cb_requested == 1)
gen byte cb_sreq_i = (cb_requested == 1 & ok)
gen byte cb_del_i = (cb_requested == 1 & ok & cb_delivered == 1)
foreach z in lin cont {
    gen double `z'_draw_req = `z'_requested_B if target_`z'
    gen double `z'_draw_valid = `z'_valid_draws if target_`z'
    gen byte `z'_draw_account_req = (target_`z' & ok)
    gen byte `z'_draw_reported = (target_`z' & ok & ///
        !missing(`z'_valid_draws))
    gen double `z'_draw_req_reported = `z'_requested_B ///
        if `z'_draw_reported
}
gen byte ci_draw_account_req = (target_ci & ok)
gen byte ci_draw_reported = (target_ci & ok & !missing(gridboot_min_draws))

* B3/LIN has one scientific estimand only: the null rejection probability of
* the linearity test.  Its shared FULL execution path still produces nuisance
* point estimates and diagnostics, but those are not turned into bias/RMSE or
* coefficient-coverage summaries that could be mistaken for B3 estimands.
gen double err_gamma = gamma_hat-gamma0 if study != "LIN" & ok_twostep & ///
    !missing(gamma_hat, gamma0)
gen double sq_gamma = err_gamma^2 if !missing(err_gamma)
local names lagyb xb qb consd qd lagyd xd
local ests b_Lagyb b_xb b_qb b_consd b_qd b_Lagyd b_xd
local ses se_Lagyb se_xb se_qb se_consd se_qd se_Lagyd se_xd
local trus Lagyb xb qb consd qd Lagyd xd
local bcvs cbcov_Lagyb cbcov_xb cbcov_qb cbcov_consd cbcov_qd ///
           cbcov_Lagyd cbcov_xd
forvalues j = 1/7 {
    local nm : word `j' of `names'
    local ee : word `j' of `ests'
    local ss : word `j' of `ses'
    local tt : word `j' of `trus'
    local bc : word `j' of `bcvs'
    gen double err_`nm' = `ee'-`tt' if study != "LIN" & ok_twostep & ///
        !missing(`ee', `tt')
    gen double sq_`nm' = err_`nm'^2 if !missing(err_`nm')
    gen byte wcov_`nm' = (abs(err_`nm') <= 1.96*`ss') ///
        if study != "LIN" & ok_twostep & !missing(err_`nm', `ss')
    gen byte bcov_`nm' = `bc' if study != "LIN" & ok_twostep & ///
        cb_delivered == 1 & !missing(`bc')
}
local BY run_id harness_version schema_version code_version_expected master ///
    B grid gridci refine R dgp_spec cell_id pair_id pairmode study design method N T ///
    miss missp qstatus boottype coefboot kappa has_x infmode ///
    maxlag_lo maxlag_hi trim gridtype gridsample

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
            if success & (study!="LIN") & (`gate') & !missing(`ss') & `ss'>=0
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
collapse ///
 (count) `XCNT' n_rep=rep n_gamma=err_gamma n_ci_cov=ci_cov_obs ///
    n_ci_cov_2s=ci_cov_2s_obs ///
    n_citest_cov=citest_acc_obs n_citest_cov_2s=citest_acc_2s_obs ///
    n_ci_inc_eval=ci_inc_obs n_ci_empty_eval=ci_empty_obs ///
    n_ci_disc_eval=ci_disconnected_obs n_ci_bound_eval=ci_boundary_obs ///
    n_hansen=hansen_rej n_ar1=ar1_rej n_ar2=ar2_rej ///
    n_lagyb=err_lagyb n_xb=err_xb n_qb=err_qb n_consd=err_consd ///
    n_qd=err_qd n_lagyd=err_lagyd n_xd=err_xd ///
    n_wcov_lagyb=wcov_lagyb n_wcov_xb=wcov_xb n_wcov_qb=wcov_qb ///
    n_wcov_consd=wcov_consd n_wcov_qd=wcov_qd ///
    n_wcov_lagyd=wcov_lagyd n_wcov_xd=wcov_xd ///
    n_bcov_lagyb=bcov_lagyb n_bcov_xb=bcov_xb n_bcov_qb=bcov_qb ///
    n_bcov_consd=bcov_consd n_bcov_qd=bcov_qd ///
    n_bcov_lagyd=bcov_lagyd n_bcov_xd=bcov_xd ///
 (sum) n_contract=ok n_success=ok_twostep n_fallback=one_step_fallback ///
    n_ci_req=ci_req_i n_ci_sreq=ci_sreq_i ///
    n_ci_sreq_2s=ci_sreq_2s_i n_ci_del_2s=ci_del_2s_i ///
    n_ci_del=ci_del_i n_ci_hit=ci_cov_hit n_ci_eff_req=ci_eff_req ///
    n_ci_hit_2s=ci_cov_2s_hit n_ci_eff_hit=ci_eff_hit ///
    n_ci_eff_2s_hit=ci_eff_2s_hit ///
    n_citest_req=citest_req_i n_citest_sreq=citest_sreq_i ///
    n_citest_returned=citest_ret_i n_citest_eval=citest_eval_i ///
    n_citest_eval_2s=citest_eval_2s_i n_citest_hit=citest_acc_hit ///
    n_citest_hit_2s=citest_acc_2s_hit ///
    n_citest_eff_hit=citest_eff_hit n_citest_eff_2s_hit=citest_eff_2s_hit ///
    n_lin_req=lin_req_i n_lin_sreq=lin_sreq_i ///
    n_lin_sreq_2s=lin_sreq_2s_i n_lin_del_2s=lin_del_2s_i ///
    n_lin_del=lin_del_i n_lin_rej=lin_rej_obs ///
    n_lin_rej_2s=lin_rej_2s_obs n_lin_eff=lin_eff_hit ///
    n_lin_eff_2s=lin_eff_2s_hit ///
    n_cont_req=cont_req_i n_cont_sreq=cont_sreq_i ///
    n_cont_sreq_2s=cont_sreq_2s_i n_cont_del=cont_del_i ///
    n_cont_del_2s=cont_del_2s_i ///
    n_cont_rej=cont_rej_obs n_cont_rej_2s=cont_rej_2s_obs ///
    n_cont_eff=cont_eff_hit n_cont_eff_2s=cont_eff_2s_hit ///
    n_cb_req=cb_req_i ///
    n_cb_sreq=cb_sreq_i n_cb_del=cb_del_i ///
    n_lin_draw_account_req=lin_draw_account_req ///
    n_lin_draw_reported=lin_draw_reported ///
    n_cont_draw_account_req=cont_draw_account_req ///
    n_cont_draw_reported=cont_draw_reported ///
    n_ci_draw_account_req=ci_draw_account_req ///
    n_ci_draw_reported=ci_draw_reported ///
 (sum) lin_draws_requested=lin_draw_req ///
    lin_draws_valid=lin_draw_valid ///
    lin_draws_requested_reported=lin_draw_req_reported ///
    cont_draws_requested=cont_draw_req ///
    cont_draws_valid=cont_draw_valid ///
    cont_draws_requested_reported=cont_draw_req_reported ///
 (max) target_ci target_citest target_lin target_cont ///
 (mean) `XMEAN' ci_incomplete_rate=ci_inc_obs ci_empty_rate=ci_empty_obs ///
    ci_disconnected_rate=ci_disconnected_obs ///
    ci_boundary_rate=ci_boundary_obs mean_set_length=ci_len_obs ///
    mean_lin_valid=lin_valid_draws ///
    mean_cont_valid=cont_valid_draws mean_hansen_p=hansen_p ///
    citest_unresolved_rate=citest_unresolved_obs mean_citest_p=citest_p ///
    mean_ar1_p=ar1_p mean_ar2_p=ar2_p mean_N_used=N_used ///
    mean_N_trans=N_trans mean_N_iv=N_iv mean_N_units=N_units ///
    mean_refine_complete=refine_complete ///
    mean_refine_iterations=refine_iterations mean_refine_added=refine_added ///
    mean_refine_remaining=refine_remaining ///
    mean_refine_neigh_unevaluated=refine_neigh_unevaluated ///
    mean_grid_effective=grid_effective mean_grid_admitted=grid_admitted ///
    mean_grid_structural=grid_structural ///
    mean_grid_twostep_admitted=grid_twostep_admitted ///
    mean_gridci_effective=gridci_effective ///
    mean_gridci_admitted=gridci_admitted ///
    mean_gridci_evaluated=gridci_evaluated ///
    mean_ci_unresolved=ci_unresolved ///
    mean_gridboot_min_draws=gridboot_min_draws ///
    mean_continuity_common_grid=continuity_common_grid ///
    mean_units_realized=units_realized ///
    mean_units_dropped=units_dropped mean_missing_rate=missing_rate ///
    mean_gap_events=gap_events mean_gap_periods=gap_periods ///
    mean_elapsed_s=elapsed_s bias_gamma=err_gamma mse_gamma=sq_gamma ///
    bias_lagyb=err_lagyb mse_lagyb=sq_lagyb bias_xb=err_xb mse_xb=sq_xb ///
    bias_qb=err_qb mse_qb=sq_qb bias_consd=err_consd mse_consd=sq_consd ///
    bias_qd=err_qd mse_qd=sq_qd bias_lagyd=err_lagyd mse_lagyd=sq_lagyd ///
    bias_xd=err_xd mse_xd=sq_xd wald_cov_lagyb=wcov_lagyb ///
    wald_cov_xb=wcov_xb wald_cov_qb=wcov_qb wald_cov_consd=wcov_consd ///
    wald_cov_qd=wcov_qd wald_cov_lagyd=wcov_lagyd wald_cov_xd=wcov_xd ///
    boot_cov_lagyb=bcov_lagyb boot_cov_xb=bcov_xb boot_cov_qb=bcov_qb ///
    boot_cov_consd=bcov_consd boot_cov_qd=bcov_qd ///
    boot_cov_lagyd=bcov_lagyd boot_cov_xd=bcov_xd ///
    hansen_reject5=hansen_rej ar1_reject5=ar1_rej ar2_reject5=ar2_rej ///
 (sd) `XSD' sd_set_length=ci_len_obs ///
 (p50) median_set_length=ci_len_obs ///
 (p95) p95_set_length=ci_len_obs, ///
    by(`BY')

capture assert n_ci_cov == n_ci_del & n_ci_cov_2s == n_ci_del_2s & ///
    n_ci_inc_eval == n_ci_sreq & n_ci_bound_eval == n_ci_sreq & ///
    n_ci_disc_eval == n_ci_del
if _rc {
    di as err "merge_study_b: collapsed CI denominators are inconsistent"
    exit 459
}

gen double contract_rate=n_contract/n_rep
gen double contract_mcse=sqrt(contract_rate*(1-contract_rate)/n_rep)
gen double success_rate=n_success/n_rep
gen double success_mcse=sqrt(success_rate*(1-success_rate)/n_rep)
gen double fallback_rate_uncond=n_fallback/n_rep
gen double fallback_rate_cond=n_fallback/n_contract
gen double lin_draw_valid_rate=lin_draws_valid/lin_draws_requested
gen double cont_draw_valid_rate=cont_draws_valid/cont_draws_requested
gen double lin_draw_accounting_rate=n_lin_draw_reported/n_lin_draw_account_req
gen double cont_draw_accounting_rate=n_cont_draw_reported/n_cont_draw_account_req
gen double ci_draw_accounting_rate=n_ci_draw_reported/n_ci_draw_account_req
gen double lin_draw_valid_rate_reported=lin_draws_valid/lin_draws_requested_reported
gen double cont_draw_valid_rate_reported=cont_draws_valid/cont_draws_requested_reported
gen double test_nominal_finite_B=(ceil(.05*(B+1))-1)/(B+1)
foreach z in ci lin cont cb {
    gen double `z'_delivery_uncond=n_`z'_del/n_`z'_req
    gen double `z'_delivery_uncond_mcse=sqrt(`z'_delivery_uncond*(1-`z'_delivery_uncond)/n_`z'_req)
    gen double `z'_delivery_cond=n_`z'_del/n_`z'_sreq
    gen double `z'_delivery_cond_mcse=sqrt(`z'_delivery_cond*(1-`z'_delivery_cond)/n_`z'_sreq)
}
foreach z in ci lin cont {
    gen double `z'_delivery_twostep=n_`z'_del_2s/n_`z'_sreq_2s
    gen double `z'_delivery_twostep_mcse=sqrt(`z'_delivery_twostep* ///
        (1-`z'_delivery_twostep)/n_`z'_sreq_2s)
}
gen double ci_coverage=n_ci_hit/n_ci_cov
gen double ci_coverage_mcse=sqrt(ci_coverage*(1-ci_coverage)/n_ci_cov)
gen double ci_coverage_twostep=n_ci_hit_2s/n_ci_cov_2s
gen double ci_coverage_twostep_mcse=sqrt(ci_coverage_twostep* ///
    (1-ci_coverage_twostep)/n_ci_cov_2s)
gen double ci_effective_coverage=n_ci_eff_hit/n_ci_eff_req
gen double ci_effective_coverage_mcse=sqrt(ci_effective_coverage*(1-ci_effective_coverage)/n_ci_eff_req)
gen double ci_effective_coverage_twostep=n_ci_eff_2s_hit/n_ci_eff_req
gen double ci_eff_cov_2s_mcse=sqrt( ///
    ci_effective_coverage_twostep*(1-ci_effective_coverage_twostep)/n_ci_eff_req)
gen double citest_delivery=n_citest_eval/n_citest_sreq
gen double citest_delivery_mcse=sqrt(citest_delivery*(1-citest_delivery)/n_citest_sreq)
gen double citest_coverage=n_citest_hit/n_citest_cov
gen double citest_coverage_mcse=sqrt(citest_coverage*(1-citest_coverage)/n_citest_cov)
gen double citest_coverage_twostep=n_citest_hit_2s/n_citest_cov_2s
gen double citest_cov_2s_mcse=sqrt(citest_coverage_twostep* ///
    (1-citest_coverage_twostep)/n_citest_cov_2s)
gen double citest_effective_coverage=n_citest_eff_hit/n_citest_req
gen double citest_eff_coverage_mcse=sqrt(citest_effective_coverage* ///
    (1-citest_effective_coverage)/n_citest_req)
gen double citest_eff_coverage_2s=n_citest_eff_2s_hit/n_citest_req
gen double citest_eff_cov_2s_mcse=sqrt(citest_eff_coverage_2s* ///
    (1-citest_eff_coverage_2s)/n_citest_req)

* Pre-specified Gong-Seo Table 1 Grid-B reference values for the FD-only B6
* calibration bridge.  Their outer M is 2,000.  Effective coverage is the
* primary comparison because it retains every requested outer replication;
* conditional coverage is retained as a delivery diagnostic.  These are
* calibration differences, not an exact-replication acceptance test.
gen double gscal_benchmark_coverage = .
replace gscal_benchmark_coverage = .992 if study=="GSCAL" & N==400 & kappa==0
replace gscal_benchmark_coverage = .966 if study=="GSCAL" & N==400 & kappa==1
replace gscal_benchmark_coverage = .986 if study=="GSCAL" & N==800 & kappa==0
replace gscal_benchmark_coverage = .955 if study=="GSCAL" & N==800 & kappa==1
quietly count if study=="GSCAL" & missing(gscal_benchmark_coverage)
if r(N) {
    di as err "merge_study_b: a B6 summary row has no pre-specified benchmark"
    exit 459
}
gen double gscal_benchmark_mcse = sqrt(gscal_benchmark_coverage* ///
    (1-gscal_benchmark_coverage)/2000)
gen double gscal_diff_cond = ci_coverage-gscal_benchmark_coverage ///
    if study=="GSCAL"
gen double gscal_diff_effective = ///
    ci_effective_coverage-gscal_benchmark_coverage if study=="GSCAL"
gen double gscal_diff_citest = ///
    citest_coverage-gscal_benchmark_coverage if study=="GSCAL"
gen double gscal_diff_citest_effective = ///
    citest_effective_coverage-gscal_benchmark_coverage if study=="GSCAL"
gen double mean_set_length_mcse=sd_set_length/sqrt(n_ci_del)
foreach z in lin cont {
    gen double `z'_reject5_cond=n_`z'_rej/n_`z'_del
    gen double `z'_reject5_cond_mcse=sqrt(`z'_reject5_cond*(1-`z'_reject5_cond)/n_`z'_del)
    gen double `z'_reject5_twostep=n_`z'_rej_2s/n_`z'_del_2s
    gen double `z'_reject5_twostep_mcse=sqrt(`z'_reject5_twostep* ///
        (1-`z'_reject5_twostep)/n_`z'_del_2s)
    gen double `z'_reject5_effective=n_`z'_eff/n_`z'_req
    gen double `z'_reject5_effective_mcse=sqrt(`z'_reject5_effective*(1-`z'_reject5_effective)/n_`z'_req)
    gen double `z'_reject5_effective_twostep=n_`z'_eff_2s/n_`z'_req
    gen double `z'_rej5_eff_2s_mcse=sqrt( ///
        `z'_reject5_effective_twostep*(1-`z'_reject5_effective_twostep)/n_`z'_req)
}
foreach z in hansen ar1 ar2 {
    gen double `z'_reject5_mcse=sqrt(`z'_reject5*(1-`z'_reject5)/n_`z')
}
foreach z in gamma lagyb xb qb consd qd lagyd xd {
    gen double rmse_`z'=sqrt(mse_`z')
    drop mse_`z'
}
foreach z in lagyb xb qb consd qd lagyd xd {
    gen double wald_cov_`z'_mcse=sqrt(wald_cov_`z'*(1-wald_cov_`z')/n_wcov_`z')
    gen double boot_cov_`z'_mcse=sqrt(boot_cov_`z'*(1-boot_cov_`z')/n_bcov_`z')
}
order `BY' target_ci target_citest target_lin target_cont n_rep n_contract contract_rate ///
    contract_mcse n_success success_rate success_mcse n_fallback ///
    fallback_rate_uncond fallback_rate_cond n_ci_req n_ci_sreq ///
    n_ci_sreq_2s n_ci_del n_ci_del_2s ci_delivery_uncond ///
    ci_delivery_uncond_mcse ci_delivery_twostep ///
    ci_delivery_twostep_mcse ///
    ci_delivery_cond ci_delivery_cond_mcse n_ci_cov n_ci_hit ci_coverage ///
    ci_coverage_mcse n_ci_eff_req n_ci_eff_hit ci_effective_coverage ///
    ci_effective_coverage_mcse n_citest_req n_citest_sreq ///
    n_citest_returned n_citest_eval citest_delivery citest_delivery_mcse ///
    n_citest_cov n_citest_hit citest_coverage citest_coverage_mcse ///
    citest_effective_coverage citest_eff_coverage_mcse ///
    gscal_benchmark_coverage gscal_benchmark_mcse ///
    gscal_diff_citest gscal_diff_citest_effective ///
    gscal_diff_effective gscal_diff_cond ///
    mean_set_length mean_set_length_mcse ///
    median_set_length p95_set_length n_ci_inc_eval ci_incomplete_rate ///
    n_ci_empty_eval ci_empty_rate n_ci_disc_eval ci_disconnected_rate ///
    n_ci_bound_eval ci_boundary_rate
sort study design N T miss missp kappa method boottype coefboot
compress

foreach nm in lagyb qb consd qd lagyd {
    foreach kind in j f c {
        gen double se_sd_`kind'_`nm' = meanse_`kind'_`nm'/esd_`kind'_`nm' ///
            if study!="LIN" & esd_`kind'_`nm'>0 & !missing(esd_`kind'_`nm')
        gen double cov_mcse_`kind'_`nm' = ///
            sqrt(cov_`kind'_`nm'*(1-cov_`kind'_`nm')/n_`kind'_`nm') ///
            if n_`kind'_`nm'>0
        * Effective joint delivery+coverage, NOT coverage conditional on delivery.
        gen double cov_eff_`kind'_`nm' = ///
            cond(n_`kind'_`nm'>0,cov_`kind'_`nm'*n_`kind'_`nm',0)/n_rep if study!="LIN"
    }
}
save `OUTSUMMARY', replace
export delimited using "`OUTSUMMARYCSV'", replace

* Publish into a nonce-specific generation directory.  The PowerShell verifier
* hashes this complete set and atomically advances a JSON pointer only after all
* eight files exist, so no mixed old/new physical output set is ever certified.
local OUTPUT_GENERATION "_merge_stage_`MERGE_NONCE'"
capture mkdir "`OUTPUT_GENERATION'"
if _rc {
    di as err "merge_study_b: cannot create output generation `OUTPUT_GENERATION'"
    exit 603
}
copy "`OUTALL'" "`OUTPUT_GENERATION'/study_b_all.dta"
copy "`OUTALLCSV'" "`OUTPUT_GENERATION'/study_b_all.csv"
copy "`OUTPAIR'" "`OUTPUT_GENERATION'/study_b_paired_fd_fod.dta"
copy "`OUTPAIRCSV'" "`OUTPUT_GENERATION'/study_b_paired_fd_fod.csv"
copy "`OUTPAIRSUM'" "`OUTPUT_GENERATION'/study_b_paired_summary.dta"
copy "`OUTPAIRSUMCSV'" "`OUTPUT_GENERATION'/study_b_paired_summary.csv"
copy "`OUTSUMMARY'" "`OUTPUT_GENERATION'/study_b_summary.dta"
copy "`OUTSUMMARYCSV'" "`OUTPUT_GENERATION'/study_b_summary.csv"

tempname AH
file open `AH' using "`MERGE_MARKER'", write replace text
file write `AH' "run_id,manifest_sha256,output_generation,rows,cells,nshard,exit_code" _n
file write `AH' "`ACT_RUN',`MANIFEST_SHA',`OUTPUT_GENERATION',`ACT_ROWS',`M_ExpectedCells',`ACT_NSHARD',0" _n
file close `AH'

di as res "MERGE COMPLETE AND VALIDATED: `ACT_RUN'"
di as txt "52 cells; `ACT_NSHARD' shards; no duplicate, missing, extra, or mixed rows"
di as txt "Staged one complete output generation: `OUTPUT_GENERATION'"
exit 0
