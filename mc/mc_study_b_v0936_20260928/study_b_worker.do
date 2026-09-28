*!==========================================================================
*! study_b_worker.do  --  frozen Study-B inference worker for xtdpthresh
*!--------------------------------------------------------------------------
*! Gong-Seo (2026) benchmark DGP: dynamic y, persistent q that is BOTH the
*! threshold variable and a regressor, with the exogeneity of q matching how
*! it is declared in estimation.
*!   QSTATUS = pred  : q weakly exogenous / predetermined (Gong-Seo benchmark);
*!                     corr(e_t, u_{t+1})=RHOEU ; estimate qx(q) predetermined(q).
*!   QSTATUS = endog : q contemporaneously endogenous (stress test);
*!                     corr(e_t, u_t)=RHOEU     ; estimate qx(q) endogenous(q).
*! Study-B additions: inference delivery/valid-draw accounting, complete-set
*! coverage via e(ci_segments), and explicit two-step/fallback diagnostics.
*!
*! Args: STUDY DESIGN METHOD N T MISS MISSP QSTATUS R BB GRID GRIDCI REFINE ///
*!       MASTER REP_START REP_END SHARD OUTSTEM BOOTTYPE COEFBOOT ///
*!       RUN_ID HARNESS_VERSION NSHARD CODE_VERSION_EXPECTED RNG_KIND ///
*!       DGP_SPEC KAPPA CELL_ID PAIR_ID INFMODE ///
*!       MAXLAG_LO MAXLAG_HI TRIM GRIDTYPE GRIDSAMPLE PAIRMODE
*!   DESIGN in: kink jump none weakkink    METHOD in: fd fod
*!   MISS   in: balanced mcar attrition    MISSP = mcar deletion probability
*!       or attrition's target expected calendar-period loss fraction
*! KAPPA is the actual discontinuity at gamma.  For any threshold design,
*! cons_d = KAPPA - q_d*gamma.  Thus KAPPA=0 is continuous.
*!==========================================================================
version 15.0
* `clear all` (NOT bare `clear`). This worker is do'd once per cell from the
* _cell program in study_b_shard.do. Stata protects the running _cell program from
* clear all's `program drop _all`, so _cell survives (empirically: mc_a ran 126
* cells this way). clear all resets matrices/Mata/results each cell, preventing
* the cross-cell resource accumulation that a bare `clear` caused (a hard crash
* around cell 41). The per-cell Mata recompile it triggers is amortized over the
* hundreds/thousands of reps in each real cell.
clear all
set more off
set varabbrev off

args STUDY DESIGN METHOD N T MISS MISSP QSTATUS R BB GRID GRIDCI REFINE ///
     MASTER REP_START REP_END SHARD OUTSTEM BOOTTYPE COEFBOOT ///
     RUN_ID HARNESS_VERSION NSHARD CODE_VERSION_EXPECTED RNG_KIND ///
     DGP_SPEC KAPPA CELL_ID PAIR_ID INFMODE ///
     MAXLAG_LO MAXLAG_HI TRIM GRIDTYPE GRIDSAMPLE PAIRMODE

* The defaults keep old smoke invocations usable.  Formal runs must pass the
* provenance arguments explicitly from study_b_shard.do.
if "`INFMODE'"              == "" local INFMODE FULL
if "`BOOTTYPE'"             == "" local BOOTTYPE wild
if "`QSTATUS'"              == "" local QSTATUS pred
if "`COEFBOOT'"             == "" local COEFBOOT none
if "`RUN_ID'"               == "" local RUN_ID legacy
if "`HARNESS_VERSION'"      == "" local HARNESS_VERSION study_b_type4_v0936
if "`NSHARD'"               == "" local NSHARD 1
if "`CODE_VERSION_EXPECTED'"== "" local CODE_VERSION_EXPECTED 0.9.36
if "`RNG_KIND'"             == "" local RNG_KIND mt64
if "`DGP_SPEC'"             == "" local DGP_SPEC gs26_base_v2
if "`MAXLAG_LO'"            == "" local MAXLAG_LO 1
if "`MAXLAG_HI'"            == "" local MAXLAG_HI 3
if "`TRIM'"                 == "" local TRIM .15
if "`GRIDTYPE'"             == "" local GRIDTYPE uniform
if "`GRIDSAMPLE'"           == "" local GRIDSAMPLE effective
if "`PAIRMODE'"             == "" local PAIRMODE paired

local SCHEMA_VERSION study_b_worker_v0936
local SEED_SCHEME    hash131_mod3_v2

*----- fail-fast argument contract ----------------------------------------
* Enumerated arguments are whitelisted.  Identifier-like provenance fields
* are restricted to a safe CSV/file-key alphabet, so a malformed launcher
* cannot silently shift positional arguments or corrupt the output schema.
local STUDY_OK = inlist("`STUDY'", "GSBAL", "GSGAP", "LIN", "BOOT", ///
        "ATTR", "TROB", "ENDOG", "GRID") | ///
        inlist("`STUDY'", "COEF", "GSCAL", "C1", "C2", "C3", "C4", ///
        "C5", "C6") | inlist("`STUDY'", "C7", "C8", "C9")
if !`STUDY_OK' {
    di as err "study_b_worker: STUDY is not in the registered design whitelist"
    exit 198
}
if !inlist("`DESIGN'", "kink", "jump", "weakkink", "none") {
    di as err "study_b_worker: DESIGN must be kink, jump, weakkink, or none"
    exit 198
}
if !inlist("`INFMODE'", "POINT", "COVERAGE", "FULL") {
    di as err "study_b_worker: INFMODE must be POINT, COVERAGE, or FULL"
    exit 198
}
if !inlist("`METHOD'", "fd", "fod") {
    di as err "study_b_worker: METHOD must be fd or fod"
    exit 198
}
if !inlist("`MISS'", "balanced", "mcar", "attrition") {
    di as err "study_b_worker: MISS must be balanced, mcar, or attrition"
    exit 198
}
if !inlist("`QSTATUS'", "pred", "endog") {
    di as err "study_b_worker: QSTATUS must be pred or endog"
    exit 198
}
if !inlist("`BOOTTYPE'", "wild", "unit") {
    di as err "study_b_worker: BOOTTYPE must be wild or unit"
    exit 198
}
if !inlist("`COEFBOOT'", "none", "onestep", "twostep") {
    di as err "study_b_worker: COEFBOOT must be none, onestep, or twostep"
    exit 198
}
if !inlist("`GRIDTYPE'", "uniform", "quantile") {
    di as err "study_b_worker: GRIDTYPE must be uniform or quantile"
    exit 198
}
if !inlist("`GRIDSAMPLE'", "effective", "observed") {
    di as err "study_b_worker: GRIDSAMPLE must be effective or observed"
    exit 198
}
if !inlist("`PAIRMODE'", "paired", "fdonly") | ///
        ("`PAIRMODE'" == "fdonly" & "`METHOD'" != "fd") {
    di as err "study_b_worker: PAIRMODE must be paired, or fdonly with METHOD=fd"
    exit 198
}
if !inlist("`RNG_KIND'", "mt64", "kiss32") {
    di as err "study_b_worker: RNG_KIND must be mt64 or kiss32"
    exit 198
}
if "`HARNESS_VERSION'" != "study_b_type4_v0936" | ///
   "`CODE_VERSION_EXPECTED'" != "0.9.36" | "`RNG_KIND'" != "mt64" {
    di as err "study_b_worker: unsupported frozen Study-B contract"
    exit 198
}
foreach tok in STUDY RUN_ID HARNESS_VERSION CODE_VERSION_EXPECTED DGP_SPEC {
    if !regexm("``tok''", "^[A-Za-z0-9][A-Za-z0-9_.-]*$") {
        di as err "study_b_worker: `tok' is empty or contains an unsafe character"
        exit 198
    }
}
if !inlist("`DGP_SPEC'", "gs26_base_v2", "gs26_official_v1", ///
        "gs26_endog_v2", "gs26_augmented_v2") {
    di as err "study_b_worker: unsupported DGP_SPEC"
    exit 198
}
if "`OUTSTEM'" == "" | strpos("`OUTSTEM'", ",") | ///
        strpos("`OUTSTEM'", char(34)) | strpos("`OUTSTEM'", char(10)) | ///
        strpos("`OUTSTEM'", char(13)) {
    di as err "study_b_worker: OUTSTEM is empty or unsafe for CSV output"
    exit 198
}
foreach z in N T R BB GRID GRIDCI REFINE MASTER REP_START REP_END SHARD NSHARD ///
             MAXLAG_LO MAXLAG_HI {
    capture confirm integer number ``z''
    if _rc {
        di as err "study_b_worker: `z' must be an integer; received [``z'']"
        exit 198
    }
}
foreach z in MISSP KAPPA TRIM {
    if "`z'" == "KAPPA" & "`KAPPA'" == "" continue
    capture confirm number ``z''
    if _rc {
        di as err "study_b_worker: `z' must be numeric; received [``z'']"
        exit 198
    }
}
if `N' < 10 | `N' > 100000 | `T' < 5 | `T' > 200 | ///
        `R' < 1 | `R' > 1000000 | `BB' < 10 | `BB' > 100000 | ///
        `GRID' < 10 | `GRID' > 10000 | `GRIDCI' < 10 | `GRIDCI' > 10000 | ///
        `REFINE' < 0 | `REFINE' > 20 | `MAXLAG_LO' < 1 | ///
        `MAXLAG_HI' < `MAXLAG_LO' | `MAXLAG_HI' > 200 | ///
        `TRIM' < .01 | `TRIM' > .45 {
    di as err "study_b_worker: N/T/R/B/grid/gridci/refine outside the supported range"
    exit 198
}
if `MASTER' < 0 | `MASTER' > 2147483000 | ///
        `REP_START' < 1 | `REP_END' < `REP_START' | `REP_END' > `R' | ///
        `NSHARD' < 1 | `SHARD' < 1 | `SHARD' > `NSHARD' {
    di as err "study_b_worker: invalid seed, replication interval, or shard contract"
    exit 198
}
if `MISSP' < 0 | `MISSP' >= 1 | ///
        ("`MISS'" == "balanced" & `MISSP' != 0) | ///
        ("`MISS'" == "attrition" & ///
         (2*`MISSP'*`T' < 1 | 2*`MISSP'*`T' > `T'-4)) {
    di as err "study_b_worker: invalid MISSP for balanced/mcar/attrition design"
    exit 198
}

* Legacy KAPPA defaults reproduce the old DGP exactly.  In particular, the
* old jump had cons_d=1 and q_d=2 at gamma=.25, hence actual jump 1.5.
local KAPPA_LEGACY = ("`KAPPA'" == "")
if `KAPPA_LEGACY' {
    if "`DESIGN'" == "jump" local KAPPA 1.5
    else                      local KAPPA 0
}
if `KAPPA' < -5 | `KAPPA' > 5 {
    di as err "study_b_worker: KAPPA must lie in [-5,5]"
    exit 198
}
if inlist("`DESIGN'", "kink", "weakkink", "none") & `KAPPA' != 0 {
    di as err "study_b_worker: kink, weakkink, and none require KAPPA=0"
    exit 198
}
if "`DESIGN'" == "jump" & !`KAPPA_LEGACY' & ///
        !inlist(`KAPPA', 0.1, 0.2, 0.5, 1) {
    di as err "study_b_worker: explicit jump KAPPA must be 0.1, 0.2, 0.5, or 1"
    exit 198
}
* gs26_official_v1 is retained as a manual legacy alias.  Formal Study B uses
* gs26_base_v2, the same label and exact benchmark DGP used by Study A.
local OFFICIAL_WE = inlist("`DGP_SPEC'", "gs26_base_v2", ///
    "gs26_official_v1")
local OFFICIAL_EN = ("`DGP_SPEC'" == "gs26_endog_v2")
local OFFICIAL = (`OFFICIAL_WE' | `OFFICIAL_EN')
if `OFFICIAL_WE' & "`QSTATUS'" != "pred" {
    di as err "study_b_worker: the standard Gong-Seo DGP requires QSTATUS=pred"
    exit 198
}
if `OFFICIAL_EN' & "`QSTATUS'" != "endog" {
    di as err "study_b_worker: gs26_endog_v2 requires QSTATUS=endog"
    exit 198
}
if inlist("`STUDY'", "GSBAL", "GSGAP", "LIN", "BOOT", "GSCAL", "ATTR") & ///
        !`OFFICIAL' {
    di as err "study_b_worker: registered benchmark studies require a Gong-Seo DGP"
    exit 198
}
if "`STUDY'" == "GSCAL" & ("`PAIRMODE'" != "fdonly" | ///
        "`METHOD'" != "fd" | "`MISS'" != "balanced") {
    di as err "study_b_worker: GSCAL must be balanced, FD-only calibration"
    exit 198
}
if ("`STUDY'" == "GSCAL" & `REFINE' != 0) | ///
   ("`STUDY'" != "GSCAL" & `REFINE' != 4) {
    di as err "study_b_worker: B1--B5 require refine(4); B6 calibration requires refine(0)"
    exit 198
}
local SEARCHMODE fixed

* Registry identifiers are appended args so the shard registry, worker output,
* merge audit, and paired FD/FOD comparisons share the same primary keys.
* Derived defaults are only for legacy/manual smoke invocations.
if "`CELL_ID'" == "" local CELL_ID "`STUDY'_`DESIGN'_`METHOD'_N`N'_T`T'_`MISS'_`MISSP'_`QSTATUS'_K`KAPPA'_B`BB'_G`GRID'_GC`GRIDCI'_R`REFINE'_`BOOTTYPE'_`COEFBOOT'"
if "`PAIR_ID'" == "" local PAIR_ID "`STUDY'_`DESIGN'_N`N'_T`T'_`MISS'_`MISSP'_`QSTATUS'_K`KAPPA'_B`BB'_G`GRID'_GC`GRIDCI'_R`REFINE'_`BOOTTYPE'_`COEFBOOT'"
foreach tok in CELL_ID PAIR_ID {
    if !regexm("``tok''", "^[A-Za-z0-9][A-Za-z0-9_.-]*$") {
        di as err "study_b_worker: `tok' is empty or contains an unsafe character"
        exit 198
    }
}

capture set rng `RNG_KIND'
if _rc {
    di as err "study_b_worker: RNG_KIND=`RNG_KIND' is unavailable in this Stata"
    exit 198
}
local RNG_ACTUAL "`c(rng_current)'"
local STATA_VERSION "`c(stata_version)'"
capture set processors 1
adopath ++ "."

* --- DGP constants (Gong-Seo 2026 published MC) ---
local RHOY   = 0.6
local BX     = 0.5      // slope on exogenous x (x_b)
local BQ     = 1.0      // slope on q, lower regime (q_b)
local D2     = 2.0      // slope-on-q shift at gamma (q_d), kink/jump
local D2W    = 0.5      // weak-kink slope shift
local GAMMA  = 0.25     // true threshold
local RHOQ   = 0.7
* Gong-Seo (2026) eq.(12): (e_it, u_i,t+1)' ~ N(0, [[1,rho_eu],[rho_eu,1]]) with
* rho_eu = 0.5.  The published correlation is between e_t and the NEXT period's
* q-innovation, which is the weakly-exogenous lag structure used below.  This
* value is shared with study_a_worker.do so the two studies describe one DGP.
local RHOEU  = 0.5
local SIGE   = 0.5      // sd of the idiosyncratic error
local SIGETA = cond(`OFFICIAL', 0, 1)
local TBURN  = cond(`OFFICIAL', 20, 30)
local HAS_X  = !`OFFICIAL'

* True coefficient values in the command's [cons_d + q_d*q]I(q>gamma)
* parameterisation.  KAPPA, not cons_d, is the economically meaningful jump:
*       KAPPA = cons_d + q_d*gamma.
if "`DESIGN'" == "none" {
    local T_consd = 0
    local T_qd    = 0
    local GAMMA0  = .
}
else {
    local T_qd = cond("`DESIGN'" == "weakkink", `D2W', `D2')
    local T_consd = `KAPPA' - `T_qd'*`GAMMA'
    local GAMMA0  = `GAMMA'
}
local T_kappa = `KAPPA'
local T_Lagyb = `RHOY'
local T_xb    = cond(`HAS_X', `BX', .)
local T_qb    = `BQ'
local T_Lagyd = 0            // persistence does NOT shift at gamma (spurious-regime check)
local T_xd    = cond(`HAS_X', 0, .)
local TMAX    = `T' + `TBURN'

* Numeric namespaces are needed both when generating a fresh replication and
* when recomputing the frozen seed contract for an existing resume row.
local DCODE = cond("`DESIGN'"=="kink",1,cond("`DESIGN'"=="jump",2,cond("`DESIGN'"=="weakkink",3,4)))
local MCODE = cond("`MISS'"=="balanced",1,cond("`MISS'"=="mcar",2,3))
local QCODE = cond("`QSTATUS'"=="pred",1,2)
local SPCODE = cond(inlist("`DGP_SPEC'","gs26_base_v2","gs26_official_v1"),1, ///
    cond("`DGP_SPEC'"=="gs26_endog_v2",2,3))
local BTCODE = cond("`BOOTTYPE'"=="wild",1,2)
local CBCODE = cond("`COEFBOOT'"=="none",1,cond("`COEFBOOT'"=="onestep",2,3))
local KCODE = round((`KAPPA' + 5)*1000)
local PCODE = round(`MISSP'*1000000)

* --- output/resume contract (one file per run/shard) ----------------------
* CELL_ID is the registry/resume key; PAIR_ID deliberately matches paired
* FD/FOD cells.  CONFIG_KEY independently records every cell-level choice and
* detects accidental reuse of a registry ID with different arguments.
local CONFIG_KEY "`STUDY'|`DESIGN'|`METHOD'|N=`N'|T=`T'|`MISS'=`MISSP'|`QSTATUS'|R=`R'|B=`BB'|G=`GRID'|GC=`GRIDCI'|RF=`REFINE'|`BOOTTYPE'|`COEFBOOT'|K=`KAPPA'|`DGP_SPEC'|INF=`INFMODE'|ML=`MAXLAG_LO':`MAXLAG_HI'|TR=`TRIM'|GT=`GRIDTYPE'|GS=`GRIDSAMPLE'|PM=`PAIRMODE'"
local OUT "`OUTSTEM'_SH`SHARD'.csv"
capture confirm file "`OUT'"
local newfile = _rc

local EXPECT_HEADER "run_id,harness_version,schema_version,code_version_expected,code_version_actual,stata_version,rng_kind,rng_actual,seed_scheme,dgp_spec,has_x,infmode,cell_id,pair_id,pairmode,config_key,full_key,study,design,method,N,T,miss,missp,qstatus,boottype,coefboot,kappa,R,B,grid,gridci,refine,master,shard,nshard,rep,dgp_seed,missing_seed,boot_seed,seed_threshold,seed_linearity,seed_continuity,seed_coefficient,rc,fit_rc0,cmd_ok,version_ok,contract_ok,success,estimator_twostep,gamma0,gamma_hat,g_bias,g_sqerr,obj,Lagyb,xb,qb,consd,qd,Lagyd,xd,b_Lagyb,b_xb,b_qb,b_consd,b_qd,b_Lagyd,b_xd,se_Lagyb,se_xb,se_qb,se_consd,se_qd,se_Lagyd,se_xd,ci_requested,ci_delivered,ci_incomplete,ci_empty,n_seg,boundary_warn,covered,set_length,lin_requested,lin_requested_B,lin_valid_draws,lin_delivered,p_lin,lin_reject5_cond,lin_reject5_uncond,cont_requested,cont_requested_B,cont_valid_draws,cont_delivered,p_cont,cont_reject5_cond,cont_reject5_uncond,hansen_p,ar1_p,ar2_p,N_used,N_trans,N_iv,N_units,refine_complete,refine_iterations,refine_added,refine_remaining,refine_neigh_unevaluated,searchmode,searchtol,search_converged,search_incomplete,search_hit_max,search_W2_builds,search_stage1_level,search_stage2_level,search_stage1_points,search_stage2_points,search_stage1_same_split,search_stage2_same_split,search_stage1_rel_gain,search_stage2_rel_gain,grid_max_requested,grid_requested,grid_effective,grid_admitted,grid_structural,grid_twostep_admitted,gridci_requested,gridci_effective,gridci_admitted,gridci_evaluated,ci_unresolved,gridboot_min_draws,threshold_requested_B,continuity_common_grid,ci_bootstrap_certified,units_target,units_realized,units_dropped,analysis_potential,analysis_observed,missing_n,missing_rate,gap_events,gap_periods,fd_pair_rows_potential,fod_rows_potential,cb_requested,cb_delivered,cbcov_Lagyb,cbcov_xb,cbcov_qb,cbcov_consd,cbcov_qd,cbcov_Lagyd,cbcov_xd,rhoy_true,bx_true,bq_true,rhoq,rhoeu_effective,sige,sigeta,tburn,maxlag_lo,maxlag_hi,trim,history,gridtype,gridsample,level,vce,elapsed_s,joint_vce,ar_joint,vce_applied,bwscale,gamma_bw,q_nvals_bw,N_iv_dep,N_iv_dep_near,iv_dep_res,ar1_cond,ar2_cond,ar1_p_cond,ar2_p_cond,se_cond_Lagyb,se_cond_qb,se_cond_consd,se_cond_qd,se_cond_Lagyd,se_delivered,ci_criterion_code,citest_requested,citest_returned,citest_evaluable,citest_gamma,citest_accept,citest_p,citest_D,citest_crit,citest_status,citest_draws,seed_citest,gammahat_in_set"
local EXPECT_COMMAS = length("`EXPECT_HEADER'") - ///
    length(subinstr("`EXPECT_HEADER'", ",", "", .))

* Resume only a byte-compatible v4 file.  A crash can interrupt the final
* append, so first repair exactly one demonstrably torn FINAL physical line.
* Keep elapsed_s as the final (non-scientific) field: a cut inside that last
* field can retain all commas, whereas every earlier-field cut loses a comma.
* Any malformed non-final line, duplicate key, or semantic mismatch remains a
* hard failure: completed MC failures are outcomes and must never be retried.
local done_reps
if !`newfile' {
    tempname RH
    file open `RH' using "`OUT'", read text
    file read `RH' firstline
    if `"`firstline'"' != `"`EXPECT_HEADER'"' {
        file close `RH'
        di as err "study_b_worker: existing `OUT' has an incompatible schema"
        di as err "Archive it or use a new RUN_ID/run directory."
        exit 459
    }

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
        local lastfield = substr(`"`thisline'"', ///
            strrpos(`"`thisline'"', ",") + 1, .)
        local final_numeric = !missing(real(strtrim(`"`lastfield'"')))
        if `ncomma' != `EXPECT_COMMAS' | !`final_numeric' {
            if `malformed_line' {
                file close `RH'
                file close `WH'
                capture erase "`REPAIRTMP'"
                di as err "study_b_worker: multiple malformed rows in resume CSV"
                exit 459
            }
            local malformed_line `physical_line'
        }
        else {
            if `malformed_line' {
                file close `RH'
                file close `WH'
                capture erase "`REPAIRTMP'"
                di as err "study_b_worker: malformed non-final row in resume CSV"
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
        di as txt "study_b_worker: removed demonstrably torn final CSV line `malformed_line'; replication will be rerun"
    }
    capture erase "`REPAIRTMP'"

    if `complete_rows' {
        capture quietly import delimited using "`OUT'", clear varnames(1) ///
            case(preserve) stringcols(_all) bindquote(strict)
        if _rc {
            di as err "study_b_worker: cannot parse existing resume file `OUT'"
            exit 459
        }

        * Freeze both field order and field types.  Importing as string first
        * prevents numeric-looking provenance (for example Stata 17.0) from
        * changing type according to the particular values in a shard.
        local STRVARS run_id harness_version schema_version ///
            code_version_expected code_version_actual stata_version rng_kind ///
            rng_actual seed_scheme dgp_spec infmode cell_id pair_id pairmode config_key ///
            full_key study design method miss qstatus boottype coefboot ///
            searchmode history gridtype gridsample vce
        local EXPECT_VARS : subinstr local EXPECT_HEADER "," " ", all
        local EXPECT_VARS : list retokenize EXPECT_VARS
        unab GOT_VARS : _all
        if `"`GOT_VARS'"' != `"`EXPECT_VARS'"' {
            di as err "study_b_worker: resume file has an incompatible or reordered schema"
            exit 459
        }
        local NUMVARS : list EXPECT_VARS - STRVARS
        capture quietly destring `NUMVARS', replace
        if _rc {
            di as err "study_b_worker: nonnumeric token in a numeric resume column"
            exit 459
        }
        foreach vv of local STRVARS {
            capture confirm string variable `vv'
            if _rc {
                di as err "study_b_worker: incompatible string column `vv'"
                exit 459
            }
            quietly count if `vv' == ""
            if r(N) {
                di as err "study_b_worker: blank mandatory string `vv' in resume CSV"
                exit 459
            }
        }
        foreach vv of local NUMVARS {
            capture confirm numeric variable `vv'
            if _rc {
                di as err "study_b_worker: incompatible numeric column `vv'"
                exit 459
            }
        }
        local CORENUM has_x N T missp kappa R B grid gridci refine master ///
            shard nshard rep dgp_seed missing_seed boot_seed rc fit_rc0 ///
            cmd_ok version_ok contract_ok success ci_requested ci_delivered ///
            lin_requested lin_requested_B lin_delivered cont_requested ///
            cont_requested_B cont_delivered cb_requested cb_delivered ///
            citest_requested citest_returned citest_evaluable ///
            units_target units_realized units_dropped analysis_potential ///
            analysis_observed missing_n missing_rate gap_events gap_periods ///
            fd_pair_rows_potential fod_rows_potential rhoy_true bq_true rhoq ///
            rhoeu_effective sige sigeta tburn maxlag_lo maxlag_hi trim level ///
            elapsed_s
        foreach vv of local CORENUM {
            quietly count if missing(`vv')
            if r(N) {
                di as err "study_b_worker: missing mandatory numeric `vv' in resume CSV"
                exit 459
            }
        }

        capture isid full_key
        if _rc {
            di as err "study_b_worker: duplicate full_key detected in existing `OUT'"
            exit 459
        }
        tempvar _expected_key
        quietly gen strL `_expected_key' = run_id + "|" + cell_id + ///
            "|rep=" + trim(string(rep, "%12.0f"))
        quietly count if full_key != `_expected_key'
        if r(N) {
            di as err "study_b_worker: full_key content is inconsistent in existing `OUT'"
            exit 459
        }
        quietly count if run_id != "`RUN_ID'" | ///
            harness_version != "`HARNESS_VERSION'" | ///
            schema_version != "`SCHEMA_VERSION'" | ///
            code_version_expected != "`CODE_VERSION_EXPECTED'" | ///
            stata_version != "`STATA_VERSION'" | rng_kind != "`RNG_KIND'" | ///
            rng_actual != "`RNG_ACTUAL'" | seed_scheme != "`SEED_SCHEME'" | ///
            master != `MASTER' | shard != `SHARD' | nshard != `NSHARD'
        if r(N) {
            di as err "study_b_worker: run-level provenance mismatch in existing `OUT'"
            exit 459
        }

        quietly count if cell_id == "`CELL_ID'" & ///
            (config_key != "`CONFIG_KEY'" | pair_id != "`PAIR_ID'" | ///
             pairmode != "`PAIRMODE'" | ///
             dgp_spec != "`DGP_SPEC'" | infmode != "`INFMODE'" | ///
             study != "`STUDY'" | design != "`DESIGN'" | ///
             method != "`METHOD'" | has_x != `HAS_X' | N != `N' | T != `T' | ///
             miss != "`MISS'" | abs(missp-`MISSP') > 1e-12 | ///
             qstatus != "`QSTATUS'" | boottype != "`BOOTTYPE'" | ///
             coefboot != "`COEFBOOT'" | abs(kappa-`KAPPA') > 1e-12 | ///
             R != `R' | B != `BB' | grid != `GRID' | gridci != `GRIDCI' | ///
             refine != `REFINE' | searchmode != "`SEARCHMODE'" | ///
             history != "panel" | ///
             gridtype != "`GRIDTYPE'" | gridsample != "`GRIDSAMPLE'" | ///
             level != 95 | vce != "robust" | units_target != `N' | ///
             analysis_potential != `N'*`T' | maxlag_lo != `MAXLAG_LO' | ///
             maxlag_hi != `MAXLAG_HI' | abs(trim-`TRIM') > 1e-12 | ///
             gamma0 != `GAMMA0' | abs(Lagyb-`T_Lagyb') > 1e-12 | ///
             xb != `T_xb' | abs(qb-`T_qb') > 1e-12 | ///
             abs(consd-`T_consd') > 1e-12 | abs(qd-`T_qd') > 1e-12 | ///
             Lagyd != `T_Lagyd' | xd != `T_xd' | ///
             abs(rhoy_true-`RHOY') > 1e-12 | bx_true != `T_xb' | ///
             abs(bq_true-`BQ') > 1e-12 | abs(rhoq-`RHOQ') > 1e-12 | ///
             abs(rhoeu_effective-`RHOEU') > 1e-12 | ///
             abs(sige-`SIGE') > 1e-12 | sigeta != `SIGETA' | ///
             tburn != `TBURN')
        if r(N) {
            di as err "study_b_worker: CELL_ID reused or completed row violates its frozen cell contract"
            exit 459
        }

        tempvar _h_expected _dseed_expected _mseed_expected _bseed_expected
        quietly gen double `_h_expected' = mod(`MASTER'*1000003 + ///
            rep*10007 + `N'*101 + `T'*1009 + `SPCODE'*1000033 + ///
            `QCODE'*1000037, 699999999) if cell_id == "`CELL_ID'"
        quietly gen double `_dseed_expected' = 1 + `_h_expected' ///
            if cell_id == "`CELL_ID'"
        quietly gen double `_mseed_expected' = 700000001 + ///
            mod(`_h_expected'*1009 + `MCODE'*10007 + `PCODE'*101 + 7919, ///
                699999999) if cell_id == "`CELL_ID'"
        quietly gen double `_bseed_expected' = 1400000001 + ///
            mod(`_h_expected'*1013 + `DCODE'*100003 + `MCODE'*10007 + ///
                `PCODE'*101 + `QCODE'*1009 + `BTCODE'*17 + `CBCODE'*29 + ///
                `KCODE'*37 + `BB'*131 + `GRID'*137 + `GRIDCI'*139 + ///
                `REFINE'*149 + `MAXLAG_LO'*151 + `MAXLAG_HI'*157 + ///
                104729, 699999999) ///
            if cell_id == "`CELL_ID'"
        quietly count if cell_id == "`CELL_ID'" & ///
            (dgp_seed != `_dseed_expected' | ///
             missing_seed != `_mseed_expected' | ///
             boot_seed != `_bseed_expected' | rep < `REP_START' | ///
             rep > `REP_END' | rep > `R' | !inlist(fit_rc0,0,1) | ///
             !inlist(cmd_ok,0,1) | !inlist(version_ok,0,1) | ///
             !inlist(contract_ok,0,1) | !inlist(success,0,1) | ///
             !inlist(ci_requested,0,1) | !inlist(ci_delivered,0,1) | ///
             !inlist(lin_requested,0,1) | !inlist(lin_delivered,0,1) | ///
             !inlist(cont_requested,0,1) | !inlist(cont_delivered,0,1) | ///
             !inlist(cb_requested,0,1) | !inlist(cb_delivered,0,1) | ///
             !inlist(citest_requested,0,1) | ///
             !inlist(citest_returned,0,1) | !inlist(citest_evaluable,0,1) | ///
             fit_rc0 != (rc == 0) | cmd_ok > fit_rc0 | ///
             version_ok > fit_rc0 | contract_ok > cmd_ok | ///
             contract_ok > version_ok | ///
             success != (contract_ok == 1 & estimator_twostep == 1) | ///
             ci_requested != (infmode != "POINT") | ///
             lin_requested != (infmode == "FULL") | ///
             cont_requested != (infmode == "FULL") | ///
             cb_requested != (coefboot != "none") | ///
             lin_requested_B != cond(lin_requested,B,0) | ///
             cont_requested_B != cond(cont_requested,B,0) | ///
             ci_delivered > ci_requested | ci_delivered > contract_ok | ///
             lin_delivered > lin_requested | lin_delivered > contract_ok | ///
             cont_delivered > cont_requested | cont_delivered > contract_ok | ///
             cb_delivered > cb_requested | cb_delivered > contract_ok | ///
             units_realized + units_dropped != units_target | ///
             analysis_observed + missing_n != analysis_potential | ///
             abs(missing_rate-missing_n/analysis_potential) > 1e-12 | ///
             units_realized < 0 | units_dropped < 0 | missing_n < 0 | ///
             gap_events < 0 | gap_periods < 0 | ///
             fd_pair_rows_potential < 0 | fod_rows_potential < 0)
        if r(N) {
            di as err "study_b_worker: completed row has invalid seed/outcome/accounting semantics"
            exit 459
        }

        quietly count if cell_id == "`CELL_ID'" & fit_rc0 == 1 & ///
            ((infmode == "POINT" & ///
                (!missing(seed_threshold) | !missing(seed_linearity) | ///
                 !missing(seed_continuity) | !missing(seed_coefficient))) | ///
             (infmode != "POINT" & seed_threshold != boot_seed) | ///
             (infmode == "FULL" & ///
                (seed_linearity != mod(boot_seed+104729,2147483648) | ///
                 seed_continuity != mod(boot_seed+224737,2147483648))) | ///
             (infmode != "FULL" & ///
                (!missing(seed_linearity) | !missing(seed_continuity))) | ///
             (coefboot == "none" & !missing(seed_coefficient)) | ///
             (coefboot != "none" & infmode != "POINT" & ///
                seed_coefficient != mod(boot_seed+350377,2147483648)) | ///
             (citest_requested == 1 & ///
                seed_citest != mod(boot_seed+477377,2147483648)) | ///
             (citest_requested == 0 & !missing(seed_citest)))
        if r(N) {
            di as err "study_b_worker: component bootstrap seeds violate their contract"
            exit 459
        }

        quietly count if cell_id == "`CELL_ID'" & contract_ok == 1 & ///
            (missing(estimator_twostep,N_used,N_trans,N_iv,N_units, ///
                     refine_complete,refine_iterations,refine_added, ///
                     refine_remaining,refine_neigh_unevaluated,obj, ///
                     search_hit_max,search_W2_builds, ///
                     search_stage1_level,search_stage2_level, ///
                     search_stage1_points,search_stage2_points, ///
                     grid_max_requested,grid_requested,grid_effective, ///
                     grid_admitted,grid_structural,gridci_requested, ///
                     gridci_effective,gridci_admitted,gridci_evaluated, ///
                     threshold_requested_B,ci_bootstrap_certified) | ///
             !inlist(estimator_twostep,0,1) | ///
             N_used < 1 | N_trans < 1 | N_iv < 1 | N_units < 1 | ///
             (study == "GSCAL" & N_iv != 28) | ///
             refine != cond(study=="GSCAL",0,4) | ///
             !inrange(refine_iterations,0,refine) | refine_added < 0 | ///
             refine_remaining < 0 | refine_neigh_unevaluated < 0 | ///
             !inlist(refine_complete,0,1) | ///
             (refine_complete == 1 & ///
                 (refine_remaining != 0 | refine_neigh_unevaluated != 0)) | ///
             searchmode != "fixed" | !missing(searchtol) | ///
             !missing(search_converged) | !missing(search_incomplete) | ///
             search_hit_max != 0 | search_W2_builds != 1 | ///
             search_stage1_level != 1 | ///
             !inlist(search_stage2_level,0,1) | ///
             search_stage1_points != grid | ///
             (search_stage2_level==0 & search_stage2_points != 0) | ///
             (search_stage2_level==1 & search_stage2_points != grid) | ///
             !missing(search_stage1_same_split) | ///
             !missing(search_stage2_same_split) | ///
             !missing(search_stage1_rel_gain) | ///
             !missing(search_stage2_rel_gain) | ///
             grid_requested != grid | gridci_requested != gridci | ///
             grid_max_requested != grid | ///
             grid_effective != grid + refine_added | ///
             grid_admitted < 2 | grid_structural < grid_admitted | ///
             grid_structural > grid_effective | ///
             (estimator_twostep == 1 & ///
                 (missing(grid_twostep_admitted) | ///
                  !inrange(grid_twostep_admitted,2,grid_effective))) | ///
             (estimator_twostep == 0 & !missing(grid_twostep_admitted) & ///
                 !inrange(grid_twostep_admitted,0,grid_effective)) | ///
             gridci_effective < 1 | gridci_admitted < 0 | ///
             gridci_evaluated < 0 | gridci_evaluated > gridci_admitted | ///
             threshold_requested_B != B | ci_bootstrap_certified != 0 | ///
             (infmode == "FULL" & ///
                 (missing(continuity_common_grid) | ///
                  !inrange(continuity_common_grid,0,grid_effective))) | ///
             (!missing(gridboot_min_draws) & ///
                 !inrange(gridboot_min_draws,0,B)))
        if r(N) {
            di as err "study_b_worker: contract-valid resume row has invalid estimator diagnostics"
            exit 459
        }

        capture quietly levelsof rep if cell_id == "`CELL_ID'", local(done_reps)
        clear
    }
}

* Open/close for each completed replication.  This makes monitor row counts
* durable, releases the Windows file lock, and confines interruption damage to
* the final physical line that the resume repair above can identify safely.
if `newfile' {
    tempname HH
    file open `HH' using "`OUT'", write replace text
    file write `HH' "`EXPECT_HEADER'" _n
    file close `HH'
}

forvalues rep = `REP_START'/`REP_END' {

    * A recorded failure is still a completed MC outcome.  Retrying only
    * failed rows under the same key conditions the simulation on success and
    * biases failure/delivery rates; use a new RUN_ID for such sensitivity.
    if strpos(" `done_reps' ", " `rep' ") {
        di as txt "  resume: skip existing `CELL_ID' rep `rep'"
        continue
    }

    * Three disjoint numeric namespaces rule out DGP/missing/bootstrap seed
    * collisions.  METHOD is absent by design, preserving paired FD/FOD
    * outer samples and synchronized bootstrap seeds.  Numerical tuning enters
    * only the bootstrap hash so calibration/sensitivity blocks do not reuse an
    * inner stream merely because their DGP parameters coincide.
    * Design/missingness are absent from the DGP hash so variants share base
    * shocks; missingness is paired across methods/designs; bootstrap streams
    * are paired only across methods for the same substantive cell.
    local _h0 = mod(`MASTER'*1000003 + `rep'*10007 + `N'*101 + ///
        `T'*1009 + `SPCODE'*1000033 + `QCODE'*1000037, 699999999)
    local dgp_seed = 1 + `_h0'
    local _hm = mod(`_h0'*1009 + `MCODE'*10007 + `PCODE'*101 + 7919, 699999999)
    local missing_seed = 700000001 + `_hm'
    local _hb = mod(`_h0'*1013 + `DCODE'*100003 + `MCODE'*10007 + ///
        `PCODE'*101 + `QCODE'*1009 + `BTCODE'*17 + `CBCODE'*29 + ///
        `KCODE'*37 + `BB'*131 + `GRID'*137 + `GRIDCI'*139 + ///
        `REFINE'*149 + `MAXLAG_LO'*151 + `MAXLAG_HI'*157 + ///
        104729, 699999999)
    local boot_seed = 1400000001 + `_hb'
    local FULL_KEY "`RUN_ID'|`CELL_ID'|rep=`rep'"

    timer clear 8
    timer on 8

    *----- generate panel -----
    clear
    quietly set seed `dgp_seed'
    quietly set obs `=`N'*`TMAX''
    quietly gen long id = ceil(_n/`TMAX')
    bysort id: gen int t = _n
    quietly xtset id t

    * Exact official mode intentionally does not generate x or eta (and hence
    * does not consume their random draws).  The augmented robustness mode
    * retains the prior worker DGP.
    if !`OFFICIAL' {
        quietly by id: gen double eta = rnormal()*`SIGETA' if _n == 1
        quietly by id: replace eta = eta[1]
        quietly gen double x = rnormal()
    }
    quietly gen double e  = rnormal()*`SIGE'
    quietly gen double es = e/`SIGE'
    quietly gen double w  = rnormal()
    quietly gen double u = .
    if `OFFICIAL_WE' {
        * Published Gong-Seo weak-exogeneity DGP: corr(u_t, e*_{t-1}) = 0.5,
        * var(u) = .25 + .75 = 1.  Identical to study_a_worker.do.
        bysort id (t): replace u = .5*es[_n-1] + sqrt(.75)*w if _n > 1
        bysort id (t): replace u = w if _n == 1
    }
    else if `OFFICIAL_EN' {
        * Gong-Seo Appendix C.3 timing: q innovation loads contemporaneously
        * on e_t; all other official benchmark constants remain unchanged.
        quietly replace u = .5*es + sqrt(.75)*w
    }
    else if "`QSTATUS'" == "endog" {
        quietly replace u = `RHOEU'*es + sqrt(1 - `RHOEU'^2)*w
    }
    else {
        bysort id (t): replace u = `RHOEU'*es[_n-1] + ///
            sqrt(1 - `RHOEU'^2)*w if _n > 1
        bysort id (t): replace u = w if _n == 1
    }
    quietly gen double q = .
    if `OFFICIAL' bysort id (t): replace q = rnormal() if _n == 1
    else          bysort id (t): replace q = u/sqrt(1 - `RHOQ'^2) if _n == 1
    bysort id (t): replace q = `RHOQ'*q[_n-1] + u if _n > 1

    * regime effect g(q): the generated jump at gamma is exactly KAPPA.
    tempvar g
    if "`DESIGN'" == "none" quietly gen double `g' = 0
    else quietly gen double `g' = (`T_consd' + `T_qd'*q)*(q > `GAMMA')

    * y: arbitrary finite start followed by the declared burn-in.
    quietly gen double y = .
    if `OFFICIAL' {
        bysort id (t): replace y = rnormal() if _n == 1
        bysort id (t): replace y = `RHOY'*y[_n-1] + `BQ'*q + `g' + e if _n > 1
    }
    else {
        bysort id (t): replace y = eta + `BX'*x + `BQ'*q + `g' + e if _n == 1
        bysort id (t): replace y = eta + `RHOY'*y[_n-1] + `BX'*x + `BQ'*q + `g' + e if _n > 1
    }

    quietly gen byte analysis = (t > `TBURN')
    * Burn-in did its job during GENERATION; now drop every presample period so
    * history(panel) cannot use the burn-in periods as free instrument history.
    * The official mode exposes exactly T observations per unit with NO extra
    * presample row, matching the published Gong-Seo sample geometry: the
    * observed panel is y_i1..y_iT and z_it=(y_i,t-2,..,y_i1, q_i,t-1,..,q_i1)
    * for t=t_0..T with t_0=3 gives 3+5+7+9=24 lag instruments at T=6.
    * xtdpthresh indexes a potential block from t_min+1 (=2); dynamic-FD
    * complete-case geometry leaves that t=2 block without a transformed row,
    * so its all-zero columns are pruned.  For the four surviving t=3..6
    * blocks, the command also retains one block constant apiece.  Thus B6's
    * maxlag(1 5) preserves q's valid lag-1 moments but reports e(N_iv)=28,
    * not 24.  This command-specific count is not an exact GS moment match.
    * Hence "T=6" counts OBSERVATIONS, not equations: the first differenced
    * equation is at t=3 and each unit contributes 4 equations.
    * study_a_worker.do uses the identical rule, so both studies describe one
    * sample geometry.
    * The augmented (non-official) robustness mode retains one y_i0 row instead.
    if `OFFICIAL' quietly drop if t <= `TBURN'
    else          quietly drop if t <  `TBURN'
    quietly xtset id t

    *----- missingness -----
    if "`MISS'" == "mcar" & `MISSP' > 0 {
        quietly set seed `missing_seed'
        bysort id (t): egen double _tmin = min(cond(analysis, t, .))
        bysort id (t): egen double _tmax = max(cond(analysis, t, .))
        quietly gen double _du = runiform() if analysis
        quietly drop if analysis & _du < `MISSP' & t != _tmin & t != _tmax
        quietly drop _du _tmin _tmax
        bysort id: egen int _na = total(analysis)
        quietly drop if _na < 4
        quietly drop _na
        quietly xtset id t
    }
    else if "`MISS'" == "attrition" {
        * Monotone attrition targeting MISSP expected calendar-period loss
        * while retaining at least four observations per unit.  Conditional on
        * attrition (probability .5), d has mean 2*MISSP*T.  A two-point draw
        * around that target gives the mean exactly without creating an
        * under-length unit.
        quietly set seed `missing_seed'
        quietly gen double _au = runiform() if t == `=`TBURN'+1'
        bysort id (t): egen double _a = max(_au)
        quietly gen double _du = runiform() if t == `=`TBURN'+1'
        bysort id (t): egen double _d0 = max(_du)
        local DTARGET = 2*`MISSP'*`T'
        local DLO = floor(`DTARGET')
        local DHI = ceil(`DTARGET')
        local DMAX = `T' - 4
        if `DLO' < 1 local DLO = 1
        if `DHI' < 1 local DHI = 1
        if `DLO' > `DMAX' local DLO = `DMAX'
        if `DHI' > `DMAX' local DHI = `DMAX'
        local PHI = cond(`DHI'==`DLO', 0, (`DTARGET'-`DLO')/(`DHI'-`DLO'))
        quietly gen int _d = `DLO' + (`DHI'-`DLO')*(_d0 < `PHI') if _a < 0.5
        * drop the last _d analysis periods (t beyond TBURN + T - _d)
        quietly drop if analysis & _a < 0.5 & t > `=`TBURN'+`T'' - _d
        quietly drop _au _a _du _d0 _d
        bysort id: egen int _na = total(analysis)
        quietly drop if _na < 4
        quietly drop _na
        quietly xtset id t
    }

    *----- realised panel and transformation-opportunity audit ------------
    * These quantities are measured before estimation.  The two *_potential
    * metrics are transparent calendar/data opportunities, not substitutes
    * for the command's exact e(N_trans), which is recorded below.
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
    bysort id (t): gen int `_gap' = max(t-t[_n-1]-1,0) ///
        if analysis & analysis[_n-1]
    quietly count if `_gap' > 0 & !missing(`_gap')
    local gap_events = r(N)
    quietly summarize `_gap' if analysis, meanonly
    local gap_periods = cond(r(N), r(sum), 0)
    bysort id (t): gen byte `_fdpot' = analysis & _n > 2 & ///
        t[_n-1] == t-1 & t[_n-2] == t-2
    quietly count if `_fdpot' == 1
    local fd_pair_rows_potential = r(N)
    bysort id (t): gen byte `_complete' = analysis & _n > 1 & t[_n-1] == t-1
    bysort id: egen int `_lastcomplete' = max(cond(`_complete',t,.))
    quietly gen byte `_fodpot' = `_complete' & t < `_lastcomplete'
    quietly count if `_fodpot' == 1
    local fod_rows_potential = r(N)

    *----- estimate: q declared with the correct exogeneity status -----
    local qopt = cond("`QSTATUS'"=="endog", "endogenous(q)", "predetermined(q)")
    local rhs = cond(`HAS_X', "x", "")

    * INFMODE controls the inference machinery on the SAME point-estimation path:
    *   FULL     = threshold CI + linearity + continuity bootstrap tests.
    *   COVERAGE = threshold CI only (-notest- skips the grid-profile tests).
    *   POINT    = -noboot- (no CI, no tests); reserved, not used in Study B.
    if "`INFMODE'" == "POINT" {
        local INFOPT "noboot"
    }
    else if "`INFMODE'" == "COVERAGE" {
        local INFOPT "boot(`BB') boottype(`BOOTTYPE') rseed(`boot_seed') notest"
    }
    else {
        local INFOPT "boot(`BB') boottype(`BOOTTYPE') rseed(`boot_seed')"
    }
    local CITESTOPT ""
    if !missing(`GAMMA0') & "`INFMODE'" != "POINT" ///
        local CITESTOPT "citest(`GAMMA0')"
    capture quietly xtdpthresh y `rhs' if analysis, ///
        qx(q) `qopt' method(`METHOD') maxlag(`MAXLAG_LO' `MAXLAG_HI') ///
        grid(`GRID') gridci(`GRIDCI') gridtype(`GRIDTYPE') ///
        gridsample(`GRIDSAMPLE') refine(`REFINE') trim(`TRIM') ///
        `INFOPT' `CITESTOPT' ///
        coefboot(`COEFBOOT') history(panel) level(95) vce(robust) nowarn
    local rc = _rc
    timer off 8
    quietly timer list 8
    local el = r(t8)

    *----- collect -----
    foreach v in gh gbias gsq obj b_Lagyb b_xb b_qb b_consd b_qd b_Lagyd b_xd ///
                 se_Lagyb se_xb se_qb se_consd se_qd se_Lagyd se_xd ///
                 cbcov_Lagyb cbcov_xb cbcov_qb cbcov_consd cbcov_qd cbcov_Lagyd cbcov_xd ///
                 inc empty nseg bwarn cov slen ///
                 pl lr pc cr hp a1 a2 nu ntr niv nunits est2 ///
                 refcomp refit refadd refrem refnrem ///
                 searchtol searchconv searchinc searchmax searchw2 ///
                 s1lev s2lev s1pts s2pts s1same s2same s1gain s2gain gmax ///
                 greq geff gadm gstruct gadm2 ///
                 gcireq gcieff gciadm gcieval ciunres gbmindraw ///
                 thr_req cont_common ci_cert ///
                 lin_valid cont_valid eseed_thr eseed_lin eseed_cont eseed_coef ///
                 citest_gamma citest_accept citest_p citest_D citest_crit ///
                 citest_status citest_draws seed_citest {
        local `v' = .
    }
    local code_actual NA
    local fit_rc0 = (`rc' == 0)
    local cmdok = 0
    local vok = 0
    local contract = 0
    local succ = 0
    local deliv = 0
    local smode "`SEARCHMODE'"
    * Requested-flag defaults must reflect what THIS INFMODE actually asked the
    * command to do, because on an estimation failure (rc!=0) there is no e()
    * to read them back from.  Defaulting to 1 would put -notest- COVERAGE cells
    * into the linearity/continuity denominators whenever a fit failed, silently
    * deflating the reported test rates.  A 100%-success smoke run never exposes
    * this; only failure-bearing production cells do.
    local ci_requested = ("`INFMODE'" != "POINT")
    local lin_requested = ("`INFMODE'" == "FULL")
    local lin_requested_B = cond(`lin_requested', `BB', 0)
    local lin_delivered = 0
    local lin_uncond = cond(`lin_requested', 0, .)
    local cont_requested = ("`INFMODE'" == "FULL")
    local cont_requested_B = cond(`cont_requested', `BB', 0)
    local cont_delivered = 0
    local cont_uncond = cond(`cont_requested', 0, .)
    local cb_requested = ("`COEFBOOT'" != "none")
    local cb_deliv = 0
    local citest_requested = (!missing(`GAMMA0') & "`INFMODE'" != "POINT")
    local citest_returned = 0
    local citest_evaluable = 0
    if `rc' == 0 {
        local code_actual "`e(cmdversion)'"
        local cmdok = ("`e(cmd)'" == "xtdpthresh")
        local vok = ("`e(cmdversion)'" == "`CODE_VERSION_EXPECTED'")
        if !`cmdok' | !`vok' {
            di as err "study_b_worker: loaded command identity/version mismatch"
            di as err "expected xtdpthresh `CODE_VERSION_EXPECTED'; received `e(cmd)' `e(cmdversion)'"
            exit 459
        }
        local contract = (`cmdok' & `vok')
        capture matrix _b_contract = e(b)
        if _rc local contract = 0
        capture assert !missing(e(N), e(N_trans), e(N_iv), e(N_units), ///
            e(gamma), e(obj), e(refine_requested), e(refine_complete), ///
            e(refine_iterations), e(refine_added), e(refine_remaining), ///
            e(refine_neigh_unevaluated), e(search_hit_max), ///
            e(search_W2_builds), e(search_stage1_level), ///
            e(search_stage2_level), e(search_stage1_points), ///
            e(search_stage2_points), e(grid_max_requested), ///
            e(grid_requested), e(grid_effective), e(grid_admitted), ///
            e(grid_structural), e(gridci_requested), e(gridci_effective), ///
            e(gridci_admitted), e(gridci_evaluated), ///
            e(boot_threshold_requested), e(ci_bootstrap_certified)) & ///
            inlist(e(estimator_twostep), 0, 1)
        if _rc local contract = 0
        foreach lab in Lag_y_b q_b cons_d q_d Lag_y_d {
            capture scalar __contract = _b[`lab']
            if _rc local contract = 0
            else if missing(__contract) local contract = 0
        }
        if `HAS_X' {
            foreach lab in x_b x_d {
                capture scalar __contract = _b[`lab']
                if _rc local contract = 0
                else if missing(__contract) local contract = 0
            }
        }
        capture local est2 = e(estimator_twostep)
        capture local nunits = e(N_units)
        capture local niv = e(N_iv)
        if "`STUDY'" == "GSCAL" & (missing(`niv') | `niv' != 28) {
            di as err "study_b_worker: GSCAL requires exactly 28 realized instruments under maxlag(1 5)"
            exit 459
        }
        capture local refcomp = e(refine_complete)
        capture local refit = e(refine_iterations)
        capture local refadd = e(refine_added)
        capture local refrem = e(refine_remaining)
        capture local refnrem = e(refine_neigh_unevaluated)
        capture local smode "`e(searchmode)'"
        capture local searchtol = e(searchtol)
        capture local searchconv = e(search_converged)
        capture local searchinc = e(search_incomplete)
        capture local searchmax = e(search_hit_max)
        capture local searchw2 = e(search_W2_builds)
        capture local s1lev = e(search_stage1_level)
        capture local s2lev = e(search_stage2_level)
        capture local s1pts = e(search_stage1_points)
        capture local s2pts = e(search_stage2_points)
        capture local s1same = e(search_stage1_same_split)
        capture local s2same = e(search_stage2_same_split)
        capture local s1gain = e(search_stage1_rel_gain)
        capture local s2gain = e(search_stage2_rel_gain)
        capture local gmax = e(grid_max_requested)
        capture local greq = e(grid_requested)
        capture local geff = e(grid_effective)
        capture local gadm = e(grid_admitted)
        capture local gstruct = e(grid_structural)
        capture local gadm2 = e(grid_twostep_admitted)
        capture local gcireq = e(gridci_requested)
        capture local gcieff = e(gridci_effective)
        capture local gciadm = e(gridci_admitted)
        capture local gcieval = e(gridci_evaluated)
        capture local ciunres = e(ci_unresolved)
        capture local gbmindraw = e(gridboot_min_draws)
        capture local thr_req = e(boot_threshold_requested)
        capture local cont_common = e(continuity_common_grid)
        capture local ci_cert = e(ci_bootstrap_certified)
        capture assert "`smode'" == "fixed" & ///
            missing(`searchtol',`searchconv',`searchinc',`s1same', ///
                    `s2same',`s1gain',`s2gain') & ///
            `searchmax' == 0 & `searchw2' == 1 & ///
            `s1lev' == 1 & inlist(`s2lev',0,1) & ///
            `s1pts' == `GRID' & ///
            `s2pts' == cond(`s2lev'==1,`GRID',0) & ///
            `gmax' == `GRID' & `greq' == `GRID' & ///
            `geff' == `GRID' + `refadd' & ///
            `gadm' >= 2 & `gadm' <= `gstruct' & `gstruct' <= `geff' & ///
            e(refine_requested) == `REFINE' & ///
            inrange(`refit',0,`REFINE') & `refadd' >= 0 & ///
            `refrem' >= 0 & `refnrem' >= 0 & ///
            inlist(`refcomp',0,1) & ///
            (`refcomp' == 0 | (`refrem' == 0 & `refnrem' == 0)) & ///
            missing(e(search_level2_points)) & ///
            missing(e(search_level3_points)) & ///
            `gcireq' == `GRIDCI' & `gcieff' >= 1 & ///
            `gciadm' >= 0 & `gcieval' >= 0 & `gcieval' <= `gciadm' & ///
            `thr_req' == `BB' & `ci_cert' == 0 & ///
            (missing(`gbmindraw') | inrange(`gbmindraw',0,`BB'))
        if _rc local contract = 0
        if "`INFMODE'" == "FULL" {
            capture assert !missing(`cont_common') & ///
                inrange(`cont_common',0,`geff')
            if _rc local contract = 0
        }
        if `est2' == 1 {
            capture assert !missing(`gadm2') & inrange(`gadm2',2,`geff')
            if _rc local contract = 0
        }
        else if `est2' == 0 {
            capture assert missing(`gadm2') | inrange(`gadm2',0,`geff')
            if _rc local contract = 0
        }
        * A contract-valid one-step fallback is a delivered estimate, but it is
        * not counted as primary two-step success.  The merger reports both.
        local succ = (`contract' & `est2' == 1)
        local gh  = e(gamma)
        local obj = e(obj)
        if !missing(`GAMMA0') & !missing(`gh') {
            local gbias = `gh' - `GAMMA0'
            local gsq   = (`gh' - `GAMMA0')^2
        }
        foreach nm in Lagyb xb qb consd qd Lagyd xd {
            local lab = cond("`nm'"=="Lagyb","Lag_y_b", ///
                        cond("`nm'"=="xb","x_b", ///
                        cond("`nm'"=="qb","q_b", ///
                        cond("`nm'"=="consd","cons_d", ///
                        cond("`nm'"=="qd","q_d", ///
                        cond("`nm'"=="Lagyd","Lag_y_d","x_d"))))))
            capture local b_`nm' = _b[`lab']
            capture local se_`nm' = _se[`lab']
        }
        * Preserve missing exactly: when inversion is incomplete, v0.9.36
        * deliberately withdraws ci_empty and the formal confidence set.
        local inc   = e(ci_incomplete)
        local empty = e(ci_empty)
        local nseg  = e(ci_nseg)
        local bwarn = e(boundary_warn)
        local pl  = e(pval_lin)
        local pc  = e(pval_cont)
        capture local lin_requested_B = e(boot_linearity_requested)
        capture local lin_valid = e(boot_linearity_valid)
        capture local cont_requested_B = e(boot_continuity_requested)
        capture local cont_valid = e(boot_continuity_valid)
        local lin_requested = (`lin_requested_B' > 0)
        local cont_requested = (`cont_requested_B' > 0)
        if !`lin_requested' local lin_uncond = .
        if !`cont_requested' local cont_uncond = .
        capture local eseed_thr = e(seed_threshold)
        capture local eseed_lin = e(seed_linearity)
        capture local eseed_cont = e(seed_continuity)
        capture local eseed_coef = e(seed_coefficient)
        capture local citest_gamma = e(citest_gamma)
        capture local citest_accept = e(citest_accept)
        capture local citest_p = e(citest_p)
        capture local citest_D = e(citest_D)
        capture local citest_crit = e(citest_crit)
        capture local citest_status = e(citest_status)
        capture local citest_draws = e(citest_draws)
        capture local seed_citest = e(seed_citest)
        local citest_returned = (`citest_requested' & !missing(`citest_status'))
        local citest_evaluable = (`citest_returned' & inlist(`citest_status',1,2) & ///
            !missing(`citest_gamma',`citest_accept',`citest_p',`citest_D'))
        if `citest_requested' {
            capture assert `citest_returned' == 1 & ///
                abs(`citest_gamma'-`GAMMA0') <= 1e-12*max(1,abs(`GAMMA0')) & ///
                inrange(`citest_status',1,6) & ///
                `seed_citest' == mod(`boot_seed'+477377,2147483648)
            if _rc {
                di as err "study_b_worker: citest() result/seed contract failed"
                exit 459
            }
            if `citest_status' == 1 {
                capture assert `citest_evaluable' == 1 & `citest_draws' == `BB' & ///
                    !missing(`citest_crit') & inlist(`citest_accept',0,1) & ///
                    inrange(`citest_p',0,1) & ///
                    `citest_accept' == (`citest_D' <= `citest_crit') & ///
                    `citest_accept' == (`citest_p' > .05)
                if _rc {
                    di as err "study_b_worker: regular citest() result is inconsistent"
                    exit 459
                }
            }
            else if `citest_status' == 2 {
                capture assert `citest_evaluable' == 1 & `citest_D' == 0 & ///
                    `citest_accept' == 1 & `citest_p' == 1 & ///
                    missing(`citest_crit',`citest_draws')
                if _rc {
                    di as err "study_b_worker: mechanical citest() result is inconsistent"
                    exit 459
                }
            }
            else if `citest_evaluable' != 0 {
                di as err "study_b_worker: unresolved citest() was marked evaluable"
                exit 459
            }
        }
        else {
            capture assert `citest_returned' == 0 & `citest_evaluable' == 0 & ///
                missing(`citest_gamma',`citest_accept',`citest_p',`citest_D', ///
                    `citest_crit',`citest_status',`citest_draws',`seed_citest')
            if _rc {
                di as err "study_b_worker: unrequested citest() returned state"
                exit 459
            }
        }
        local hp  = e(hansen_p)
        local a1  = e(ar1_p)
        local a2  = e(ar2_p)
        local nu  = e(N)
        local ntr = e(N_trans)
        local niv = e(N_iv)
        if `inc' != 1 & `empty' != 1 {
            capture matrix _cis = e(ci_segments)
            if _rc == 0 {
                local slen = 0
                local cov  = 0
                local nr = rowsof(_cis)
                local anyseg = 0
                forvalues rr = 1/`nr' {
                    local Lo = _cis[`rr', 1]
                    local Hi = _cis[`rr', 2]
                    if !missing(`Lo') & !missing(`Hi') {
                        local anyseg = 1
                        local slen = `slen' + (`Hi' - `Lo')
                        if !missing(`GAMMA0') & `Lo' <= `GAMMA0' & `GAMMA0' <= `Hi' local cov = 1
                    }
                }
                local deliv = `anyseg'
                if !`deliv' {
                    local cov = .
                    local slen = .
                }
                if missing(`GAMMA0') local cov = .
            }
        }
        if !missing(`pl') {
            local lin_delivered = 1
            local lr = (`pl' < 0.05)
            local lin_uncond = `lr'
        }
        if !missing(`pc') {
            local cont_delivered = 1
            local cr = (`pc' < 0.05)
            local cont_uncond = `cr'
        }
        * Delivered inference must carry an auditable positive valid-draw
        * denominator.  Missing counts remain legitimate only for requested
        * inference that the command could not deliver on this outer sample.
        if `lin_delivered' & (missing(`lin_valid') | `lin_valid' < 1 | ///
                `lin_valid' > `lin_requested_B' | ///
                `lin_valid' != floor(`lin_valid')) {
            di as err "study_b_worker: delivered linearity p-value lacks a valid draw count"
            exit 459
        }
        if `cont_delivered' & (missing(`cont_valid') | `cont_valid' < 1 | ///
                `cont_valid' > `cont_requested_B' | ///
                `cont_valid' != floor(`cont_valid')) {
            di as err "study_b_worker: delivered continuity p-value lacks a valid draw count"
            exit 459
        }
        if `deliv' & (missing(`gbmindraw') | `gbmindraw' < 1 | ///
                `gbmindraw' > `BB' | `gbmindraw' != floor(`gbmindraw')) {
            di as err "study_b_worker: delivered threshold set lacks grid-bootstrap draw accounting"
            exit 459
        }
        * coefficient bootstrap-CI coverage (only when coefboot() posted a matrix)
        if "`COEFBOOT'" != "none" {
            capture matrix _cb = e(b_bootci)
            local cb_deliv = (_rc == 0)      // was a coefficient bootstrap CI delivered?
            if `cb_deliv' {
                foreach nm in Lagyb xb qb consd qd Lagyd xd {
                    local lab = cond("`nm'"=="Lagyb","Lag_y_b", ///
                                cond("`nm'"=="xb","x_b", ///
                                cond("`nm'"=="qb","q_b", ///
                                cond("`nm'"=="consd","cons_d", ///
                                cond("`nm'"=="qd","q_d", ///
                                cond("`nm'"=="Lagyd","Lag_y_d","x_d"))))))
                    local cc = colnumb(_cb, "`lab'")
                    if !missing(`cc') {
                        local cblo = _cb[1, `cc']
                        local cbhi = _cb[2, `cc']
                        if !missing(`cblo') & !missing(`cbhi') ///
                            local cbcov_`nm' = (`cblo' <= `T_`nm'' & `T_`nm'' <= `cbhi')
                    }
                }
            }
        }
    }


    * 0.9.36 telemetry: missing SE is not a point-estimation failure.
    foreach zz in joint_vce ar_joint vce_applied bwscale gamma_bw q_nvals_bw N_iv_dep N_iv_dep_near iv_dep_res ar1_cond ar2_cond ar1_p_cond ar2_p_cond se_cond_Lagyb se_cond_qb se_cond_consd se_cond_qd se_cond_Lagyd se_delivered ci_criterion_code gammahat_in_set {
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

        local ci_criterion_code = cond("`e(ci_criterion)'"=="twostep",2, ///
            cond("`e(ci_criterion)'"=="onestep",1,.))
        if `deliv' {
            local gammahat_in_set = 0
            capture matrix _mc_segments = e(ci_segments)
            if !_rc {
                forvalues rr = 1/`=rowsof(_mc_segments)' {
                    if !missing(_mc_segments[`rr',1],_mc_segments[`rr',2]) & ///
                        _mc_segments[`rr',1] <= e(gamma)+1e-12*max(1,abs(e(gamma))) & ///
                        e(gamma)-1e-12*max(1,abs(e(gamma))) <= _mc_segments[`rr',2] ///
                        local gammahat_in_set = 1
                }
            }
        }
        if `contract' & `ci_criterion_code' != cond(`est2'==1,2,1) {
            di as error "Study B: unexpected confidence-set criterion"
            exit 459
        }
        if `deliv' & (`gbmindraw' != `BB' | `gammahat_in_set' != 1) {
            di as error "Study B: draws=`gbmindraw' / `BB'; gamma-hat membership=`gammahat_in_set'"
            exit 459
        }
        if (`lin_delivered' & `lin_valid' != `lin_requested_B') | ///
           (`cont_delivered' & `cont_valid' != `cont_requested_B') {
            di as error "Study B: a reported test did not use all requested draws"
            exit 459
        }
    }

    tempname FH
    file open `FH' using "`OUT'", write append text
    file write `FH' "`RUN_ID',`HARNESS_VERSION',`SCHEMA_VERSION',`CODE_VERSION_EXPECTED',`code_actual',`STATA_VERSION',`RNG_KIND',`RNG_ACTUAL',`SEED_SCHEME',`DGP_SPEC',`HAS_X',`INFMODE'," ///
        "`CELL_ID',`PAIR_ID',`PAIRMODE',`CONFIG_KEY',`FULL_KEY',`STUDY',`DESIGN',`METHOD',`N',`T',`MISS',`MISSP',`QSTATUS',`BOOTTYPE',`COEFBOOT',`KAPPA'," ///
        "`R',`BB',`GRID',`GRIDCI',`REFINE',`MASTER',`SHARD',`NSHARD',`rep'," ///
        "`dgp_seed',`missing_seed',`boot_seed',`eseed_thr',`eseed_lin',`eseed_cont',`eseed_coef'," ///
        "`rc',`fit_rc0',`cmdok',`vok',`contract',`succ',`est2'," ///
        "`GAMMA0',`gh',`gbias',`gsq',`obj'," ///
        "`T_Lagyb',`T_xb',`T_qb',`T_consd',`T_qd',`T_Lagyd',`T_xd'," ///
        "`b_Lagyb',`b_xb',`b_qb',`b_consd',`b_qd',`b_Lagyd',`b_xd'," ///
        "`se_Lagyb',`se_xb',`se_qb',`se_consd',`se_qd',`se_Lagyd',`se_xd'," ///
        "`ci_requested',`deliv',`inc',`empty',`nseg',`bwarn',`cov',`slen'," ///
        "`lin_requested',`lin_requested_B',`lin_valid',`lin_delivered',`pl',`lr',`lin_uncond'," ///
        "`cont_requested',`cont_requested_B',`cont_valid',`cont_delivered',`pc',`cr',`cont_uncond'," ///
        "`hp',`a1',`a2',`nu',`ntr',`niv',`nunits'," ///
        "`refcomp',`refit',`refadd',`refrem',`refnrem'," ///
        "`smode',`searchtol',`searchconv',`searchinc',`searchmax',`searchw2'," ///
        "`s1lev',`s2lev',`s1pts',`s2pts',`s1same',`s2same',`s1gain',`s2gain',`gmax'," ///
        "`greq',`geff',`gadm',`gstruct',`gadm2'," ///
        "`gcireq',`gcieff',`gciadm',`gcieval',`ciunres',`gbmindraw'," ///
        "`thr_req',`cont_common',`ci_cert'," ///
        "`N',`units_realized',`units_dropped',`analysis_potential',`analysis_observed',`missing_n',`missing_rate'," ///
        "`gap_events',`gap_periods',`fd_pair_rows_potential',`fod_rows_potential'," ///
        "`cb_requested',`cb_deliv',`cbcov_Lagyb',`cbcov_xb',`cbcov_qb',`cbcov_consd',`cbcov_qd',`cbcov_Lagyd',`cbcov_xd'," ///
        "`RHOY',`T_xb',`BQ',`RHOQ',`RHOEU',`SIGE',`SIGETA',`TBURN',`MAXLAG_LO',`MAXLAG_HI',`TRIM',panel,`GRIDTYPE',`GRIDSAMPLE',95,robust,`el',`joint_vce',`ar_joint',`vce_applied',`bwscale',`gamma_bw',`q_nvals_bw',`N_iv_dep',`N_iv_dep_near',`iv_dep_res',`ar1_cond',`ar2_cond',`ar1_p_cond',`ar2_p_cond',`se_cond_Lagyb',`se_cond_qb',`se_cond_consd',`se_cond_qd',`se_cond_Lagyd',`se_delivered',`ci_criterion_code',`citest_requested',`citest_returned',`citest_evaluable',`citest_gamma',`citest_accept',`citest_p',`citest_D',`citest_crit',`citest_status',`citest_draws',`seed_citest',`gammahat_in_set'" _n
    file close `FH'

    if mod(`rep', 25) == 0 di as txt "  [`STUDY'/`DESIGN'/`METHOD'/N`N'/`MISS'`MISSP'/`QSTATUS'] rep `rep' g=`gh' cov=`cov' deliv=`deliv'"
}
di as res "WORKER DONE: `STUDY' `DESIGN' `METHOD' N`N' T`T' `MISS'`MISSP' `QSTATUS' reps `REP_START'-`REP_END' shard `SHARD'"
