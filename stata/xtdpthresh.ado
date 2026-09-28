*! version 0.9.37  28sep2026
*! xtdpthresh -- dynamic panel threshold regression (Seo-Shin 2016; Gong-Seo 2026)
*! Duy Chinh Nguyen (IU VNU-HCM) & Nhat Duy Lai (SGU, corresponding). See -help xtdpthresh-.

program define xtdpthresh, eclass sortpreserve
    version 15.0

    // Standard eclass replay: -xtdpthresh- and -xtdpthresh, level()- after
    // estimation redisplay the active results instead of being reparsed as a
    // new model with a missing varlist/qx().
    if replay() {
        if "`e(cmd)'" != "xtdpthresh" error 301
        // v0.9.29: default to the level of the estimation, not c(level).
        syntax [, Level(string)]
        if `"`level'"' == "" {
            local level = e(level)
            if missing(`level') local level = c(level)
        }
        ereturn display, level(`level')
        di as text "γ̂ = " as res %9.0g e(gamma) _c
        if !missing(e(gamma_lo)) & !missing(e(gamma_hi)) {
            di as text "   threshold set hull (" as res e(level) as text "%) = [" ///
               as res %9.0g e(gamma_lo) as text ", " as res %9.0g e(gamma_hi) as text "]"
        }
        else di ""
        exit
    }

    // Capture the full command line BEFORE any parsing (the iv() sub-parser
    // below reuses local 0) so that e(cmdline) can be stored (v0.7.0).
    local cmdline `"`0'"'

    syntax varlist(min=1 numeric ts) [if] [in] ,   ///
        QX(varname numeric)                         ///
        [                                           ///
        IV(string)                                  ///
        ENDOgenous(varlist numeric ts)              ///
        PREDetermined(varlist numeric ts)           ///
        EXOgenous(varlist numeric ts)               ///
        KINK                                        ///
        STATIC                                      ///
        TD                                          ///
        COLLAPSE                                    ///
        MAXLAG(numlist max=2 min=1 integer >0)      ///
        LEVMAXLAG(string)                           ///
        METHOD(string)                              ///
        GRID(integer 100)                           ///
        SEARCHMode(string)                          ///
        SEARCHMAX(integer -1)                       ///
        SEARCHTol(real 1e-8)                        ///
        GRIDCI(integer 100)                         ///
        CITest(numlist max=1 min=1)                 ///
        GRIDType(string)                            ///
        MINREGime(integer 0)                        ///
        GRIDSample(string)                          ///
        BOOTType(string)                            ///
        HISTory(string)                             ///
        REFine(integer 0)                           ///
        BWscale(real 1.5)                           ///
        TRIM(real 0.10)                             ///
        BOOT(integer 299)                           ///
        RSEED(string)                               ///
        Level(cilevel)                              ///
        NOBOOT                                      ///
        NOWARN                                      ///
        EXPORTGMM                                   ///
        NOTEST                                      ///
        VCE(string)                                 ///
        COEFCItype(string)                          ///
        COEFBoot(string)                            ///
        NOCENTER                                    ///
        VERBOSE                                     ///
        ]
    local flag_verbose = cond("`verbose'" != "", 1, 0)
    local flag_exportgmm = cond("`exportgmm'" != "", 1, 0)
    // v0.7.13 (audit R4, C1) / v0.9.9 R25: vce(robust) is the default two-step
    // cluster-robust sandwich (fixed-weight, no small-sample adjustment);
    // vce(windmeijer) applies the Windmeijer (2005) finite-sample correction
    // to the two-step variance (computed once at the final estimate; no
    // effect on the grid search or the CI/test bootstraps, which run on the
    // fast one-step machinery; the legacy coefficient bootstrap can replay
    // the reported two-step estimator when explicitly requested.
    local vce = lower(trim("`vce'"))
    if "`vce'" == "" local vce "robust"
    // v0.9.9 R25: "uncorrected" read as model-based/nonrobust in the Stata
    // ecosystem, but the default has always been the CLUSTER-ROBUST
    // two-step sandwich -- merely without the Windmeijer small-sample
    // correction. Renamed vce(robust). No alias: vce() only ever existed
    // in unreleased builds (introduced v0.7.13, post the 0.7.12 submission
    // package), so there is no installed base to stay compatible with.
    if !inlist("`vce'", "robust", "windmeijer") {
        di as err "option vce() must be robust or windmeijer"
        exit 198
    }
    local flag_vce_wind = cond("`vce'" == "windmeijer", 1, 0)
    // v0.8.1 (audit R6): the CENTERED clustered moment covariance of
    // Seo-Shin (2016, eq. 11) / xthenreg is now the DEFAULT -- it is the
    // convention of the estimator this command implements. -nocenter-
    // restores the uncentered Arellano-Bond / xtabond2 form for
    // cross-checking Hansen/AR against xtabond2 (difference is O(1/n)).
    local flag_center = cond("`nocenter'" == "", 1, 0)
    // v0.8.1 (audit R6, #5): symmetric coefficient bootstrap CIs when the
    // legacy coefficient-bootstrap replay is explicitly requested.
    local coefcitype = lower(trim("`coefcitype'"))
    if "`coefcitype'" == "" local coefcitype "symmetric"
    if !inlist("`coefcitype'", "symmetric", "percentile") {
        di as err "option coefcitype() must be symmetric or percentile"
        exit 198
    }
    local flag_coefci_sym = cond("`coefcitype'" == "symmetric", 1, 0)
    // v0.9.26: coefficient bootstrap is outside the supported core. Keep the
    // existing implementation as an opt-in legacy compatibility path, but do
    // not run it silently as part of every threshold-CI request.
    local coefboot = lower(trim("`coefboot'"))
    if "`coefboot'" == "" local coefboot "none"
    if "`coefboot'" == "gs" {
        // v0.9.3 R19 (#5): honesty gate. The implemented scheme is a
        // threshold-search-aware cluster wild residual bootstrap -- NOT the
        // Gong-Seo coefficient bootstrap (synchronized unit resampling of
        // regressors/instruments/residuals, moment recentering, rebuilt
        // weights, shrinkage theta0*). gs stays locked until that exists.
        di as err "coefboot(gs) is not implemented: the available scheme is a"
        di as err "threshold-search-aware cluster wild residual bootstrap, not the"
        di as err "Gong-Seo coefficient bootstrap. Use coefboot(twostep|onestep|none)."
        exit 198
    }
    if !inlist("`coefboot'", "twostep", "onestep", "none") {
        di as err "option coefboot() must be twostep, onestep, or none"
        exit 198
    }
    if inlist("`coefboot'", "twostep", "onestep") & "`nowarn'" == "" {
        di as txt "note: coefboot(`coefboot') is a legacy compatibility feature of xtdpthresh;"
        di as txt "      coefficient-bootstrap inference is outside the supported core."
    }
    local flag_coefboot_2s  = cond("`coefboot'" == "twostep", 1, 0)
    local flag_coefboot_off = cond("`coefboot'" == "none", 1, 0)
    // v0.7.13 (audit R4, C3): gridtype(quantile) places the γ grid on
    // empirical quantiles of q over the effective sample (equal observation
    // counts between consecutive points — the Gong-Seo application layout);
    // duplicates from ties are collapsed, so the effective grid may hold
    // fewer points than requested. Default uniform (equally spaced values)
    // preserves the xthenreg-comparable convention.
    local gridtype = lower(trim("`gridtype'"))
    if "`gridtype'" == "" local gridtype "uniform"
    if !inlist("`gridtype'", "uniform", "quantile") {
        di as err "option gridtype() must be uniform or quantile"
        exit 198
    }
    local flag_grid_quant = cond("`gridtype'" == "quantile", 1, 0)
    // v0.8.2 (audit R9, #2): optional user floor for per-regime support
    // counts at each candidate gamma; 0 = default trim-based rule.
    if missing(`minregime') | `minregime' < 0 {
        di as err "minregime() must be >= 0"
        exit 198
    }
    // v0.8.2 R10 (#5): support used for trim bounds and quantile grids.
    //   effective (default): observations entering the effective GMM
    //                        criterion (q_t and q_{t-1} under FD, + future
    //                        equation-row q under FOD; level rows current).
    //   observed:            current-row q of the retained equation rows --
    //                        closer to the xthenreg convention (quantiles
    //                        of the observed threshold variable) for
    //                        replication on balanced FD panels.
    // v0.9.1 R17 (#4): history(panel|sample). L.y and every ts-operator
    // regressor are materialized by tsrevar on the FULL panel before if/in
    // bites, so out-of-scope rows already reach the RHS through lags.
    // history(panel) (default) makes the instrument/lag HISTORY consistent
    // with that: if/in restricts the EQUATION sample only. history(sample)
    // treats if/in as a hard boundary instead and nulls the auto L.y
    // wherever its source row falls outside the history sample.
    // v0.9.25: refine() is a FINAL, jump-only support search under the one
    // fixed W2 constructed after the global stage-1 search.  Unlike the old
    // implementation it never restarts the estimator and never rebuilds W2.
    if `bwscale' <= 0 | missing(`bwscale') {
        di as err "option bwscale() must be positive"
        exit 198
    }
    if `refine' < 0 | `refine' > 20 {
        di as err "option refine() must be an integer between 0 and 20"
        exit 198
    }
    if `refine' > 0 & "`kink'" != "" {
        di as err "refine() is supported for the jump specification only"
        di as err "(the kink regressor (q-gamma)*1(q>gamma) varies continuously in gamma)"
        exit 198
    }
    // v0.9.26: the numerical estimator is the profiled argmin over the
    // finite grid explicitly requested by grid(), followed (for jump models)
    // by the optional final fixed-W2 refine().  Retired adaptive controls are
    // parsed solely to return an informative compatibility error.
    local _legacy_searchopt = ("`searchmode'" != "" | `searchmax' != -1 | ///
        `searchtol' != 1e-8 | strpos(lower(`"`cmdline'"'), "searchtol(") > 0)
    if `_legacy_searchopt' {
        di as err "searchmode(), searchmax(), and searchtol() are no longer supported in xtdpthresh 0.9.26"
        di as err "The threshold search is a fixed finite grid selected by grid(#), with optional refine(#)."
        exit 198
    }
    local searchmode "fixed"
    local searchmax_set 0
    local search_max_level 1
    local searchmax_effective `grid'
    local searchtol .

    local history = lower(trim("`history'"))
    if "`history'" == "" local history "panel"
    if !inlist("`history'", "panel", "sample") {
        di as err "option history() must be panel or sample"
        exit 198
    }
    // v0.9.2 R18 (#4): history(sample) can null the auto L.y, but user
    // ts-operator terms (L.x, D.x, ...) are materialized on the FULL panel
    // by tsrevar and cannot be retro-restricted -- with them, "sample"
    // would be a hard boundary for instruments but not for the RHS.
    // Reject the combination instead of silently half-honoring it.
    if "`history'" == "sample" {
        local _tschk `varlist' `endogenous' `predetermined' `exogenous' `iv'
        if strpos("`_tschk'", ".") {
            di as err "history(sample) cannot retro-restrict time-series-operator terms"
            di as err "(L.x, D.x, ...): they are materialized on the full panel before"
            di as err "if/in applies. Pre-generate the lagged variables as plain"
            di as err "variables, or use history(panel)."
            exit 198
        }
    }

    // v0.9.26: boottype(wild|unit). wild is the supported fast cluster wild
    // residual bootstrap (default; a computational approximation of
    // Gong-Seo Alg. 1, see Remarks). unit = EXPERIMENTAL unit-multiplicity
    // resampling ORIENTED at Gong-Seo Alg. 1 (unrestricted-residual DGP,
    // recentering at theta-hat, fixed sample W1 and per-draw Omega/W2*) but NOT
    // certified against the paper -- hence not named "exact". In 0.9.26 it
    // is retained strictly as a reproducibility/verification path. Threshold-CI
    // inversion only (linearity/continuity/coefficient bootstraps keep the
    // wild scheme); fd without kink only (the Alg. 1 theory is developed
    // for the first-differenced jump estimator).
    local boottype = lower(trim("`boottype'"))
    if "`boottype'" == "" local boottype "wild"
    if "`boottype'" == "exact" {
        di as err "boottype(exact) has been renamed boottype(unit): the unit-resampling"
        di as err "scheme is Alg. 1-oriented but not certified as the exact algorithm"
        exit 198
    }
    if !inlist("`boottype'", "wild", "unit") {
        di as err "option boottype() must be wild or unit"
        exit 198
    }
    local flag_boot_exact = cond("`boottype'" == "unit", 1, 0)
    if "`noboot'" == "" & `flag_boot_exact' {
        local _bt_m = lower(trim("`method'"))
        if !inlist("`_bt_m'", "", "fd") {
            di as err "boottype(unit) currently supports method(fd) only (Gong-Seo Alg. 1"
            di as err "is developed for the first-differenced estimator)"
            exit 198
        }
        if "`kink'" != "" {
            di as err "boottype(unit) does not support kink (Alg. 1 targets the"
            di as err "unrestricted jump estimator)"
            exit 198
        }
    }

    local gridsample = lower(trim("`gridsample'"))
    if "`gridsample'" == "" local gridsample "effective"
    if !inlist("`gridsample'", "effective", "observed") {
        di as err "option gridsample() must be effective or observed"
        exit 198
    }
    local flag_notest = cond("`notest'" != "", 1, 0)

    // === Parse iv() with sub-options: iv(varlist [, collapse]) ===
    // v0.7.0 semantics fix: the collapse sub-option no longer overrides the
    // top-level collapse; it collapses ONLY the user-IV block. The maxlag()
    // sub-option is rejected with an explanation — it never affected user IVs
    // (they always enter as the period-t value) but silently overrode the
    // top-level maxlag() for GMM-style instruments, a serious silent trap.
    // levmaxlag() belonged exclusively to the retired System-GMM level
    // block. Check it before the nested iv() parser can clobber option locals.
    if "`levmaxlag'" != "" {
        di as err "levmaxlag() was retired with method(system); use maxlag() for FD/FOD moments"
        exit 198
    }
    local flag_ivcol_sub 0
    if "`iv'" != "" {
        // Save outer locals — the inner `syntax' call below clobbers `varlist',
        // `if', `in', and option locals such as `maxlag' / `collapse'.
        local _save_varlist    `"`varlist'"'
        local _save_if         `"`if'"'
        local _save_in         `"`in'"'
        local _outer_maxlag    `"`maxlag'"'
        local _outer_collapse  `"`collapse'"'

        // Robust parser for iv(z1 z2 [, collapse]).
        local 0 `"`iv'"'
        // v0.9.31: keep Stata's own message and return code (111 for a
        // variable that does not exist, 109 for a string variable, 198 for
        // a bad option); every error used to be reported as 198.
        cap noisily syntax varlist(numeric ts) [, MAXLAG(numlist max=2 min=1 integer >0) COLLAPSE]
        if _rc {
            local _iv_rc = _rc
            di as err "  (in iv(); use iv(varlist [, collapse]))"
            exit `_iv_rc'
        }
        local iv_vars `"`varlist'"'
        if "`maxlag'" != "" {
            di as err "iv() sub-option maxlag() is not supported (v0.7.0):"
            di as err "user-supplied IVs always enter as their period-t value, so"
            di as err "maxlag() never affected them — it only (silently) overrode the"
            di as err "top-level maxlag() controlling GMM-style instruments. Use the"
            di as err "top-level maxlag() option instead."
            exit 198
        }
        local flag_ivcol_sub = cond("`collapse'" != "", 1, 0)

        // Restore outer positional/if/in and option locals clobbered by the
        // inner syntax call. Top-level maxlag()/collapse are authoritative.
        local varlist  `"`_save_varlist'"'
        local if       `"`_save_if'"'
        local in       `"`_save_in'"'
        local maxlag   `"`_outer_maxlag'"'
        local collapse `"`_outer_collapse'"'

        local iv `iv_vars'
    }

    marksample touse, novarlist

    // === Option normalization & validation ===
    // v0.7.12: normalize case so method(FOD)/method(FD) etc. are accepted
    // rather than rejected by the case-sensitive inlist checks below.
    local method = lower(trim("`method'"))
    if "`method'" == "" local method "fd"
    if "`method'" == "system" {
        di as err "method(system) is no longer supported as of xtdpthresh 0.9.25"
        di as err "The supported estimators are method(fd) and method(fod)."
        exit 198
    }
    if !inlist("`method'", "fd", "fod") {
        di as err "option method() must be fd or fod"
        exit 198
    }
    if "`levmaxlag'" != "" {
        di as err "levmaxlag() was retired with method(system); use maxlag() for FD/FOD moments"
        exit 198
    }
    // v0.8.0 (audit R5): citype() removed. It duplicated -noboot- exactly
    // (citype(none) == noboot; grid was the only other value and the default)
    // and its reserved extension citype(asym) will never exist -- asymptotic
    // threshold CIs are invalid under continuity (Gong-Seo 2026).

    // v0.7.13/0.8.0: the grid bootstrap runs unless -noboot-. gridci()/
    // boot()/rseed() are irrelevant on the point-estimate-only path, so they
    // are neither validated nor applied there. grid() always matters (it
    // sets the gamma search grid for the point estimate too).
    local _will_boot = ("`noboot'" == "")
    if missing(`grid') | `grid' < 10 {
        di as err "grid() must be at least 10"
        exit 198
    }
    if `_will_boot' & (missing(`gridci') | `gridci' < 10) {
        di as err "gridci() must be at least 10"
        exit 198
    }
    // v0.9.36: citest(#) runs the grid-bootstrap test of H0: gamma = # alone
    // (a diagnostic: its rejection rate at the true threshold separates the
    // size of the test from the discretization of the confidence set)
    if "`citest'" != "" & !`_will_boot' {
        di as err "citest() requires the grid bootstrap (remove noboot)"
        exit 198
    }
    if `_will_boot' & `gridci' < 100 & "`nowarn'" == "" {
        di as txt "note: gridci(`gridci') uses a coarse CI-inversion grid; gridci(100) or larger is recommended for final threshold inference"
    }
    if missing(`trim') | `trim' < 0.01 | `trim' > 0.45 {
        di as err "trim must be in [0.01, 0.45]"
        exit 198
    }
    if `_will_boot' & (missing(`boot') | `boot' < 10) {
        di as err "boot should be at least 10 (99+ recommended for production)"
        exit 198
    }
    // Validate rseed() without changing the caller's RNG state. The seed is
    // applied inside Mata immediately before the first bootstrap draw, after
    // point estimation and all unit-bootstrap preflight checks have passed.
    if "`rseed'" != "" & `_will_boot' {
        cap confirm integer number `rseed'
        if _rc {
            di as err "rseed() must be an integer in [0, 2147483647]"
            exit 198
        }
        if `rseed' < 0 | `rseed' > 2147483647 {
            di as err "rseed() must be an integer in [0, 2147483647]"
            exit 198
        }
    }

    local flag_kink = cond("`kink'" != "", 1, 0)
    local flag_static = cond("`static'" != "", 1, 0)
    local flag_collapse = cond("`collapse'" != "", 1, 0)
    // v0.7.0: user-IV collapse flag — set by the iv(, collapse) sub-option;
    // top-level collapse implies it (collapsing everything includes user IVs).
    local flag_iv_collapse = cond(`flag_ivcol_sub' | `flag_collapse', 1, 0)

    // Parse maxlag(# [#]): interval form a..b; default unlimited
    local n_ml : word count `maxlag'
    if `n_ml' == 0 {
        local maxlag_lo = 1
        // A true open upper bound. The old 9999 sentinel silently omitted
        // valid calendar lags on very long/sparse delta-1 time indexes.
        // Mata caps this value to the retained history span before any loop.
        local maxlag_hi = 1e300
    }
    else if `n_ml' == 1 {
        local maxlag_lo = 1
        local maxlag_hi : word 1 of `maxlag'
    }
    else {
        local maxlag_lo : word 1 of `maxlag'
        local maxlag_hi : word 2 of `maxlag'
    }
    if `maxlag_lo' > `maxlag_hi' {
        di as err "maxlag(# #) requires min <= max"
        exit 198
    }

    // === Parse varlist (xthreg2-style self-documenting syntax) ===
    // Syntax: xtdpthresh depvar [indepvars], qx(threshold_var) [options]
    // q_var is ALWAYS via qx() — explicit and required.
    gettoken depvar indepvars : varlist
    local indepvars = trim("`indepvars'")
    local q_var "`qx'"

    // Dependent variable must be a plain variable name. Time-series operators
    // are allowed for regressors/options only and are expanded below via tsrevar.
    if strpos("`depvar'", ".") {
        di as err "dependent variable may not contain time-series operators; create the lag/difference as a separate variable if needed"
        exit 198
    }

    // Save user-facing lists before tsrevar expansion. Expanded temporary
    // variable names are used internally; these labels are used for display,
    // ereturn metadata, and coefficient names.
    local indepvars_lab `"`indepvars'"'

    // Note: Seo-Shin (2016) permits q_var to appear as a regressor (in indepvars
    // or endogenous()). In that case the slope on q changes at γ — a natural
    // specification (e.g. "does the marginal effect of debt on invest shift
    // above some debt threshold?"). No overlap check needed.

    // === Collect all regressors ===
    // Exog regressors (default: all indepvars are exogenous unless endogenous()
    //                  or predetermined() is specified)
    // Endog regressors     (option endogenous())    : instruments from t-2
    // Predet regressors    (option predetermined()) : instruments from t-1
    // None of these three groups may overlap with each other or with indepvars.
    local exog_extra : list clean exogenous
    local endog      : list clean endogenous
    local predet     : list clean predetermined
    local inst_extra : list clean iv

    // User-facing copies after option macros are normalized.
    local exog_extra_lab `"`exog_extra'"'
    local endog_lab      `"`endog'"'
    local predet_lab     `"`predet'"'
    local inst_extra_lab `"`inst_extra'"'

    // Check no overlap between endogenous and indepvars
    local overlap : list endog & indepvars
    if "`overlap'" != "" {
        di as err "endogenous() vars must not appear in indepvars: `overlap'"
        exit 198
    }
    // Check no overlap between endogenous and exogenous
    local overlap2 : list endog & exog_extra
    if "`overlap2'" != "" {
        di as err "endogenous() and exogenous() vars must not overlap: `overlap2'"
        exit 198
    }
    // Check no overlap between indepvars and exogenous (would duplicate in X)
    local overlap3 : list indepvars & exog_extra
    if "`overlap3'" != "" {
        di as err "indepvars and exogenous() vars must not overlap: `overlap3'"
        exit 198
    }
    // Check no overlap between predetermined and (indepvars / exog / endog)
    local overlap4 : list predet & indepvars
    if "`overlap4'" != "" {
        di as err "predetermined() vars must not appear in indepvars: `overlap4'"
        exit 198
    }
    local overlap5 : list predet & exog_extra
    if "`overlap5'" != "" {
        di as err "predetermined() and exogenous() vars must not overlap: `overlap5'"
        exit 198
    }
    local overlap6 : list predet & endog
    if "`overlap6'" != "" {
        di as err "predetermined() and endogenous() vars must not overlap: `overlap6'"
        exit 198
    }
    // An explicitly endogenous/predetermined regressor cannot simultaneously
    // be inserted as its own contemporaneous external IV. That contradicts
    // the declared timing status and silently creates invalid moments.
    local overlap7 : list inst_extra & endog
    local overlap8 : list inst_extra & predet
    if "`overlap7'`overlap8'" != "" {
        di as err "iv() may not contain variables declared endogenous() or predetermined()"
        di as err "  conflicting variables: `overlap7' `overlap8'"
        exit 198
    }

    local k_exog   : word count `indepvars' `exog_extra'
    local k_endog  : word count `endog'
    local k_predet : word count `predet'
    local k_inst   : word count `inst_extra'

    // v0.7.13 (audit): a static model with no regressors leaves only the
    // block-constant moment; the transformed-equation availability filter
    // then drops every row (no data-driven instrument), losing the whole
    // sample. Reject it with a clear message rather than failing opaquely.
    // (The dynamic default always has L.`depvar', so this only bites an
    // explicit -static- with an empty regressor list.)
    if `flag_static' & `k_exog' == 0 & `k_endog' == 0 & `k_predet' == 0 & `k_inst' == 0 {
        di as err "a static model requires at least one regressor or external iv()"
        di as err "  add regressors/an external iv(), or drop -static- to use the dynamic L.`depvar' model"
        exit 198
    }

    // v0.7.12: under kink the level term is [X, (q-gamma)*1(q>gamma)]. For a
    // two-sided kink (a baseline slope on q below gamma plus a slope change
    // above), q must ALSO be a regressor. If qx() is absent from the RHS the
    // model reduces to a one-sided hinge (flat in q below gamma); warn unless
    // nowarn. Checked pre-tsrevar on the user's original names. A base term
    // equal to q -- possibly written with a no-op ts operator like L0.q --
    // supplies the baseline slope; a genuine lag (L.q) does not. Match on the
    // contemporaneous level form: strip a single #0. no-op operator, then
    // compare to q_var, so e.g. L0.q does NOT trigger a false warning.
    if `flag_kink' & "`nowarn'" == "" {
        local _qrhs 0
        foreach _t in `indepvars' `exog_extra' `endog' `predet' {
            local _tc = regexr("`_t'", "^[LFDSlfds]0\.", "")
            if "`_tc'" == "`q_var'" local _qrhs 1
        }
        if !`_qrhs' {
            di as text "{err}Warning:{txt} " as res "kink" as text " specified but " ///
               as res "`q_var'" as text " is not a base regressor (indepvars, "       ///
               as res "exogenous()" as text ", " as res "endogenous()" as text ", or " ///
               as res "predetermined()" as text ")."
            di as text "  This estimates a one-sided hinge (q-γ)·1(q>γ), not a" ///
               as text " two-sided slope-kink in " as res "`q_var'" as text "."
            di as text "  Add " as res "`q_var'" as text " to the RHS for the standard kink" ///
               as text " model (" as res "nowarn" as text " hides this)."
            di ""
        }
    }

    // A maxlag() interval ending at 1 supplies no internal lagged-level
    // moments for L.y/endogenous regressors. Do not reject it outright:
    // transformed exogenous/predetermined moments or external iv() variables
    // can still identify the model, and the downstream rank/conditioning
    // gates fail closed when they do not.
    if `maxlag_hi' < 2 & (!`flag_static' | `k_endog' > 0) {
        if "`nowarn'" == "" {
            di as text "Note: maxlag() supplies no internal lagged-level moments for"
            di as text "L.`depvar' or endogenous() regressors; identification must come"
            di as text "from exogenous/predetermined moments or external iv() variables."
        }
    }

    // === xtset check ===
    capture xtset
    if _rc {
        di as err "must xtset panelid timevar before using xtdpthresh"
        exit 459
    }
    local panelvar = r(panelvar)
    local timevar  = r(timevar)
    local tdelta = r(tdelta)
    // v0.7.13 (audit): panel-only -xtset id- leaves r(timevar) as "." (a
    // literal dot), not "", and r(tdelta) missing; the old code then fell
    // through to the delta!=1 message, which misdiagnosed the problem. Catch
    // both the empty and dot forms.
    if "`timevar'" == "" | "`timevar'" == "." {
        di as err "must xtset panelid timevar (a time variable is required, e.g. xtset `panelvar' year)"
        exit 459
    }
    // Time delta validation follows.
    if r(tdelta) != 1 {
        di as err "xtdpthresh currently requires an xtset time delta of 1"
        exit 459
    }
    local is_balanced = ("`r(balanced)'" == "strongly balanced")

    // === Expand time-series operators for Mata st_data() =======================
    // Stata's parser can accept numeric ts varlists, but Mata's st_data() cannot
    // reliably read expressions such as L.x, D.x, or L(1/2).x directly. tsrevar
    // materializes them as temporary variables. All downstream data handling uses
    // the expanded names, while *_lab locals preserve the user's original syntax.
    if `"`indepvars'"' != "" {
        local _old_type "`c(type)'"
        quietly set type double
        cap tsrevar `indepvars'
        local _ts_rc = _rc
        if !`_ts_rc' local _expanded `"`r(varlist)'"'
        quietly set type `_old_type'
        if `_ts_rc' {
            di as err "could not expand indepvars with time-series operators"
            exit `_ts_rc'
        }
        local _n_user : word count `indepvars_lab'
        local _n_exp  : word count `_expanded'
        if `_n_user' != `_n_exp' {
            di as err "range time-series operators such as L(1/2).x are not supported here; spell them out as separate terms"
            exit 198
        }
        local indepvars `"`_expanded'"'
        _xdpt_tsdouble, expanded(`indepvars') labels(`indepvars_lab')
    }
    if `"`exog_extra'"' != "" {
        local _old_type "`c(type)'"
        quietly set type double
        cap tsrevar `exog_extra'
        local _ts_rc = _rc
        if !`_ts_rc' local _expanded `"`r(varlist)'"'
        quietly set type `_old_type'
        if `_ts_rc' {
            di as err "could not expand exogenous() with time-series operators"
            exit `_ts_rc'
        }
        local _n_user : word count `exog_extra_lab'
        local _n_exp  : word count `_expanded'
        if `_n_user' != `_n_exp' {
            di as err "range time-series operators such as L(1/2).x are not supported here; spell them out as separate terms"
            exit 198
        }
        local exog_extra `"`_expanded'"'
        _xdpt_tsdouble, expanded(`exog_extra') labels(`exog_extra_lab')
    }
    if `"`endog'"' != "" {
        local _old_type "`c(type)'"
        quietly set type double
        cap tsrevar `endog'
        local _ts_rc = _rc
        if !`_ts_rc' local _expanded `"`r(varlist)'"'
        quietly set type `_old_type'
        if `_ts_rc' {
            di as err "could not expand endogenous() with time-series operators"
            exit `_ts_rc'
        }
        local _n_user : word count `endog_lab'
        local _n_exp  : word count `_expanded'
        if `_n_user' != `_n_exp' {
            di as err "range time-series operators such as L(1/2).x are not supported here; spell them out as separate terms"
            exit 198
        }
        local endog `"`_expanded'"'
        _xdpt_tsdouble, expanded(`endog') labels(`endog_lab')
    }
    if `"`predet'"' != "" {
        local _old_type "`c(type)'"
        quietly set type double
        cap tsrevar `predet'
        local _ts_rc = _rc
        if !`_ts_rc' local _expanded `"`r(varlist)'"'
        quietly set type `_old_type'
        if `_ts_rc' {
            di as err "could not expand predetermined() with time-series operators"
            exit `_ts_rc'
        }
        local _n_user : word count `predet_lab'
        local _n_exp  : word count `_expanded'
        if `_n_user' != `_n_exp' {
            di as err "range time-series operators such as L(1/2).x are not supported here; spell them out as separate terms"
            exit 198
        }
        local predet `"`_expanded'"'
        _xdpt_tsdouble, expanded(`predet') labels(`predet_lab')
    }
    if `"`inst_extra'"' != "" {
        local _old_type "`c(type)'"
        quietly set type double
        cap tsrevar `inst_extra'
        local _ts_rc = _rc
        if !`_ts_rc' local _expanded `"`r(varlist)'"'
        quietly set type `_old_type'
        if `_ts_rc' {
            di as err "could not expand iv() variables with time-series operators"
            exit `_ts_rc'
        }
        local _n_user : word count `inst_extra_lab'
        local _n_exp  : word count `_expanded'
        if `_n_user' != `_n_exp' {
            di as err "range time-series operators such as L(1/2).x are not supported here; spell them out as separate terms"
            exit 198
        }
        local inst_extra `"`_expanded'"'
        _xdpt_tsdouble, expanded(`inst_extra') labels(`inst_extra_lab')
    }

    // Re-check overlap AFTER tsrevar expansion. Textually different terms can
    // resolve to the same variable (for example x and L0.x); allowing them in
    // different status groups silently duplicates columns and IV status.
    local overlap : list endog & indepvars
    local overlap2 : list endog & exog_extra
    local overlap3 : list indepvars & exog_extra
    local overlap4 : list predet & indepvars
    local overlap5 : list predet & exog_extra
    local overlap6 : list predet & endog
    local overlap7 : list inst_extra & endog
    local overlap8 : list inst_extra & predet
    if "`overlap'`overlap2'`overlap3'`overlap4'`overlap5'`overlap6'" != "" {
        di as err "regressor groups overlap after time-series expansion"
        di as err "  use each expanded variable in exactly one regressor group"
        exit 198
    }
    if "`overlap7'`overlap8'" != "" {
        di as err "iv() conflicts with endogenous()/predetermined() after time-series expansion"
        di as err "  use a variable in iv() only when its contemporaneous value is exogenous"
        exit 198
    }

    // The automatic continuity comparison is nested only when q enters the
    // unrestricted jump model contemporaneously: q*r and r can then impose
    // (q-gamma)*r. A lag of q is not a substitute. Expansion canonicalizes
    // no-op operators such as L0.q to q, so this test is exact.
    local _rhs_expanded "`indepvars' `exog_extra' `endog' `predet'"
    local _q_rhs : list q_var in _rhs_expanded
    local flag_cont_test = cond(!`flag_kink' & `_q_rhs', 1, 0)
    if !`flag_kink' & !`_q_rhs' & `_will_boot' & !`flag_notest' & "`nowarn'" == "" {
        di as text "Note: continuity test omitted because " as res "`q_var'" ///
            as text " is not a contemporaneous RHS regressor;"
        di as text "the kink model would not be nested in the estimated jump model."
    }

    // v0.7.13 (audit): duplicates WITHIN one group also survive -syntax-
    // (e.g. "L.x l1.x" canonicalize to the same term; "iv(z z)") and resolve
    // to the same expanded variable, silently producing collinear columns
    // (a coefficient of 0 with no "omitted" note) and inflating the
    // instrument count behind the Hansen J df. Reject them.
    local _grp_names `""indepvars" "exogenous()" "endogenous()" "predetermined()" "iv()""'
    local _gi = 0
    foreach _g in indepvars exog_extra endog predet inst_extra {
        local ++_gi
        local _dups : list dups `_g'
        if "`_dups'" != "" {
            local _gn : word `_gi' of `_grp_names'
            di as err "duplicate variables in `_gn' after time-series expansion"
            di as err "  each variable may appear only once within a list"
            exit 198
        }
    }

    // v0.7.13 (audit): the dependent variable may not be its own regressor,
    // instrument, or threshold. In particular endogenous(`depvar') is a
    // plausible misreading ("declare the depvar endogenous") that would put
    // y_t on its own RHS and produce degenerate GMM estimates silently.
    local _all_rhs "`indepvars' `exog_extra' `endog' `predet' `inst_extra'"
    local _dv_hit : list depvar in _all_rhs
    if `_dv_hit' {
        di as err "the dependent variable may not appear in indepvars, exogenous(), endogenous(), predetermined(), or iv()"
        di as err "  (the dynamic model adds L.`depvar' automatically; use static to suppress it)"
        exit 198
    }
    if "`q_var'" == "`depvar'" {
        di as err "qx() may not be the dependent variable"
        di as err "  for a self-exciting threshold, create the lag first: gen Lq = L.`depvar'"
        exit 198
    }

    // Dynamic models already add L.depvar with lagged-level instruments.
    // Adding it again creates two identical regressor columns.
    if !`flag_static' {
        local _user_terms "`indepvars_lab' `endog_lab' `predet_lab' `exog_extra_lab'"
        foreach _ul of local _user_terms {
            // v0.7.13 (audit): compare the operator case-insensitively but
            // the variable name case-SENSITIVELY. Stata variable names are
            // case-sensitive: in a dataset holding both Y and y, "L.y" is a
            // different variable from the auto-added "L.Y" and must not be
            // rejected. (The old code lowercased both sides.)
            local _dot = strpos("`_ul'", ".")
            if `_dot' {
                local _op  = lower(substr("`_ul'", 1, `_dot'-1))
                local _bas = substr("`_ul'", `_dot'+1, .)
                if "`_bas'" == "`depvar'" {
                    // v0.9.19: reduce every pure L/F operator chain to its
                    // net shift. Stata accepts equivalent spellings such as
                    // FL2.y and L2F.y; both equal L.y and used to evade the
                    // four-literal gate, silently duplicating auto L.y.
                    local _rest "`_op'"
                    local _netlag = 0
                    local _purelf = 1
                    while "`_rest'" != "" & `_purelf' {
                        local _tok ""
                        local _tn = .
                        if regexm("`_rest'", "^([lf])\(([0-9]+)/([0-9]+)\)") {
                            local _tok = regexs(1)
                            local _n1  = real(regexs(2))
                            local _n2  = real(regexs(3))
                            local _hit = regexs(0)
                            if `_n1' != `_n2' local _purelf = 0
                            else local _tn = `_n1'
                        }
                        else if regexm("`_rest'", "^([lf])\(([0-9]+)\)") {
                            local _tok = regexs(1)
                            local _tn  = real(regexs(2))
                            local _hit = regexs(0)
                        }
                        else if regexm("`_rest'", "^([lf])([0-9]+)") {
                            local _tok = regexs(1)
                            local _tn  = real(regexs(2))
                            local _hit = regexs(0)
                        }
                        else if regexm("`_rest'", "^([lf])") {
                            local _tok = regexs(1)
                            local _tn  = 1
                            local _hit = regexs(0)
                        }
                        else local _purelf = 0
                        if `_purelf' {
                            if "`_tok'" == "l" local _netlag = `_netlag' + `_tn'
                            else                    local _netlag = `_netlag' - `_tn'
                            local _rest = substr("`_rest'", strlen("`_hit'") + 1, .)
                        }
                    }
                    if `_purelf' & `_netlag' == 1 {
                        di as err "L.`depvar' is added automatically in the dynamic model"
                        di as err "  `_ul' is algebraically the same lag; remove it, or specify static"
                        exit 198
                    }
                }
            }
        }
    }

    // v0.9.28/0.9.29: classification guards for time-series-operated terms.
    // (a) The dependent variable. In indepvars/exogenous() every operator on
    //     depvar is rejected: a lag depends on earlier errors, so it is not
    //     strictly exogenous (it belongs in predetermined()); D., S., and F.
    //     terms contain the current or a future value of depvar. Under
    //     strict exogeneity such a term would be instrumented by its own
    //     transformed value, which is invalid, with no warning. In
    //     endogenous(), predetermined(), and iv() only pure lags (net lag of at
    //     least 1, no D., S., or F.) are allowed.
    // (b) A strictly exogenous term may not be a lag (or the same date) of a
    //     variable declared in endogenous()/predetermined(): if x_{t-k} may
    //     respond to past errors, so may x_{t-j} for j >= k, and those lags
    //     are not strictly exogenous. The reverse -- x strictly exogenous and
    //     a further lag such as L.x declared predetermined() -- is only
    //     conservative and is allowed. Terms with D., S., or F. are treated
    //     conservatively as containing the contemporaneous value.
    // Terms are decomposed by _xdpt_tsterm (defined after this program).
    local _ep_bases ""
    local _ep_minlag ""
    local _ep_kind ""
    foreach _ul in `endog_lab' {
        _xdpt_tsterm `_ul'
        local _ep_bases  "`_ep_bases' `r(base)'"
        local _ep_minlag "`_ep_minlag' `=cond(r(hasdsf), min(r(netlag), 0), r(netlag))'"
        local _ep_kind   "`_ep_kind' en"
    }
    foreach _ul in `predet_lab' {
        _xdpt_tsterm `_ul'
        local _ep_bases  "`_ep_bases' `r(base)'"
        local _ep_minlag "`_ep_minlag' `=cond(r(hasdsf), min(r(netlag), 0), r(netlag))'"
        local _ep_kind   "`_ep_kind' pr"
    }
    local _n_ep : word count `_ep_bases'
    foreach _grp in sx ep iv {
        if "`_grp'" == "sx" local _terms "`indepvars_lab' `exog_extra_lab'"
        if "`_grp'" == "ep" local _terms "`endog_lab' `predet_lab'"
        if "`_grp'" == "iv" local _terms "`inst_extra_lab'"
        foreach _ul of local _terms {
            _xdpt_tsterm `_ul'
            local _b   "`r(base)'"
            local _isop = r(isop)
            local _nl   = r(netlag)
            local _dsf  = r(hasdsf)
            if "`_b'" == "`depvar'" & `_isop' {
                if "`_grp'" == "sx" {
                    if !`_dsf' & `_nl' >= 1 {
                        di as err "`_ul' is a lag of the dependent variable and cannot be strictly exogenous"
                        di as err "  declare it in predetermined() instead of indepvars or exogenous()"
                    }
                    else {
                        di as err "`_ul' contains the current or a future value of the dependent variable"
                        di as err "  differences, seasonal differences, and leads of `depvar' cannot be regressors"
                    }
                    exit 198
                }
                if `_dsf' | `_nl' < 1 {
                    di as err "`_ul' contains the current or a future value of the dependent variable"
                    if "`_grp'" == "iv" {
                        di as err "  it cannot be an instrument; only lags of `depvar' may appear in iv()"
                    }
                    else {
                        di as err "  it cannot be a regressor; lags of `depvar' belong in predetermined()"
                    }
                    exit 198
                }
            }
            // v0.9.34: an iv() term must be uncorrelated with the transformed
            // error: under FD with e(t) and e(t-1), under FOD with e(t) and
            // later errors. A lag of order j of depvar or of an endogenous
            // variable contains e(t-j), and a predetermined variable may
            // respond to the previous error. So under FD the net lag must be
            // at least 2 for depvar and endogenous bases and 1 for
            // predetermined bases, counted from the lag at which the base is
            // declared; under FOD, 1 and 0. iv(L.depvar) under FD, accepted
            // before, biased the estimates without a Hansen-test signal.
            if "`_grp'" == "iv" {
                local _need = .
                if "`_b'" == "`depvar'" local _need = cond("`method'" == "fd", 2, 1)
                forvalues _k = 1/`_n_ep' {
                    local _eb : word `_k' of `_ep_bases'
                    if "`_eb'" == "`_b'" {
                        local _em : word `_k' of `_ep_minlag'
                        local _ek : word `_k' of `_ep_kind'
                        local _rq = `_em' + cond("`_ek'" == "en", ///
                            cond("`method'" == "fd", 2, 1), cond("`method'" == "fd", 1, 0))
                        if `_need' >= . | `_rq' > `_need' local _need = `_rq'
                    }
                }
                local _nle = cond(`_dsf', min(`_nl', 0), `_nl')
                if `_need' < . & `_nle' < `_need' {
                    di as err "`_ul' is not a valid instrument under method(`method'): it is"
                    di as err "  correlated with the transformed error; instruments built on"
                    di as err "  `_b' need a lag of at least `_need' here"
                    exit 198
                }
            }
            if "`_grp'" == "sx" & `_n_ep' > 0 {
                local _hit = 0
                forvalues _k = 1/`_n_ep' {
                    local _eb : word `_k' of `_ep_bases'
                    local _em : word `_k' of `_ep_minlag'
                    if "`_eb'" == "`_b'" {
                        if `_dsf' local _hit = 1
                        else if `_em' <= `_nl' local _hit = 1
                    }
                }
                if `_hit' {
                    di as err "`_ul' is declared strictly exogenous, but it is a lag (or the same"
                    di as err "  date) of `_b', which is declared in endogenous() or predetermined();"
                    di as err "  declare `_ul' in predetermined() or endogenous() instead"
                    exit 198
                }
            }
        }
    }

    // === Time-effect treatment (td) ============================================
    // v0.7.13 (audit R4, C2): two treatments.
    //   td      -> FWL-CORRECT: common-across-regime time dummies are
    //              partialled out of the FINAL transformed system. dY, every
    //              column of dW(γ) — including 1(q>γ) and the interactions —
    //              and every column of Z are cross-sectionally demeaned
    //              within each time cell AFTER stacking, per γ for dW (dY/Z
    //              are γ-invariant). Algebraically identical to including
    //              the dummies in both W and Z and applying FWL: M_t[x·1(q>γ)],
    //              not M_t(x)·1(q>γ). Exact under FD (Δλ_t is common at each
    //              t by construction) AND under FOD on any panel: the time
    //              dummies receive each unit's own FOD operator and are
    //              partialled out by projection (v0.8.0 #6 fix).
    // Both leave q untouched so γ retains its interpretation.
    // v0.8.0 (audit R5): tdpurge (the legacy pre-demeaning construction)
    // REMOVED from the public surface. It implemented exactly the
    // M_t(x)*1(q>gamma) construction the review identified as not equivalent
    // to time dummies; the package is pre-release, so no user results depend
    // on it. td (FWL-correct) is the only time-effects treatment.
    local flag_td_fwl   = cond("`td'" != "", 1, 0)
    local flag_td = `flag_td_fwl'
    // v0.9.21 R42: td is FWL-correct on the estimation sample, but the
    // projection itself changes when units are resampled. The experimental
    // unit threshold bootstrap currently reweights the already-projected
    // cache and therefore cannot reproduce that draw-specific projection.
    // Likewise, the wild TWO-step coefficient replay would need to project
    // each draw's residual again before constructing its cluster Omega*.
    // The wild threshold/test bootstraps and coefboot(onestep) use only
    // global moments with the projected Z, so those combinations remain
    // algebraically valid. Fail fast rather than report mislabelled draws;
    // -noboot- remains a point-estimation-only escape hatch.
    if `_will_boot' & `flag_td' & `flag_boot_exact' {
        di as err "boottype(unit) is not available with td: unit resampling changes"
        di as err "the time-effects FWL projection in every draw, which is not replayed."
        di as err "Use boottype(wild), or remove td."
        exit 198
    }
    if `_will_boot' & `flag_td' & "`coefboot'" == "twostep" {
        di as err "coefboot(twostep) is not available with td: the two-step replay"
        di as err "requires recomputing the time-effects FWL projection in every draw."
        di as err "Use coefboot(onestep) or coefboot(none)."
        exit 198
    }
    // v0.8.0 (audit R5 #6): td is now EXACT under method(fod) on unbalanced
    // panels too -- the stacked system is partialled out on the per-unit
    // FOD-transformed time dummies (the exact operator the data received),
    // not on naive within-time means. No approximation note needed.

    // === Auto-add L.y as first regressor if dynamic (non-static) ===
    // This matches xthenreg: dynamic model automatically includes L.y
    tempvar Ly
    if !`flag_static' {
        // v0.9.7: an untyped -generate-
        // creates a FLOAT -- the auto lag then loses precision for large-
        // magnitude depvars (|y| ~ 1e9+: GDP, VND-denominated series) and
        // contaminates every FD/FOD difference built from it. double it.
        qui gen double `Ly' = L.`depvar'
    }

    // === Sample handling: v0.8.1 (audit R6, #1) SPLIT SAMPLES ===
    // EQUATION sample (touse): strict complete-case rows, INCLUDING the
    // auto L.depvar for dynamic models -- only these rows form GMM
    // equations. HISTORY sample (touse_hist): every in-scope row, kept so
    // that lagged LEVELS remain available as instrument sources. This
    // matches xthenreg, whose Mata reshapes the complete in-scope y matrix
    // and constructs L.y internally, so y_i1 IS an instrument (the previous
    // keep-if-touse deleted each unit's first row, losing y_i1 and the
    // earliest first-difference equation). Instrument reads are
    // value-guarded, so partially-missing history rows are safe.
    // if/in defines BOTH the equation sample and the available history:
    // v0.9.1 R17 (#4): under history(sample), observations excluded by
    // if/in cannot serve as lagged instrument sources; under the default
    // history(panel) the full panel is the history (consistent with how
    // the materialized lag regressors are built).
    // v0.8.3 R12 (#7): -marksample ..., novarlist- leaves the panel/time
    // keys unmarked, and they are the join keys for unit construction, lag
    // histories, and predict's residual merge. xtset usually guarantees
    // them nonmissing, but a production command should not rely on that
    // implicitly. Done BEFORE the history copy so both samples inherit it.
    markout `touse' `panelvar' `timevar'
    tempvar touse_hist
    if "`history'" == "panel" {
        // Design A: the equation sample honors if/in; lag histories and
        // instrument sources come from the full panel, matching the
        // pre-restriction materialization of L.y / ts-operator regressors.
        qui gen byte `touse_hist' = 1
        markout `touse_hist' `panelvar' `timevar'
    }
    else {
        // Design B: if/in is a hard sample boundary.
        qui gen byte `touse_hist' = `touse'
        // The auto lag was materialized BEFORE if/in bit, so a lag whose
        // source row lies outside the history sample must be made missing
        // -- otherwise out-of-scope data enters the RHS while being banned
        // from the instrument history (inconsistent semantics). User-typed
        // ts-operator terms are materialized on the full panel and CANNOT
        // be retro-restricted here; the help documents this and recommends
        // history(panel) when if/in is combined with such terms.
        if !`flag_static' {
            qui replace `Ly' = . if L.`touse_hist' != 1
        }
    }
    // v0.9.29: iv() variables are instruments only and no longer define the
    // equation sample. A missing value contributes a zero instrument in that
    // row (the Z builder reads iv() values under a nonmissing guard), exactly
    // as for the internal lag instruments. Before, a missing value removed
    // the equation, and under FD also the next period's equation.
    markout `touse' `depvar' `q_var' `indepvars' `endog' `predet' `exog_extra'
    if !`flag_static' {
        markout `touse' `Ly'
    }
    quietly count if `touse'
    if r(N) == 0 {
        di as err "no usable complete-case observations in the requested sample"
        exit 2000
    }

    // v0.9.29: regressors that the transformation (or td) removes exactly.
    // FD turned such a column into exact zeros and every candidate threshold
    // failed with a generic message; FOD left rounding noise that the
    // estimator then fitted, giving absurd coefficients and silently wrong
    // estimates of the other parameters. Checked on the equation-eligible
    // rows with exact comparisons, so no admissible specification is affected.
    local _cw_vars "`indepvars' `exog_extra' `endog' `predet'"
    local _cw_labs "`indepvars_lab' `exog_extra_lab' `endog_lab' `predet_lab'"
    // v0.9.34: with no regressor besides the automatic L.depvar the lists
    // are empty and the concatenation is blanks; without -list clean- the
    // test below passed and st_data() failed (r(3598)) on every such model.
    local _cw_vars : list clean _cw_vars
    local _cw_labs : list clean _cw_labs
    if `"`_cw_vars'"' != "" {
        local _cw_flags ""
        mata: st_local("_cw_flags", invtokens(strofreal( ///
            xdpt2_const_within("`_cw_vars'", "`panelvar'", "`touse'"))))
        local _j = 0
        foreach _f of local _cw_flags {
            local ++_j
            if `_f' {
                local _lab : word `_j' of `_cw_labs'
                di as err "`_lab' does not vary over time within any unit of the estimation sample;"
                di as err "  the `method' transformation removes it, so its coefficient is not identified"
                di as err "  (remove it, or model it through the unit effects)"
                exit 498
            }
        }
        if `flag_td' {
            mata: st_local("_cw_flags", invtokens(strofreal( ///
                xdpt2_const_within("`_cw_vars'", "`timevar'", "`touse'"))))
            local _j = 0
            foreach _f of local _cw_flags {
                local ++_j
                if `_f' {
                    local _lab : word `_j' of `_cw_labs'
                    di as err "`_lab' takes the same value for all units in every period;"
                    di as err "  td removes it, so its coefficient is not identified"
                    exit 498
                }
            }
        }
    }
    // The threshold variable: if it never changes within a unit, the regime
    // indicator is time-invariant and the regime intercept (and, under kink,
    // the kink regressor) is removed by the transformation; under td, a q
    // common to all units in every period makes the indicator a time effect.
    local _cw_flags ""
    mata: st_local("_cw_flags", strofreal( ///
        xdpt2_const_within("`q_var'", "`panelvar'", "`touse'")))
    if `_cw_flags' {
        di as err "the threshold variable `q_var' does not vary over time within any unit;"
        di as err "  the regime indicator is removed by the `method' transformation, so the"
        di as err "  regime shift is not identified"
        exit 498
    }
    if `flag_td' {
        mata: st_local("_cw_flags", strofreal( ///
            xdpt2_const_within("`q_var'", "`timevar'", "`touse'")))
        if `_cw_flags' {
            di as err "the threshold variable `q_var' takes the same value for all units in"
            di as err "  every period; td removes the regime indicator, so the regime shift is"
            di as err "  not identified"
            exit 498
        }
    }

    // v0.9.30: regressors and iv() variables that take one value for every
    // unit in each period (macro variables, trends). Without td they are
    // identified through the per-period constants in Z; the Z builder drops
    // their own instrument columns, which only repeat those constants, and
    // the output note names them. _cc_vars/_cc_labs are also read by the
    // Mata rank check to name nearly (not exactly) common variables.
    local _cc_vars "`_cw_vars' `inst_extra'"
    local _cc_vars : list clean _cc_vars
    local _cc_labs "`_cw_labs' `inst_extra_lab'"
    local _cc_labs : list clean _cc_labs
    local iv_common_vars ""
    if `"`_cc_vars'"' != "" {
        local _cw_flags ""
        mata: st_local("_cw_flags", invtokens(strofreal( ///
            xdpt2_const_within("`_cc_vars'", "`timevar'", "`touse'"))))
        local _j = 0
        foreach _f of local _cw_flags {
            local ++_j
            if `_f' {
                local _lab : word `_j' of `_cc_labs'
                local iv_common_vars "`iv_common_vars' `_lab'"
            }
        }
        local iv_common_vars : list clean iv_common_vars
    }

    tempvar eqflag
    qui gen byte `eqflag' = `touse'   // equation-eligible marker for Mata

    // Immutable threshold copy. Under td, q may also be a regressor and must
    // be demeaned in that role without changing regime membership/gamma scale.
    tempvar q_threshold
    qui gen double `q_threshold' = `q_var' if `touse_hist'

    // === Display header ===
    local method_lab "`method'"
    if "`method'" == "fod" local method_lab "FOD (Arellano-Bover 1995)"
    if "`method'" == "fd"  local method_lab "FD (Arellano-Bond 1991)"

    di ""
    di as text "{hline 78}"
    di as text "Dynamic Panel Threshold Model (Seo-Shin 2016, Gong-Seo 2026)"
    di as text "{hline 78}"
    di as text "Transformation: " as res "`method_lab'" ///
       as text "   Panel: " as res "`panelvar'" ///
       as text "   Time: " as res "`timevar'"
    local restr_lab = cond(`flag_kink', "   Restriction: kink", "")
    di as text "Dep. var: " as res "`depvar'" ///
       as text "   Threshold (q): " as res "`q_var'" ///
       as text "`restr_lab'"
    local reglist "`indepvars_lab' `exog_extra_lab'"
    local reglist : list clean reglist
    if !`flag_static' local reglist "L.`depvar' (auto) `reglist'"
    if `flag_td' local reglist "`reglist' (time-demeaned)"
    di as text "Regressors: " as res "`reglist'"
    if `k_endog'  > 0 di as text "Endogenous:    " as res "`endog_lab'"
    if `k_predet' > 0 di as text "Predetermined: " as res "`predet_lab'"
    if `k_inst'   > 0 di as text "Extra IVs:     " as res "`inst_extra_lab'"
    // Run-time note on the theoretical status of non-FD methods (no ranking
    // language: the choice of moments depends on the data, not a fixed tier).
    if "`method'" != "fd" & "`nowarn'" == "" {
        di as text "Note: for " as res "method(`method')" as text " the analytic VCE, Hansen J, AR"
        di as text "diagnostics, and the threshold CI rely on Monte Carlo evidence rather than"
        di as text "on the cited first-difference theorems of Seo-Shin (2016) / Gong-Seo (2026)."
    }
    di ""

    // === Build regressor list for Mata: order matters ===
    // Column layout in X_mat:
    //   [L.y (if dynamic), exog_regressors, endog_regressors, predet_regressors]
    // Lagged y is column 1 when dynamic (handled via var_type).
    local all_exog "`indepvars' `exog_extra'"
    local all_exog : list clean all_exog
    local all_exog_lab "`indepvars_lab' `exog_extra_lab'"
    local all_exog_lab : list clean all_exog_lab

    // === Compute trim range for γ grid from q_var distribution ===
    // Match xthenreg convention: trim(0.2) = trim 0.1 each tail → p10 to p90
    tempname q_lo q_hi
    local trim_lo = (`trim' / 2) * 100
    local trim_hi = 100 - `trim_lo'
    capture quietly _pctile `q_var' if `touse_hist', percentiles(`trim_lo' `trim_hi')
    if _rc {
        di as err "qx() has no usable values in the marked history sample"
        exit 498
    }
    scalar `q_lo' = r(r1)
    scalar `q_hi' = r(r2)
    if missing(`q_lo') | missing(`q_hi') {
        di as err "qx() has no usable values in the marked sample"
        exit 498
    }
    // These are initialization bounds only; Mata recomputes percentiles on
    // the effective GMM stack. If raw-sample quantiles tie, use the raw range
    // so the effective-sample check—not an irrelevant boundary row—decides.
    if `q_lo' >= `q_hi' {
        qui summarize `q_var' if `touse', meanonly
        scalar `q_lo' = r(min)
        scalar `q_hi' = r(max)
        if missing(`q_lo') | missing(`q_hi') | `q_lo' >= `q_hi' {
            di as err "qx() has insufficient variation; threshold grid is empty"
            exit 498
        }
    }

    // === Dispatch to Mata ===
    // v0.7.12: -noboot- (renamed from the misleading -nosearch-, which never
    // skipped the gamma grid search) turns off ALL bootstrap inference -- the
    // grid CI and the linearity/continuity tests -- leaving the point estimate.
    local do_grid_ci = cond("`noboot'" == "", 1, 0)
    // v0.9.2: unit resampling is far slower than wild -- say so upfront.
    if `do_grid_ci' & `flag_boot_exact' & "`nowarn'" == "" {
        di as txt "boottype(unit): VERIFICATION ONLY; not certified as exact Gong-Seo Algorithm 1."
        di as txt "It uses unit-multiplicity resampling with fixed sample W1 and per-draw"
        di as txt "recentered Omega/W2*; expect >= 50x the wild-bootstrap runtime."
    }

    // Unique per-fit token generated by Stata (independent of Mata state and
    // the statistical RNG). It survives in e() and prevents a restarted Mata
    // serial counter from ever aliasing a different cached fit.
    tempfile _p_cache_token

    preserve
    qui keep if `touse_hist'
    sort `panelvar' `timevar'

    // v0.8.0: legacy tdpurge pre-demeaning block removed (see td parse note).

    tempname b V gam obj nused gam_lo gam_hi pval_lin pval_cont
    tempname n_raw n_trans n_iv n_units
    tempname hansen hansen_df hansen_p ar1 ar1_p ar2 ar2_p
    tempname ci_empty ci_nseg
    local eqvar `eqflag'
    mata: xtdpthresh_run("`depvar'", "`Ly'", "`all_exog'", "`endog'",    ///
                          "`predet'", "`inst_extra'", "`q_threshold'",     ///
                          "`panelvar'", "`timevar'",                       ///
                          "`method'", `flag_static', `flag_kink',           ///
                          `flag_collapse', `maxlag_lo', `maxlag_hi',         ///
                          `grid', `gridci', `trim', `=`q_lo'', `=`q_hi'',   ///
                          `do_grid_ci', `boot', `=(100-`level')/100',       ///
                          `flag_iv_collapse', `flag_exportgmm', `flag_notest', ///
                          `flag_cont_test')

    // === Retrieve results from r() ===
    // v0.7.0 (A2 fix): retrieval moved BEFORE -restore-. Stata does not
    // guarantee that r() survives -restore-, so reading r() afterwards was
    // version-fragile.
    matrix `b'       = r(xdpt2_theta)
    matrix `V'       = r(xdpt2_V)
    // v0.9.32: under kink, e(V) is the slope block of the joint variance of
    // the slopes and gamma-hat (1), or the conditional variance when that
    // could not be computed (0); missing for the jump model.
    local kink_joint = r(xdpt2_kink_joint)
    tempname V_cond
    if "`kink_joint'" == "1" matrix `V_cond' = r(xdpt2_V_cond)
    scalar `gam'     = r(xdpt2_gamma)
    scalar `obj'     = r(xdpt2_obj)
    scalar `nused'   = r(xdpt2_nused)
    scalar `gam_lo'  = r(xdpt2_gam_lo)
    scalar `gam_hi'  = r(xdpt2_gam_hi)
    scalar `pval_lin' = r(xdpt2_pval_lin)
    scalar `pval_cont' = r(xdpt2_pval_cont)
    local lin_valid = r(xdpt2_lin_valid)
    local cont_valid = r(xdpt2_cont_valid)
    local cont_common = r(xdpt2_cont_common)
    scalar `ci_empty' = r(xdpt2_ci_empty)
    scalar `ci_nseg'  = r(xdpt2_ci_nseg)
    scalar `n_raw'    = r(xdpt2_n_raw)
    scalar `n_trans'  = r(xdpt2_n_trans)
    scalar `n_iv'     = r(xdpt2_n_iv)
    scalar `n_units'  = r(xdpt2_n_units)
    local n_switch = r(xdpt2_n_switch)
    // v0.9.35: kernel bandwidth of the jump model's joint variance
    local gamma_bw = r(xdpt2_gamma_bw)
    // v0.9.34 (C3): 1 if the kink refinement moved gamma-hat off the grid
    local kink_refined = r(xdpt2_kink_refined)
    // v0.9.29: 0 = Arellano-Bond MA(1) first-step weight as documented (or
    // FOD's 2SLS weight); 1 = (Z'Z)^-1 substituted; 2 = identity substituted.
    local w1_fallback = r(xdpt2_w1_fallback)
    if missing(`w1_fallback') local w1_fallback = 0
    // v0.9.30: instrument columns dropped as multiples of the per-period
    // constants (see xdpt2_drop_cellconst); 0 on every design without a
    // variable common to all units in a period.
    local iv_common = r(xdpt2_iv_common)
    if missing(`iv_common') local iv_common = 0
    local iv_dep = r(xdpt2_iv_dep)
    if missing(`iv_dep') local iv_dep = 0
    local iv_dep_res = r(xdpt2_iv_dep_res)
    local iv_dep_near = r(xdpt2_iv_dep_near)
    if missing(`iv_dep_near') local iv_dep_near = 0
    local balanced_eff = r(xdpt2_balanced_eff)
    if missing(`balanced_eff') local balanced_eff = 0
    scalar `hansen'    = r(xdpt2_hansen)
    scalar `hansen_df' = r(xdpt2_hansen_df)
    scalar `hansen_p'  = r(xdpt2_hansen_p)
    scalar `ar1'       = r(xdpt2_ar1)
    scalar `ar1_p'     = r(xdpt2_ar1_p)
    scalar `ar2'       = r(xdpt2_ar2)
    scalar `ar2_p'     = r(xdpt2_ar2_p)
    // v0.9.35: AR with gamma-hat treated as known; 1 if the AR statistics
    // include the estimation of gamma-hat
    tempname ar1_cond ar2_cond
    scalar `ar1_cond'    = r(xdpt2_ar1_cond)
    scalar `ar2_cond'    = r(xdpt2_ar2_cond)
    local ar_joint = r(xdpt2_ar_joint)
    if missing(`ar_joint') local ar_joint = 0
    local q_nvals_bw = r(xdpt2_q_nvals_bw)
    tempname ar1_b0 ar1_T1 ar1_TT ar2_b0 ar2_T1 ar2_TT
    scalar `ar1_b0' = r(xdpt2_ar1_b0)
    scalar `ar1_T1' = r(xdpt2_ar1_T1)
    scalar `ar1_TT' = r(xdpt2_ar1_TT)
    scalar `ar2_b0' = r(xdpt2_ar2_b0)
    scalar `ar2_T1' = r(xdpt2_ar2_T1)
    scalar `ar2_TT' = r(xdpt2_ar2_TT)
    tempname p_serial p_sig ar1_np ar2_np
    scalar `p_serial' = r(xdpt2_p_serial)
    scalar `p_sig'    = r(xdpt2_p_sig)
    scalar `ar1_np'   = r(xdpt2_ar1_np)
    scalar `ar2_np'   = r(xdpt2_ar2_np)
    local ar1_nclust = r(xdpt2_ar1_nclust)
    local ar2_nclust = r(xdpt2_ar2_nclust)
    if missing(`ar1_nclust') local ar1_nclust = 0
    if missing(`ar2_nclust') local ar2_nclust = 0
    scalar `q_lo'     = r(xdpt2_q_lo)
    scalar `q_hi'     = r(xdpt2_q_hi)
    local seed_threshold   = r(xdpt2_seed_threshold)
    local seed_linearity   = r(xdpt2_seed_linearity)
    local seed_continuity  = r(xdpt2_seed_continuity)
    local seed_coefficient = r(xdpt2_seed_coefficient)
    local wind_applied = r(xdpt2_wind_applied)
    if missing(`wind_applied') local wind_applied = 0
    local same_threshold = r(xdpt2_same_threshold)
    // v0.8.0 (#2): coefficient percentile bootstrap (threshold-search aware)
    tempname bci_mat
    local bci_B = r(xdpt2_bci_B)
    if missing(`bci_B') local bci_B = 0
    local bci_2s = r(xdpt2_bci_2s)
    if missing(`bci_2s') local bci_2s = 0
    local bci_fb = r(xdpt2_bci_fb)
    if missing(`bci_fb') local bci_fb = 0
    local bci_skip = r(xdpt2_bci_skip)
    if missing(`bci_skip') local bci_skip = 0
    // v0.8.3 R12 (#2): replay search-space sizes (threaded out of
    // coefboot as output args; in-helper r() writes are wiped by the
    // st_rclear() that precedes xdpt2_run's export block)
    local bci_g1 = r(xdpt2_bci_g1)
    if missing(`bci_g1') local bci_g1 = 0
    local bci_g2 = r(xdpt2_bci_g2)
    if missing(`bci_g2') local bci_g2 = 0
    // v0.8.3 R12 (#1/#5): CI validity and attempt accounting. bci_B > 0
    // does NOT imply r(xdpt2_bci) exists (1-9 successful draws return no
    // matrix), so every consumer of the matrix gates on bci_valid.
    local bci_valid = r(xdpt2_bci_valid)
    if missing(`bci_valid') local bci_valid = 0
    local bci_att = r(xdpt2_bci_att)
    if missing(`bci_att') local bci_att = 0
    // No estimator mixture exists: failed fixed-B two-step draws are
    // discarded, not replaced by one-step estimates or replacement draws.
    // The helper requires >=90% and >=10 valid draws; bci_fb/rate expose
    // every failed draw.
    local bci_rate = cond(`bci_att' > 0, `bci_fb'/`bci_att', 0)
    local est_2s = r(xdpt2_twostep)
    if missing(`est_2s') local est_2s = 0
    local grid_req = r(xdpt2_grid_req)
    local grid_eff = r(xdpt2_grid_eff)
    local grid_adm = r(xdpt2_grid_adm)
    tempname grid_lo grid_hi grid2_lo grid2_hi
    scalar `grid_lo' = r(xdpt2_grid_lo)
    scalar `grid_hi' = r(xdpt2_grid_hi)
    scalar `grid2_lo' = r(xdpt2_grid2_lo)
    scalar `grid2_hi' = r(xdpt2_grid2_hi)
    // v0.8.2 R11 (#2/#3): finer admission counts + CI-grid span
    local grid_struct = r(xdpt2_grid_struct)
    local search_l1_n = r(xdpt2_search_l1_n)
    local search_l2_n = r(xdpt2_search_l2_n)
    local search_l3_n = r(xdpt2_search_l3_n)
    local search_s1_level = r(xdpt2_search_s1_level)
    local search_s2_level = r(xdpt2_search_s2_level)
    local search_s1_n = r(xdpt2_search_s1_n)
    local search_s2_n = r(xdpt2_search_s2_n)
    local search_s1_same = r(xdpt2_search_s1_same)
    local search_s2_same = r(xdpt2_search_s2_same)
    // v0.9.31: real-valued diagnostic kept in a scalar (a local rounds it)
    tempname search_s1_gain
    scalar `search_s1_gain' = r(xdpt2_search_s1_gain)
    // v0.9.31: real-valued diagnostic kept in a scalar (a local rounds it)
    tempname search_s2_gain
    scalar `search_s2_gain' = r(xdpt2_search_s2_gain)
    local search_s1_conv = r(xdpt2_search_s1_conv)
    local search_s2_conv = r(xdpt2_search_s2_conv)
    local search_hit_max = r(xdpt2_search_hit_max)
    local search_W2_builds = r(xdpt2_search_W2_builds)
    // v0.9.31: real-valued diagnostic kept in a scalar (a local rounds it)
    tempname search_g1
    scalar `search_g1' = r(xdpt2_search_g1)
    // v0.9.31: real-valued diagnostic kept in a scalar (a local rounds it)
    tempname search_o1
    scalar `search_o1' = r(xdpt2_search_o1)
    // v0.9.31: real-valued diagnostic kept in a scalar (a local rounds it)
    tempname search_g2_global
    scalar `search_g2_global' = r(xdpt2_search_g2_global)
    // v0.9.31: real-valued diagnostic kept in a scalar (a local rounds it)
    tempname search_o2_global
    scalar `search_o2_global' = r(xdpt2_search_o2_global)
    local search_converged = .
    local search_incomplete = .
    if "`searchmode'" == "adaptive" {
        local search_converged = (`search_s1_conv' == 1 & ///
            `search_s2_conv' == 1 & `est_2s' == 1)
        local search_incomplete = (`search_converged' == 0)
        if `search_incomplete' & "`nowarn'" == "" {
            di as text "Note: adaptive threshold search did not certify stability within its requested cap (" ///
                as res `searchmax_effective' as text " points);"
            di as text "see e(search_incomplete), e(search_stage1_converged), and e(search_stage2_converged)."
        }
    }
    local ref_it = r(xdpt2_ref_it)
    if missing(`ref_it') local ref_it = 0
    local ref_add = r(xdpt2_ref_add)
    if missing(`ref_add') local ref_add = 0
    local ref_pool = r(xdpt2_ref_pool)
    local ref_rem = r(xdpt2_ref_rem)
    local ref_exh = r(xdpt2_ref_exh)
    tempname ref_lo ref_hi
    scalar `ref_lo' = r(xdpt2_ref_lo)
    scalar `ref_hi' = r(xdpt2_ref_hi)
    local ref_inb = r(xdpt2_ref_inb)
    local ref_nrem = r(xdpt2_ref_nrem)
    local ref_comp = r(xdpt2_ref_comp)
    // v0.9.31: real-valued diagnostic kept in a scalar (a local rounds it)
    tempname ref_obj_gain
    scalar `ref_obj_gain' = r(xdpt2_ref_obj_gain)
    if missing(`grid_struct') local grid_struct = 0
    local grid_adm2 = r(xdpt2_grid_adm2)
    local gci_eff = r(xdpt2_gci_eff)
    local gci_adm = r(xdpt2_gci_adm)
    local gci_eval = r(xdpt2_gci_eval)
    tempname gci_lo gci_hi
    scalar `gci_lo' = r(xdpt2_gci_lo)
    scalar `gci_hi' = r(xdpt2_gci_hi)
    local ci_minB = r(xdpt2_ci_minB)
    tempname ci_tab_m ci_seg_m
    cap matrix `ci_tab_m' = r(xdpt2_ci_grid)
    cap matrix `ci_seg_m' = r(xdpt2_ci_segments)
    local ci_unres = r(xdpt2_ci_unres)
    tempname citest_m
    local _has_citest 0
    cap matrix `citest_m' = r(xdpt2_citest)
    if !_rc {
        if colsof(`citest_m') == 7 local _has_citest 1
    }
    local seed_citest = r(xdpt2_seed_citest)
    // v0.8.2 R11 (#7): grid/floor reproducibility metadata
    local minreg_def = r(xdpt2_minreg_def)
    local minreg_app = r(xdpt2_minreg_app)
    if `bci_valid' {
        matrix `bci_mat' = r(xdpt2_bci)
    }

    restore

    // v0.9.6 R22 (#5): sign the SOURCE COLUMNS of every cached predict
    // series (keys, depvar, threshold, and the base variables behind all
    // regressors/instruments -- ts-operator terms reduce to their base).
    // predict verifies this signature: cached residuals/fits are
    // estimation-time values and are invalid once the data change or a
    // different dataset with coincident keys is loaded.
    local _dsvars "`panelvar' `timevar' `depvar' `q_var'"
    local _dsterms "`indepvars_lab' `exog_extra_lab' `endog_lab' `predet_lab' `inst_extra_lab'"
    if `"`_dsterms'"' != "" {
        // Let Stata parse nested operators (for example L.D.x) rather than
        // stripping only the first prefix and accidentally signing D.x.
        capture quietly tsrevar `_dsterms', list
        if _rc {
            di as err "could not resolve source variables for the predict data signature"
            exit 498
        }
        local _dsvars "`_dsvars' `r(varlist)'"
    }
    local _dsvars : list uniq _dsvars
    // v0.9.29: _datasignature checksums each column on its own, so values
    // moved between rows (or relabelled panel ids) left it unchanged and
    // predict served stale residuals. xdpt2_rowsig2 (v0.9.33) ties every value
    // to its (panel, time) key; e(p_dsig_type) tells predict which method to use.
    local _dsig ""
    capture mata: st_local("_dsig", xdpt2_rowsig2("`_dsvars'"))
    if _rc | `"`_dsig'"' == "" {
        di as err "could not construct the predict data-integrity signature"
        di as err "cached fitted values would be unsafe; estimation results were not posted"
        exit 498
    }

    // e(sample) marks the dated row (period t) of each transformed equation,
    // so e(N) = e(N_stack) = the number of equations; rows that enter only
    // as lags or as later FOD rows are not marked.
    tempvar _es_value _esample_actual
    qui gen double `_es_value' = . if `touse'
    mata: xdpt2_p_fill("`panelvar'", "`timevar'", "`_es_value'", ///
                        "`touse'", 2, 1, st_numscalar("`p_serial'"), ///
                        st_numscalar("`p_sig'"), st_local("_p_cache_token"))
    qui gen byte `_esample_actual' = !missing(`_es_value')

    // === Coefficient labels (match xthenreg convention) ===
    local cnames ""
    if !`flag_static' local cnames "Lag_y_b"
    foreach v of local all_exog_lab {
        local _cv = subinstr("`v'", ".", "_", .)
        local _cv = subinstr("`_cv'", "/", "_", .)
        local _cv = subinstr("`_cv'", "(", "", .)
        local _cv = subinstr("`_cv'", ")", "", .)
        local cnames "`cnames' `_cv'_b"
    }
    foreach v of local endog_lab {
        local _cv = subinstr("`v'", ".", "_", .)
        local _cv = subinstr("`_cv'", "/", "_", .)
        local _cv = subinstr("`_cv'", "(", "", .)
        local _cv = subinstr("`_cv'", ")", "", .)
        local cnames "`cnames' `_cv'_b"
    }
    foreach v of local predet_lab {
        local _cv = subinstr("`v'", ".", "_", .)
        local _cv = subinstr("`_cv'", "/", "_", .)
        local _cv = subinstr("`_cv'", "(", "", .)
        local _cv = subinstr("`_cv'", ")", "", .)
        local cnames "`cnames' `_cv'_b"
    }
    if !`flag_kink' {
        local cnames "`cnames' cons_d"
        if !`flag_static' local cnames "`cnames' Lag_y_d"
        foreach v of local all_exog_lab {
            local _cv = subinstr("`v'", ".", "_", .)
            local _cv = subinstr("`_cv'", "/", "_", .)
            local _cv = subinstr("`_cv'", "(", "", .)
            local _cv = subinstr("`_cv'", ")", "", .)
            local cnames "`cnames' `_cv'_d"
        }
        foreach v of local endog_lab {
            local _cv = subinstr("`v'", ".", "_", .)
            local _cv = subinstr("`_cv'", "/", "_", .)
            local _cv = subinstr("`_cv'", "(", "", .)
            local _cv = subinstr("`_cv'", ")", "", .)
            local cnames "`cnames' `_cv'_d"
        }
        foreach v of local predet_lab {
            local _cv = subinstr("`v'", ".", "_", .)
            local _cv = subinstr("`_cv'", "/", "_", .)
            local _cv = subinstr("`_cv'", "(", "", .)
            local _cv = subinstr("`_cv'", ")", "", .)
            local cnames "`cnames' `_cv'_d"
        }
    }
    else {
        local cnames "`cnames' kink_slope"
    }
    // v0.7.0 (D3 fix): de-duplicate sanitized names — e.g. L.x and a variable
    // literally named L_x both map to "L_x_b"; append "_" until unique.
    local _cn_out ""
    foreach _c of local cnames {
        if strlen("`_c'") > 32 {
            local _tail = substr("`_c'", strlen("`_c'") - 1, 2)
            local _c2 = substr("`_c'", 1, 30) + "`_tail'"
        }
        else local _c2 "`_c'"
        local _dupno = 1
        while `: list _c2 in _cn_out' {
            local _suf "_`_dupno'"
            local _keep = 32 - strlen("`_suf'")
            local _c2 = substr("`_c'", 1, `_keep') + "`_suf'"
            local ++_dupno
        }
        local _cn_out "`_cn_out' `_c2'"
    }
    local cnames : list clean _cn_out

    matrix colnames `b' = `cnames'
    matrix rownames `b' = y1
    matrix colnames `V' = `cnames'
    matrix rownames `V' = `cnames'
    if "`kink_joint'" == "1" {
        matrix colnames `V_cond' = `cnames'
        matrix rownames `V_cond' = `cnames'
    }

    // === Compact final report (xthreg2-style) ===
    di as text "{hline 78}"
    // Compute boundary-pin flag (used both for display and e(boundary_warn))
    //   0 = neither bound pins  |  1 = lower pins  |  2 = upper pins  |  3 = both pin
    // v0.8.2 R11 (#3): pinning is judged against the CI grid's OWN
    // admitted span (gci_lo/gci_hi), not the estimation grid's -- the two
    // grids can admit different ranges, and the CI is inverted on the
    // former.
    // v0.9.2 R18 (user): surface silent per-point draw loss
    if `do_grid_ci' & "`nowarn'" == "" & !missing(`ci_minB') & `ci_minB' < `boot' {
        di as text "Note: some bootstrap draws were not valid (a statistic was not finite);"
        di as text "the candidate thresholds concerned are unresolved. Smallest per-point"
        di as text "valid count = " as res `ci_minB' as text " of " ///
            as res `boot' as text " (see e(gridboot_min_draws))."
    }
    // v0.9.5 R21: the old "may be understated" warning is superseded --
    // the inversion summary is withdrawn outright when incomplete (the main
    // display branch explains it), and _bwarn stays 0 automatically since
    // gam_lo/gam_hi are missing.
    local _bwarn = 0
    if `do_grid_ci' {
        // v0.8.2 R10 (#2): pinning is judged against the ADMITTED grid
        // span, not the nominal trim bounds -- minregime/ties/rank pruning
        // can make an interior-looking endpoint the true edge of the
        // search space.
        local _range_bnd = (`=`gci_hi'') - (`=`gci_lo'')
        if !missing(`=`gci_lo'') & !missing(`=`gci_hi'') & `_range_bnd' > 0 ///
            & !missing(`=`gam_lo'') & !missing(`=`gam_hi'') {
            local _eps_bnd = 1e-4 * `_range_bnd'
            local _lo_pin = (abs((`=`gam_lo'') - (`=`gci_lo'')) < `_eps_bnd')
            local _hi_pin = (abs((`=`gam_hi'') - (`=`gci_hi'')) < `_eps_bnd')
            if `_lo_pin' & `_hi_pin' local _bwarn = 3
            else if `_lo_pin'        local _bwarn = 1
            else if `_hi_pin'        local _bwarn = 2
        }
    }

    if `do_grid_ci' & !missing(`ci_unres') & `ci_unres' > 0 {
        // v0.9.5 R21 (blocker): incomplete inversion -- no reported set.
        di as text "Threshold estimate:"
        di as text "   γ̂ = " as res %9.0g `gam' ///
           as text "   GMM obj = " as res %7.3f `obj'
        di ""
        di as text "   {err}Grid-bootstrap inversion INCOMPLETE:{txt} " as res `ci_unres' ///
            as text " gamma point(s) could not be"
        di as text "   evaluated (status 4-6 in e(ci_grid)). Unevaluated is NOT rejected, so"
        di as text "   no complete bootstrap inversion set is reported: e(gamma_lo)/e(gamma_hi)/"
        di as text "   e(ci_empty)/e(ci_nseg) are missing and e(ci_incomplete) = 1."
        di as text "   Acceptance runs over the evaluated points only are stored in"
        di as text "   e(ci_segments_evaluated); inspect e(ci_grid) for the point-level detail."
        di as text "   A larger boot() or a different gridci()/trim() may resolve the points."
    }
    else if `do_grid_ci' & `=`ci_empty'' == 1 {
        // v0.7.0 (B3 fix): an empty acceptance set is reported, not hidden
        // behind a degenerate point CI.
        di as text "Threshold estimate:"
        di as text "   γ̂ = " as res %9.0g `gam' ///
           as text "   GMM obj = " as res %7.3f `obj'
        di ""
        di as text "   {err}Warning:{txt} grid bootstrap rejected ALL candidate γ at the " ///
           as res "`level'%" as text " level."
        di as text "   No bootstrap interval summary is reported: e(gamma_lo)/e(gamma_hi) are missing"
        di as text "   and e(ci_empty) = 1. This usually signals weak identification of γ or"
        di as text "   an unsuitable grid; re-run with different trim()/gridci() to check."
    }
    else if `do_grid_ci' {
        local _ci_kind = cond(`flag_boot_exact', "verification-only unit-bootstrap", ///
            "approximate wild-bootstrap")
        di as text "Threshold estimate (" as res "`level'% `_ci_kind' interval summary" as text "):"
        di as text "   γ̂ = " as res %9.0g `gam' ///
           as text "   interval = [" as res %9.0g `gam_lo' ", " %9.0g `gam_hi' "]" ///
           as text "   GMM obj = " as res %7.3f `obj'

        // v0.7.0 (B3 fix): disconnected acceptance regions are flagged; the
        // hull is still reported for continuity with earlier versions.
        // v0.9.29: shown even under nowarn -- it describes the result (the
        // hull covers rejected values), it is not a cautionary warning.
        if !missing(`=`ci_nseg'') & `=`ci_nseg'' > 1 {
            di ""
            di as text "   {err}Note:{txt} bootstrap acceptance region is disconnected (" ///
               as res `=`ci_nseg'' as text " segments);"
            di as text "   the CI shown is its convex hull -- a SUMMARY that also covers"
            di as text "   rejected or inadmissible gamma between segments. e(ci_segments)"
            di as text "   lists the accepted segments; e(ci_grid) holds the full inversion table."
        }

        // Display warning unless nowarn set. e(boundary_warn) flag is always
        // ereturn'd below regardless of display suppression.
        if `_bwarn' > 0 & "`nowarn'" == "" {
            di ""
            di as text "   {err}Warning:{txt} CI " _c
            // v0.8.2 R11 (#3): report the CI-grid admitted span -- the
            // frame the pin was detected in -- not the nominal trim bounds.
            if `_bwarn' == 3 {
                di as text "BOTH bounds pin to the CI-grid edges [" ///
                   as res %9.0g `=`gci_lo'' as text ", " ///
                   as res %9.0g `=`gci_hi'' as text "]"
            }
            else if `_bwarn' == 1 {
                di as text "lower bound pins to the CI-grid lower edge (" ///
                   as res %9.0g `=`gci_lo'' as text ")"
            }
            else {
                di as text "upper bound pins to the CI-grid upper edge (" ///
                   as res %9.0g `=`gci_hi'' as text ")"
            }
            di as text "   Possible causes: weak identification in the affected regime,"
            di as text "   or grid edge at trim(" as res %4.2f `trim' as text ") cuts close to γ̂."
            di as text "   Re-run with a different trim to check γ̂ / CI stability."
            di as text "   Alternatives: a smaller trim() widens the grid; trim(0.20) with"
            di as text "   gridtype(quantile) matches the p10-p90 grid of the Gong-Seo"
            di as text "   application; trim(0.40) is the xthenreg default. Suppress this"
            di as text "   warning with the " as res "nowarn" as text " option."
        }
    }
    else {
        di as text "Threshold estimate:"
        di as text "   γ̂ = " as res %9.0g `gam' ///
           as text "   GMM obj = " as res %7.3f `obj'
    }
    di ""
    di as text "Specification tests:"
    // v0.7.13 (audit): under -notest- the bootstrap tests are skipped and
    // their p-values are missing; printing "p = ." implied the tests ran and
    // failed. Gate the lines on !notest as well.
    if `do_grid_ci' & !`flag_notest' {
        // v0.9.4 R20 (#6) / v0.9.5 R21: a bare "p = ." hides WHY. Early
        // returns in the test (sample statistic not computable) leave the
        // valid count missing -- that case gets its own explanation.
        if missing(`=`pval_lin'') {
            if missing(`lin_valid') {
                di as text "   Linearity test not reported: the sample statistic could not be computed."
            }
            else {
                di as text "   Linearity test not reported: only " as res `lin_valid' ///
                    as text " of " as res `boot' as text " bootstrap draws were valid."
            }
        }
        else {
            di as text "   Linearity (H0: δ=0)      p = " as res %6.4f `pval_lin'
        }
        if !`flag_kink' & `flag_cont_test' {
            if missing(`=`pval_cont'') {
                if !missing(`cont_common') & `cont_common' < 2 {
                    di as text "   Continuity test not reported: fewer than two gamma points"
                    di as text "   were jointly feasible for the nested kink/jump comparison."
                }
                else if missing(`cont_valid') {
                    di as text "   Continuity test not reported: the sample statistic could not be computed."
                }
                else {
                    di as text "   Continuity test not reported: only " as res `cont_valid' ///
                        as text " of " as res `boot' as text " bootstrap draws were valid."
                }
            }
            else {
                di as text "   Continuity (H0: kink)    p = " as res %6.4f `pval_cont'
            }
        }
    }
    // v0.9.36: pointwise threshold test of citest(#)
    if `_has_citest' {
        local _ct_st = `citest_m'[1, 6]
        if inlist(`_ct_st', 1, 2) {
            di as text "   Threshold test (H0: γ = " as res %9.0g `citest_m'[1, 1] ///
               as text ")  D = " as res %7.3f `citest_m'[1, 2] ///
               as text "  p = " as res %6.4f `citest_m'[1, 7] ///
               as text "  " cond(`citest_m'[1, 4] == 1, "not rejected", "rejected") ///
               as text " at the " as res "`=100-`level''%" as text " level"
        }
        else {
            di as text "   Threshold test (H0: γ = " as res %9.0g `citest_m'[1, 1] ///
               as text ") not evaluated (status " as res `_ct_st' ///
               as text "; 3 = not admissible, 4-6 = unresolved)"
        }
    }
    // v0.9.10 R28: gamma-hat is grid-SELECTED and can be irregular under
    // the null, so the chi-square reference for J is a conditional
    // DIAGNOSTIC, not a fully standard specification test -- say so on the
    // line itself (the last place the output still read like plain GMM).
    di as text "   Diagnostic Hansen J (conditional on γ̂) = " as res %6.3f `hansen' ///
       as text "  (df=" as res %2.0f `hansen_df' ///
       as text ")  p = " as res %6.4f `hansen_p'
    // v0.9.29: say why J is missing instead of printing a bare ".".
    if missing(`=`hansen'') {
        if `est_2s' == 0 {
            di as text "   Hansen J not available: the reported estimate is one-step GMM."
        }
        else {
            di as text "   Hansen J not available: no overidentifying restrictions (df <= 0)."
        }
    }
    if "`nowarn'" == "" {
        di as text "   J is a diagnostic; its chi-square reference assumes regular identification."
    }
    di as text "   AR(1): m = " as res %6.3f `ar1' ///
       as text "  p = " as res %6.4f `ar1_p' ///
       as text "    AR(2): m = " as res %6.3f `ar2' ///
       as text "  p = " as res %6.4f `ar2_p'
    if "`nowarn'" == "" {
        if missing(`=`ar1'') | missing(`=`ar2'') {
            di as text "   Note: a missing AR statistic means too few lag pairs or an"
            di as text "   unavailable/nonpositive full Arellano-Bond variance; no simplified"
            di as text "   T1-only p-value is substituted."
        }
        if "`ar_joint'" == "1" {
            di as text "   AR includes the estimation of γ̂; e(ar1_cond) and e(ar2_cond)"
            di as text "   treat γ̂ as known."
        }
        else di as text "   AR p-values are conditional on the selected threshold."
    }
    di ""
    // v0.8.0 (audit R5): display counts that match e(): equations = e(N)
    // (= e(N_stack): one row of e(sample) per transformed equation);
    // complete-case rows = e(N_raw).
    // v0.9.29: relabelled -- "obs used" suggested that complete-case rows
    // outside e(sample) were unused, although they enter as t-1 (FD) or
    // forward-mean (FOD) values.
    qui count if `_esample_actual'
    local n_used_raw = r(N)
    // v0.9.34: the admitted count is stage 1 on the initial grid, so show it
    // against the initial grid and list refine() points apart
    local _grid_s1 = cond(missing(`search_s1_n'), `grid_eff', `search_s1_n')
    local _grid_ref = `grid_eff' - `_grid_s1' - ("`kink_refined'" == "1")
    di as text "Sample: " ///
       as text "equations = " as res %5.0f `n_used_raw' ///
       as text "  complete-case rows = " as res %5.0f `n_raw'  ///
       as text "  units = " as res %3.0f `n_units' ///
       as text "  #IV = "   as res %4.0f `n_iv' ///
       as text "  grid = " as res "`grid_adm'" as text "/" as res "`_grid_s1'" ///
       as text " (req " as res "`grid_req'" as text ")" ///
       as text cond(`_grid_ref' > 0, " + `_grid_ref' refine", "")
    di as text "{hline 78}"

    // v0.9.30: shown even under nowarn, because the instrument set differs
    // from the declared one (as xtabond2's "dropped due to collinearity").
    if `iv_common' > 0 {
        di as text "Note: " as res `iv_common' as text " instrument column(s) take one value for every unit in their"
        di as text "      period and are linear combinations of the constant instrument"
        di as text "      columns (one per period, or one under collapse); they add no moment"
        di as text "      condition and were dropped (e(N_iv_common)); #IV counts the rest."
        if `"`iv_common_vars'"' != "" {
            di as text "      Variables common to all units in each period: " as res "`iv_common_vars'"
        }
    }
    // v0.9.34: shown even under nowarn, for the same reason
    // v0.9.35: "linear combinations" only when every dropped column is one up
    // to rounding on Z itself (e(N_iv_dep_near) = 0)
    if `iv_dep' > 0 & `iv_dep_near' == 0 {
        di as text "Note: " as res `iv_dep' as text " instrument column(s) are linear combinations of other columns"
        di as text "      (for example in a period with fewer units than instruments) and were"
        di as text "      dropped (e(N_iv_dep)), as xtabond2's generalized inverse does."
    }
    else if `iv_dep' > 0 {
        di as text "{err}Warning:{txt} " as res `iv_dep' as text " instrument column(s) dropped as numerically dependent"
        di as text "      (e(N_iv_dep)); " as res `iv_dep_near' as text " are only nearly dependent (e(N_iv_dep_near))."
        di as text "      Their directions are lost (largest relative residual " ///
            as res strtrim(string(`iv_dep_res', "%9.2e")) as text ";"
        di as text "      e(iv_dep_res)); estimates depend on how the instruments are written."
    }

    // v0.7.12: post-estimation diagnostic warnings (respect nowarn).
    if "`nowarn'" == "" {
        // v0.9.29: >= rather than >. The centered moment covariance has rank
        // at most #units - 1, so from #instruments = #units the second-step
        // weight is singular and the command falls back to one-step GMM.
        if `=`n_iv'' >= `=`n_units'' {
            di as text "{err}Warning:{txt} #instruments (" as res `=`n_iv''      ///
               as text ") is not below #units (" as res `=`n_units'' as text ")." ///
               as text " The second-step weight is unreliable; use " as res "collapse" ///
               as text " or a tighter " as res "maxlag()" as text "."
        }
        if `est_2s' == 0 {
            di as text "{err}Warning:{txt} the two-step search could not be completed (the"
            di as text "second-step weight is singular or ill-conditioned, or no candidate"
            di as text "threshold gives a usable two-step fit), so the reported estimates are"
            di as text "ONE-STEP GMM with the first-step weight (e(estimator_twostep) = 0);"
            di as text "Hansen J is not reported."
        }
        if `w1_fallback' > 0 {
            di as text "{err}Warning:{txt} the Arellano-Bond first-step weight matrix is singular or"
            di as text "ill-conditioned; " ///
                cond(`w1_fallback' == 1, "(Z'Z)^-1", "the identity matrix") ///
                " was used instead (e(W1_fallback) = " `w1_fallback' ")."
        }
        if `=`n_units'' < 30 {
            di as text "Note: few cross-sectional units (" as res `=`n_units''     ///
               as text "); GMM estimates and diagnostics may be unreliable."
        }
        // v0.9.34: the change in the intercept is identified by the units
        // whose regime changes within their equations (v0.9.35: only it; the
        // slope changes also use units that stay in the upper regime)
        if !missing(`n_switch') & `n_switch' < 10 {
            di as text "Note: only " as res `n_switch' as text                        ///
               " units change regime at the estimated threshold (e(N_switch));"
            di as text "  the change in the intercept rests on them."
        }
        if `do_grid_ci' & `boot' < 999 {
            di as text "Note: " as res "boot(`boot')" as text                       ///
               " is fine for exploration; " as res "999+" as text " reduces Monte Carlo"
            di as text "error but does not certify the bootstrap design."
        }
        // v0.9.34: with (1 - level)(B + 1) < 1 no draw can reject
        if `do_grid_ci' & (1 - `level'/100)*(`boot' + 1) < 1 {
            di as text "Note: " as res "boot(`boot')" as text " is too small for level(" ///
               as res "`level'" as text "): no candidate threshold can be rejected."
        }
        if `flag_vce_wind' & `wind_applied' == 0 {
            // v0.9.29: distinguish the two causes.
            if `est_2s' == 0 {
                di as text "Note: " as res "vce(windmeijer)" as text " requested but the estimator used the"
                di as text "one-step fallback (no two-step weight); the cluster-robust paired"
                di as text "sandwich is reported. See " as res "e(vce_applied)" as text "."
            }
            else {
                di as text "Note: " as res "vce(windmeijer)" as text " requested but the corrected variance could"
                di as text "not be computed or was not positive semidefinite; the cluster-robust"
                di as text "sandwich is reported. See " as res "e(vce_applied)" as text "."
            }
        }
        if `wind_applied' & "`same_threshold'" == "0" {
            di as text "Note: the Windmeijer correction is derived for a second-step weight built"
            di as text "at the reported threshold; here W2 was built at the stage-1 threshold"
            di as text "(e(gamma_stage1) != e(gamma)), so the corrected SEs are approximate"
            di as text "(e(wind_same_threshold) = 0)."
        }
        // v0.8.0 (audit R5): make the conditioning of the analytic SEs explicit.
        // v0.9.29: the static-panel convention does not justify them here.
        // v0.9.32: under kink the SEs include the estimation error of gamma-hat.
        // v0.9.35: so do those of the jump model (kernel derivative).
        if "`kink_joint'" == "1" {
            if "`kink'" != "" {
                di as text "Note: slope SEs include the estimation error of gamma-hat (joint variance);"
                di as text "e(V_cond) is conditional on gamma-hat."
            }
            else {
                di as text "Note: slope SEs include the estimation error of gamma-hat (joint variance;"
                di as text "bandwidth " as res strtrim(string(`gamma_bw', "%9.4g")) ///
                    as text ", e(gamma_bw)); e(V_cond) is conditional on gamma-hat."
            }
            // v0.9.35: a discrete q has no density for the kernel derivative
            if "`kink'" == "" & !missing(`q_nvals_bw') & `q_nvals_bw' < 10 {
                di as text "{err}Warning:{txt} q takes only " as res `q_nvals_bw' as text " distinct " cond(`q_nvals_bw' == 1, "value", "values") " within two bandwidths of gamma-hat"
                di as text "(e(q_nvals_bw)); the joint variance and AR assume a continuous q (Seo and Shin"
                di as text "2016, Assumption 2); e(V_cond), e(ar1_cond), and e(ar2_cond) do not."
            }
        }
        else {
            di as text "Note: analytic slope SEs are CONDITIONAL on the estimated threshold"
            di as text "(gamma-hat treated as fixed). The GMM threshold estimator converges at"
            di as text "rate sqrt(n) or slower, so they omit a first-order source of uncertainty."
            if "`kink_joint'" == "0" {
                di as text "The joint variance of the slopes and gamma-hat could not be computed"
                di as text "(e(joint_vce) = 0)."
            }
        }
        if `do_grid_ci' {
            di as text "The threshold confidence set above is approximate and not certified as"
            di as text "the exact Gong-Seo procedure; see e(ci_bootstrap_certified)."
        }
        if `wind_applied' & "`kink_joint'" != "1" {
            di as text "      vce(windmeijer): Windmeijer-type correction CONDITIONAL on the"
            di as text "      selected threshold (it does not add gamma-search variability)."
        }
        if `bci_valid' {
            di as text "      Threshold-search-aware bootstrap CIs for the slopes are in"
            di as text "      " as res "e(b_bootci)" as text " (cluster wild residual bootstrap, B=" as res `bci_B' as text " valid draws;"
            di as text "      continuity robustness is not theoretically established)."
            // Failed fixed-B draws are exposed and never replaced by a
            // different estimator or a replacement random draw.
            if `bci_fb' > 0 {
                di as text "Note: " as res `bci_fb' as text " fixed-B draw(s) failed numerically (" ///
                    as res `bci_att' as text " attempted for " as res `bci_B' ///
                    as text " valid); see e(boot_coef_failed)."
            }
            if `bci_B' < 100 {
                di as text "Note: fewer than 100 valid draws -- bootstrap quantiles are noisy;"
                di as text "consider a larger boot()."
            }
            if `bci_skip' > 0 {
                di as text "Note: the gamma* re-search skipped " as res `bci_skip' ///
                    as text " estimation-grid point(s) lacking the fast path"
                di as text "(see e(boot_grid_skipped))."
            }
        }
    }
    // A requested inference object that failed is not a cosmetic warning:
    // report it even under nowarn so a missing e(b_bootci) cannot pass silently.
    if `do_grid_ci' & "`coefboot'" != "none" & !`bci_valid' {
        local _cbmode = cond(`est_2s' == 1 & "`coefboot'" == "twostep", "two-step", "one-step")
        if `bci_B' >= 10 & `bci_B' >= ceil(.9*`boot') {
            di as text "{err}Warning:{txt} successful coefficient-bootstrap replays did not yield"
            di as text "finite interval bounds; e(b_bootci) not stored."
        }
        else if `bci_att' > 0 {
            di as text "{err}Warning:{txt} coefficient bootstrap could not reach 90% valid `_cbmode'"
            di as text "draws (valid " as res `bci_B' as text " of " as res `boot' as text " requested; " ///
                as res `bci_att' as text " attempted); e(b_bootci) not stored."
        }
        else {
            di as text "{err}Warning:{txt} coefficient bootstrap could not start a valid `_cbmode' replay;"
            di as text "no coefficient interval was delivered (e(b_bootci) not stored)."
        }
    }

    // v0.7.0 (A1 fix): esample() marks the estimation sample so that
    // post-estimation commands relying on e(sample) work correctly.
    // e(N) counts raw panel-time observations in the estimation sample;
    // e(N_stack) counts transformed-equation rows used by GMM.
    ereturn post `b' `V', obs(`n_used_raw') esample(`_esample_actual') depname("`depvar'")
    ereturn scalar N_stack   = `=`nused''
    ereturn local predict    "xtdpthresh_p"
    ereturn local cmdline    `"xtdpthresh `cmdline'"'
    ereturn local cmdversion "0.9.37"
    ereturn local searchmode "`searchmode'"
    ereturn scalar searchtol = `searchtol'
    ereturn scalar searchmax = `searchmax_effective'
    ereturn scalar searchmax_specified = `searchmax_set'
    ereturn scalar search_max_level = `search_max_level'
    ereturn local boottype "`boottype'"
    if `flag_boot_exact' ereturn local boottype_status "verification-only"
    else ereturn local boottype_status "supported"
    ereturn local history "`history'"
    // v0.8.0 (audit R5): explicit bootstrap metadata so users/scripts can
    // see WHICH bootstrap produced the CI without reading the help.
    // v0.9.2 R18 (#3): metadata SPLIT per inference object -- a single
    // bootstrap_method field misdescribed the threshold CI under
    // boottype(unit) and invited "everything is exact" readings.
    if `do_grid_ci' {
        if `flag_boot_exact' {
            ereturn local threshold_bootstrap   "verification-only unit-multiplicity resampling (Gong-Seo Alg. 1-oriented: unrestricted-residual DGP, fixed sample W1, per-draw recentered Omega/W2*; NOT certified as Algorithm 1)"
            ereturn local threshold_resampling  "panel unit (iid with replacement, multiplicity weights)"
            ereturn local threshold_recentering "explicit (sample moment at the reported theta-hat subtracted from every bootstrap moment)"
        }
        else {
            // v0.9.34 (C1): the criterion of the reported estimator, its
            // weight held at the sample value; one set of Mammen weights
            // serves every candidate (O1)
            local _tbw = cond(`est_2s' == 1, "two-step criterion, second-step weight W2", ///
                                             "one-step criterion, first-step weight W1")
            ereturn local threshold_bootstrap   "cluster wild residual (`_tbw' fixed; the same draws at every candidate; approximation of Gong-Seo Alg. 1)"
            ereturn local threshold_resampling  "panel unit (multiplicative Mammen weights)"
            ereturn local threshold_recentering "implicit (E*[eta]=0 centers wild moments at the restricted fit)"
        }
        local _ccrit = cond(`flag_boot_exact' | `est_2s' == 1, "twostep", "onestep")
        ereturn local ci_criterion "`_ccrit'"
    }
    else {
        ereturn local threshold_bootstrap   "none"
        ereturn local threshold_resampling  "none"
        ereturn local threshold_recentering "none"
    }
    if `do_grid_ci' & !`flag_notest' {
        if missing(`lin_valid') ereturn local linearity_bootstrap "requested cluster wild residual; sample statistic unavailable"
        else ereturn local linearity_bootstrap "cluster wild residual"
    }
    else ereturn local linearity_bootstrap "none"
    if `do_grid_ci' & !`flag_notest' & `flag_cont_test' {
        if missing(`cont_valid') ereturn local continuity_bootstrap "requested cluster wild residual; sample statistic unavailable"
        else ereturn local continuity_bootstrap "cluster wild residual"
    }
    else ereturn local continuity_bootstrap "none"
    if !`do_grid_ci' | "`coefboot'" == "none" {
        ereturn local coefficient_bootstrap "none"
    }
    else if !`bci_valid' {
        ereturn local coefficient_bootstrap "requested cluster wild residual replay; no coefficient interval delivered"
    }
    else {
        ereturn local coefficient_bootstrap "fixed-B cluster wild residual (threshold-search-aware replay; NOT the Gong-Seo coefficient bootstrap)"
    }

    ereturn local depvar     "`depvar'"
    ereturn local q_var      "`q_var'"
    ereturn local indepvars  "`indepvars_lab'"
    ereturn local endog      "`endog_lab'"
    ereturn local predet     "`predet_lab'"
    ereturn local exog_extra "`exog_extra_lab'"
    ereturn local inst       "`inst_extra_lab'"
    ereturn local method     "`method'"
    // v0.7.13 (C1): e(vce) records the request; e(vce_applied)=1 only when
    // the Windmeijer correction actually replaced the reported V (two-step
    // path). On the one-step fallback the correction is undefined and the
    // paired robust sandwich is reported unchanged.
    // v0.8.0 (audit R5): e(vce) reports the VCE actually delivered;
    // e(vce_requested) preserves the request (they differ only on the
    // one-step fallback where the correction is undefined). e(vcetype)
    // says whether the slope SEs include the estimation error of gamma-hat
    // (joint variance; since 0.9.35 in both models) or treat it as fixed
    // (only when the joint variance could not be computed). Neither is
    // continuity-robust; threshold inference runs through the grid CI.
    ereturn local vce_requested "`vce'"
    ereturn local vce        = cond(`wind_applied', "windmeijer", "robust")
    ereturn scalar vce_applied = `wind_applied'
    if "`kink_joint'" == "1" {
        ereturn local vcetype "Joint with estimated threshold"
        ereturn matrix V_cond = `V_cond'
    }
    else ereturn local vcetype "Conditional on estimated threshold"
    // v0.9.35: the AR statistics include gamma-hat when e(ar_joint) = 1
    if "`ar_joint'" == "1" ereturn local ar_vcetype "Joint with estimated threshold"
    else ereturn local ar_vcetype "Conditional on estimated threshold"
    // v0.9.35: e(joint_vce) for both models; e(kink_joint_vce) kept for kink
    if "`kink_joint'" == "1" | "`kink_joint'" == "0" {
        ereturn scalar joint_vce = `kink_joint'
    }
    else ereturn scalar joint_vce = .
    if "`kink'" != "" & ("`kink_joint'" == "1" | "`kink_joint'" == "0") {
        ereturn scalar kink_joint_vce = `kink_joint'
    }
    else ereturn scalar kink_joint_vce = .
    if "`kink'" == "" & "`gamma_bw'" != "" ereturn scalar gamma_bw = `gamma_bw'
    else ereturn scalar gamma_bw = .
    ereturn scalar bwscale = `bwscale'
    // v0.9.35: distinct values of q within two bandwidths of gamma-hat
    ereturn scalar q_nvals_bw = `q_nvals_bw'
    if "`kink'" != "" & "`kink_refined'" != "" ereturn scalar kink_refined = `kink_refined'
    else ereturn scalar kink_refined = .
    ereturn local hansen_reference "diagnostic chi-square; regular joint threshold rank not certified"
    ereturn scalar gamma_regular_rank_certified = 0
    if `refine' > 0 {
        ereturn local threshold_search "fixed global profile grid; one W2 constructed after stage 1; final support refinement under that W2"
    }
    else if "`kink'" != "" {
        ereturn local threshold_search "fixed discrete profile grid; one W2 constructed after the stage-1 pass; continuous minimization between the grid neighbours of the grid argmin"
    }
    else {
        ereturn local threshold_search "fixed discrete profile grid; one W2 constructed after the stage-1 pass"
    }
    ereturn local continuity_kink_search "initial grid, as for the jump alternative"
    if `do_grid_ci' & !`flag_notest' {
        ereturn local linearity_statistic "profile GMM-distance (wild bootstrap; not Seo-Shin sup-Wald)"
    }
    else ereturn local linearity_statistic "not run"
    if `flag_cont_test' {
        if `do_grid_ci' & !`flag_notest' {
            if !missing(`cont_common') & `cont_common' < 2 {
                ereturn local continuity_test "unavailable: fewer than two jointly feasible kink/jump grid points"
            }
            else if missing(`=`pval_cont'') ereturn local continuity_test "nested common-grid comparison attempted; no p-value delivered"
            else ereturn local continuity_test "nested profile GMM-distance on the common feasible grid (wild bootstrap; heuristic)"
        }
        else ereturn local continuity_test "nested comparison; not run"
    }
    else if `flag_kink' ereturn local continuity_test "not run; the fitted model is the kink model"
    else ereturn local continuity_test "not run; q is not a contemporaneous regressor, so the kink model is not nested"
    if `do_grid_ci' {
        ereturn local threshold_bootstrap_conditioning "valid fixed-B solves only; unresolved points are withdrawn under the validity rule"
    }
    else ereturn local threshold_bootstrap_conditioning "not applicable"
    if `do_grid_ci' & "`coefboot'" != "none" & `bci_att' > 0 {
        if !`bci_valid' & `bci_B' >= 10 & `bci_B' >= ceil(.9*`boot') {
            ereturn local coef_bootstrap_conditioning "successful replay draws; interval construction failed finite-value checks"
        }
        else ereturn local coef_bootstrap_conditioning = cond(`bci_fb' > 0, ///
            "quantiles condition on successful fixed-B solves; see e(boot_coef_fail_rate)", ///
            "no numerical failures observed")
    }
    else ereturn local coef_bootstrap_conditioning "not applicable"
    // Neither implemented resampling scheme is certified as the exact
    // Gong-Seo algorithm; this flag prevents numerical completeness from
    // being mistaken for theoretical certification.
    ereturn scalar ci_bootstrap_certified = 0
    // v0.8.0/0.8.1 (#2): threshold-search-aware cluster wild bootstrap CIs
    // for the slopes (rows lo/hi, columns follow e(b)). They complement the
    // conditional analytic SEs by adding gamma-search variability; their
    // continuity robustness is NOT theoretically established (the wild
    // scheme is a computational approximation of the Gong-Seo bootstrap).
    // Fixed-B draws that fail numerically are not replaced. A coefficient
    // interval is posted only when at least 90% (and at least 10) succeed;
    // attempted/success/failed expose the complete accounting.
    if `bci_valid' {
        matrix colnames `bci_mat' = `cnames'
        matrix rownames `bci_mat' = lo hi
        ereturn matrix b_bootci = `bci_mat'
    }
    ereturn scalar boot_coef_B = `bci_B'
    // v0.8.3 R12 (#5): attempted = draws actually attempted (999 attempts
    // with 0 successes is NOT "attempted 0"); success = successful draws
    // (B_eff, informative even when < 10 and no CI is produced); valid =
    // a CI matrix exists. boot_coef_B kept as the legacy success count.
    // v0.8.5 R14: under -noboot- nothing was requested of the coefficient
    // bootstrap (the boot() default would otherwise masquerade as a request)
    ereturn scalar boot_coef_requested = cond(`do_grid_ci' & "`coefboot'" != "none", `boot', 0)
    ereturn scalar boot_coef_attempted = `bci_att'
    ereturn scalar boot_coef_success   = `bci_B'
    ereturn scalar boot_coef_valid     = `bci_valid'
    // Failed fixed-B draws are discarded, never replaced or mixed in.
    ereturn scalar boot_coef_failed    = `bci_fb'
    ereturn scalar boot_coef_fail_rate = `bci_rate'
    // v0.8.1 R7 (#4.2/#4.3): composition of the bootstrap distribution and
    // grid coverage, so users can see when the CI mixes estimators or the
    // gamma* re-search ran on a strict subset of the estimation grid.
    ereturn scalar boot_coef_twostep  = `bci_2s'
    // (boot_coef_fallback removed v0.9.3: no fallback draws exist)
    ereturn scalar boot_grid_skipped  = `bci_skip'
    // v0.8.2 R11 (#4): sizes of the two replay search spaces (stage 2 may
    // exceed stage 1 when points fail under W1 but solve under W2).
    ereturn scalar boot_grid_stage1 = `bci_g1'
    ereturn scalar boot_grid_stage2 = `bci_g2'
    ereturn local coefcitype "`coefcitype'"
    // v0.8.1 R8 (#5): e(coefboot) reports what the bootstrap ACTUALLY
    // replayed -- on the one-step sample fallback a twostep request is
    // downgraded (sample-level, not a draw-level failure). The request is
    // preserved separately, and the sample estimator's own mode is stored.
    ereturn local coefboot_requested "`coefboot'"
    // v0.9.4 R20 (#1): under coefboot(none) nothing replayed anything --
    // reporting "onestep" contradicted e(coefficient_bootstrap)="none".
    if !`do_grid_ci' | "`coefboot'" == "none" | `bci_att' == 0 {
        ereturn local coefboot "none"
    }
    else {
        ereturn local coefboot = cond(`est_2s' == 1 & "`coefboot'" == "twostep", "twostep", "onestep")
    }
    ereturn scalar estimator_twostep = `est_2s'
    // v0.9.29: first-step weight actually used (0 = as documented) and, for
    // vce(windmeijer), whether W2 was built at the reported threshold.
    ereturn scalar W1_fallback = `w1_fallback'
    ereturn scalar wind_same_threshold = cond(`wind_applied', `same_threshold', .)
    // v0.9.30: instrument columns dropped as multiples of the constant
    // instrument columns, and the variables common to all units in each period.
    ereturn scalar N_iv_common = `iv_common'
    ereturn scalar N_iv_dep = `iv_dep'
    // v0.9.35: dropped columns that are not linear combinations of the kept
    // ones (relative residual on Z above 1e-10), and the largest residual
    ereturn scalar N_iv_dep_near = `iv_dep_near'
    ereturn scalar iv_dep_res = `iv_dep_res'
    ereturn local iv_common = cond(`iv_common' > 0, "`iv_common_vars'", "")
    ereturn local td_mode    = cond(`flag_td_fwl', "fwl", "")
    ereturn local panelvar   "`panelvar'"
    ereturn local timevar    "`timevar'"
    ereturn scalar gamma     = `gam'
    ereturn scalar obj       = `obj'
    ereturn scalar gamma_lo  = `gam_lo'
    ereturn scalar gamma_hi  = `gam_hi'
    ereturn scalar pval_lin  = `pval_lin'
    ereturn scalar pval_cont = `pval_cont'
    // v0.9.4 R20 (#6): per-test bootstrap accounting
    ereturn scalar boot_threshold_requested = cond(`do_grid_ci', `boot', 0)
    ereturn scalar boot_linearity_requested  = cond(`do_grid_ci' & !`flag_notest', `boot', 0)
    ereturn scalar boot_linearity_valid      = `lin_valid'
    ereturn scalar boot_continuity_requested = cond(`do_grid_ci' & !`flag_notest' & `flag_cont_test', `boot', 0)
    ereturn scalar boot_continuity_valid     = `cont_valid'
    ereturn scalar continuity_common_grid    = `cont_common'
    if `do_grid_ci' & "`rseed'" != "" {
        ereturn scalar rseed = `rseed'
    }
    else ereturn scalar rseed = .
    ereturn local rng "`c(rng_current)'"
    ereturn scalar seed_threshold   = `seed_threshold'
    ereturn scalar seed_linearity   = `seed_linearity'
    ereturn scalar seed_continuity  = `seed_continuity'
    ereturn scalar seed_coefficient = `seed_coefficient'
    ereturn scalar ci_empty  = `ci_empty'
    ereturn scalar ci_nseg   = `ci_nseg'
    // v0.7.11: effective-sample trim bounds (the gamma grid domain), for
    // scripts and tests asserting gamma-hat and the CI lie inside them.
    ereturn scalar q_lo      = `q_lo'
    ereturn scalar q_hi      = `q_hi'
    ereturn scalar N_raw     = `n_raw'
    ereturn scalar N_trans   = `n_trans'
    ereturn scalar N_iv      = `n_iv'
    ereturn scalar N_units   = `n_units'
    ereturn scalar N_switch  = `n_switch'
    ereturn scalar hansen    = `hansen'
    ereturn scalar hansen_df = `hansen_df'
    ereturn scalar hansen_p  = `hansen_p'
    ereturn scalar ar1       = `ar1'
    ereturn scalar ar1_p     = `ar1_p'
    ereturn scalar ar2       = `ar2'
    ereturn scalar ar2_p     = `ar2_p'
    ereturn scalar ar_joint  = `ar_joint'
    ereturn scalar ar1_cond  = `ar1_cond'
    ereturn scalar ar2_cond  = `ar2_cond'
    ereturn scalar ar1_b0    = `ar1_b0'
    ereturn scalar ar1_T1    = `ar1_T1'
    ereturn scalar ar1_TT    = `ar1_TT'
    ereturn scalar ar2_b0    = `ar2_b0'
    ereturn scalar ar2_T1    = `ar2_T1'
    ereturn scalar ar2_TT    = `ar2_TT'
    ereturn scalar p_serial  = `p_serial'
    ereturn scalar p_cache_sig = `p_sig'
    ereturn local p_cache_token `"`_p_cache_token'"'
    // v0.9.6 R22 (#5): predict-cache data signature
    ereturn local p_dsig `"`_dsig'"'
    ereturn local p_dsig_vars "`_dsvars'"
    ereturn local p_dsig_type "rowsig2"

    // Signed AR pair counts: positive means the full AB (1991, eq. 8)
    // variance was used; negative means that variance was unavailable and
    // the statistic/p-value were deliberately left missing (no T1 fallback).
    ereturn scalar ar1_np    = `ar1_np'
    ereturn scalar ar2_np    = `ar2_np'
    ereturn scalar ar1_N_clust = `ar1_nclust'
    ereturn scalar ar2_N_clust = `ar2_nclust'
    ereturn scalar k_exog    = `k_exog'
    ereturn scalar k_endog   = `k_endog'
    ereturn scalar k_predet  = `k_predet'
    ereturn scalar k_inst    = `k_inst'
    ereturn scalar flag_kink = `flag_kink'
    ereturn scalar flag_static = `flag_static'
    // e(balanced) describes the panel that actually enters the final GMM
    // stack. Keep the pre-estimation xtset result separately for diagnostics.
    ereturn scalar panel_balanced = `is_balanced'
    ereturn scalar balanced       = `balanced_eff'
    ereturn local balanced_definition "same transformed-equation time support across contributing units"
    ereturn scalar flag_td   = `flag_td'
    ereturn scalar boundary_warn = `_bwarn'
    // v0.7.13 (audit): store the confidence level like standard estimation
    // commands, so a later -ereturn display- (e.g. after estimates restore)
    // reproduces the same CI level instead of reverting to c(level).
    ereturn scalar level     = `level'
    // v0.8.2 R10 (#2/#4): the admitted search space (what the boundary
    // warning and any endpoint interpretation must reference), plus the
    // grid bookkeeping and the regime floor actually applied.
    // v0.8.3 R12 (#4): stage-labeled spans. grid1 = the one-step search
    // space (ok & fast_ok); grid2 = points solvable under W_n_2 (the
    // two-step search space; missing on one-step-only paths). A reported
    // two-step gamma-hat is selected over grid2, which can extend beyond
    // grid1 -- so neither span alone bounds both estimators.
    ereturn scalar gamma_grid1_lo  = `grid_lo'
    ereturn scalar gamma_grid1_hi  = `grid_hi'
    ereturn scalar gamma_grid2_lo  = `grid2_lo'
    ereturn scalar gamma_grid2_hi  = `grid2_hi'
    ereturn scalar grid_requested  = `grid_req'
    ereturn scalar grid_effective  = `grid_eff'
    ereturn scalar grid_admitted   = `grid_adm'
    ereturn scalar grid_max_requested = `searchmax_effective'
    ereturn scalar search_level1_points = `search_l1_n'
    ereturn scalar search_level2_points = `search_l2_n'
    ereturn scalar search_level3_points = `search_l3_n'
    ereturn scalar search_stage1_level = `search_s1_level'
    ereturn scalar search_stage2_level = `search_s2_level'
    ereturn scalar search_stage1_points = `search_s1_n'
    ereturn scalar search_stage2_points = `search_s2_n'
    ereturn scalar search_stage1_same_split = `search_s1_same'
    ereturn scalar search_stage2_same_split = `search_s2_same'
    ereturn scalar search_stage1_rel_gain = `search_s1_gain'
    ereturn scalar search_stage2_rel_gain = `search_s2_gain'
    ereturn scalar search_stage1_converged = `search_s1_conv'
    ereturn scalar search_stage2_converged = `search_s2_conv'
    ereturn scalar search_converged = `search_converged'
    // Legacy adaptive-search fields are retained for one compatibility
    // release. They are non-applicable under the v0.9.26 fixed-grid contract.
    ereturn scalar search_cap_exhausted = `search_hit_max'
    ereturn scalar search_hit_max = `search_hit_max'
    ereturn local search_hit_max_definition ///
        "retired adaptive-search field; always 0 under the fixed-grid estimator"
    ereturn scalar search_incomplete = `search_incomplete'
    ereturn scalar search_W2_builds = `search_W2_builds'
    ereturn scalar gamma_stage1 = `search_g1'
    ereturn scalar obj_stage1 = `search_o1'
    ereturn scalar gamma_stage2_global = `search_g2_global'
    ereturn scalar obj_stage2_global = `search_o2_global'
    // v0.8.2 R11 (#2): structural = ok-only (sample-size/rank) admission;
    // twostep_admitted = solvable under W_n_2 (missing on one-step paths)
    ereturn scalar grid_structural = `grid_struct'
    ereturn scalar grid_twostep_admitted = `grid_adm2'
    // v0.8.2 R11 (#3): CI-grid bookkeeping (the boundary warning's frame)
    ereturn scalar gridci_requested = `gridci'
    ereturn scalar gridci_effective = `gci_eff'
    // "admitted" is retained for compatibility and means the sample-side
    // solve succeeded; statuses 5/6 may still be bootstrap-unresolved.
    ereturn scalar gridci_admitted        = `gci_adm'
    ereturn scalar gridci_sample_admitted = `gci_adm'
    ereturn scalar gridci_evaluated       = `gci_eval'
    ereturn scalar gamma_ci_grid_lo = `gci_lo'
    ereturn scalar gamma_ci_grid_hi = `gci_hi'
    // v0.9.2 R18 (user): smallest per-gamma-point count of valid bootstrap
    // replications behind the CI inversion (missing if no CI ran)
    ereturn scalar gridboot_min_draws = `ci_minB'
    // v0.9.4 R20 (#3): unresolved gamma points (status 4-6) and the
    // incompleteness flag -- unresolved is NOT rejected.
    ereturn scalar ci_unresolved = `ci_unres'
    ereturn scalar ci_incomplete = cond(missing(`ci_unres'), ., cond(`ci_unres' > 0, 1, 0))
    // v0.9.36: citest(#) -- the pointwise grid-bootstrap test at gamma = #
    if `_has_citest' {
        ereturn scalar citest_gamma  = `citest_m'[1, 1]
        ereturn scalar citest_D      = `citest_m'[1, 2]
        ereturn scalar citest_crit   = `citest_m'[1, 3]
        ereturn scalar citest_accept = `citest_m'[1, 4]
        ereturn scalar citest_draws  = `citest_m'[1, 5]
        ereturn scalar citest_status = `citest_m'[1, 6]
        ereturn scalar citest_p      = `citest_m'[1, 7]
        ereturn scalar seed_citest   = `seed_citest'
    }
    // v0.9.3 R19 (#8): the confidence SET, not just its hull
    // v0.9.3 hotfix: -matrix X = r(name)- with a nonexistent r() matrix
    // silently creates a 1x1 missing matrix (the scalar-expression reading;
    // same trap as R12 #1 on the bci matrix) -- gate on the column count,
    // not on mere existence.
    cap confirm matrix `ci_seg_m'
    if !_rc {
        if colsof(`ci_seg_m') == 2 {
            matrix colnames `ci_seg_m' = lower upper
            // v0.9.5 R21 (blocker): only a COMPLETE inversion yields a
            // complete inversion summary; otherwise the acceptance runs cover
            // the evaluated points only and are stored under an explicitly
            // non-formal name.
            if !missing(`ci_unres') & `ci_unres' > 0 {
                ereturn matrix ci_segments_evaluated = `ci_seg_m'
            }
            else {
                ereturn matrix ci_segments = `ci_seg_m'
            }
        }
    }
    cap confirm matrix `ci_tab_m'
    if !_rc {
        if colsof(`ci_tab_m') == 6 {
            matrix colnames `ci_tab_m' = gamma D_stat crit accepted B_valid status
            ereturn matrix ci_grid = `ci_tab_m'
        }
    }
    ereturn scalar minregime_requested = `minregime'
    // v0.9.10 R27: refinement bookkeeping
    ereturn scalar refine_requested  = `refine'
    ereturn scalar refine_iterations = `ref_it'
    ereturn scalar refine_added      = `ref_add'
    // Basin-pool accounting. Refined points enlarge the estimation grid and
    // the unrestricted bootstrap searches; the separate CI candidate grid is
    // still formed below from gridci() plus the relevant zero point(s).
    ereturn scalar refine_pool      = `ref_pool'
    ereturn scalar refine_remaining = `ref_rem'
    ereturn scalar refine_exhausted = `ref_exh'
    // v0.9.16 R35: these are the CONVEX HULL of the pooled coverage --
    // with expand-only unions the covered basins need not be contiguous.
    ereturn scalar refine_hull_lo = `ref_lo'
    ereturn scalar refine_hull_hi = `ref_hi'
    // v0.9.15 R34 (#1): completeness = pool consumed AND the final
    // gamma-hat's neighbourhood fully evaluated -- exhausted alone says
    // nothing about the second when gamma-hat migrates between basins.
    ereturn scalar refine_final_in_initial_basin = `ref_inb'
    ereturn scalar refine_neigh_unevaluated = `ref_nrem'
    ereturn scalar refine_complete = `ref_comp'
    ereturn scalar refine_obj_rel_gain = `ref_obj_gain'
    ereturn local refine_stage = cond(`refine' > 0, "stage2_fixed_W2", "none")
    // v0.8.2 R11 (#7): grid-config reproducibility without e(cmdline) parsing
    ereturn scalar minregime_default = `minreg_def'
    ereturn scalar minregime_applied = `minreg_app'
    ereturn scalar trim = `trim'
    ereturn local gridtype "`gridtype'"
    ereturn local gridsample "`gridsample'"

    // Stamp e(cmd) last so an error while posting metadata cannot leave a
    // partial result set that falsely identifies itself as a valid fit.
    ereturn local cmd "xtdpthresh"
    ereturn display, level(`level')
end


// v0.9.31: -tsrevar- ignores -set type-: a temporary takes the type of its
// source variable, so an operator on a float variable (D.x, LD.x, S.x) was
// evaluated in float precision. Each operator temporary is recomputed here in
// double precision; for double sources the values are unchanged.
program define _xdpt_tsdouble
    version 15.0
    syntax , EXPanded(string) LABels(string)
    local _n : word count `expanded'
    forvalues _i = 1/`_n' {
        local _tv : word `_i' of `expanded'
        local _ut : word `_i' of `labels'
        if strpos(`"`_ut'"', ".") & `"`_tv'"' != `"`_ut'"' {
            capture confirm variable `_tv', exact
            if !_rc {
                quietly recast double `_tv'
                quietly replace `_tv' = `_ut'
            }
        }
    }
end


// v0.9.29: decompose a (possibly) time-series-operated term, in the canonical
// spelling produced by -syntax varlist(ts)- (e.g. L2.x, LD.y, S12.x, F.y),
// into its base variable, its net lag (L minus F), and whether it contains a
// D., S., or F. operator, i.e. the current or a future value of the base.
// Unrecognized operator syntax is treated conservatively (hasdsf = 1).
program define _xdpt_tsterm, rclass
    args term
    local dot = strpos(`"`term'"', ".")
    if !`dot' {
        return local base `"`term'"'
        return scalar isop   = 0
        return scalar netlag = 0
        return scalar hasdsf = 0
        exit
    }
    local op   = lower(substr(`"`term'"', 1, `dot' - 1))
    local base = substr(`"`term'"', `dot' + 1, .)
    local netlag = 0
    local hasdsf = 0
    local rest "`op'"
    while "`rest'" != "" {
        local hit ""
        local tok ""
        local num ""
        if regexm("`rest'", "^([lfds])\(([0-9]+)\)") {
            local tok = regexs(1)
            local num = regexs(2)
            local hit = regexs(0)
        }
        else if regexm("`rest'", "^([lfds])([0-9]*)") {
            local tok = regexs(1)
            local num = regexs(2)
            local hit = regexs(0)
        }
        if "`hit'" == "" {
            local hasdsf = 1
            local rest ""
            continue
        }
        if "`num'" == "" local num 1
        if "`tok'" == "l" local netlag = `netlag' + `num'
        else if "`tok'" == "f" {
            local netlag = `netlag' - `num'
            local hasdsf = 1
        }
        else local hasdsf = 1
        local rest = substr("`rest'", strlen("`hit'") + 1, .)
    }
    return local base `"`base'"'
    return scalar isop   = 1
    return scalar netlag = `netlag'
    return scalar hasdsf = `hasdsf'
end


// ==================================================================
// MATA BACKEND — point estimation (Stage B.1)
// ==================================================================
// v0.9.3 R19 (#9): explicit top-level -version- for the Mata source (the
// program-level version statement does not govern this block). matastrict
// stays OFF deliberately: the backend uses implicit locals throughout and
// enabling strict would require a full declaration audit (future work).
version 15.0
mata:
mata set matastrict off

// Package-level Mata scalars (xdpt_collapse, xdpt_lag_lo, xdpt_lag_hi,
// xdpt_verbose) are declared "external" inside
// each function that uses them (Mata does not allow file-scope declarations
// at the top of a mata: block). xtdpthresh_run() assigns the values once
// per invocation; helpers read them via "external real scalar ..." locals.
// v0.9.30: never declare an external in a function that runs once per unit or
// per bootstrap draw -- Mata binds externals at every call, at a cost that
// grows with the number of live objects; read them once in the caller and
// pass them as arguments (see xdpt2_transform_unit, xdpt2_unit_cfg).

// Built-in-safe replacement for rangen(): n equally spaced points from a to b.
// This avoids relying on version-specific Mata helpers.
// v0.7.13 (audit R4, C3): n grid points on empirical quantiles of q between
// probabilities p_lo and p_hi (inclusive); duplicate quantiles from ties are
// collapsed. Endpoints equal the trim-bound quantiles by construction.
real colvector xdpt2_quantile_grid(real colvector qv, real scalar p_lo,
                                    real scalar p_hi, real scalar n)
{
    real colvector g
    real scalar i
    if (n <= 1) return(J(1, 1, xdpt2_quantile(qv, p_lo)))
    g = J(n, 1, .)
    for (i = 1; i <= n; i++) {
        g[i] = xdpt2_quantile(qv, p_lo + (p_hi - p_lo) * (i - 1) / (n - 1))
    }
    return(uniqrows(g))
}

real colvector xdpt2_rangen(real scalar a, real scalar b, real scalar n)
{
    real scalar i
    real colvector out
    if (n <= 1) return(J(1, 1, a))
    out = J(n, 1, .)
    for (i = 1; i <= n; i++) {
        out[i] = a + (b - a) * (i - 1) / (n - 1)
    }
    // v0.9.34: the last point is b exactly (a + (b - a) could round below b)
    out[n] = b
    return(out)
}

// Scale-equivariant tolerance for comparing GMM objective values.  A fixed
// unit floor (for example max(1, |a|, |b|)) makes search/tie decisions depend
// on the units of the outcome whenever both objectives are below one.
real scalar xdpt2_objtol(real scalar a, real scalar b, real scalar rel)
{
    real scalar s
    s = max((abs(a), abs(b)))
    if (s == 0) return(0)
    return(rel * s)
}

// Scale-equivariant inverse/admission check for symmetric normal and moment
// matrices.  Admission is always decided after Jacobi equilibration, so rank
// and conditioning do not depend on column units.  A truly rank-deficient
// matrix is rejected by the 1e12 relative gate.  v0.9.35 (external review,
// R1): the inverse is always that of the equilibrated matrix, mapped back.
// The raw path invsym(A), kept before when the raw matrix also passed, drops
// pivots below an absolute tolerance: a matrix of tiny scale (1e-20 times a
// matrix with condition number 3) came back as a zero "inverse" that was
// accepted, and the estimates depended on the units of y. The result is now
// checked -- no dropped pivot, and C * inv(C) within 10 k eps cond(C) of the
// identity -- and ok = 0 otherwise.
void xdpt2_syminv(real matrix A, real scalar ok, real matrix Ainv)
{
    real scalar cA, tol
    real colvector d
    real matrix C, L, Ci

    ok = 0
    Ainv = J(rows(A), cols(A), .)
    if (rows(A) == 0 | rows(A) != cols(A) | hasmissing(A)) return

    d = sqrt(diagonal(A))
    if (hasmissing(d) | any(d :<= 0)) return
    C = A :/ (d * d')
    C = (C + C') / 2
    // cond() is based on singular values and therefore does not distinguish
    // a positive-definite matrix from an invertible indefinite one. Every
    // live caller supplies a covariance, weight, or GMM normal matrix, all of
    // which must be positive definite. Reject indefinite matrices before
    // invsym() can turn them into negative weights/objectives.
    L = cholesky(C)
    if (hasmissing(L)) return
    cA = cond(C)
    if (cA >= . | cA > 1e12) return

    Ci = invsym(C)
    if (hasmissing(Ci) | diag0cnt(Ci) > 0) return
    tol = max((1e-8, 10 * rows(C) * epsilon(1) * cA))
    if (max(abs(C * Ci - I(rows(C)))) > tol) return
    Ainv = Ci :/ (d * d')
    if (hasmissing(Ainv)) return
    ok = 1
}

// v0.9.35 (external review, R2): 1 if V is finite, symmetric, and positive
// semidefinite up to a relative tolerance (the smallest eigenvalue at least
// -1e-10 times the largest in absolute value). Checks every direction, not
// only the diagonal: a contrast can have negative variance although every
// coefficient's variance is positive.
real scalar xdpt2_psd_ok(real matrix V)
{
    real rowvector ev
    real scalar mx
    if (rows(V) == 0 | rows(V) != cols(V) | hasmissing(V)) return(0)
    if (max(abs(V - V')) > 1e-8 * max(abs(V))) return(0)
    ev = symeigenvalues((V + V') / 2)
    if (hasmissing(ev)) return(0)
    mx = max(abs(ev))
    if (mx >= . | mx <= 0) return(0)
    return(min(ev) >= -1e-10 * mx)
}

// Per-unit data structure (unbalanced-aware)
struct xdpt2_unit {
    real scalar    id
    real scalar    t0       // first observed time (v0.7.0: O(1) lookup anchor)
    real colvector tpos     // tpos[t - t0 + 1] = row index in u.t; 0 = absent
                            // (empty => xdpt2_find_t falls back to linear scan)
    real colvector t        // observed times, sorted
    real colvector y        // y at observed times
    real matrix    X        // regressors at observed times (n_i × K)
                            // col 1: L.y (if dynamic); then exog; then endog;
                            //        then predetermined
    real colvector q        // q at observed times
    real colvector eq       // v0.8.1: 1 = equation-eligible (complete) row;
                            // 0 = history-only row (instrument source)
    real rowvector var_type // per col: 1=lag_y, 2=exog, 3=endog, 4=predet
    real scalar    k_endog_start  // col index where endog starts (0 if none)
    real matrix    X_inst   // user-supplied instrument values (n_i × k_inst)
}

// Build per-unit structs from long-form data
struct xdpt2_unit rowvector xdpt2_build_units(
    real colvector y, real matrix Ly, real matrix X_exog,
    real matrix X_endog, real matrix X_predet,
    real matrix X_inst, real colvector q,
    real colvector pid, real colvector tid,
    real colvector eqf,
    real scalar flag_static)
{
    struct xdpt2_unit rowvector U
    struct xdpt2_unit scalar u
    real colvector idx, ord, var_type
    real matrix pinfo
    real scalar n_units, i, K, k_ex, k_en, k_pd, k_in, n_ok, k_u

    // The ado sorts (panel,time) before entering Mata. panelsetup() therefore
    // yields all unit runs in O(N), avoiding one full pid scan per unit.
    pinfo = panelsetup(pid, 1)
    n_units = rows(pinfo)
    // v0.9.34 (SPEEDUP): the units that enter are counted first and the
    // vector is allocated once. Appending (U = U, u) copied the whole vector
    // at every unit, O(N^2): 16 of 22 seconds of a fit with N = 3000.
    n_ok = 0
    for (i = 1; i <= n_units; i++) {
        idx = (pinfo[i, 1]::pinfo[i, 2])
        if (sum(eqf[idx]) >= 2) n_ok++
    }
    U = xdpt2_unit(n_ok)
    k_u = 0

    k_ex = cols(X_exog)
    k_en = cols(X_endog)
    k_pd = cols(X_predet)
    k_in = cols(X_inst)
    real scalar lag_y_present
    lag_y_present = (flag_static ? 0 : 1)
    K = lag_y_present + k_ex + k_en + k_pd

    var_type = J(1, 0, 0)
    if (!flag_static) var_type = var_type, 1
    if (k_ex > 0)     var_type = var_type, J(1, k_ex, 2)
    if (k_en > 0)     var_type = var_type, J(1, k_en, 3)
    if (k_pd > 0)     var_type = var_type, J(1, k_pd, 4)

    for (i = 1; i <= n_units; i++) {
        idx = (pinfo[i, 1]::pinfo[i, 2])
        // The hard prerequisite for FD/FOD is two
        // equation-eligible rows (no FD/FOD pair can form otherwise) --
        // never a dynamic length prefilter, which would silently discard
        // short-but-valid panels and select units by panel length.
        if (sum(eqf[idx]) < 2) continue

        u.id = pid[idx[1]]
        u.t = tid[idx]
        u.y = y[idx]
        u.q = q[idx]

        ord = order(u.t, 1)
        u.t = u.t[ord]
        u.y = u.y[ord]
        u.q = u.q[ord]
        u.eq = eqf[idx][ord]

        u.X = J(rows(idx), 0, 0)
        if (!flag_static) u.X = u.X, Ly[idx][ord]
        if (k_ex > 0)     u.X = u.X, X_exog[idx, .][ord, .]
        if (k_en > 0)     u.X = u.X, X_endog[idx, .][ord, .]
        if (k_pd > 0)     u.X = u.X, X_predet[idx, .][ord, .]

        u.X_inst = J(rows(idx), 0, 0)
        if (k_in > 0) u.X_inst = X_inst[idx, .][ord, .]

        u.var_type = var_type
        u.k_endog_start = (k_en + k_pd > 0 ? K - (k_en + k_pd) + 1 : 0)

        // O(1) time-position lookup for dense integer calendars. Its storage
        // is capped in proportion to observed rows; sparse calendars use the
        // lower-bound search in xdpt2_find_t() and never allocate O(span).
        real scalar _span, _jj
        u.t0 = u.t[1]
        _span = u.t[rows(u.t)] - u.t0 + 1
        if (_span >= rows(u.t) & _span <= 100000 &
            _span <= 4 * rows(u.t) & min(u.t :== floor(u.t)) == 1) {
            u.tpos = J(_span, 1, 0)
            for (_jj = 1; _jj <= rows(u.t); _jj++) {
                u.tpos[u.t[_jj] - u.t0 + 1] = _jj
            }
        }
        else u.tpos = J(0, 1, 0)

        k_u++
        U[k_u] = u
    }
    return(U)
}

// Helper: find index in u.t where u.t[j] == target; 0 if not found.
// O(1) via the dense tpos table when available; otherwise O(log n_i) via
// lower-bound search on the sorted time vector.
// v0.9.1 R17 (#3): lower-bound binary search on a sorted colvector.
// Returns the first index i with tv[i] >= t, or rows(tv)+1 when every
// element is smaller. Used to rank calendar times within the global
// observed-equation-time vector xdpt_teq, so that instrument blocks and
// td-fod dummy columns never require span-sized allocations.
real scalar xdpt2_tpos(real colvector tv, real scalar t)
{
    real scalar lo, hi, mid
    lo = 1
    hi = rows(tv) + 1
    while (lo < hi) {
        mid = floor((lo + hi) / 2)
        if (tv[mid] < t) lo = mid + 1
        else hi = mid
    }
    return(lo)
}

real scalar xdpt2_find_t(struct xdpt2_unit scalar u, real scalar target)
{
    real scalar j, off
    if (rows(u.tpos) > 0) {
        off = target - u.t0 + 1
        if (off != trunc(off)) return(0)
        if (off < 1 | off > rows(u.tpos)) return(0)
        return(u.tpos[off])
    }
    j = xdpt2_tpos(u.t, target)
    if (j > rows(u.t)) return(0)
    if (u.t[j] == target) return(j)
    return(0)
}

// Map stacked (unit-index,time) rows back to the immutable threshold value.
// v0.8.2 R10: live again -- gridsample(observed) uses it for the
// xthenreg-style current-row support.
real colvector xdpt2_q_at_rows(struct xdpt2_unit rowvector units,
                                real colvector times,
                                real colvector uid)
{
    real scalar r, p
    real colvector out
    out = J(rows(times), 1, .)
    for (r = 1; r <= rows(times); r++) {
        p = xdpt2_find_t(units[uid[r]], times[r])
        if (p > 0) out[r] = units[uid[r]].q[p]
    }
    return(out)
}

// v0.8.2 R9 (audit): support of the EFFECTIVE criterion, deduplicated by
// LEVEL-OBSERVATION KEY via per-unit markers. Each transformed row
// contributes the level observations whose indicators enter it: FD rows at
// (i,t) -> {(i,t),(i,t-1)}; FOD transformed rows -> {(i,t)} + future
// equation-row keys in the forward mean.
// R9 (#4): markers replace the R8 key multiset -- that allocated ~O(N*T^2)
// rows before uniqrows (each FOD row reserved the unit's whole history) and
// deep-copied the unit struct once per row; this is O(N*T) memory with
// direct field access. The resulting support SET is identical.
real colvector xdpt2_q_support(struct xdpt2_unit rowvector units,
                                real colvector times, real colvector uid,
                                string scalar method)
{
    real colvector out
    real scalar r, j, jp, u_i, m, total, nu
    pointer(real colvector) rowvector pused
    nu = length(units)
    pused = J(1, nu, NULL)
    for (r = 1; r <= rows(times); r++) {
        u_i = uid[r]
        if (pused[u_i] == NULL) {
            pused[u_i] = &(J(rows(units[u_i].t), 1, 0))
        }
        j = xdpt2_find_t(units[u_i], times[r])
        if (j == 0) continue
        (*pused[u_i])[j] = 1
        if (method == "fd") {
            jp = xdpt2_find_t(units[u_i], times[r] - 1)
            if (jp > 0) (*pused[u_i])[jp] = 1
        }
        else {
            for (jp = j + 1; jp <= rows(units[u_i].t); jp++) {
                if (units[u_i].eq[jp]) (*pused[u_i])[jp] = 1
            }
        }
    }
    total = 0
    for (u_i = 1; u_i <= nu; u_i++) {
        if (pused[u_i] != NULL) total = total + sum(*pused[u_i])
    }
    if (total == 0) return(J(0, 1, .))
    out = J(total, 1, .)
    m = 0
    for (u_i = 1; u_i <= nu; u_i++) {
        if (pused[u_i] == NULL) continue
        for (j = 1; j <= rows(*pused[u_i]); j++) {
            if ((*pused[u_i])[j]) {
                m = m + 1
                out[m] = units[u_i].q[j]
            }
        }
    }
    return(select(out, out :< .))
}

// Helper: true if any element of a row vector / matrix block is missing.
real scalar xdpt2_hasmiss(real matrix A)
{
    if (rows(A) == 0 | cols(A) == 0) return(0)
    return(sum(A :>= .) > 0)
}

// Transform one unit: FD or FOD for both y and X.
// Returns (dy, dW(γ), Z, retained_times) for this unit.
// W includes regime regressors: W = [X_trans, r, r·y_lag, r·X_exog, r·X_endog, r·X_predet]
// For non-kink model. For kink, W has fewer cols (see separate function).
void xdpt2_transform_unit(struct xdpt2_unit scalar u,
                           real scalar gamma, string scalar method,
                           real scalar flag_static, real scalar flag_kink,
                           real scalar t_min_global, real scalar t_max_global,
                           real matrix dy_out, real matrix dW_out,
                           real matrix Z_out, real colvector times_out,
                           real scalar xdpt_lag_lo, real scalar xdpt_lag_hi,
                           real scalar xdpt_collapse, real scalar xdpt_iv_collapse,
                           real colvector xdpt_teq)
{
    // v0.9.30 (SPEEDUP, bit-for-bit): the five settings below used to be
    // declared -external- here. Mata binds externals at every call, at a cost
    // that grows with the number of live Mata objects, and this function runs
    // once per unit per grid point: the binding took about 90% of a point
    // estimate and made it O(N^2) (0.48 ms per call at N = 400, 0.85 ms at
    // N = 800, against 0.05 ms of work). xdpt2_stack_at_gamma now reads them
    // once (xdpt2_unit_cfg) and passes the values under the same names, so
    // the body below is unchanged.
    real scalar n, K, j, t, Tf, c, lag_max, b, base_col, block_K, n_blocks
    real scalar block_start
    real scalar y_lag_t, x_lag_t, v, lag_needed, n_iv_cols
    real colvector r, dy_list, times_list, iv_row
    real matrix W_lvl, dW_list, w_row, Z_list
    real rowvector fut_mean_w

    n = rows(u.y)
    K = cols(u.X)

    // Build level regressors w_it(γ):
    //   Jump (non-kink): [X_it, r_it, r_it·X_it]  → 2K+1 cols
    //   Kink:            [X_it, (q_it-γ)·r_it]    → K+1 cols
    //                     where the kink term has coefficient δ_3 (slope change)
    real scalar k_W_cols
    real colvector kink_var
    r = (u.q :> gamma)
    if (flag_kink) {
        kink_var = (u.q :- gamma) :* r
        W_lvl = u.X, kink_var       // (n × (K+1))
        k_W_cols = K + 1
    }
    else {
        W_lvl = u.X, r, u.X :* r    // (n × (2K+1))
        k_W_cols = 2*K + 1
    }

    // Per-unit transformation loop
    // v0.7.8 (SPEEDUP, bit-for-bit): preallocate the per-unit lists and fill
    // by row counter instead of growing them with the \ operator (each append
    // copies the whole accumulated matrix -> O(T^2 k) per unit per gamma).
    // The values written are IDENTICAL; only the memory pattern changes.
    real scalar m_tr
    dy_list = J(n, 1, .)
    dW_list = J(n, k_W_cols, .)
    times_list = J(n, 1, .)
    m_tr = 0

    // Missing internal lags leave their IV cells at zero; they must not drop
    // an otherwise instrumented equation. External IVs or valid exogenous
    // moments can identify an early row. The structural iv_avail filter below
    // is the single source of truth for row-level instrument availability.

    if (method == "fd") {
        // v0.8.1 (R6 #1): equations form only on EQUATION rows (both t and
        // t-1 complete); history-only rows serve as instrument sources via
        // xdpt2_find_t below, exactly as in xthenreg's full-matrix build.
        real scalar jp
        for (j = 2; j <= n; j++) {
            if (!u.eq[j]) continue
            jp = xdpt2_find_t(u, u.t[j] - 1)
            if (jp == 0) continue
            if (!u.eq[jp]) continue
            // Skip candidate rows whose transformed regressor would contain missing values
            if (u.y[j] >= . | u.y[jp] >= . | u.q[j] >= . | u.q[jp] >= .) continue
            if (xdpt2_hasmiss(W_lvl[j, .]) | xdpt2_hasmiss(W_lvl[jp, .])) continue
            m_tr = m_tr + 1
            dy_list[m_tr]    = u.y[j] - u.y[jp]
            dW_list[m_tr, .] = W_lvl[j, .] - W_lvl[jp, .]
            times_list[m_tr] = u.t[j]
        }
    }
    else {  // fod
        // v0.8.1 (R6 #1): FOD equations and their forward means are defined
        // over EQUATION rows only (identical to the pre-split estimator when
        // no history-only rows exist); history rows feed instruments only.
        real colvector ei, fut
        real scalar m_e, ne
        ei = selectindex(u.eq)
        ne = rows(ei)
        for (m_e = 1; m_e <= ne - 1; m_e++) {
            j = ei[m_e]
            // Future equation rows define the forward mean. Instrument
            // availability is assessed later after every IV source is built.
            fut = ei[|m_e + 1 \ ne|]
            Tf = rows(fut)
            if (Tf < 1) continue
            // Skip candidate rows whose FOD-transformed regressor would contain missing values
            if (u.y[j] >= . | u.q[j] >= . | xdpt2_hasmiss(u.y[fut]) | xdpt2_hasmiss(u.q[fut])) continue
            if (xdpt2_hasmiss(W_lvl[j, .]) | xdpt2_hasmiss(W_lvl[fut, .])) continue
            c = sqrt(Tf / (Tf + 1))
            m_tr = m_tr + 1
            dy_list[m_tr]    = c * (u.y[j] - mean(u.y[fut]))
            fut_mean_w = mean(W_lvl[fut, .])
            dW_list[m_tr, .] = c * (W_lvl[j, .] - fut_mean_w)
            times_list[m_tr] = u.t[j]
        }
    }

    // Truncate to the filled rows (empty -> 0-row matrices, exactly as the
    // old append version produced)
    if (m_tr == 0) {
        dy_list    = J(0, 1, 0)
        dW_list    = J(0, k_W_cols, 0)
        times_list = J(0, 1, 0)
    }
    else if (m_tr < n) {
        dy_list    = dy_list[|1 \ m_tr|]
        dW_list    = dW_list[|1, 1 \ m_tr, k_W_cols|]
        times_list = times_list[|1 \ m_tr|]
    }

    // === Build block-diagonal Z matrix per xthenreg moment structure ===
    // Per time block t ∈ [t_min+2, t_max]:
    //   col 1:        constant (1)
    //   cols 2..:     y lags (y_{t-2}, y_{t-3}, ..., up to lag_max = t-t_min)
    //   next cols:    for each exog x: Δx_t (1 IV per t)
    //   next cols:    for each endog x: x_lags (x_{t-2}, x_{t-3}, ...)
    //   next cols:    for each predet x: x_lags (x_{t-1}, x_{t-2}, ...)
    //   next cols:    user-supplied external IVs (1 per var per block)
    // Block-diagonal across t: row at time t has nonzero only in block b(t).

    // Determine compact per-block lag widths. Allocate only the requested
    // interval that can exist in the retained calendar span: L.y/endogenous
    // start at lag 2, predetermined variables at lag 1. The old allocation
    // started every block at lag 1 and discarded the leading all-zero columns
    // later, which could be enormous for maxlag(lo hi) with a large lo.
    lag_max = t_max_global - t_min_global
    if (xdpt_lag_hi < lag_max) lag_max = xdpt_lag_hi
    real scalar lag_lo_y, lag_lo_p, n_lag_y, n_lag_p
    lag_lo_y = (xdpt_lag_lo > 2 ? xdpt_lag_lo : 2)
    lag_lo_p = (xdpt_lag_lo > 1 ? xdpt_lag_lo : 1)
    n_lag_y = (lag_max >= lag_lo_y ? lag_max - lag_lo_y + 1 : 0)
    n_lag_p = (lag_max >= lag_lo_p ? lag_max - lag_lo_p + 1 : 0)
    // v0.9.27: FOD lags count from t+1 (xtabond2 convention, Roodman 2009); FD unchanged.
    real scalar fod_shift
    fod_shift = (method == "fod")
    // Per block: constant + L.y lags + transformed exogenous moments +
    // separate endogenous/predetermined lag intervals. User IVs live in the
    // tail region below and retain their collapse semantics.
    real scalar k_exog, k_endog, k_predet, k_inst, iv_width
    real scalar core_cols, inst_cols, inst_base
    k_exog   = sum(u.var_type :== 2)
    k_endog  = sum(u.var_type :== 3)
    k_predet = sum(u.var_type :== 4)
    k_inst = cols(u.X_inst)
    // v0.7.0 (C1): user IVs moved OUT of the per-block core into a tail
    // region, so iv(..., collapse) can collapse ONLY the user-IV block.
    // GMM estimates are invariant to this column permutation.
    iv_width = 1                     // constant
    if (!flag_static) iv_width = iv_width + n_lag_y   // lag y
    iv_width = iv_width + k_exog                      // Δx per exog
    iv_width = iv_width + k_endog*n_lag_y + k_predet*n_lag_p

    block_start = t_min_global
    if (method == "fd") block_start = t_min_global + 1
    // v0.9.1 R17 (#3): one block per OBSERVED equation time >= block_start
    // (rank-indexed via xdpt_teq), not per calendar integer in the span --
    // a sparse delta-1 index would otherwise allocate instrument columns
    // for thousands of never-observed periods (the zero columns were
    // dropped later, but the RAM was already spent). Gap-free index:
    // rank == t - block_start + 1, so results are bit-for-bit unchanged.
    real scalar tq_off
    tq_off = xdpt2_tpos(xdpt_teq, block_start)
    n_blocks = rows(xdpt_teq) - tq_off + 1
    if (n_blocks < 1) n_blocks = 1
    block_K = iv_width
    // Collapsed: single shared block across all t; else block-diagonal by t
    if (xdpt_collapse) core_cols = block_K
    else               core_cols = n_blocks * block_K
    if (k_inst > 0)    inst_cols = (xdpt_iv_collapse ? k_inst : n_blocks * k_inst)
    else               inst_cols = 0
    n_iv_cols = core_cols + inst_cols

    Z_list = J(rows(times_list), n_iv_cols, 0)

    real scalar i, col_off, lag_idx, pos
    real colvector cons_pos, iv_avail
    cons_pos = J(rows(times_list), 1, 0)
    // v0.7.13 (audit): structural-availability mask. iv_avail[i] = 1 iff the
    // row has at least one data-driven instrument column that STRUCTURALLY
    // exists (a lag row that is present, an exogenous transform, or a user
    // IV) — independent of that instrument's numeric value. Replaces the old
    // rowsum(abs(Z)) > 1e-12 test, which treated a genuinely zero-valued or
    // tiny-scaled valid instrument as absent, wrongly dropping the row (and
    // thereby changing e(sample), the AR test, and the bootstrap).
    iv_avail = J(rows(times_list), 1, 0)
    for (i = 1; i <= rows(times_list); i++) {
        t = times_list[i]
        b = xdpt2_tpos(xdpt_teq, t)
        if (b > rows(xdpt_teq)) continue
        if (xdpt_teq[b] != t) continue
        b = b - tq_off + 1
        if (b < 1 | b > n_blocks) continue
        if (xdpt_collapse) base_col = 0
        else               base_col = (b - 1) * block_K

        // Col 1: constant — v0.7.13 (audit, B5): position recorded here but
        // WRITTEN only after the zero-IV row filter below. Writing it up-front made every
        // rowsum >= 1, so the filter was dead and rows with no data-driven
        // instrument survived on the constant alone.
        cons_pos[i] = base_col + 1
        col_off = 1

        // Lagged y (if dynamic), compactly indexed over the effective
        // maxlag() interval.
        if (!flag_static) {
            for (lag_idx = lag_lo_y; lag_idx <= lag_max; lag_idx++) {
                pos = xdpt2_find_t(u, t + fod_shift - lag_idx)
                if (pos > 0) {
                    if (u.y[pos] < .) {
                        iv_avail[i] = 1   // v0.8.1: value must exist too
                        Z_list[i, base_col + col_off +
                            (lag_idx - lag_lo_y + 1)] = u.y[pos]
                    }
                }
            }
            col_off = col_off + n_lag_y
        }

        // Strictly exogenous x instruments itself after the chosen transform:
        // Δx for FD and the forward-deviation of x for FOD.
        // v0.7.11: pos_tm1 lookup removed -- dead since the v0.7.10
        // exog-IV change (the instrument now comes from dW_list);
        // pos_t is still needed by the user-inst block below.
        real scalar vt, vi, pos_t
        pos_t = xdpt2_find_t(u, t)
        vi = 0
        for (vt = 1; vt <= cols(u.X); vt++) {
            if (u.var_type[vt] == 2) {
                vi = vi + 1
                iv_avail[i] = 1   // Δx / FOD-x instrument is structurally the
                                  // row's own transform, always present
                Z_list[i, base_col + col_off + vi] = dW_list[i, vt]
            }
        }
        col_off = col_off + k_exog

        // Endogenous and predetermined variables use different admissible
        // lower lags, so their compact column regions have separate widths.
        real scalar vi_en, vi_pr
        vi_en = 0
        vi_pr = 0
        for (vt = 1; vt <= cols(u.X); vt++) {
            if (u.var_type[vt] == 3) {
                for (lag_idx = lag_lo_y; lag_idx <= lag_max; lag_idx++) {
                    pos = xdpt2_find_t(u, t + fod_shift - lag_idx)
                    if (pos > 0) {
                        if (u.X[pos, vt] < .) {
                            iv_avail[i] = 1   // v0.8.1: value must exist too
                            Z_list[i, base_col + col_off + vi_en*n_lag_y +
                                (lag_idx - lag_lo_y + 1)] = u.X[pos, vt]
                        }
                    }
                }
                vi_en = vi_en + 1
            }
            else if (u.var_type[vt] == 4) {
                for (lag_idx = lag_lo_p; lag_idx <= lag_max; lag_idx++) {
                    pos = xdpt2_find_t(u, t + fod_shift - lag_idx)
                    if (pos > 0) {
                        if (u.X[pos, vt] < .) {
                            iv_avail[i] = 1
                            Z_list[i, base_col + col_off + k_endog*n_lag_y +
                                vi_pr*n_lag_p + (lag_idx - lag_lo_p + 1)] = u.X[pos, vt]
                        }
                    }
                }
                vi_pr = vi_pr + 1
            }
        }
        col_off = col_off + k_endog*n_lag_y + k_predet*n_lag_p

        // User-supplied instruments (inst): value at time t, one IV per inst
        // var, in the tail region (v0.7.0: collapsible independently of the
        // GMM-style core via iv(..., collapse)).
        if (k_inst > 0 & pos_t > 0) {
            real scalar ii
            inst_base = core_cols + (xdpt_iv_collapse ? 0 : (b - 1) * k_inst)
            for (ii = 1; ii <= k_inst; ii++) {
                if (u.X_inst[pos_t, ii] < .) {
                    iv_avail[i] = 1   // user IV present at time t
                    Z_list[i, inst_base + ii] = u.X_inst[pos_t, ii]
                }
            }
        }
    }

    // Drop rows without a data-driven instrument. Their only moment is the
    // period constant (valid, and kept by xtabond2); they are dropped by
    // design so that their residuals do not enter the wild-bootstrap pool
    // and the AR tests (documented under Instrument set in the help).
    // The filter depends only on Z, which is γ-invariant, so per-γ row counts
    // stay matched across the bootstrap caches.
    // v0.7.13 (audit, B5): the filter drops rows whose STRUCTURAL instrument
    // availability mask is 0 (no lag row / exog transform / user IV present),
    // the constant deferred so it does not mask the test. Using iv_avail
    // instead of rowsum(abs(Z)) makes the drop independent of instrument
    // magnitude — a valid instrument equal to 0 no longer looks absent.
    if (rows(times_list) > 0) {
        real colvector _keep
        real scalar _ki
        _keep = selectindex(iv_avail :!= 0)
        if (length(_keep) == 0) {
            dy_list    = J(0, 1, 0)
            dW_list    = J(0, k_W_cols, 0)
            Z_list     = J(0, n_iv_cols, 0)
            times_list = J(0, 1, 0)
        }
        else {
            for (_ki = 1; _ki <= length(_keep); _ki++) {
                Z_list[_keep[_ki], cons_pos[_keep[_ki]]] = 1
            }
            if (length(_keep) < rows(times_list)) {
                dy_list    = dy_list[_keep]
                dW_list    = dW_list[_keep, .]
                Z_list     = Z_list[_keep, .]
                times_list = times_list[_keep]
            }
        }
    }

    dy_out = dy_list
    dW_out = dW_list
    Z_out = Z_list
    times_out = times_list
}

// Helper: stack (dY, dW, Z, times, unit_id) across all units at given γ

// v0.7.13 (audit R4, C2): in-place cross-sectional demeaning within each
// time cell — the FWL partialling of common-across-regime time dummies out
// of the stacked system. A singleton cell demeans to exactly zero.
void xdpt2_demean_bytime(real matrix M, real colvector times)
{
    real colvector ut, idx
    real scalar ti
    if (rows(M) == 0) return
    ut = uniqrows(times)
    for (ti = 1; ti <= rows(ut); ti++) {
        idx = selectindex(times :== ut[ti])
        if (rows(idx) > 1) {
            M[idx, .] = M[idx, .] :- mean(M[idx, .])
        }
        else {
            M[idx, .] = J(1, cols(M), 0)
        }
    }
}

// v0.9.34: indices of the instrument columns kept when numerically
// dependent columns are dropped. Cholesky in column order on the Gram matrix
// scaled to a unit diagonal: column j is kept when its squared relative
// residual after the kept earlier columns is at least 1e-13 (a relative
// residual of 3.2e-7). A Z whose scaled Gram matrix passes the rank gate of
// xdpt2_syminv (condition number <= 1e12, so every pivot >= 1e-12) keeps
// every column. Exactly dependent columns -- a period with fewer units than
// its block has columns, or period constants made collinear by the FOD
// partialling of time effects -- are dropped; xtabond2 handles them with a
// generalized inverse, which gives the same estimates. v0.9.35: rounding in
// the Gram matrix is of order 1e-8 in relative residual, so this rule also
// drops columns that are only nearly dependent (within 3.2e-7) and cannot
// tell them apart. When xdpt_ivc_diag = 1 (the final stack at gamma-hat),
// the relative residual of each dropped column is measured on Z itself
// against the kept columns (Householder QR, rounding of order 1e-12 or
// below): xdpt_ivc_dep_res is the largest, and xdpt_ivc_dep_near counts
// those above 1e-10, which are not linear combinations of the kept columns.
// The decision is the same with or without the flag.
real rowvector xdpt2_indep_cols(real matrix Z)
{
    real matrix G, L, A, B, R1
    real colvector d, l, g
    real rowvector keep, mk, drp, tau, r
    real scalar k, j, m, p, nk, nd
    external real scalar xdpt_ivc_diag, xdpt_ivc_dep_res, xdpt_ivc_dep_near
    if (xdpt_ivc_diag == 1) {
        xdpt_ivc_dep_res = .
        xdpt_ivc_dep_near = 0
    }
    k = cols(Z)
    if (k < 2 | rows(Z) == 0) return(1..k)
    G = cross(Z, Z)
    d = sqrt(diagonal(G))
    if (hasmissing(d) | any(d :<= 0)) return(1..k)
    G = G :/ (d * d')
    L = J(k, k, 0)
    keep = J(1, 0, .)
    m = 0
    for (j = 1; j <= k; j++) {
        if (m == 0) {
            l = J(0, 1, .)
            p = G[j, j]
        }
        else {
            g = G[keep', j]
            l = solvelower(L[|1, 1 \ m, m|], g)
            p = G[j, j] - l' * l
        }
        if (p >= 1e-13) {
            m = m + 1
            keep = keep, j
            if (m > 1) L[|m, 1 \ m, m - 1|] = l'
            L[m, m] = sqrt(p)
        }
    }
    nk = cols(keep)
    nd = k - nk
    if (xdpt_ivc_diag == 1 & nd > 0) {
        mk = J(1, k, 1)
        mk[keep] = J(1, nk, 0)
        drp = selectindex(mk)
        if (nk >= rows(Z)) r = J(1, nd, 0)
        else {
            A = Z[., keep]
            A = A :/ sqrt(colsum(A :^ 2))
            _hqrd(A, tau, R1)
            B = Z[., drp]
            B = hqrdmultq(A, tau, B :/ sqrt(colsum(B :^ 2)), 1)
            r = sqrt(colsum(B[|nk + 1, 1 \ rows(Z), nd|] :^ 2))
        }
        xdpt_ivc_dep_res = max(r)
        xdpt_ivc_dep_near = sum(r :> 1e-10)
    }
    return(keep)
}

// v0.9.28: indices of the instrument columns to keep when exact duplicates
// are dropped (the first occurrence is kept). Two declared variables that are
// lags of one another -- x and L.x, or L2.depvar next to the automatic
// L.depvar -- have overlapping lag windows and generate identical (variable,
// date) columns. A duplicate adds no moment but makes Z'Z singular, which
// used to reject every candidate threshold. Columns are compared exactly;
// two cheap signatures restrict the element-by-element comparisons to ties.
// Any fit that succeeded before has no duplicate (its Z'Z was nonsingular),
// so for it the returned index is 1..cols(Z) and nothing changes.
real rowvector xdpt2_nodup_cols(real matrix Z)
{
    real scalar k, n, m, g0, i, j
    real colvector ord, dup
    real matrix S

    k = cols(Z)
    n = rows(Z)
    if (k < 2 | n == 0) return(1..k)
    S = (colsum(Z)', colsum(Z :* (1::n))', (1::k))
    ord = order(S, (1, 2, 3))
    dup = J(k, 1, 0)
    m = 1
    while (m <= k) {
        g0 = m
        while (m < k) {
            if (S[ord[m + 1], 1] != S[ord[g0], 1] |
                S[ord[m + 1], 2] != S[ord[g0], 2]) break
            m++
        }
        // ord[g0..m] share both signatures and are in column order
        for (j = g0 + 1; j <= m; j++) {
            for (i = g0; i < j; i++) {
                if (dup[ord[i]]) continue
                if (Z[., ord[i]] == Z[., ord[j]]) {
                    dup[ord[j]] = 1
                    break
                }
            }
        }
        m++
    }
    return(selectindex(!dup)')
}

// v0.9.30: zero the instrument columns that are exact linear combinations of
// the constant columns kept before them, and return how many were zeroed (the
// all-zero drop in xdpt2_stack_at_gamma then removes them). A regressor or
// iv() variable that takes one value for every unit in a period -- a macro
// variable, a trend -- gives such columns: in the per-period layout its
// column in block t is that value times block t's constant, and under
// collapse a column that is constant over all rows repeats the constant.
// These columns add no moment, but they made Z'Z singular, so the command
// stopped at the rank check; xtabond2 and xthenreg discard them silently
// through a generalized inverse, and the estimates here equal theirs.
// Only columns that take one value within every time cell are candidates. A
// cell with a single row cannot show that, so a column that is nonzero in
// such a cell is never a candidate and a period observed for one unit keeps
// its old behavior. A candidate is dropped when the constants kept so far
// already span it: single-cell indicators span their cell, a kept column
// that is constant over all rows spans the all-ones direction, and any other
// kept candidate is not used to span later ones. Constancy within a cell is
// tested exactly (units in a period share the same inputs, so their values
// are bitwise equal). Equality ACROSS cells -- a column equal to one value in
// every cell, i.e. a multiple of the all-ones column -- uses a relative
// tolerance of 1e-10: a trend such as (year - 2000)/10 has first differences
// that differ from year to year in the last bits. A column that is within
// 1e-10 of the span has a condition number far above the 1e12 gate of
// xdpt2_syminv, so a Z that passed the rank check before contains no column
// this rule drops, and such fits are unchanged.
real scalar xdpt2_drop_cellconst(real matrix Z, real colvector times)
{
    real scalar k, n, nc, c, j, r0, r1, n_unc, all_ind, n_drop, drop
    real colvector ord, csize, cellid, kept, rows_nz, vbc, u_nz, uval
    real matrix info, sub
    real rowvector cc, mx, mn, nzc, newf, nzcount, first_cell, first_val
    real rowvector allsame

    k = cols(Z)
    n = rows(Z)
    if (k < 2 | n < 2) return(0)
    ord   = order(times, 1)
    info  = panelsetup(times[ord], 1)
    nc    = rows(info)
    csize = info[., 2] - info[., 1] :+ 1
    cellid = J(n, 1, 0)
    cc         = J(1, k, 1)
    nzcount    = J(1, k, 0)
    first_cell = J(1, k, 0)
    first_val  = J(1, k, 0)
    allsame    = J(1, k, 1)
    for (c = 1; c <= nc; c++) {
        r0 = info[c, 1]
        r1 = info[c, 2]
        cellid[ord[|r0 \ r1|]] = J(r1 - r0 + 1, 1, c)
        sub = Z[ord[|r0 \ r1|], .]
        if (csize[c] == 1) {
            cc = cc :& (sub :== 0)
            continue
        }
        mx = colmax(sub)
        mn = colmin(sub)
        cc = cc :& (mx :== mn)
        nzc = cc :& (mx :!= 0)
        newf = nzc :& (first_cell :== 0)
        first_cell = first_cell + newf :* c
        first_val  = first_val + newf :* mx
        allsame = allsame :& (!nzc :| (abs(mx :- first_val) :<=
                              1e-10 :* abs(first_val)))
        nzcount = nzcount + nzc
    }
    if (!any(cc :& (nzcount :> 0))) return(0)

    kept = J(nc, 1, 0)
    n_unc = nc
    all_ind = 0
    n_drop = 0
    for (j = 1; j <= k; j++) {
        if (!cc[j] | nzcount[j] == 0) continue
        if (nzcount[j] == 1) {
            // c * (indicator of one cell)
            c = first_cell[j]
            drop = (kept[c] | (all_ind & n_unc == 1))
            if (!drop) {
                kept[c] = 1
                n_unc--
            }
        }
        else if (nzcount[j] == nc & allsame[j]) {
            // c * (all ones)
            drop = (n_unc == 0 | all_ind)
            if (!drop) all_ind = 1
        }
        else if (n_unc == 0) {
            drop = 1
        }
        else {
            // one value per cell: spanned when it vanishes on every cell not
            // yet spanned, or, with the all-ones direction kept, when it is
            // one common value on all of those cells
            rows_nz = selectindex(Z[., j] :!= 0)
            vbc = J(nc, 1, 0)
            vbc[cellid[rows_nz]] = Z[rows_nz, j]
            u_nz = selectindex((vbc :!= 0) :& !kept)
            if (rows(u_nz) == 0) drop = 1
            else if (all_ind & rows(u_nz) == n_unc) {
                uval = vbc[u_nz]
                drop = all(abs(uval :- uval[1]) :<= 1e-10 * max(abs(uval)))
            }
            else drop = 0
        }
        if (drop) {
            Z[., j] = J(n, 1, 0)
            n_drop++
        }
    }
    return(n_drop)
}

// v0.9.30: user-facing labels of the variables in -vars- that take nearly,
// but not exactly, one value for every unit in each period (within-period
// range at most 1e-6 of the variable's largest absolute value, and nonzero
// in at least one period). xdpt2_drop_cellconst recognizes only exactly
// common variables; a nearly common one (for example a macro series merged
// with rounding differences) can leave Z ill-conditioned. Used only to
// explain a rank failure.
string scalar xdpt2_near_common(string scalar vars, string scalar labs,
                                string scalar grp, string scalar touse)
{
    real matrix X, info
    real colvector g, ord, x
    real scalar j, i, r0, r1, scale, rng, exact, near
    string rowvector lb
    string scalar out

    out = ""
    if (vars == "") return(out)
    lb = tokens(labs)
    X = st_data(., vars, touse)
    g = st_data(., grp, touse)
    if (rows(X) < 2 | cols(lb) != cols(X)) return(out)
    ord = order(g, 1)
    X = X[ord, .]
    g = g[ord]
    info = panelsetup(g, 1)
    for (j = 1; j <= cols(X); j++) {
        x = X[., j]
        scale = max(abs(x))
        if (scale >= . | scale == 0) continue
        exact = 1
        near = 1
        for (i = 1; i <= rows(info); i++) {
            r0 = info[i, 1]
            r1 = info[i, 2]
            if (r1 <= r0) continue
            rng = max(x[|r0 \ r1|]) - min(x[|r0 \ r1|])
            if (rng >= .) continue
            if (rng != 0) exact = 0
            if (rng > 1e-6 * scale) {
                near = 0
                break
            }
        }
        if (near & !exact) out = out + (out == "" ? "" : " ") + lb[j]
    }
    return(out)
}

// v0.9.30: the settings xdpt2_transform_unit needs, read once per stack
// (see the note there on the cost of -external- in a hot function).
void xdpt2_unit_cfg(real rowvector cfg, real colvector teq)
{
    external real scalar xdpt_lag_lo, xdpt_lag_hi, xdpt_collapse, xdpt_iv_collapse
    external real colvector xdpt_teq
    cfg = (xdpt_lag_lo, xdpt_lag_hi, xdpt_collapse, xdpt_iv_collapse)
    teq = xdpt_teq
}

void xdpt2_stack_at_gamma(struct xdpt2_unit rowvector units,
                           real scalar gamma, string scalar method,
                           real scalar flag_static, real scalar flag_kink,
                           real scalar t_min, real scalar t_max,
                           real matrix dY_out, real matrix dW_out,
                           real matrix Z_out, real colvector times_out,
                           real colvector unit_id_out)
{
    real scalar i, K, n_units, n_rows_i, k_W_cols
    real colvector dy_i, time_i
    real matrix dW_i, Z_i

    n_units = length(units)
    K = cols(units[1].X)
    k_W_cols = (flag_kink ? K + 1 : 2*K + 1)

    // === FOD/FD equation rows ===
    // v0.7.8 (SPEEDUP, bit-for-bit): two-pass stacking. Pass 1 calls the
    // per-unit transform ONCE per unit and parks the results behind pointers
    // (as copies -- the locals are overwritten by the next call); pass 2
    // allocates each stacked matrix once and fills it by row ranges. The old
    // one-pass version grew the stacks with the \ operator, which copies the
    // ENTIRE accumulated matrix on every append -> O(N_units^2) copying per
    // gamma point; with thousands of units this dominated cache building.
    // Row values and row order are IDENTICAL to the append version, so the
    // stacks -- and everything downstream -- are bit-for-bit unchanged.

    pointer() rowvector pY_s, pW_s, pZ_s, pT_s
    real colvector nr_s
    real scalar n_tot_s, r0_s, n_iv_s

    pY_s = J(1, n_units, NULL)
    pW_s = J(1, n_units, NULL)
    pZ_s = J(1, n_units, NULL)
    pT_s = J(1, n_units, NULL)
    nr_s = J(n_units, 1, 0)
    n_tot_s = 0
    n_iv_s = 0   // set from the first unit with rows, like the old code
    // v0.9.30: settings for xdpt2_transform_unit, read once per stack.
    real rowvector _ucfg
    real colvector _uteq
    xdpt2_unit_cfg(_ucfg, _uteq)

    for (i = 1; i <= n_units; i++) {
        xdpt2_transform_unit(units[i], gamma, method,
                              flag_static, flag_kink,
                              t_min, t_max,
                              dy_i, dW_i, Z_i, time_i,
                              _ucfg[1], _ucfg[2], _ucfg[3], _ucfg[4], _uteq)
        n_rows_i = rows(dy_i)
        if (n_rows_i == 0) continue
        nr_s[i] = n_rows_i
        if (n_iv_s == 0) n_iv_s = cols(Z_i)
        // &(x[.,.]) parks a COPY of x behind the pointer
        pY_s[i] = &(dy_i[., .])
        pW_s[i] = &(dW_i[., .])
        pZ_s[i] = &(Z_i[., .])
        pT_s[i] = &(time_i[., .])
        n_tot_s = n_tot_s + n_rows_i
    }

    if (n_tot_s == 0) {
        dY_out = J(0, 1, 0)
        dW_out = J(0, k_W_cols, 0)
        Z_out = J(0, 0, 0)
        times_out = J(0, 1, 0)
        unit_id_out = J(0, 1, 0)
    }
    else {
        dY_out      = J(n_tot_s, 1, .)
        dW_out      = J(n_tot_s, k_W_cols, .)
        Z_out       = J(n_tot_s, n_iv_s, .)
        times_out   = J(n_tot_s, 1, .)
        unit_id_out = J(n_tot_s, 1, .)
        r0_s = 1
        for (i = 1; i <= n_units; i++) {
            if (nr_s[i] == 0) continue
            dY_out[|r0_s \ r0_s + nr_s[i] - 1|]                = *pY_s[i]
            dW_out[|r0_s, 1 \ r0_s + nr_s[i] - 1, k_W_cols|]   = *pW_s[i]
            Z_out[|r0_s, 1 \ r0_s + nr_s[i] - 1, n_iv_s|]      = *pZ_s[i]
            times_out[|r0_s \ r0_s + nr_s[i] - 1|]             = *pT_s[i]
            unit_id_out[|r0_s \ r0_s + nr_s[i] - 1|]           = J(nr_s[i], 1, i)
            r0_s = r0_s + nr_s[i]
            // free the parked copy early to cap peak memory
            pY_s[i] = NULL
            pW_s[i] = NULL
            pZ_s[i] = NULL
            pT_s[i] = NULL
        }
    }

    // v0.9.30: instrument columns spanned by the per-period constants (a
    // regressor or iv() variable common to all units in a period) are zeroed
    // here and removed by the all-zero drop below. Z is gamma-invariant, so
    // the same columns go at every grid point; xdpt_ivc_drop keeps the count
    // of the last call for the output note. Placed before the td partialling,
    // which only acts column by column, so no other column is affected.
    // The mask is computed once per distinct Z: Z does not depend on gamma,
    // so later stacks of the same run reuse it under an exact bitwise guard
    // (the reference copy is kept only while Z has at most 2e7 cells).
    external real scalar xdpt_ivc_drop
    external real matrix xdpt_ivc_Zref
    external real rowvector xdpt_ivc_mask
    external real colvector xdpt_ivc_tref
    real scalar _ivc_hit, _ivc_j
    real matrix _ivc_Z0
    xdpt_ivc_drop = 0
    if (rows(Z_out) > 1 & cols(Z_out) > 1) {
        _ivc_hit = 0
        if (rows(xdpt_ivc_Zref) == rows(Z_out) &
            cols(xdpt_ivc_Zref) == cols(Z_out) &
            rows(xdpt_ivc_tref) == rows(times_out)) {
            _ivc_hit = (xdpt_ivc_tref == times_out)
            if (_ivc_hit) _ivc_hit = (xdpt_ivc_Zref == Z_out)
        }
        if (_ivc_hit) {
            for (_ivc_j = 1; _ivc_j <= cols(Z_out); _ivc_j++) {
                if (xdpt_ivc_mask[_ivc_j]) {
                    Z_out[., _ivc_j] = J(rows(Z_out), 1, 0)
                    xdpt_ivc_drop = xdpt_ivc_drop + 1
                }
            }
        }
        else {
            if (rows(Z_out) * cols(Z_out) <= 2e7) _ivc_Z0 = Z_out
            xdpt_ivc_drop = xdpt2_drop_cellconst(Z_out, times_out)
            if (rows(Z_out) * cols(Z_out) <= 2e7) {
                xdpt_ivc_Zref = _ivc_Z0
                xdpt_ivc_tref = times_out
                xdpt_ivc_mask = (colsum(abs(_ivc_Z0)) :> 0) :& (colsum(abs(Z_out)) :== 0)
            }
            else {
                xdpt_ivc_Zref = J(0, 0, .)
                xdpt_ivc_tref = J(0, 1, .)
                xdpt_ivc_mask = J(1, 0, .)
            }
        }
    }

    // v0.9.33 (F4): when the cache builder asks for it (xdpt_tpl_rec = 1),
    // record the stacked rows before the time-effect partialling, dW at that
    // point, and the gamma-invariant pieces of the partialling, so that
    // xdpt2_tpl_dW can rebuild dW at another gamma without restacking.
    external real scalar xdpt_tpl_rec, xdpt_tpl_td
    external real colvector xdpt_tpl_upre, xdpt_tpl_tpre, xdpt_tpl_keep
    external real matrix xdpt_tpl_dWpre, xdpt_tpl_D, xdpt_tpl_DtD
    if (xdpt_tpl_rec == 1) {
        xdpt_tpl_upre  = unit_id_out
        xdpt_tpl_tpre  = times_out
        xdpt_tpl_dWpre = dW_out
        xdpt_tpl_td    = 0
        xdpt_tpl_keep  = J(0, 1, .)
        xdpt_tpl_D     = J(0, 0, .)
        xdpt_tpl_DtD   = J(0, 0, .)
    }

    // v0.7.13 (audit R4, C2): FWL time-dummy partialling — demean dY, dW(γ),
    // and Z within each time cell. Placed BEFORE the zero-column drop so
    // instrument columns that are constant within every time cell (which
    // demean to numerical dust) can be snapped to exact zero — RELATIVE to
    // their pre-demeaning scale, so the test is scale-free — and dropped.
    external real scalar xdpt_td_fwl
    // v0.9.31: base regressor columns (the first K of dW, which do not depend
    // on gamma) that the time-effect partialling removes: a regressor of the
    // form a_i + g_t (firm age = year - founding year) survives the unit-level
    // and the period-level constancy checks, but FD/FOD removes a_i and td
    // removes g_t. Recorded here, reported by xtdpthresh_run (error 498).
    external real rowvector xdpt_td_gone
    xdpt_td_gone = J(1, 0, .)
    if (xdpt_td_fwl == 1 & rows(Z_out) > 0) {
        real rowvector _preZ, _preW
        real scalar _cj
        _preZ = colsum(abs(Z_out))
        _preW = colsum(abs(dW_out))
        if (method == "fd") {
            // FD: the transformed time effect (delta-lambda_t) is common to
            // every unit at each t, so within-time demeaning is EXACT.
            xdpt2_demean_bytime(dY_out,    times_out)
            xdpt2_demean_bytime(dW_out,    times_out)
            xdpt2_demean_bytime(Z_out,     times_out)
            // v0.8.0 (audit R5): singleton time cells demean to all-zero
            // rows -- no moment information, but they inflate n_rows, the
            // residual pools, and e(N_stack). Drop them (IDs/times aligned).
            // The drop set depends only on times (gamma-invariant), so
            // per-gamma row counts stay matched across bootstrap caches.
            real colvector _ut2, _idx2, _keepr
            real scalar _ti2
            _ut2 = uniqrows(times_out)
            _keepr = J(rows(times_out), 1, 1)
            for (_ti2 = 1; _ti2 <= rows(_ut2); _ti2++) {
                _idx2 = selectindex(times_out :== _ut2[_ti2])
                if (rows(_idx2) == 1) _keepr[_idx2] = 0
            }
            if (xdpt_tpl_rec == 1) {
                xdpt_tpl_td   = 1
                xdpt_tpl_keep = selectindex(_keepr)
            }
            if (sum(_keepr) < rows(times_out)) {
                _idx2 = selectindex(_keepr)
                dY_out      = dY_out[_idx2]
                dW_out      = dW_out[_idx2, .]
                Z_out       = Z_out[_idx2, .]
                times_out   = times_out[_idx2]
                unit_id_out = unit_id_out[_idx2]
            }
        }
        else {
            // v0.8.0 (audit R5 #6 FIX, FOD): on unbalanced panels the
            // FOD-transformed time effect c_it(lambda_t - mean of future
            // lambda_s) varies with each unit's future-observation set, so
            // within-time demeaning is NOT exact. Build the FOD transform of
            // each time dummy row-by-row (same operator the data received:
            // +c at t, -c/Tf at each future date) and partial the stacked
            // system out on that matrix. The projection uses invsym, which
            // handles the structural rank deficiency (the FOD transform of
            // the all-ones column is zero, so the dummy columns sum to zero).
            real scalar _rr, _uu, _jj, _kk, _nf, _cc, _tp
            real matrix _D, _DtD, _Dcoef
            real colvector _dcols
            struct xdpt2_unit scalar _up
            external real colvector xdpt_teq
            // v0.9.1 R17 (#3): dummy columns indexed by rank in the global
            // observed-equation-time vector -- NO span-sized allocation of
            // any kind (R16 still used a span-length marker vector, which a
            // millisecond-scale delta-1 index could blow up). Never-touched
            // columns stay all-zero and are dropped by the _dcols filter
            // exactly as before; column order follows time order, so the
            // projection is bit-for-bit unchanged on any panel.
            _D = J(rows(times_out), max((rows(xdpt_teq), 1)), 0)
            for (_rr = 1; _rr <= rows(times_out); _rr++) {
                _uu = unit_id_out[_rr]
                _up = units[_uu]
                _jj = xdpt2_find_t(_up, times_out[_rr])
                if (_jj == 0) continue
                // v0.8.1: the FOD operator applied to the data averages over
                // future EQUATION rows only -- the dummies must receive the
                // identical operator.
                _nf = 0
                for (_kk = _jj + 1; _kk <= rows(_up.t); _kk++) {
                    if (_up.eq[_kk]) _nf = _nf + 1
                }
                if (_nf < 1) continue
                _cc = sqrt(_nf / (_nf + 1))
                _tp = xdpt2_tpos(xdpt_teq, times_out[_rr])
                if (_tp <= rows(xdpt_teq)) {
                    if (xdpt_teq[_tp] == times_out[_rr]) _D[_rr, _tp] = _cc
                }
                for (_kk = _jj + 1; _kk <= rows(_up.t); _kk++) {
                    if (!_up.eq[_kk]) continue
                    _tp = xdpt2_tpos(xdpt_teq, _up.t[_kk])
                    if (_tp > rows(xdpt_teq)) continue
                    if (xdpt_teq[_tp] != _up.t[_kk]) continue
                    _D[_rr, _tp] = _D[_rr, _tp] - _cc / _nf
                }
            }
            _dcols = selectindex(colsum(abs(_D))' :> 0)
            if (rows(_dcols) > 0) {
                _D = _D[., _dcols']
                _DtD = invsym(_D' * _D)
                _Dcoef = _DtD * (_D' * dY_out)
                dY_out = dY_out - _D * _Dcoef
                _Dcoef = _DtD * (_D' * dW_out)
                dW_out = dW_out - _D * _Dcoef
                _Dcoef = _DtD * (_D' * Z_out)
                Z_out  = Z_out  - _D * _Dcoef
                // v0.8.5 R14 (#3): rows the projection SATURATES (leverage
                // ~= 1: e.g. a time cell whose FOD-transformed dummy row is
                // unique on an unbalanced panel) are annihilated by
                // M = I - P. They carry no moment information but would
                // still inflate e(N_stack), the moment-covariance
                // normalization, the residual pools, and the AR/bootstrap
                // diagnostics. Drop them with a gamma-INVARIANT mask built
                // from _D alone (never from dW_out, which varies with gamma
                // and would fragment the estimation sample across grid
                // points). Mirrors the fd-td singleton-cell drop above.
                real colvector _Hd, _keepr2
                _Hd = rowsum((_D * _DtD) :* _D)
                _keepr2 = selectindex((1 :- _Hd) :> 1e-10)
                if (xdpt_tpl_rec == 1) {
                    xdpt_tpl_td   = 2
                    xdpt_tpl_D    = _D
                    xdpt_tpl_DtD  = _DtD
                    xdpt_tpl_keep = _keepr2
                }
                if (rows(_keepr2) < rows(dY_out)) {
                    dY_out      = dY_out[_keepr2, .]
                    dW_out      = dW_out[_keepr2, .]
                    Z_out       = Z_out[_keepr2, .]
                    times_out   = times_out[_keepr2]
                    unit_id_out = unit_id_out[_keepr2]
                }
            }
        }
        for (_cj = 1; _cj <= cols(Z_out); _cj++) {
            if (_preZ[_cj] > 0) {
                if (colsum(abs(Z_out[., _cj])) < 1e-10 * _preZ[_cj]) {
                    Z_out[., _cj] = J(rows(Z_out), 1, 0)
                }
            }
        }
        real rowvector _postW
        _postW = colsum(abs(dW_out))
        for (_cj = 1; _cj <= K & _cj <= cols(dW_out); _cj++) {
            if (_preW[_cj] > 0 & _postW[_cj] <= 1e-10 * _preW[_cj]) {
                xdpt_td_gone = xdpt_td_gone, _cj
            }
        }
    }

    // Drop all-zero columns of transformed Z
    if (rows(Z_out) > 0) {
        real rowvector col_sum, keep_idx
        col_sum = colsum(abs(Z_out))
        // v0.7.13 (audit): drop only columns that are EXACTLY all-zero
        // (structurally never populated — a lag depth no unit reaches). A
        // valid instrument on a tiny numeric scale sums to a small-but-
        // positive value and is retained; the old 1e-12 floor could delete it.
        keep_idx = selectindex(col_sum :> 0)
        if (length(keep_idx) > 0 & length(keep_idx) < cols(Z_out)) {
            Z_out = Z_out[., keep_idx]
        }
    }

    // v0.9.28: drop exact duplicate instrument columns (first kept). Z does
    // not depend on gamma, so the same columns are dropped at every grid
    // point. When no column is duplicated, Z_out is left untouched.
    if (rows(Z_out) > 0 & cols(Z_out) > 1) {
        real rowvector keep_nd
        keep_nd = xdpt2_nodup_cols(Z_out)
        if (length(keep_nd) < cols(Z_out)) Z_out = Z_out[., keep_nd]
    }

    // v0.9.34: drop instrument columns that are linear combinations of
    // earlier ones (xdpt2_indep_cols); xdpt_ivc_dep counts them for the
    // output note. Z does not depend on gamma, so the same columns go at
    // every grid point, and a Z that passes the rank gate is untouched.
    external real scalar xdpt_ivc_dep
    xdpt_ivc_dep = 0
    if (rows(Z_out) > 0 & cols(Z_out) > 1) {
        real rowvector keep_dep
        keep_dep = xdpt2_indep_cols(Z_out)
        if (length(keep_dep) < cols(Z_out)) {
            xdpt_ivc_dep = cols(Z_out) - length(keep_dep)
            Z_out = Z_out[., keep_dep]
        }
    }

}


// Helper: build MA(1)-aware first-stage weight matrix for FD GMM.
// Based on xthenreg's GMM_W_n_con (Seo-Shin 2016).
// Under iid ε with σ²=1: Var(Δε_t) = 2, Cov(Δε_t, Δε_{t-1}) = -1.
// Weight W_n = [Var(g)]^{-1} ∝ (2·W2 - W1 - W1')^{-1}
// where W2 has diagonal blocks (Σ_t Z_t' Z_t / N) and W1 has off-diagonal
// blocks (Σ_t Z_{t-1}' Z_t / N) for consecutive time pairs.
real matrix xdpt2_build_W_ma1(real matrix Z, real colvector times,
                                real colvector unit_id)
{
    real scalar k_iv, n_u, i, j, n_rows, inv_ok
    real matrix W2, W1, W_mat, ma1_struct
    real colvector rows_i, times_i
    real matrix Z_i

    k_iv = cols(Z)
    n_u = max(unit_id)
    n_rows = rows(Z)

    // W2: diagonal blocks (sum across all rows with block structure preserved
    // by Z's block-diagonal IV layout; thanks to that, Z'Z already has proper
    // block-diagonal form with nonzero diag blocks)
    W2 = Z' * Z / n_u

    // W1: within-unit consecutive-time cross products
    // v0.7.8 (SPEEDUP, bit-for-bit): the consecutive pairs are located with
    // ONE vectorized pass over the stack instead of one selectindex() scan of
    // the full unit_id vector per unit (O(N_units x n_rows) comparisons, plus
    // a per-unit copy of the unit's Z block). The stack is unit-ordered with
    // time ascending within unit, so the pair list enumerates in EXACTLY the
    // order of the old double loop; the rank-1 accumulation below adds the
    // identical terms in the identical order -> W1 is bit-for-bit unchanged.
    // (Do NOT replace the loop with the single matmul
    //  Z[pair_j :- 1, .]' * Z[pair_j, .] without re-certifying: it is
    //  algebraically identical but BLAS accumulation order may differ.)
    W1 = J(k_iv, k_iv, 0)
    if (n_rows >= 2) {
        real colvector pair_j
        pair_j = selectindex(
              (unit_id[|2 \ n_rows|] :== unit_id[|1 \ n_rows - 1|])
           :& ((times[|2 \ n_rows|] :- times[|1 \ n_rows - 1|]) :== 1))
        if (rows(pair_j) > 0) {
            pair_j = pair_j :+ 1   // shift to index of the LATER row of the pair
            for (j = 1; j <= rows(pair_j); j++) {
                // Consecutive pair: add z_{t-1}' z_t
                W1 = W1 + Z[pair_j[j] - 1, .]' * Z[pair_j[j], .]
            }
        }
    }
    W1 = W1 / n_u

    // MA(1) structure: 2·W2 - W1 - W1'
    ma1_struct = 2 * W2 - W1 - W1'
    xdpt2_syminv(ma1_struct, inv_ok, W_mat)
    if (!inv_ok) {
        // v0.7.0 (D6): fallback chain. invsym() on a singular W2 would
        // silently zero out rows/columns; use the identity weight (valid,
        // merely inefficient) when even Z'Z is ill-conditioned.
        // v0.9.29: record the fallback (1 = (Z'Z)^-1, 2 = identity) so that
        // it is reported in e(W1_fallback) instead of passing silently.
        external real scalar xdpt_w1_fallback
        xdpt2_syminv(W2, inv_ok, W_mat)
        if (inv_ok) {
            if (xdpt_w1_fallback < 1) xdpt_w1_fallback = 1
            return(W_mat)
        }
        xdpt_w1_fallback = 2
        return(I(k_iv))
    }
    return(W_mat)
}

// v0.7.13 (audit R4, C4): xdpt2_solve_gmm() DELETED. It had no live callers
// after the bootstrap 2-step fallbacks were removed for criterion symmetry;
// the two-stage estimation path uses xdpt2_solve_gmm_1step(_pre) +
// xdpt2_build_cluster_omega + xdpt2_gmm_sandwich directly. See git/_archive
// history for the last version of the standalone solver.

// Helper: 1-step GMM with fixed weight matrix. Returns obj and theta only.
void xdpt2_solve_gmm_1step(real colvector Y, real matrix W_reg, real matrix Z,
                            real matrix W_wt,
                            real scalar ok, real colvector theta,
                            real scalar obj, real matrix V)
{
    real scalar n_rows, inv_ok
    real matrix ZW, A, Ainv
    real colvector ZY, r, g_bar

    ok = 0
    if (rows(Y) < 20) return
    n_rows = rows(Y)

    ZW = Z' * W_reg / n_rows
    ZY = Z' * Y / n_rows

    A = ZW' * W_wt * ZW
    xdpt2_syminv(A, inv_ok, Ainv)
    if (!inv_ok) return
    theta = Ainv * ZW' * W_wt * ZY
    if (hasmissing(theta)) return
    r = Y - W_reg * theta
    if (hasmissing(r)) return
    g_bar = Z' * r / n_rows
    if (hasmissing(g_bar)) return
    obj = n_rows * (g_bar' * W_wt * g_bar)
    V = Ainv / n_rows
    if (obj >= . | hasmissing(V)) return
    ok = 1
}

// v0.7.9 (C): variant of xdpt2_solve_gmm_1step taking the two cross
// products ZW = Z'W_reg/n and ZY = Z'Y/n precomputed. The cache stores them
// computed by the very same expressions on the very same matrices this
// function would use, so A, theta, and V are bitwise identical; the
// residual r and moment g_bar are kept VERBATIM on the n-row data to
// preserve floating-point association, so obj is bitwise identical too.
void xdpt2_solve_gmm_1step_pre(real colvector Y, real matrix W_reg,
                                real matrix Z,
                                real matrix ZW, real colvector ZY,
                                real matrix W_wt,
                                real scalar ok, real colvector theta,
                                real scalar obj, real matrix V)
{
    real scalar n_rows, inv_ok
    real matrix A, Ainv
    real colvector r, g_bar

    ok = 0
    if (rows(Y) < 20) return
    n_rows = rows(Y)

    A = ZW' * W_wt * ZW
    xdpt2_syminv(A, inv_ok, Ainv)
    if (!inv_ok) return
    theta = Ainv * ZW' * W_wt * ZY
    if (hasmissing(theta)) return
    r = Y - W_reg * theta
    if (hasmissing(r)) return
    g_bar = Z' * r / n_rows
    if (hasmissing(g_bar)) return
    obj = n_rows * (g_bar' * W_wt * g_bar)
    V = Ainv / n_rows
    if (obj >= . | hasmissing(V)) return
    ok = 1
}

// v0.7.8 (SPEEDUP): per-unit sum of moment rows, run-based.
// Replaces the row-by-row accumulation loop (one interpreted iteration plus
// two row-vector extract/store copies PER ROW) with one colsum() per
// contiguous same-unit run. The FD/FOD stack is unit-contiguous by
// construction (one run per contributing unit), and
// run subtotals are added into g_per_unit in that same order. colsum()
// accumulates top-down in double precision (Mata keeps quad accumulation in
// the separate quadcolsum()), so the result is expected bit-for-bit
// identical to the scalar loop -- certify with verify_speedup before
// relying; the scalar loop is preserved below as the reference fallback.
real matrix xdpt2_gsum_by_unit(real matrix Ze, real colvector unit_id,
                                real scalar n_units)
{
    real scalar n_rows, k, rr, u
    real matrix g_per_unit
    real colvector bnd, s_idx, e_idx

    n_rows = rows(Ze)
    k = cols(Ze)
    // v0.9.35 (external review, R3): colsum() below skips missing values; a
    // missing input (an overflow) makes every sum missing instead
    if (hasmissing(Ze)) return(J(n_units, k, .))
    g_per_unit = J(n_units, k, 0)
    if (n_rows == 0) return(g_per_unit)
    if (n_rows == 1) {
        g_per_unit[unit_id[1], .] = Ze[1, .]
        return(g_per_unit)
    }

    // Run boundaries: row rr ends a run iff unit_id[rr+1] != unit_id[rr]
    bnd = selectindex(unit_id[|2 \ n_rows|] :!= unit_id[|1 \ n_rows - 1|])
    s_idx = 1 \ (bnd :+ 1)      // run starts
    e_idx = bnd \ n_rows        // run ends
    for (rr = 1; rr <= rows(s_idx); rr++) {
        u = unit_id[s_idx[rr]]
        if (e_idx[rr] > s_idx[rr]) {
            g_per_unit[u, .] = g_per_unit[u, .] + colsum(Ze[|s_idx[rr], 1 \ e_idx[rr], k|])
        }
        else {
            g_per_unit[u, .] = g_per_unit[u, .] + Ze[s_idx[rr], .]
        }
    }
    return(g_per_unit)

    // --- original scalar loop (bit-for-bit reference fallback) ---
    // g_per_unit = J(n_units, cols(Ze), 0)
    // for (i = 1; i <= rows(Ze); i++) {
    //     u = unit_id[i]
    //     g_per_unit[u, .] = g_per_unit[u, .] + Ze[i, .]
    // }
    // return(g_per_unit)
}

// Helper: compute cluster-robust Ω from residuals (unit-level clustering)
real matrix xdpt2_build_cluster_omega(real matrix Z, real colvector r,
                                        real colvector unit_id)
{
    external real scalar xdpt_center
    return(xdpt2_build_cluster_omega_c(Z, r, unit_id, xdpt_center))
}

// v0.9.30: the same computation with the centering flag as an argument, so
// that a per-draw caller (the coefficient bootstrap) does not bind the
// -external- at every call (see xdpt2_transform_unit). Bit-for-bit the body
// of the former xdpt2_build_cluster_omega.
real matrix xdpt2_build_cluster_omega_c(real matrix Z, real colvector r,
                                          real colvector unit_id,
                                          real scalar xdpt_center)
{
    real scalar n_rows, n_units
    real matrix Ze, g_per_unit, Omega

    n_rows = rows(Z)
    n_units = max(unit_id)
    Ze = Z :* r
    // v0.7.8: run-based per-unit aggregation (see xdpt2_gsum_by_unit)
    g_per_unit = xdpt2_gsum_by_unit(Ze, unit_id, n_units)
    Omega = g_per_unit' * g_per_unit / n_rows
    // v0.8.0 (audit R5, -center-): subtract the mean-moment outer product
    // (Seo-Shin 2016 eq. 11 / xthenreg convention), expressed consistently
    // with this function's per-n_rows scaling: with s = sum_i g_i,
    // sum_i (g_i - s/n_u)(g_i - s/n_u)' = sum_i g_i g_i' - s s'/n_u.
    if (xdpt_center == 1) {
        real rowvector s_c
        real scalar n_contrib
        // v0.8.1 (audit R6): unit_id keeps ORIGINAL unit indices, so a unit
        // with no transformed rows leaves a gap and max(unit_id) overcounts
        // the clusters (g_per_unit then holds all-zero ghost rows -- they do
        // not change Omega or s_c, but they would inflate this denominator).
        n_contrib = rows(uniqrows(unit_id))
        s_c = colsum(g_per_unit)
        Omega = Omega - (s_c' * s_c) / (n_contrib * n_rows)
    }
    return(Omega)
}

// Cluster-robust sandwich variance for the estimator that actually used A:
// V = B (D' A Omega A D) B / n, B=(D' A D)^-1, D=Z'X/n.
// This remains valid when A is a preliminary fixed weight and therefore must
// be used instead of pretending that a subsequently recomputed Omega^-1 was
// paired with theta_hat.
real matrix xdpt2_gmm_sandwich(real matrix ZW, real matrix A,
                                real matrix Omega, real scalar n_rows)
{
    real scalar inv_ok
    real matrix bread0, B, meat, Vout
    if (rows(A) != rows(Omega) | cols(A) != cols(Omega)) return(J(0, 0, .))
    bread0 = ZW' * A * ZW
    xdpt2_syminv(bread0, inv_ok, B)
    if (!inv_ok) return(J(0, 0, .))
    meat = ZW' * A * Omega * A * ZW
    if (hasmissing(meat)) return(J(0, 0, .))
    Vout = B * meat * B / n_rows
    if (hasmissing(Vout)) return(J(0, 0, .))
    return(Vout)
}

// v0.9.33 (F4, SPEEDUP): rebuild dW at a new gamma without restacking. Of
// the stacked system only dW depends on gamma, and within dW only the regime
// columns: the row selection (missing values, instrument availability), dY,
// Z, the base columns of dW and the pieces of the time-effect partialling do
// not. xdpt2_tpl_build takes them from one full xdpt2_stack_at_gamma call
// (which records the rows before the partialling) and finds the level rows
// behind every stacked row. xdpt2_tpl_dW then forms the regime columns with
// the arithmetic of xdpt2_transform_unit -- FD: W(t) - W(t-1); FOD:
// c*(W(t) - mean of the later equation rows), the mean as a quad-precision
// sum over those rows, in increasing order, divided by their number, which
// is how mean() computes it (quadcross) -- and applies the recorded
// partialling. v0.9.34: the sums run forward, as in mean(); 0.9.33 summed
// them in reverse order, which differs in the last bit when a unit's values
// span more than about 2^35 and cancel exactly. xdpt2_build_gamma_cache
// compares the result with a full stack at its second gamma and restacks at
// every gamma if any bit differs.
struct xdpt2_stack_tpl {
    real scalar    ok
    real scalar    fd, kink, td, K
    real colvector qlev        // q at the level rows of all units, stacked
    real matrix    Xlev        // X at the level rows
    real colvector J, JP       // level rows of t and (FD) t-1, per stacked row
    real colvector c, Tf, iF   // FOD: scale, number of later rows, row of F
    real colvector EQF         // FOD: level rows of each unit's equation
                               // rows, in increasing order, units stacked
    real matrix    RU          // FOD: (first, last) row of each unit in EQF
    real matrix    FG          // FOD: units grouped by their number v >= 2
                               // of equation rows: (v, first, last) in FP
    real colvector FP          // FOD: positions in EQF, by group, unit blocks
    real matrix    dXpre       // base columns of dW before the partialling
    real colvector tpre, keep  // partialling: times, rows kept after it
    real matrix    D, DtD      // FOD partialling: dummies and inverse
    real colvector dY, times, uid
    real matrix    Z
}

struct xdpt2_stack_tpl scalar xdpt2_tpl_build(
    struct xdpt2_unit rowvector units,
    string scalar method,
    real scalar flag_kink,
    real colvector dY, real matrix Z,
    real colvector times, real colvector uid)
{
    struct xdpt2_stack_tpl scalar tp
    external real scalar xdpt_tpl_td
    external real colvector xdpt_tpl_upre, xdpt_tpl_tpre, xdpt_tpl_keep
    external real matrix xdpt_tpl_dWpre, xdpt_tpl_D, xdpt_tpl_DtD
    real scalar nu, n, i, r, u, j, jp, m, ne, a, gi, v, nuv, p0
    real colvector nlev, off, ei, neq, qpos, vv, uu

    tp.ok = 0
    nu = length(units)
    n = rows(xdpt_tpl_upre)
    if (nu == 0 | n == 0) return(tp)
    tp.K = cols(units[1].X)
    tp.fd = (method == "fd")
    tp.kink = flag_kink
    // v0.9.35 (external review, R4): K = 0 -- a static model without base
    // regressors, whose only regime term is the intercept (or the hinge
    // under kink) -- is allowed; the joint variance and the AR statistics
    // with gamma-hat use this template
    if (rows(xdpt_tpl_tpre) != n | rows(xdpt_tpl_dWpre) != n) return(tp)
    if (cols(xdpt_tpl_dWpre) != (flag_kink ? tp.K + 1 : 2 * tp.K + 1)) return(tp)

    // level data of all units, stacked in unit order
    nlev = J(nu, 1, 0)
    for (i = 1; i <= nu; i++) nlev[i] = rows(units[i].q)
    off = J(nu, 1, 0)
    for (i = 2; i <= nu; i++) off[i] = off[i - 1] + nlev[i - 1]
    tp.qlev = J(sum(nlev), 1, .)
    tp.Xlev = J(sum(nlev), tp.K, .)
    for (i = 1; i <= nu; i++) {
        if (nlev[i] == 0) continue
        if (rows(units[i].X) != nlev[i] | cols(units[i].X) != tp.K) return(tp)
        tp.qlev[|off[i] + 1 \ off[i] + nlev[i]|] = units[i].q
        if (tp.K > 0) {
            tp.Xlev[|off[i] + 1, 1 \ off[i] + nlev[i], tp.K|] = units[i].X
        }
    }

    // FOD: every unit's equation rows, in order, stacked (for later sums)
    if (!tp.fd) {
        neq = J(nu, 1, 0)
        for (i = 1; i <= nu; i++) neq[i] = sum(units[i].eq :!= 0)
        tp.RU = J(nu, 2, 0)
        tp.EQF = J(sum(neq), 1, .)
        a = 0
        for (i = 1; i <= nu; i++) {
            ne = neq[i]
            tp.RU[i, .] = (a + 1, a + ne)
            if (ne > 0) {
                ei = selectindex(units[i].eq)
                if (rows(ei) != ne) return(tp)
                tp.EQF[|a + 1 \ a + ne|] = off[i] :+ ei
            }
            a = a + ne
        }
        // units with the same number v >= 2 of equation rows form a group;
        // FP lists the EQF positions of its units, one block of v per unit
        vv = select(neq, neq :>= 2)
        if (rows(vv) > 0) vv = uniqrows(vv)
        tp.FG = J(rows(vv), 3, .)
        tp.FP = J(0, 1, .)
        p0 = 0
        for (gi = 1; gi <= rows(vv); gi++) {
            v = vv[gi]
            uu = selectindex(neq :== v)
            nuv = rows(uu)
            tp.FP = tp.FP \ vec(J(v, 1, 1) * tp.RU[uu, 1]' :+
                                (0::v - 1) * J(1, nuv, 1))
            tp.FG[gi, .] = (v, p0 + 1, p0 + v * nuv)
            p0 = p0 + v * nuv
        }
        tp.c  = J(n, 1, .)
        tp.Tf = J(n, 1, .)
        tp.iF = J(n, 1, .)
    }

    // the level rows behind every stacked row
    tp.J = J(n, 1, .)
    if (tp.fd) tp.JP = J(n, 1, .)
    for (r = 1; r <= n; r++) {
        u = xdpt_tpl_upre[r]
        if (u < 1 | u > nu | u != trunc(u)) return(tp)
        j = xdpt2_find_t(units[u], xdpt_tpl_tpre[r])
        if (j == 0) return(tp)
        tp.J[r] = off[u] + j
        if (tp.fd) {
            jp = xdpt2_find_t(units[u], xdpt_tpl_tpre[r] - 1)
            if (jp == 0) return(tp)
            tp.JP[r] = off[u] + jp
        }
        else {
            // position m of row j among the unit's equation rows; the later
            // equation rows are those after it
            ne = tp.RU[u, 2] - tp.RU[u, 1] + 1
            if (ne < 2) return(tp)
            qpos = selectindex(tp.EQF[|tp.RU[u, 1] \ tp.RU[u, 2]|] :== off[u] + j)
            if (rows(qpos) != 1) return(tp)
            m = qpos[1]
            tp.Tf[r] = ne - m
            if (tp.Tf[r] < 1) return(tp)
            tp.c[r]  = sqrt(tp.Tf[r] / (tp.Tf[r] + 1))
            tp.iF[r] = tp.RU[u, 1] - 1 + m
        }
    }

    if (tp.K > 0) tp.dXpre = xdpt_tpl_dWpre[., 1..tp.K]
    else          tp.dXpre = J(n, 0, .)
    tp.td    = xdpt_tpl_td
    tp.tpre  = xdpt_tpl_tpre
    tp.keep  = xdpt_tpl_keep
    tp.D     = xdpt_tpl_D
    tp.DtD   = xdpt_tpl_DtD
    if (tp.td != 0 & tp.td != 1 & tp.td != 2) return(tp)
    tp.dY    = dY
    tp.Z     = Z
    tp.times = times
    tp.uid   = uid
    tp.ok    = 1
    return(tp)
}

real matrix xdpt2_tpl_dW(struct xdpt2_stack_tpl scalar tp, real scalar gamma)
{
    real colvector r, ix
    real matrix Wreg, dreg, dW, F, WE, X, S
    real scalar i, j, v
    r = (tp.qlev :> gamma)
    if (tp.kink) Wreg = (tp.qlev :- gamma) :* r
    else         Wreg = r, tp.Xlev :* r
    if (tp.fd) {
        dreg = Wreg[tp.J, .] - Wreg[tp.JP, .]
    }
    else {
        // row m of a unit's block of F: the sum of its equation rows after
        // row m, accumulated forward in quad precision (quadcross, the sum
        // mean() forms; the zero terms of S leave the sums unchanged). The
        // units of a group (same v) are summed by one quadcross over their
        // v x (units) blocks: every element is the same forward quad sum.
        WE = Wreg[tp.EQF, .]
        F = J(rows(tp.EQF), cols(Wreg), .)
        for (i = 1; i <= rows(tp.FG); i++) {
            v  = tp.FG[i, 1]
            ix = tp.FP[|tp.FG[i, 2] \ tp.FG[i, 3]|]
            S  = lowertriangle(J(v, v, 1), 0)
            for (j = 1; j <= cols(Wreg); j++) {
                X = colshape(WE[ix, j], v)'
                F[ix, j] = vec(quadcross(S, X))
            }
        }
        dreg = tp.c :* (Wreg[tp.J, .] - F[tp.iF, .] :/ tp.Tf)
    }
    dW = tp.dXpre, dreg
    if (tp.td == 1) {
        xdpt2_demean_bytime(dW, tp.tpre)
        dW = dW[tp.keep, .]
    }
    else if (tp.td == 2) {
        dW = dW - tp.D * (tp.DtD * (tp.D' * dW))
        dW = dW[tp.keep, .]
    }
    return(dW)
}

// v0.9.34: number of units whose regime 1(q > gamma) changes within their
// transformed equations -- under FD between t-1 and t, under FOD between t
// and a later equation period. Only these units identify the regime terms
// (the regime moments and the clusters of their weight); rows are the final
// stack at gamma (uid, times).
real scalar xdpt2_n_switch(struct xdpt2_unit rowvector units,
                           real scalar gamma, string scalar method,
                           real colvector uid, real colvector times)
{
    real scalar r, u, j, jp, k
    real colvector sw, ei
    sw = J(length(units), 1, 0)
    for (r = 1; r <= rows(uid); r++) {
        u = uid[r]
        if (u < 1 | u > length(units)) continue
        if (sw[u]) continue
        j = xdpt2_find_t(units[u], times[r])
        if (j == 0) continue
        if (method == "fd") {
            jp = xdpt2_find_t(units[u], times[r] - 1)
            if (jp == 0) continue
            if ((units[u].q[j] > gamma) != (units[u].q[jp] > gamma)) sw[u] = 1
        }
        else {
            ei = selectindex(units[u].eq)
            for (k = 1; k <= rows(ei); k++) {
                if (ei[k] <= j) continue
                if ((units[u].q[ei[k]] > gamma) != (units[u].q[j] > gamma)) {
                    sw[u] = 1
                    break
                }
            }
        }
    }
    return(sum(sw))
}

// v0.9.34: release the template records of xdpt2_stack_at_gamma (under
// FOD + td they include an n x T matrix); called at the start and the end
// of every run.
void xdpt2_tpl_release()
{
    external real scalar xdpt_tpl_rec, xdpt_tpl_td
    external real colvector xdpt_tpl_upre, xdpt_tpl_tpre, xdpt_tpl_keep
    external real matrix xdpt_tpl_dWpre, xdpt_tpl_D, xdpt_tpl_DtD
    xdpt_tpl_rec   = 0
    xdpt_tpl_td    = 0
    xdpt_tpl_upre  = J(0, 1, .)
    xdpt_tpl_tpre  = J(0, 1, .)
    xdpt_tpl_keep  = J(0, 1, .)
    xdpt_tpl_dWpre = J(0, 0, .)
    xdpt_tpl_D     = J(0, 0, .)
    xdpt_tpl_DtD   = J(0, 0, .)
}

// Per-gamma cache: avoid repeated xdpt2_stack_at_gamma calls.
// Used by both main grid search and bootstrap loops.
struct xdpt2_gamma_cache {
    real scalar    ok
    real scalar    gamma
    real colvector dY
    real matrix    dW
    // v0.7.13 (audit R4, C5): Z and W_first are γ-INVARIANT and were the
    // dominant per-entry memory (Z is n×L; W_first L×L). Mata struct
    // assignment REAL-copies matrices (verified empirically: 50 assignments
    // of an 80MB matrix peaked at ~4GB), so storing them by value duplicated
    // them across every grid/CI/kink cache entry. They are now heap objects
    // shared via pointers: entries that pass the bitwise reuse guard copy
    // the POINTER (8 bytes), not the data. Deref with (*e.pZ) / (*e.pW1).
    pointer(real matrix) scalar pZ
    pointer(real matrix) scalar pW1   // first-stage weight (MA(1) fd, ZZ_inv else)
    real colvector times
    real colvector uid
    // NEW: precomputed for the fast 1-step wild-bootstrap GMM (xthenreg-style)
    //   C_g = invsym(ZW' W_first ZW) * ZW' * W_first    // (k_W × n_iv)
    //   Given Y_boot, θ̂ = C_g · (Z'Y_boot/n); single matmul, no cluster-Ω loop.
    real matrix    C_g
    real scalar    n_rows     // cached rows(dY)
    real scalar    fast_ok    // 1 if C_g valid (non-singular A)
    // v0.7.9 (C): cross products stored at build time so the grid search
    // does not recompute them (they were rebuilt twice per gamma: stage 1
    // and stage 2). Set for every ok entry.
    real matrix    ZW         // Z'dW/n  (k_iv x k_W)
    real colvector ZY         // Z'dY/n  (k_iv x 1)
    // v0.9.34 (C1): the fixed-W2 solve of the confidence set, filled by
    // xdpt2_cache_w2 when the reported estimator is two-step
    real matrix    C_g2       // invsym(ZW' W2 ZW) ZW' W2
    real scalar    fast2_ok   // 1 if C_g2 valid
}

struct xdpt2_gamma_cache rowvector xdpt2_build_gamma_cache(
    struct xdpt2_unit rowvector units,
    real colvector gamma_grid,
    string scalar method,
    real scalar flag_static,
    real scalar flag_kink,
    real scalar t_min,
    real scalar t_max,
    real colvector q_supp,
    real scalar min_user)
{
    struct xdpt2_stack_tpl scalar tpl
    real scalar tpl_state
    tpl_state = 0
    return(xdpt2_build_gamma_cache_t(units, gamma_grid, method, flag_static,
                                     flag_kink, t_min, t_max, q_supp,
                                     min_user, tpl, tpl_state))
}

// v0.9.34: the template and its state are passed in and out, so that later
// builds of the same model (the refinements, the confidence-set grid) use the
// template of the main build and do not restack. With refc (the main cache),
// its first admitted entry is the reference of the exact reuse guards below,
// so its Z, W1 and ZY are shared instead of rebuilt.
struct xdpt2_gamma_cache rowvector xdpt2_build_gamma_cache_t(
    struct xdpt2_unit rowvector units,
    real colvector gamma_grid,
    string scalar method,
    real scalar flag_static,
    real scalar flag_kink,
    real scalar t_min,
    real scalar t_max,
    real colvector q_supp,
    real scalar min_user,
    struct xdpt2_stack_tpl scalar tpl,
    real scalar tpl_state,
    | struct xdpt2_gamma_cache rowvector refc)
{
    real scalar g, G, n_rows
    struct xdpt2_gamma_cache rowvector cache
    real colvector dY_cur, times_cur, uid_cur
    real matrix dW_cur, Z_cur, ZZ, ZZ_inv, W_first
    real matrix ZW_cur, A_cur, Ainv_cur
    real scalar reuse_w, reuse_y, reuse_tu, inv_ok
    // reference of the reuse guards; the comparisons of the template pieces
    // (the same objects at every template-built entry) are made once
    pointer(real matrix) scalar R_pZ, R_pW1
    real colvector R_dY, R_ZY, R_times, R_uid
    real scalar have_R, R_ext, k, from_tpl, tc_known, tc_tu, tc_w, tc_y
    real colvector ZY_cur
    real scalar min_reg
    external real scalar xdpt_trim_rate
    external real scalar xdpt_minreg_def
    // v0.9.33 (F4): tpl_state 0 = no template yet; 1 = template taken from
    // the first stack; 2 = its dW matched a second full stack bit for bit,
    // so dW is rebuilt by xdpt2_tpl_dW; -1 = full stack at every gamma.
    external real scalar xdpt_tpl_rec
    real matrix dW_chk

    G = rows(gamma_grid)
    cache = xdpt2_gamma_cache(1, G)
    have_R = 0
    R_ext = 0
    tc_known = 0
    if (args() == 12) {
        for (k = 1; k <= cols(refc); k++) {
            if (!refc[k].ok) continue
            R_pZ    = refc[k].pZ
            R_pW1   = refc[k].pW1
            R_dY    = refc[k].dY
            R_ZY    = refc[k].ZY
            R_times = refc[k].times
            R_uid   = refc[k].uid
            have_R  = 1
            R_ext   = 1
            break
        }
    }

    for (g = 1; g <= G; g++) {
        cache[g].ok = 0
        cache[g].fast_ok = 0
        cache[g].gamma = gamma_grid[g]

        from_tpl = (tpl_state == 2)
        if (tpl_state == 2) {
            dY_cur    = tpl.dY
            dW_cur    = xdpt2_tpl_dW(tpl, gamma_grid[g])
            Z_cur     = tpl.Z
            times_cur = tpl.times
            uid_cur   = tpl.uid
        }
        else {
            xdpt_tpl_rec = (tpl_state == 0)
            xdpt2_stack_at_gamma(units, gamma_grid[g], method,
                                  flag_static, flag_kink, t_min, t_max,
                                  dY_cur, dW_cur, Z_cur, times_cur, uid_cur)
            xdpt_tpl_rec = 0
            if (tpl_state == 0) {
                tpl = xdpt2_tpl_build(units, method, flag_kink, dY_cur,
                                      Z_cur, times_cur, uid_cur)
                tpl_state = (tpl.ok ? 1 : -1)
            }
            else if (tpl_state == 1) {
                tpl_state = -1
                if (rows(dY_cur) == rows(tpl.dY) &
                    rows(times_cur) == rows(tpl.times) &
                    rows(uid_cur) == rows(tpl.uid) &
                    rows(Z_cur) == rows(tpl.Z) & cols(Z_cur) == cols(tpl.Z)) {
                    dW_chk = xdpt2_tpl_dW(tpl, gamma_grid[g])
                    if (rows(dW_chk) == rows(dW_cur) &
                        cols(dW_chk) == cols(dW_cur)) {
                        if (dW_chk == dW_cur & dY_cur == tpl.dY &
                            Z_cur == tpl.Z & times_cur == tpl.times &
                            uid_cur == tpl.uid) tpl_state = 2
                    }
                }
            }
        }

        if (rows(dY_cur) == 0) continue
        n_rows = rows(dY_cur)
        // Gamma is an additional nonlinear parameter. With fewer than
        // k_W+1 moments, the linear coefficients can fit every fixed gamma
        // exactly and the threshold is not identified.
        if (cols(Z_cur) < cols(dW_cur) + 1) continue
        // v0.8.2 R9 (#2): the guard is a TRIMMING rule on the deduplicated
        // effective support (passed as an argument, R9 #7), NOT a parameter-
        // count rule. The old cols(dW)+1 per-side floor demanded 2K+2
        // observations on EACH side of every candidate gamma -- with K=6 it
        // silently shrank a trim(.10) search range to roughly [14%,86%],
        // excluding valid thresholds BEFORE the objective was evaluated
        // (identification is a rank condition on the moment matrices, and
        // the FD design moves through 1(q_t>g) AND 1(q_{t-1}>g); the kink
        // design has no second coefficient set at all). Rank/conditioning
        // is enforced where it belongs: cond(ZZ), cond(A_cur)/fast_ok, and
        // the stage guards. minregime(#) overrides the floor if set.
        real scalar n_supp
        n_supp = rows(q_supp)
        // v0.8.2 R10 (#3): minregime() is a FLOOR (max with the default
        // trim rule), not an override -- an option named "minimum" must
        // strengthen the safeguard, never weaken it below the default.
        // v0.9.34: the default floor set by xtdpthresh_run (the smaller
        // trimmed tail); the old rule if it is not set
        min_reg = xdpt_minreg_def
        if (min_reg >= .) min_reg = ceil(xdpt_trim_rate * n_supp / 2)
        if (min_reg < 2) min_reg = 2
        if (min_user > 0 & min_user > min_reg) min_reg = min_user
        if (sum(q_supp :<= gamma_grid[g]) < min_reg | ///
            sum(q_supp :>  gamma_grid[g]) < min_reg) continue
        // v0.7.9 (B): exact-guarded reuse of the gamma-invariant weight
        // pieces. Z, times, uid (and dY) do not depend on gamma by
        // construction, so ZZ, the cond(ZZ) admissibility decision, and
        // W_first (including the expensive MA(1) build under method(fd))
        // are identical across the grid. Rather than trusting that
        // invariant, each entry is compared BITWISE (matrix == matrix is a
        // scalar test in Mata) against the first admitted entry; on any
        // mismatch the entry takes the original fresh-compute path below,
        // so the stored values are bit-for-bit unchanged either way. The
        // fresh path also stops computing a dead invsym(ZZ) under
        // method(fd) (it was computed, then discarded for the MA(1)
        // weight) -- a pure dead-value elimination.
        // v0.9.34: times/uid, then Z, then dY; for template-built entries
        // the three results are those of the first such entry.
        reuse_tu = 0
        reuse_w = 0
        reuse_y = 0
        if (from_tpl & tc_known) {
            reuse_tu = tc_tu
            reuse_w  = tc_w
            reuse_y  = tc_y
        }
        else if (have_R) {
            reuse_tu = (times_cur == R_times & uid_cur == R_uid)
            if (reuse_tu) {
                reuse_w = (Z_cur == *R_pZ)
                if (reuse_w) reuse_y = (dY_cur == R_dY)
            }
            if (from_tpl) {
                tc_known = 1
                tc_tu = reuse_tu
                tc_w  = reuse_w
                tc_y  = reuse_y
            }
        }
        if (reuse_w) {
            // v0.7.13 (C5): share the reference entry's heap objects — a
            // pointer copy, not a data copy (this is the memory fix).
            cache[g].pZ  = R_pZ
            cache[g].pW1 = R_pW1
        }
        else {
            ZZ = Z_cur' * Z_cur / n_rows
            xdpt2_syminv(ZZ, inv_ok, ZZ_inv)
            if (!inv_ok) continue
            if (method == "fd") {
                W_first = xdpt2_build_W_ma1(Z_cur, times_cur, uid_cur)
            }
            else {
                W_first = ZZ_inv
            }
            // v0.7.13 (C5): (expr :+ 0) creates UNNAMED heap objects — a
            // pointer to a named local would dangle after this function
            // returns; a pointed-to unnamed object lives while referenced.
            cache[g].pZ  = &(Z_cur :+ 0)
            cache[g].pW1 = &(W_first :+ 0)
        }

        // v0.7.9 (C): ZY = Z'dY/n, reused from the reference entry under
        // the exact guard above (both Z and dY bitwise equal), computed by
        // the identical expression otherwise.
        if (reuse_y) ZY_cur = R_ZY
        else         ZY_cur = Z_cur' * dY_cur / n_rows

        cache[g].dY      = dY_cur
        cache[g].dW      = dW_cur
        cache[g].times   = times_cur
        cache[g].uid     = uid_cur
        cache[g].n_rows  = n_rows
        cache[g].ok      = 1
        cache[g].ZY      = ZY_cur
        if (!have_R | (R_ext & !reuse_w)) {
            // first admitted entry = bitwise reference (also when the
            // entries do not match the reference of refc)
            R_pZ    = cache[g].pZ
            R_pW1   = cache[g].pW1
            R_dY    = dY_cur
            R_ZY    = ZY_cur
            R_times = times_cur
            R_uid   = uid_cur
            have_R  = 1
            R_ext   = 0
            tc_known = 0
        }

        // Precompute C_g for fast bootstrap (1-step GMM with fixed W_first)
        //   θ(Y) = invsym(ZW' W_first ZW) · ZW' · W_first · (Z'Y/n)
        //        = C_g · (Z'Y/n)
        ZW_cur = Z_cur' * dW_cur / n_rows
        cache[g].ZW = ZW_cur   // v0.7.9 (C): stored for the grid search
        A_cur = ZW_cur' * (*cache[g].pW1) * ZW_cur
        xdpt2_syminv(A_cur, inv_ok, Ainv_cur)
        if (inv_ok) {
            cache[g].C_g = Ainv_cur * ZW_cur' * (*cache[g].pW1)
            if (!hasmissing(cache[g].C_g)) cache[g].fast_ok = 1
        }
    }
    return(cache)
}

// Fast bootstrap 1-step GMM: uses precomputed C_g from cache entry.
// O(n_rows·n_iv + n_iv²) per call vs O(n_rows·n_iv² + n_iv³ + k_W³) for
// xdpt2_solve_gmm. Drops cluster-Ω loop (bootstrap doesn't need 2-step).
// v0.7.13 (audit R4) label correction: this is the xthenreg-style FAST
// cluster wild residual bootstrap (unit-level Mammen weights, fixed W_first,
// 1-step GMM per draw) — NOT the exact Gong-Seo (2026) Algorithm 1, which
// resamples (x, z, resid) jointly at the unit level and recenters the
// bootstrap moments. Gong-Seo validity is proved for the exact algorithm;
// this scheme is supported by the Monte Carlo evidence in the paper.
// The fixed W_first shared between sample and bootstrap sides keeps the
// two statistics on the same criterion.
void xdpt2_fast_gmm_boot(real colvector Y_boot,
                          struct xdpt2_gamma_cache scalar gc,
                          real scalar ok, real colvector theta,
                          real scalar obj)
{
    real colvector ZY, r, g
    ok = 0
    if (gc.ok == 0 | gc.fast_ok == 0) return
    if (rows(Y_boot) != gc.n_rows) return

    ZY = (*gc.pZ)' * Y_boot / gc.n_rows
    theta = gc.C_g * ZY
    if (hasmissing(theta)) return
    r = Y_boot - gc.dW * theta
    if (hasmissing(r)) return
    g = (*gc.pZ)' * r / gc.n_rows
    if (hasmissing(g)) return
    obj = gc.n_rows * (g' * (*gc.pW1) * g)
    if (obj >= .) return
    ok = 1
}

// v0.9.31 (SPEEDUP): one-step objectives of the bootstrap samples
// Y* = F + E for the cache entries idx (in that order), one row per entry.
// F (n x 1) is the fit that generates every draw and E the reweighted
// residuals, E[., b] = R :* ETA[ud, b]: ud gives the cluster of each row and
// ETA (clusters x B) the cluster weights. The moment vector is linear in Y*,
// g(Y*) = g(F) + g(E), where
//   g(F) = Z'(F - dW theta_F)/n, theta_F = C Z'F/n   (from residuals, as before)
//   g(E) = Z'E/n - (Z'dW/n)(C Z'E/n)                  (from cross products)
// with Z'E/n and Z'F/n formed once per shared Z and Z'dW/n = ZW taken from the
// cache entry, next to the dW and C it belongs to. g(F) is O(n k) per entry
// and g(E) O(k^2 B), instead of O(n k B) for the n-row residual matrix and
// its second product. Forming g(F) from residuals matters: in a pure
// cross-product form, g(Y*) = Z'Y*/n - ZW theta, the fitted component of Y*
// cancels against ZW theta and that form is 10-50 times less accurate. Against
// a 40-digit reference computed from the same double inputs, the split is as
// accurate as the residual form (largest error about 1e-13 of the objective
// scale in both; _dev_0931 audit). Rows equal the residual form to rounding.
// v0.9.33 (SPEEDUP): all rows of a cluster share one weight, so
// Z'E/n = S'ETA/n, where S sums z_i*r_i within each cluster; E itself is
// never formed. The residuals of F at the entries of a block that shares Z
// are stacked as columns and multiplied by Z' at once, instead of one
// matrix-vector product per entry; a block holds at most 4e6 elements. The
// arithmetic is the same (the residual form for g(F)); only the order of
// summation changes, so rows equal the 0.9.32 values to rounding.
// v0.9.34 (C1): with the optional W2, every entry is solved and weighted with
// the fixed second-step weight (C_g2, W2) instead of its first-step pair.
real matrix xdpt2_fast_obj_split_list(real colvector F, real colvector R,
                                      real colvector ud, real matrix ETA,
                                      struct xdpt2_gamma_cache rowvector cache,
                                      real colvector idx, | real matrix W2)
{
    real matrix OUT, ZE, GE, G, WG, RF, GF, P
    real rowvector bad
    real colvector ZF
    real scalar j, j0, j1, e, n_ref, m, nb, w2m
    pointer(real matrix) scalar pZ_ref
    w2m = 0
    if (args() == 7) w2m = (rows(W2) > 0)
    m = rows(idx)
    OUT = J(m, cols(ETA), .)
    pZ_ref = NULL
    n_ref = .
    nb = 1
    j0 = 1
    while (j0 <= m) {
        e = idx[j0]
        if (cache[e].pZ != pZ_ref | cache[e].n_rows != n_ref) {
            pZ_ref = cache[e].pZ
            n_ref  = cache[e].n_rows
            ZE     = xdpt2_gsum_by_unit((*pZ_ref) :* R, ud, rows(ETA))' * ETA / n_ref
            ZF     = (*pZ_ref)' * F / n_ref
            nb     = max((1, floor(4e6 / n_ref)))
        }
        // block j0..j1: consecutive entries on this Z, at most nb of them
        j1 = j0
        while (j1 < m & j1 - j0 + 1 < nb) {
            e = idx[j1 + 1]
            if (cache[e].pZ != pZ_ref | cache[e].n_rows != n_ref) break
            j1++
        }
        RF = J(n_ref, j1 - j0 + 1, .)
        for (j = j0; j <= j1; j++) {
            e = idx[j]
            if (w2m) RF[., j - j0 + 1] = F - cache[e].dW * (cache[e].C_g2 * ZF)
            else     RF[., j - j0 + 1] = F - cache[e].dW * (cache[e].C_g * ZF)
        }
        GF = (*pZ_ref)' * RF / n_ref
        for (j = j0; j <= j1; j++) {
            e = idx[j]
            if (w2m) {
                GE = ZE - cache[e].ZW * (cache[e].C_g2 * ZE)
                G  = GF[., j - j0 + 1] :+ GE
                WG = W2 * G
            }
            else {
                GE = ZE - cache[e].ZW * (cache[e].C_g * ZE)
                G  = GF[., j - j0 + 1] :+ GE
                WG = (*cache[e].pW1) * G
            }
            // v0.9.35 (external review, R3): colsum() skips missing values,
            // so a draw whose terms overflowed (or whose weight is missing)
            // had a partial sum -- zero when the other terms were zero --
            // that passed as a valid objective; it is now missing, as in the
            // scalar path (xdpt2_fast_gmm_boot)
            P = G :* WG
            OUT[j, .] = n_ref :* colsum(P)
            if (hasmissing(P)) {
                bad = selectindex(colmissing(P) :> 0)
                OUT[j, bad] = J(1, cols(bad), .)
            }
        }
        j0 = j1 + 1
    }
    return(OUT)
}

// Dense draw indices for the clusters that actually contribute rows. This
// keeps seeded bootstrap results invariant to adding a pruned/ghost panel.
real colvector xdpt2_dense_uid(real colvector uid)
{
    real colvector ord, s, dense, out
    real scalar j
    if (rows(uid) == 0) return(J(0, 1, .))
    if (hasmissing(uid)) {
        errprintf("xtdpthresh internal error: bootstrap cluster ID is missing\n")
        exit(498)
    }
    ord = order(uid, 1)
    s = uid[ord]
    dense = J(rows(uid), 1, 1)
    for (j = 2; j <= rows(uid); j++) {
        dense[j] = dense[j-1] + (s[j] != s[j-1])
    }
    out = J(rows(uid), 1, .)
    out[ord] = dense
    return(out)
}

real scalar xdpt2_is_strongly_balanced(real colvector uid,
                                        real colvector times)
{
    real matrix keys, pinfo
    real colvector tref, ti
    real scalar j
    if (rows(uid) == 0 | rows(times) != rows(uid)) return(0)
    if (hasmissing(uid) | hasmissing(times)) return(0)
    keys = uniqrows((uid, times))
    pinfo = panelsetup(keys, 1)
    if (rows(pinfo) <= 1) return(1)
    tref = keys[|pinfo[1, 1], 2 \ pinfo[1, 2], 2|]
    for (j = 2; j <= rows(pinfo); j++) {
        ti = keys[|pinfo[j, 1], 2 \ pinfo[j, 2], 2|]
        if (rows(ti) != rows(tref)) return(0)
        if (ti != tref) return(0)
    }
    return(1)
}

// Vectorized Mammen 2-point draws for n_u contributing units.
// P[η = -φ]   = (√5+1)/(2√5)   where φ = (√5-1)/2
// P[η = 1/φ]  = (√5-1)/(2√5)   where 1/φ = (√5+1)/2
real colvector xdpt2_mammen_draw(real scalar n_u)
{
    real scalar phi_mam, prob_mam
    real colvector d
    phi_mam = (sqrt(5) - 1) / 2
    prob_mam = (sqrt(5) + 1) / (2 * sqrt(5))
    d = (runiform(n_u, 1) :< prob_mam)
    return(-phi_mam :* d :+ (1/phi_mam) :* (1 :- d))
}

// Deterministic component-specific seeds when rseed() is set. This decouples
// inference objects; it does not claim mathematically nonoverlapping streams.
// Without rseed(), leave the caller's sequential RNG behavior unchanged.
real scalar xdpt2_component_seed(real scalar offset)
{
    real scalar s
    if (st_local("rseed") == "") return(.)
    s = mod(strtoreal(st_local("rseed")) + offset, 2147483648)
    stata("quietly set seed " + strofreal(s, "%21.0f"))
    return(s)
}

// Quantile helper that avoids dependency on moremata's mm_quantile().
// Uses Hyndman-Fan type 7 interpolation and ignores missing values.
real scalar xdpt2_quantile(real colvector x, real scalar p)
{
    real colvector xs
    real scalar n, h, j, g

    xs = select(x, x :< .)
    n = rows(xs)
    if (n == 0) return(.)
    xs = sort(xs, 1)
    if (p <= 0) return(xs[1])
    if (p >= 1) return(xs[n])

    h = 1 + (n - 1) * p
    j = floor(h)
    g = h - j
    if (j >= n) return(xs[n])
    return((1 - g) * xs[j] + g * xs[j + 1])
}

// v0.9.32: critical value of the grid-bootstrap test, an order statistic of
// the valid bootstrap statistics; the type-7 interpolated quantile used up to
// 0.9.31 was liberal by about 0.9/(n+1).
// v0.9.34: the ceil(p*(n+1))-th order statistic (was ceil(p*n), Gong and Seo
// 2026, eq. 7). "Accept iff D <= crit" is then exactly the add-one rule of the
// linearity and continuity tests, reject iff (1 + #{D* >= D})/(1 + n) <=
// 1 - p, whose size is at most 1 - p for every n under exchangeability. The
// old rank was one lower when (1-p)*(n+1) is not an integer, a liberal test
// (size 5.9% at n = 100, 9.5% at n = 20). The two agree when (1-p)*(n+1) is an
// integer (n = 299, 499, 999 at p = .95). If the rank exceeds n, no draw can
// reject and maxdouble() is returned (every candidate is accepted). The
// product is compared with a tolerance so that representation error cannot
// move the rank up by one.
real scalar xdpt2_crit_orderstat(real colvector x, real scalar p)
{
    real colvector xs
    real scalar n, k

    xs = select(x, x :< .)
    n = rows(xs)
    if (n == 0) return(.)
    xs = sort(xs, 1)
    k = ceil(p * (n + 1) - 1e-9)
    if (k < 1) k = 1
    if (k > n) return(maxdouble())
    return(xs[k])
}

// 2-stage grid search (xthenreg-style):
//   Stage 1: grid search with W_first (MA(1) for FD, ZZ_inv for FOD)
//   Stage 2: compute W_n_2 from Stage-1 residuals, grid search again with W_n_2 fixed
// v0.7.0 (A4): the 2-stage fixed-weight path now applies to ALL methods.
// Previously FOD called 2-step solve_gmm per grid point, re-estimating
// Ω(γ) at every γ — objective values were then not comparable across the grid
// and the argmin could be distorted. W_first (MA(1) for FD; Z'Z-inverse for
// FOD) is γ-invariant because Z does not depend on γ, so Stage 1 is a
// proper fixed-weight search; Stage 2 fixes the cluster Ω from the Stage-1
// optimum across the whole grid, exactly as xthenreg does.
// v0.7.0 (D1): the per-γ cache is built once by the caller and passed in.
// v0.7.13 (audit R4, C1): Windmeijer (2005) finite-sample correction for the
// two-step GMM variance, robust variant:
//   V_c = V2 + D V2 + V2 D' + D V1r D'
// where V2 = (G'W2 G)^{-1}/n (efficient-form two-step variance), V1r is the
// robust sandwich of the STAGE-1 estimator, and column j of D is
// [v0.9.35 (external review, R2): that form holds when the stage-1 and
// stage-2 Jacobians are equal. They differ when gamma-hat_1 != gamma-hat_2
// and, in the joint variance, through the derivative columns built at the two
// sets of estimates; its cross terms then used V2 in place of the cross
// covariance L1 Omega1 L2'/n and could give negative variances. The
// correction is now the variance of the linearized two-step estimator,
//   V_c = (L2 + D L1) Omega1 (L2 + D L1)' / n,
//   L1 = (G1'W1 G1)^{-1} G1'W1,  L2 = (G2'W2 G2)^{-1} G2'W2,
// which equals the form above when G1 = G2 (W2 = Omega1^{-1}) and is
// positive semidefinite by construction.]
//   D_j = -B2 G'W2 [dOmega/dtheta_j] W2 gbar2,   B2 = (G'W2 G)^{-1},
// with the Omega derivative evaluated at the stage-1 residuals that BUILT
// W2 (linear model: dOmega/dtheta_j = -(Hj'Gm + Gm'Hj)/n, Hj/Gm the per-unit
// moment rows of Z:*x_j and Z:*r1). Per-n scaling conventions match the
// callers (ZW = Z'X/n, gbar = Z'r/n, Omega = per-unit outer sums / n).
// gamma-wrinkle: with gamma_hat_2 != gamma_hat_1, W2 is built from the
// stage-1 design at idx_1, so Hj/Gm/V1r use idx_1 pieces while G/gbar2 use
// idx_2 pieces; Z, dY, uid are gamma-invariant so the moment space matches.
// Returns J(0,0,.) on numerical failure (caller keeps the uncorrected V).
real matrix xdpt2_windmeijer(real matrix ZW1, real matrix X1, real matrix Z,
                              real colvector uid, real colvector r1,
                              real matrix Omega1, real matrix W1,
                              real matrix W2,
                              real matrix ZW2, real colvector gbar2,
                              real scalar n_rows)
{
    real scalar k, j, n_u, inv_ok
    real matrix B2, B2inv, B1, B1inv, Gm, Hj, dOm, Dmat, Vc, M
    k = cols(ZW2)
    if (cols(X1) != k | cols(ZW1) != k) return(J(0,0,.))
    B2 = ZW2' * W2 * ZW2
    xdpt2_syminv(B2, inv_ok, B2inv)
    if (!inv_ok) return(J(0,0,.))
    B2 = B2inv
    B1 = ZW1' * W1 * ZW1
    xdpt2_syminv(B1, inv_ok, B1inv)
    if (!inv_ok) return(J(0,0,.))
    n_u = max(uid)
    Gm = xdpt2_gsum_by_unit(Z :* r1, uid, n_u)
    // v0.8.2 R11 (#5): the derivative of Omega must match the Omega
    // estimator actually used. Under -center- (the default), Omega
    // subtracts s s'/(n_contrib*n_rows) with s = sum_i g_i, so
    // dOmega/dtheta_j gains +(h_j' s + s' h_j)/(n_contrib*n_rows) with
    // h_j = colsum(Hj) (the derivative of s is -h_j; the two minus signs
    // cancel). Ghost rows from gsum_by_unit are all-zero and do not affect
    // the column sums; n_contrib counts contributing units only, exactly
    // as in xdpt2_build_cluster_omega.
    external real scalar xdpt_center
    real scalar n_contrib_w
    real rowvector gsum_w, hsum_w
    n_contrib_w = rows(uniqrows(uid))
    gsum_w = colsum(Gm)
    Dmat = J(k, k, 0)
    for (j = 1; j <= k; j++) {
        Hj  = xdpt2_gsum_by_unit(Z :* X1[., j], uid, n_u)
        dOm = -(Hj' * Gm + Gm' * Hj) / n_rows
        if (xdpt_center == 1) {
            hsum_w = colsum(Hj)
            dOm = dOm + (hsum_w' * gsum_w + gsum_w' * hsum_w) / (n_contrib_w * n_rows)
        }
        Dmat[., j] = -B2 * (ZW2' * (W2 * (dOm * (W2 * gbar2))))
    }
    // v0.9.35 (R2): the variance of the linearized two-step estimator
    M = B2 * ZW2' * W2 + Dmat * (B1inv * ZW1' * W1)
    Vc = M * Omega1 * M' / n_rows
    Vc = (Vc + Vc') / 2
    if (hasmissing(Vc)) return(J(0,0,.))
    return(Vc)
}

// v0.9.32: joint variance of (theta, gamma) in the kink model. The kink term
// delta_k*(q - gamma)*1(q > gamma) is continuous in gamma, with derivative
// -delta_k*1(q > gamma), so (theta, gamma) is estimated as in a regular GMM
// problem whose Jacobian has one more column, Z'x_g/n, where
// x_g = -delta_k*T(1(q > gamma-hat)); the continuity-restricted estimator is
// asymptotically normal (Gong and Seo 2026, sec. 1 and 3). The slope
// variance is the theta block of the joint variance, computed by the same
// sandwich and Windmeijer code as the linear parameters, applied to the
// design augmented by x_g. T(1(q > gamma)) is column K+1 of the jump design
// at the same gamma: the same transformation, rows, and td partialling as the
// regressors. The jump model is untouched.
//
// x_g aligned with the rows of the kink design (dY_ref, times_ref, uid_ref);
// J(0,1,.) if the jump design does not have the same rows.
real colvector xdpt2_kink_xg(struct xdpt2_unit rowvector units,
                             real scalar gamma, string scalar method,
                             real scalar flag_static, real scalar t_min,
                             real scalar t_max, real scalar delta_k,
                             real colvector dY_ref, real colvector times_ref,
                             real colvector uid_ref)
{
    real colvector dY_j, times_j, uid_j
    real matrix dW_j, Z_j
    real scalar K, drop_save
    real rowvector gone_save
    external real scalar xdpt_ivc_drop
    external real rowvector xdpt_td_gone

    if (delta_k >= . | gamma >= .) return(J(0, 1, .))
    // The stacker resets these run-level diagnostics. The jump design has
    // the same Z and base columns, but the kink run's values are restored.
    drop_save = xdpt_ivc_drop
    gone_save = xdpt_td_gone
    xdpt2_stack_at_gamma(units, gamma, method, flag_static, 0, t_min, t_max,
                          dY_j, dW_j, Z_j, times_j, uid_j)
    xdpt_ivc_drop = drop_save
    xdpt_td_gone = gone_save
    K = cols(units[1].X)
    if (rows(dY_j) != rows(dY_ref) | cols(dW_j) < K + 1) return(J(0, 1, .))
    if (rows(dY_j) == 0) return(J(0, 1, .))
    if (any(times_j :!= times_ref) | any(uid_j :!= uid_ref)) return(J(0, 1, .))
    if (any(dY_j :!= dY_ref)) return(J(0, 1, .))
    return(-delta_k :* dW_j[., K + 1])
}

// v0.9.35: the transformation of a template applied to any level columns:
// FD or FOD with the arithmetic of xdpt2_tpl_dW, then the td partialling.
real matrix xdpt2_tpl_trans(struct xdpt2_stack_tpl scalar tp, real matrix Wlev)
{
    real colvector ix
    real matrix dreg, F, WE, X, S
    real scalar i, j, v
    if (tp.fd) {
        dreg = Wlev[tp.J, .] - Wlev[tp.JP, .]
    }
    else {
        WE = Wlev[tp.EQF, .]
        F = J(rows(tp.EQF), cols(Wlev), .)
        for (i = 1; i <= rows(tp.FG); i++) {
            v  = tp.FG[i, 1]
            ix = tp.FP[|tp.FG[i, 2] \ tp.FG[i, 3]|]
            S  = lowertriangle(J(v, v, 1), 0)
            for (j = 1; j <= cols(Wlev); j++) {
                X = colshape(WE[ix, j], v)'
                F[ix, j] = vec(quadcross(S, X))
            }
        }
        dreg = tp.c :* (Wlev[tp.J, .] - F[tp.iF, .] :/ tp.Tf)
    }
    if (tp.td == 1) {
        xdpt2_demean_bytime(dreg, tp.tpre)
        dreg = dreg[tp.keep, .]
    }
    else if (tp.td == 2) {
        dreg = dreg - tp.D * (tp.DtD * (tp.D' * dreg))
        dreg = dreg[tp.keep, .]
    }
    return(dreg)
}

// v0.9.35: joint variance of (theta, gamma) in the jump model. The sample
// moments are step functions of gamma, but their expectation is smooth:
// d/dgamma E[z 1(q > gamma) w] = -E[z w f(gamma | .)]. As in Seo and Shin
// (2016, sec. 3) and xthenreg (Seo, Kim, and Kim 2019), the derivative column
// is x_g = -T(phi_h(q - gamma-hat) * (1, x')delta-hat), phi_h(u) = phi(u/h)/h
// the Gaussian kernel, h = bwscale * 1.06 * s_q * n^(-1/5) with s_q the
// standard deviation of q over the panel rows and n the number of units (the
// form of xthenreg's bandwidth rule, with its default multiplier 1.5).
// T is applied by a template (the transformation,
// rows and td partialling of the regressors): the main build's template when
// it was verified and has the rows of the reference design, else one recorded
// at gamma-hat. It must reproduce the regime columns of the reference design
// dW_ref; otherwise, or if the rows differ, J(0,1,.) (the conditional
// variance is then kept). The slope block is invariant to the scale of x_g.
// the column from a template tp whose rows are those of the reference design
real colvector xdpt2_jump_xg_tp(struct xdpt2_stack_tpl scalar tp,
                                real scalar gamma, real colvector theta,
                                real scalar K, real matrix dW_ref,
                                real colvector uid_ref)
{
    real colvector qv, w, delta
    real matrix Wreg, dchk, dreg
    real scalar h, nu, scale
    external real scalar xdpt_bwscale, xdpt_gamma_bw

    // the template reproduces the regime columns of the design at gamma-hat
    if (cols(dW_ref) != 2 * K + 1) return(J(0, 1, .))
    Wreg = (tp.qlev :> gamma)
    Wreg = Wreg, tp.Xlev :* Wreg
    dchk = xdpt2_tpl_trans(tp, Wreg)
    dreg = dW_ref[|1, K + 1 \ rows(dW_ref), 2 * K + 1|]
    if (rows(dchk) != rows(dreg) | cols(dchk) != cols(dreg)) return(J(0, 1, .))
    scale = max((1, max(abs(dreg))))
    if (max(abs(dchk - dreg)) > 1e-9 * scale) return(J(0, 1, .))
    // bandwidth rule of xthenreg, with multiplier bwscale (its default h_0 is 1.5)
    qv = select(tp.qlev, tp.qlev :< .)
    nu = rows(uniqrows(uid_ref))
    if (rows(qv) < 2 | nu < 1) return(J(0, 1, .))
    h = xdpt_bwscale * 1.06 * sqrt(variance(qv)) * nu^(-0.2)
    if (h >= . | h <= 0) return(J(0, 1, .))
    xdpt_gamma_bw = h
    delta = theta[|K + 1 \ 2 * K + 1|]
    w = normalden((tp.qlev :- gamma) :/ h) :/ h :*
        ((J(rows(tp.qlev), 1, 1), tp.Xlev) * delta)
    return(-xdpt2_tpl_trans(tp, w))
}

real colvector xdpt2_jump_xg(struct xdpt2_unit rowvector units,
                             real scalar gamma, real colvector theta,
                             string scalar method, real scalar flag_static,
                             real scalar t_min, real scalar t_max,
                             real colvector dY_ref, real matrix dW_ref,
                             real colvector times_ref, real colvector uid_ref,
                             struct xdpt2_stack_tpl scalar tpl_in,
                             real scalar tpl_in_st)
{
    struct xdpt2_stack_tpl scalar tp
    real colvector dY_j, times_j, uid_j
    real matrix dW_j, Z_j
    real scalar K, drop_save, dep_save
    real rowvector gone_save
    external real scalar xdpt_ivc_drop, xdpt_ivc_dep, xdpt_tpl_rec
    external real rowvector xdpt_td_gone

    K = cols(units[1].X)
    if (gamma >= . | rows(theta) != 2 * K + 1) return(J(0, 1, .))
    if (hasmissing(theta) | rows(dY_ref) == 0) return(J(0, 1, .))
    // the main build's template: gamma-invariant and checked bit for bit
    if (tpl_in_st == 2) {
        if (rows(tpl_in.dY) == rows(dY_ref)) {
            if (tpl_in.dY == dY_ref & tpl_in.times == times_ref &
                tpl_in.uid == uid_ref) {
                return(xdpt2_jump_xg_tp(tpl_in, gamma, theta, K, dW_ref,
                                        uid_ref))
            }
        }
    }
    // otherwise a template recorded at gamma-hat; the stacker resets these
    // run-level diagnostics, which are restored
    drop_save = xdpt_ivc_drop
    dep_save = xdpt_ivc_dep
    gone_save = xdpt_td_gone
    xdpt_tpl_rec = 1
    xdpt2_stack_at_gamma(units, gamma, method, flag_static, 0, t_min, t_max,
                          dY_j, dW_j, Z_j, times_j, uid_j)
    xdpt_tpl_rec = 0
    tp = xdpt2_tpl_build(units, method, 0, dY_j, Z_j, times_j, uid_j)
    xdpt_ivc_drop = drop_save
    xdpt_ivc_dep = dep_save
    xdpt_td_gone = gone_save
    if (!tp.ok) return(J(0, 1, .))
    if (rows(dY_j) != rows(dY_ref) | cols(dW_j) != 2 * K + 1) return(J(0, 1, .))
    if (any(times_j :!= times_ref) | any(uid_j :!= uid_ref)) return(J(0, 1, .))
    if (any(dY_j :!= dY_ref)) return(J(0, 1, .))
    return(xdpt2_jump_xg_tp(tp, gamma, theta, K, dW_j, uid_ref))
}

// v0.9.35: number of distinct values of q (the level rows of all units, as
// for the bandwidth) within two bandwidths of gamma. Few values mean a
// discrete q, for which the kernel derivative of the jump model's joint
// variance has no density to estimate (Seo and Shin 2016, Assumption 2).
real scalar xdpt2_q_nvals(struct xdpt2_unit rowvector units, real scalar gamma,
                          real scalar h)
{
    real colvector qv
    real scalar i, n, a

    if (gamma >= . | h >= . | h <= 0) return(.)
    n = 0
    for (i = 1; i <= length(units); i++) n = n + rows(units[i].q)
    if (n == 0) return(.)
    qv = J(n, 1, .)
    a = 0
    for (i = 1; i <= length(units); i++) {
        if (rows(units[i].q) == 0) continue
        qv[|a + 1 \ a + rows(units[i].q)|] = units[i].q
        a = a + rows(units[i].q)
    }
    qv = select(qv, (qv :< .) :& (abs(qv :- gamma) :<= 2 * h))
    if (rows(qv) == 0) return(0)
    return(rows(uniqrows(qv)))
}

// v0.9.35: derivative in gamma of the residuals of a stack built with the
// template tp (the estimation stack, or under FOD the FD stack of the AR
// test): the column of the joint variance, x_g = -T(phi_h(q - gamma)
// (1, x')delta) in the jump model (bandwidth h, as used there) and
// x_g = -delta_k T(1(q > gamma)) under kink. The template must reproduce the
// regime columns of dW_ref; otherwise J(0,1,.).
real colvector xdpt2_ar_xg(struct xdpt2_stack_tpl scalar tp, real scalar gamma,
                           real colvector theta, real scalar K,
                           real scalar flag_kink, real matrix dW_ref,
                           real scalar h)
{
    real colvector ind, w, delta, xg
    real matrix Wreg, dchk, dreg
    real scalar scale

    if (!tp.ok | gamma >= . | hasmissing(theta)) return(J(0, 1, .))
    ind = (tp.qlev :> gamma)
    if (flag_kink) {
        if (cols(dW_ref) != K + 1 | rows(theta) != K + 1) return(J(0, 1, .))
        Wreg = (tp.qlev :- gamma) :* ind
        dreg = dW_ref[., K + 1]
    }
    else {
        if (cols(dW_ref) != 2 * K + 1 | rows(theta) != 2 * K + 1) return(J(0, 1, .))
        Wreg = ind, tp.Xlev :* ind
        dreg = dW_ref[|1, K + 1 \ rows(dW_ref), 2 * K + 1|]
    }
    dchk = xdpt2_tpl_trans(tp, Wreg)
    if (rows(dchk) != rows(dreg) | cols(dchk) != cols(dreg)) return(J(0, 1, .))
    scale = max((1, max(abs(dreg))))
    if (max(abs(dchk - dreg)) > 1e-9 * scale) return(J(0, 1, .))
    if (flag_kink) {
        xg = -theta[K + 1] :* xdpt2_tpl_trans(tp, ind)
    }
    else {
        if (h >= . | h <= 0) return(J(0, 1, .))
        delta = theta[|K + 1 \ 2 * K + 1|]
        w = normalden((tp.qlev :- gamma) :/ h) :/ h :*
            ((J(rows(tp.qlev), 1, 1), tp.Xlev) * delta)
        xg = -xdpt2_tpl_trans(tp, w)
    }
    if (hasmissing(xg)) return(J(0, 1, .))
    return(xg)
}

// Theta block of the joint sandwich variance (estimator with weight A).
// v0.9.35: the full matrix, gamma included, is kept in xdpt_V_joint_full for
// the AR statistics.
real matrix xdpt2_kink_joint_V(real matrix ZW, real colvector xg,
                                real matrix Z, real matrix A,
                                real matrix Omega, real scalar n_rows)
{
    real matrix Va
    real scalar k
    external real matrix xdpt_V_joint_full

    xdpt_V_joint_full = J(0, 0, .)
    if (rows(xg) == 0 | rows(xg) != rows(Z)) return(J(0, 0, .))
    k = cols(ZW)
    Va = xdpt2_gmm_sandwich((ZW, (Z' * xg) / n_rows), A, Omega, n_rows)
    if (rows(Va) == 0) return(J(0, 0, .))
    xdpt_V_joint_full = Va
    return(Va[|1, 1 \ k, k|])
}

// Theta block of the joint Windmeijer-corrected variance: xdpt2_windmeijer
// applied to the stage-1 and stage-2 designs augmented by their x_g.
real matrix xdpt2_kink_joint_wind(real matrix ZW1, real matrix X1,
                                   real colvector xg1, real matrix Z,
                                   real colvector uid, real colvector r1,
                                   real matrix Omega1, real matrix W1,
                                   real matrix W2, real matrix ZW2,
                                   real colvector xg2, real colvector gbar2,
                                   real scalar n_rows)
{
    real matrix Va
    real scalar k
    external real matrix xdpt_V_joint_full

    xdpt_V_joint_full = J(0, 0, .)
    if (rows(xg1) == 0 | rows(xg2) == 0) return(J(0, 0, .))
    if (rows(xg1) != rows(Z) | rows(xg2) != rows(Z)) return(J(0, 0, .))
    k = cols(ZW2)
    Va = xdpt2_windmeijer((ZW1, (Z' * xg1) / n_rows), (X1, xg1), Z, uid, r1,
                          Omega1, W1, W2, (ZW2, (Z' * xg2) / n_rows), gbar2,
                          n_rows)
    if (rows(Va) == 0) return(J(0, 0, .))
    xdpt_V_joint_full = Va
    return(Va[|1, 1 \ k, k|])
}

// Replace the conditional variance V by the joint theta block V_joint when it
// is available; the conditional V is kept (e(V_cond)) and gives the AR
// statistics that treat gamma-hat as known.
void xdpt2_kink_apply_joint(real matrix V, real matrix V_joint)
{
    external real matrix xdpt_V_cond_ar, xdpt_V_joint_full
    external real scalar xdpt_kink_joint

    xdpt_V_cond_ar = V
    if (rows(V_joint) == rows(V) & cols(V_joint) == cols(V) &
        !hasmissing(V_joint)) {
        V = V_joint
        xdpt_kink_joint = 1
    }
    else {
        xdpt_kink_joint = 0
        xdpt_V_joint_full = J(0, 0, .)
    }
}

// v0.9.25: nested-grid search helpers. The search engine below is the only
// point-estimation engine for searchmode(adaptive) and searchmode(fixed).
// Midpoint insertion gives g -> 2g-1 -> 4g-3 (100 -> 199 -> 397).
real colvector xdpt2_nested_midpoints(real colvector gamma_grid)
{
    real colvector g, left, right, mid, keep
    g = sort(uniqrows(gamma_grid), 1)
    if (rows(g) < 2) return(J(0, 1, .))
    left  = g[|1 \ rows(g)-1|]
    right = g[|2 \ rows(g)|]
    mid = left :+ (right :- left) :/ 2
    keep = (mid :> left) :& (mid :< right)
    if (sum(keep) == 0) return(J(0, 1, .))
    return(select(mid, keep))
}


// Append only never-before-seen candidates and build only their cache
// entries. Exact guards then rewire gamma-invariant Z/W1 pointers to an
// existing entry, avoiding retained duplicate copies of those large matrices.
// This is shared by global midpoint escalation and final local refinement.
void xdpt2_append_candidates(
    struct xdpt2_unit rowvector units,
    real colvector candidates,
    string scalar method,
    real scalar flag_static,
    real scalar flag_kink,
    real scalar t_min,
    real scalar t_max,
    real colvector q_supp,
    real scalar min_user,
    real colvector gamma_grid,
    struct xdpt2_gamma_cache rowvector cache,
    real scalar first_new,
    real scalar last_new,
    struct xdpt2_stack_tpl scalar tpl,
    real scalar tpl_st)
{
    real colvector new_g, keep
    struct xdpt2_gamma_cache rowvector cache_new
    real scalar j, ref, same_tu, same_z, same_y

    first_new = 0
    last_new = 0
    new_g = candidates
    if (rows(new_g) > 0) new_g = select(new_g, new_g :< .)
    if (rows(new_g) == 0) return
    new_g = uniqrows(sort(new_g, 1))
    keep = J(rows(new_g), 1, 1)
    for (j = 1; j <= rows(new_g); j++) {
        if (sum(gamma_grid :== new_g[j]) > 0) keep[j] = 0
    }
    if (sum(keep) == 0) return
    new_g = select(new_g, keep)
    if (rows(new_g) == 0) return

    cache_new = xdpt2_build_gamma_cache_t(units, new_g, method,
                                          flag_static, flag_kink,
                                          t_min, t_max, q_supp, min_user,
                                          tpl, tpl_st, cache)

    ref = 0
    for (j = 1; j <= cols(cache); j++) {
        if (cache[j].ok) {
            ref = j
            break
        }
    }
    if (ref > 0) {
        for (j = 1; j <= cols(cache_new); j++) {
            if (!cache_new[j].ok) continue
            // already shared by the build (v0.9.34)
            if (cache_new[j].pZ == cache[ref].pZ &
                cache_new[j].pW1 == cache[ref].pW1) continue
            same_tu = 0
            if (rows(cache_new[j].times) == rows(cache[ref].times) &
                rows(cache_new[j].uid)   == rows(cache[ref].uid)) {
                same_tu = (cache_new[j].times == cache[ref].times &
                           cache_new[j].uid   == cache[ref].uid)
            }
            if (!same_tu) continue
            same_z = 0
            if (rows(*cache_new[j].pZ) == rows(*cache[ref].pZ) &
                cols(*cache_new[j].pZ) == cols(*cache[ref].pZ)) {
                same_z = (*cache_new[j].pZ == *cache[ref].pZ)
            }
            if (!same_z) continue
            cache_new[j].pZ  = cache[ref].pZ
            cache_new[j].pW1 = cache[ref].pW1
            same_y = 0
            if (rows(cache_new[j].dY) == rows(cache[ref].dY)) {
                same_y = (cache_new[j].dY == cache[ref].dY)
            }
            if (same_y) cache_new[j].ZY = cache[ref].ZY
        }
    }

    first_new = rows(gamma_grid) + 1
    gamma_grid = gamma_grid \ new_g
    cache = cache, cache_new
    last_new = rows(gamma_grid)
}


// Append the next globally nested level: g -> 2g-1 -> 4g-3.
void xdpt2_append_nested_level(
    struct xdpt2_unit rowvector units,
    string scalar method,
    real scalar flag_static,
    real scalar flag_kink,
    real scalar t_min,
    real scalar t_max,
    real colvector q_supp,
    real scalar min_user,
    real colvector gamma_grid,
    struct xdpt2_gamma_cache rowvector cache,
    real scalar first_new,
    real scalar last_new,
    struct xdpt2_stack_tpl scalar tpl,
    real scalar tpl_st)
{
    real colvector new_g
    new_g = xdpt2_nested_midpoints(gamma_grid)
    xdpt2_append_candidates(units, new_g, method, flag_static, flag_kink,
                             t_min, t_max, q_supp, min_user,
                             gamma_grid, cache, first_new, last_new,
                             tpl, tpl_st)
}


// Evaluate one never-before-profiled cache slice and update a running argmin.
// Stage 1 uses each entry's shared W1; stage 2 uses one common fixed W2.
void xdpt2_profile_update(
    real colvector gamma_grid,
    struct xdpt2_gamma_cache rowvector cache,
    real scalar first,
    real scalar last,
    real scalar fixed_weight,
    real matrix A_fixed,
    real scalar best_idx,
    real scalar best_obj,
    real scalar best_gamma,
    real colvector best_theta,
    real matrix best_V,
    real scalar obj_hi,
    real scalar n_adm,
    real colvector gamma_admitted)
{
    real scalar gl, lo, hi, ok, obj_cur, tol
    real colvector theta_cur
    real matrix V_cur

    lo = max((1, first))
    hi = min((rows(gamma_grid), last))
    if (hi < lo) return

    for (gl = lo; gl <= hi; gl++) {
        if (!cache[gl].ok) continue
        if (rows(cache[gl].dY) < 20) continue
        if (fixed_weight) {
            if (rows(A_fixed) != cols(*cache[gl].pZ) |
                cols(A_fixed) != cols(*cache[gl].pZ)) continue
            xdpt2_solve_gmm_1step_pre(cache[gl].dY, cache[gl].dW,
                                       *cache[gl].pZ, cache[gl].ZW,
                                       cache[gl].ZY, A_fixed,
                                       ok, theta_cur, obj_cur, V_cur)
        }
        else {
            xdpt2_solve_gmm_1step_pre(cache[gl].dY, cache[gl].dW,
                                       *cache[gl].pZ, cache[gl].ZW,
                                       cache[gl].ZY, *cache[gl].pW1,
                                       ok, theta_cur, obj_cur, V_cur)
        }
        if (!ok) continue

        n_adm = n_adm + 1
        gamma_admitted = gamma_admitted \ gamma_grid[gl]
        if (obj_hi >= . | obj_cur > obj_hi) obj_hi = obj_cur

        if (best_idx == 0) {
            best_idx = gl
            best_obj = obj_cur
            best_gamma = gamma_grid[gl]
            best_theta = theta_cur
            best_V = V_cur
        }
        else {
            tol = xdpt2_objtol(obj_cur, best_obj, 1e-12)
            if (obj_cur < best_obj - tol |
                (abs(obj_cur - best_obj) <= tol &
                 gamma_grid[gl] < best_gamma)) {
                best_idx = gl
                best_obj = obj_cur
                best_gamma = gamma_grid[gl]
                best_theta = theta_cur
                best_V = V_cur
            }
        }
    }
}


real scalar xdpt2_profile_is_flat(real scalar best_obj,
                                   real scalar obj_hi,
                                   real scalar n_adm)
{
    real scalar scale
    if (best_obj >= . | obj_hi >= . | n_adm < 2) return(1)
    scale = max((abs(best_obj), abs(obj_hi)))
    if (scale == 0) return(1)
    return(abs(obj_hi - best_obj) <= 1e-12 * scale)
}


// Stop only when the realized regime split is unchanged and the relative
// objective improvement is negligible.
void xdpt2_profile_stability(
    real scalar previous_idx,
    real scalar current_idx,
    real scalar previous_obj,
    real scalar current_obj,
    real colvector gamma_grid,
    real colvector q_split_supp,
    real scalar search_tol,
    real scalar same_split,
    real scalar rel_gain,
    real scalar converged)
{
    real scalar scale, nleft_old, nleft_new
    same_split = 0
    rel_gain = .
    converged = 0
    if (previous_idx <= 0 | current_idx <= 0 |
        previous_obj >= . | current_obj >= .) return

    if (rows(q_split_supp) > 0) {
        nleft_old = sum(q_split_supp :<= gamma_grid[previous_idx])
        nleft_new = sum(q_split_supp :<= gamma_grid[current_idx])
        same_split = (nleft_old == nleft_new)
    }
    else same_split = (gamma_grid[previous_idx] == gamma_grid[current_idx])

    scale = max((abs(previous_obj), abs(current_obj)))
    if (scale == 0) rel_gain = 0
    else {
        rel_gain = (previous_obj - current_obj) / scale
        if (rel_gain < 0) rel_gain = 0
    }
    converged = (same_split & rel_gain <= search_tol)
}


// v0.9.34 (C3): kink only. The kink regressor (q - gamma)1(q > gamma), and
// with it the criterion, is continuous in gamma, and the joint variance of
// (theta, gamma) presumes the exact minimizer; the grid argmin can be half a
// grid step away from it. The criterion of the reported estimator (the fixed
// W2, or each entry's W1 when fixed_weight = 0) is minimized between the
// admitted neighbours of the grid argmin: 40 equally spaced points, then
// rounds of 10 spanning the previous spacing on either side of the best point
// so far, until the spacing is below 1e-6 of the bracket. Only the best point
// is added to the estimation grid, and only if it improves on the grid argmin
// (the rule of xdpt2_profile_update); the admitted-grid counts are unchanged.
void xdpt2_kink_refine(struct xdpt2_unit rowvector units,
                       string scalar method,
                       real scalar flag_static,
                       real scalar t_min,
                       real scalar t_max,
                       real colvector q_supp,
                       real scalar min_user,
                       real colvector anchors,
                       real scalar fixed_weight,
                       real matrix A_fixed,
                       real colvector gamma_grid,
                       struct xdpt2_gamma_cache rowvector cache,
                       real scalar best_idx,
                       real scalar best_obj,
                       real scalar best_gamma,
                       real colvector best_theta,
                       real matrix best_V,
                       struct xdpt2_stack_tpl scalar tpl,
                       real scalar tpl_state)
{
    external real scalar xdpt_kref
    real colvector a, tg, t_theta, t_adm
    real scalar j, pos, lo, hi, L, U, sp, m, r, tol, tl, have
    real scalar r_obj, r_gam, t_idx, t_obj, t_gam, t_hi, t_n
    real scalar d_hi, d_n
    real colvector d_adm
    real matrix t_V
    struct xdpt2_gamma_cache rowvector tc
    struct xdpt2_gamma_cache scalar keep

    xdpt_kref = 0
    if (best_idx == 0 | best_gamma >= . | rows(anchors) == 0) return
    a = select(anchors, anchors :< .)
    if (rows(a) == 0) return
    a = uniqrows(a)
    pos = 0
    for (j = 1; j <= rows(a); j++) {
        if (a[j] == best_gamma) pos = j
    }
    if (pos == 0) return
    lo = (pos > 1 ? a[pos - 1] : min(gamma_grid))
    hi = (pos < rows(a) ? a[pos + 1] : max(gamma_grid))
    if (!(hi > lo)) return

    tol = 1e-6 * (hi - lo)
    r_obj = best_obj
    r_gam = best_gamma
    have = 0
    L = lo
    U = hi
    m = 40
    for (r = 1; r <= 8; r++) {
        sp = (U - L) / (m + 1)
        tg = L :+ sp :* (1::m)
        tc = xdpt2_build_gamma_cache_t(units, tg, method, flag_static, 1,
                                        t_min, t_max, q_supp, min_user,
                                        tpl, tpl_state, cache)
        t_idx = 0
        t_obj = .
        t_gam = .
        t_theta = J(0, 1, .)
        t_V = J(0, 0, .)
        t_hi = .
        t_n = 0
        t_adm = J(0, 1, .)
        xdpt2_profile_update(tg, tc, 1, rows(tg), fixed_weight, A_fixed,
                              t_idx, t_obj, t_gam, t_theta, t_V, t_hi, t_n,
                              t_adm)
        if (t_idx > 0) {
            tl = xdpt2_objtol(t_obj, r_obj, 1e-12)
            if (t_obj < r_obj - tl |
                (abs(t_obj - r_obj) <= tl & t_gam < r_gam)) {
                r_obj = t_obj
                r_gam = t_gam
                keep = tc[t_idx]
                have = 1
            }
        }
        if (sp <= tol) break
        L = max((lo, r_gam - sp))
        U = min((hi, r_gam + sp))
        m = 10
    }
    if (!have) return

    // Share Z and W1 with the grid argmin's entry when bitwise equal, as
    // xdpt2_append_candidates does.
    if (keep.n_rows == cache[best_idx].n_rows) {
        if (keep.times == cache[best_idx].times &
            keep.uid == cache[best_idx].uid &
            rows(*keep.pZ) == rows(*cache[best_idx].pZ) &
            cols(*keep.pZ) == cols(*cache[best_idx].pZ)) {
            if (*keep.pZ == *cache[best_idx].pZ) {
                keep.pZ  = cache[best_idx].pZ
                keep.pW1 = cache[best_idx].pW1
            }
        }
    }
    gamma_grid = gamma_grid \ r_gam
    cache = cache, keep
    j = rows(gamma_grid)
    d_hi = .
    d_n = 0
    d_adm = J(0, 1, .)
    xdpt2_profile_update(gamma_grid, cache, j, j, fixed_weight, A_fixed,
                          best_idx, best_obj, best_gamma, best_theta, best_V,
                          d_hi, d_n, d_adm)
    xdpt_kref = (best_idx == j)
}


// Unified point-estimation search. Adaptive mode profiles nested levels under
// one W1, constructs W2 once, then profiles nested levels under that fixed W2.
// Fixed mode runs one level through the same engine (used by calibration B6).
void xdpt2_grid_search(
    struct xdpt2_unit rowvector units,
    real colvector gamma_grid,
    struct xdpt2_gamma_cache rowvector cache,
    string scalar method,
    real scalar flag_static,
    real scalar flag_kink,
    real scalar t_min,
    real scalar t_max,
    real colvector q_supp,
    real colvector q_split_supp,
    real scalar min_user,
    string scalar search_mode,
    real scalar search_tol,
    real scalar search_max_level,
    real scalar n_refine,
    real scalar best_gamma,
    real scalar best_obj,
    real colvector best_theta,
    real matrix best_V,
    real matrix best_V_influence,
    real matrix best_A,
    real scalar best_twostep,
    real scalar n_adm2,
    real scalar gamma_adm2_lo,
    real scalar gamma_adm2_hi,
    real colvector gamma_admitted,
    real scalar level1_points,
    real scalar level2_points,
    real scalar level3_points,
    real scalar stage1_level,
    real scalar stage2_level,
    real scalar stage1_points,
    real scalar stage2_points,
    real scalar stage1_same_split,
    real scalar stage2_same_split,
    real scalar stage1_rel_gain,
    real scalar stage2_rel_gain,
    real scalar stage1_converged,
    real scalar stage2_converged,
    real scalar search_hit_max,
    real scalar W2_builds,
    real scalar stage1_gamma,
    real scalar stage1_obj,
    real scalar stage2_global_gamma,
    real scalar stage2_global_obj,
    real scalar ref_it,
    real scalar ref_added,
    real scalar ref_pool_n,
    real scalar ref_remaining,
    real scalar ref_exhausted,
    real scalar ref_alo,
    real scalar ref_ahi,
    real scalar ref_in_basin,
    real scalar ref_neigh_rem,
    real scalar ref_complete,
    real scalar ref_obj_gain,
    struct xdpt2_stack_tpl scalar tpl_main,
    real scalar tpl_main_st)
{
    real scalar gl, k_W, first_new, last_new, n1, n2, n3
    real scalar best_idx_1, best_obj_1, best_gamma_1, obj_hi1, n_prof1
    real scalar best_idx_2, best_obj_2, best_gamma_2, obj_hi2, n_prof2
    real scalar prev_idx, prev_obj, inv_ok
    real colvector best_theta_1, best_theta_2, r_1
    real matrix best_V_1, best_V_2, Omega, W_n_2
    real colvector gamma_adm2

    best_gamma = .
    best_obj = .
    best_A = J(0, 0, .)
    best_twostep = 0
    n_adm2 = .
    gamma_adm2_lo = .
    gamma_adm2_hi = .
    gamma_admitted = J(0, 1, .)

    level1_points = rows(gamma_grid)
    level2_points = .
    level3_points = .
    stage1_level = .
    stage2_level = .
    stage1_points = .
    stage2_points = .
    stage1_same_split = .
    stage2_same_split = .
    stage1_rel_gain = .
    stage2_rel_gain = .
    stage1_converged = .
    stage2_converged = .
    search_hit_max = 0
    W2_builds = 0
    stage1_gamma = .
    stage1_obj = .
    stage2_global_gamma = .
    stage2_global_obj = .
    ref_it = 0
    ref_added = 0
    ref_pool_n = 0
    ref_remaining = 0
    ref_exhausted = (n_refine > 0 ? 0 : 1)
    ref_alo = .
    ref_ahi = .
    ref_in_basin = .
    ref_neigh_rem = (n_refine > 0 ? . : 0)
    ref_complete = (n_refine > 0 ? . : 1)
    ref_obj_gain = (n_refine > 0 ? . : 0)

    k_W = .
    for (gl = 1; gl <= cols(cache); gl++) {
        if (cache[gl].ok) {
            k_W = cols(cache[gl].dW)
            break
        }
    }
    if (k_W == .) {
        k_W = (flag_kink ? cols(units[1].X) + 1 :
                            2 * cols(units[1].X) + 1)
    }
    best_theta = J(k_W, 1, 0)
    best_V = J(k_W, k_W, 0)
    best_V_influence = best_V

    n1 = rows(gamma_grid)
    n2 = n1
    n3 = .
    first_new = 0
    last_new = 0
    best_idx_1 = 0
    best_obj_1 = .
    best_gamma_1 = .
    best_theta_1 = J(k_W, 1, 0)
    best_V_1 = J(k_W, k_W, 0)
    obj_hi1 = .
    n_prof1 = 0

    // Stage 1: fixed W1.
    xdpt2_profile_update(gamma_grid, cache, 1, n1, 0, J(0,0,.),
                          best_idx_1, best_obj_1, best_gamma_1,
                          best_theta_1, best_V_1, obj_hi1, n_prof1,
                          gamma_admitted)

    if (search_mode == "adaptive") {
        prev_idx = best_idx_1
        prev_obj = best_obj_1
        xdpt2_append_nested_level(units, method, flag_static, flag_kink,
                                   t_min, t_max, q_supp, min_user,
                                   gamma_grid, cache, first_new, last_new,
                                   tpl_main, tpl_main_st)
        if (first_new > 0) {
            n2 = last_new
            xdpt2_profile_update(gamma_grid, cache, first_new, last_new,
                                  0, J(0,0,.), best_idx_1, best_obj_1,
                                  best_gamma_1, best_theta_1, best_V_1,
                                  obj_hi1, n_prof1, gamma_admitted)
        }
        level2_points = n2
        stage1_level = 2
        stage1_points = n2
        xdpt2_profile_stability(prev_idx, best_idx_1, prev_obj, best_obj_1,
                                 gamma_grid, q_split_supp, search_tol,
                                 stage1_same_split, stage1_rel_gain,
                                 stage1_converged)

        // This flag means CAP EXHAUSTED WITHOUT STABILITY, not merely that
        // the current level happens to equal the configured cap.
        if (search_max_level == 2 & !stage1_converged) search_hit_max = 1

        if (!stage1_converged & search_max_level >= 3) {
            prev_idx = best_idx_1
            prev_obj = best_obj_1
            xdpt2_append_nested_level(units, method, flag_static, flag_kink,
                                       t_min, t_max, q_supp, min_user,
                                       gamma_grid, cache, first_new, last_new,
                                       tpl_main, tpl_main_st)
            if (first_new > 0) {
                n3 = last_new
                xdpt2_profile_update(gamma_grid, cache, first_new, last_new,
                                      0, J(0,0,.), best_idx_1, best_obj_1,
                                      best_gamma_1, best_theta_1, best_V_1,
                                      obj_hi1, n_prof1, gamma_admitted)
            }
            else n3 = n2
            level3_points = n3
            stage1_level = 3
            stage1_points = n3
            xdpt2_profile_stability(prev_idx, best_idx_1, prev_obj,
                                     best_obj_1, gamma_grid, q_split_supp,
                                     search_tol, stage1_same_split,
                                     stage1_rel_gain, stage1_converged)
            if (!stage1_converged) search_hit_max = 1
        }
    }
    else {
        stage1_level = 1
        stage1_points = n1
    }

    if (best_idx_1 == 0 |
        xdpt2_profile_is_flat(best_obj_1, obj_hi1, n_prof1)) return

    stage1_gamma = best_gamma_1
    stage1_obj = best_obj_1

    // W2 is constructed exactly once after stage 1 has stopped.
    r_1 = cache[best_idx_1].dY -
          cache[best_idx_1].dW * best_theta_1
    Omega = xdpt2_build_cluster_omega(*cache[best_idx_1].pZ, r_1,
                                       cache[best_idx_1].uid)
    W2_builds = 1
    xdpt2_syminv(Omega, inv_ok, W_n_2)
    if (!inv_ok) {
        // v0.9.34 (C3): kink -- the one-step criterion is minimized between
        // the grid neighbours; Omega is then taken at the reported estimate.
        if (flag_kink) {
            xdpt2_kink_refine(units, method, flag_static, t_min, t_max,
                               q_supp, min_user, gamma_admitted, 0,
                               J(0, 0, .), gamma_grid, cache, best_idx_1,
                               best_obj_1, best_gamma_1, best_theta_1,
                               best_V_1, tpl_main, tpl_main_st)
            r_1 = cache[best_idx_1].dY -
                  cache[best_idx_1].dW * best_theta_1
            Omega = xdpt2_build_cluster_omega(*cache[best_idx_1].pZ, r_1,
                                               cache[best_idx_1].uid)
        }
        best_gamma = best_gamma_1
        best_obj = best_obj_1
        best_theta = best_theta_1
        best_A = *cache[best_idx_1].pW1
        best_V = xdpt2_gmm_sandwich(cache[best_idx_1].ZW, best_A, Omega,
                                     cache[best_idx_1].n_rows)
        if (rows(best_V) == 0) {
            errprintf("xtdpthresh: cluster-robust variance failed on the one-step fallback\n")
            exit(498)
        }
        best_V_influence = best_V
        // v0.9.32: kink model, joint (theta, gamma) variance; v0.9.35: the
        // jump model too, with the kernel derivative column.
        if (flag_kink) {
            xdpt2_kink_apply_joint(best_V,
                xdpt2_kink_joint_V(cache[best_idx_1].ZW,
                    xdpt2_kink_xg(units, best_gamma_1, method, flag_static,
                                  t_min, t_max,
                                  best_theta_1[rows(best_theta_1)],
                                  cache[best_idx_1].dY,
                                  cache[best_idx_1].times,
                                  cache[best_idx_1].uid),
                    *cache[best_idx_1].pZ, best_A, Omega,
                    cache[best_idx_1].n_rows))
        }
        else {
            xdpt2_kink_apply_joint(best_V,
                xdpt2_kink_joint_V(cache[best_idx_1].ZW,
                    xdpt2_jump_xg(units, best_gamma_1, best_theta_1, method,
                                  flag_static, t_min, t_max,
                                  cache[best_idx_1].dY, cache[best_idx_1].dW,
                                  cache[best_idx_1].times,
                                  cache[best_idx_1].uid, tpl_main,
                                  tpl_main_st),
                    *cache[best_idx_1].pZ, best_A, Omega,
                    cache[best_idx_1].n_rows))
        }
        stage2_level = 0
        stage2_points = 0
        stage2_same_split = 0
        stage2_rel_gain = 0
        stage2_converged = 0
        if (n_refine > 0) {
            ref_exhausted = 0
            ref_complete = 0
        }
        return
    }

    // Stage 2: one fixed W2.
    best_idx_2 = 0
    best_obj_2 = .
    best_gamma_2 = .
    best_theta_2 = J(k_W, 1, 0)
    best_V_2 = J(k_W, k_W, 0)
    obj_hi2 = .
    n_prof2 = 0
    gamma_adm2 = J(0, 1, .)

    xdpt2_profile_update(gamma_grid, cache, 1, n1, 1, W_n_2,
                          best_idx_2, best_obj_2, best_gamma_2,
                          best_theta_2, best_V_2, obj_hi2, n_prof2,
                          gamma_adm2)

    if (search_mode == "adaptive") {
        prev_idx = best_idx_2
        prev_obj = best_obj_2
        xdpt2_profile_update(gamma_grid, cache, n1 + 1, n2, 1, W_n_2,
                              best_idx_2, best_obj_2, best_gamma_2,
                              best_theta_2, best_V_2, obj_hi2, n_prof2,
                              gamma_adm2)
        stage2_level = 2
        stage2_points = n2
        xdpt2_profile_stability(prev_idx, best_idx_2, prev_obj, best_obj_2,
                                 gamma_grid, q_split_supp, search_tol,
                                 stage2_same_split, stage2_rel_gain,
                                 stage2_converged)

        // Stage 2 may need a denser grid than stage 1 because its basin can
        // differ. It may not stop on a grid coarser than the one that selected
        // the residuals used to construct W2. This also guarantees that every
        // cached point retained for inference was evaluated under fixed W2.
        if (search_max_level == 2 & !stage2_converged) search_hit_max = 1

        if (search_max_level >= 3 &
            (!stage2_converged | stage1_level == 3)) {
            if (n3 >= .) {
                xdpt2_append_nested_level(units, method, flag_static,
                                           flag_kink, t_min, t_max, q_supp,
                                           min_user, gamma_grid, cache,
                                           first_new, last_new,
                                           tpl_main, tpl_main_st)
                if (first_new > 0) n3 = last_new
                else n3 = n2
                level3_points = n3
            }
            prev_idx = best_idx_2
            prev_obj = best_obj_2
            xdpt2_profile_update(gamma_grid, cache, n2 + 1, n3, 1, W_n_2,
                                  best_idx_2, best_obj_2, best_gamma_2,
                                  best_theta_2, best_V_2, obj_hi2, n_prof2,
                                  gamma_adm2)
            stage2_level = 3
            stage2_points = n3
            xdpt2_profile_stability(prev_idx, best_idx_2, prev_obj,
                                     best_obj_2, gamma_grid, q_split_supp,
                                     search_tol, stage2_same_split,
                                     stage2_rel_gain, stage2_converged)
            if (!stage2_converged) search_hit_max = 1
        }
    }
    else {
        stage2_level = 1
        stage2_points = n1
    }

    n_adm2 = n_prof2
    if (rows(gamma_adm2) > 0) {
        gamma_adm2_lo = min(gamma_adm2)
        gamma_adm2_hi = max(gamma_adm2)
    }
    if (best_idx_2 == 0 |
        xdpt2_profile_is_flat(best_obj_2, obj_hi2, n_prof2)) {
        // v0.9.34 (C3): as above for the one-step fallback
        if (flag_kink) {
            xdpt2_kink_refine(units, method, flag_static, t_min, t_max,
                               q_supp, min_user, gamma_admitted, 0,
                               J(0, 0, .), gamma_grid, cache, best_idx_1,
                               best_obj_1, best_gamma_1, best_theta_1,
                               best_V_1, tpl_main, tpl_main_st)
            r_1 = cache[best_idx_1].dY -
                  cache[best_idx_1].dW * best_theta_1
            Omega = xdpt2_build_cluster_omega(*cache[best_idx_1].pZ, r_1,
                                               cache[best_idx_1].uid)
        }
        best_gamma = best_gamma_1
        best_obj = best_obj_1
        best_theta = best_theta_1
        best_A = *cache[best_idx_1].pW1
        best_V = xdpt2_gmm_sandwich(cache[best_idx_1].ZW, best_A, Omega,
                                     cache[best_idx_1].n_rows)
        if (rows(best_V) == 0) {
            errprintf("xtdpthresh: cluster-robust variance failed on the one-step fallback\n")
            exit(498)
        }
        best_V_influence = best_V
        // v0.9.32: kink model, joint (theta, gamma) variance; v0.9.35: the
        // jump model too, with the kernel derivative column.
        if (flag_kink) {
            xdpt2_kink_apply_joint(best_V,
                xdpt2_kink_joint_V(cache[best_idx_1].ZW,
                    xdpt2_kink_xg(units, best_gamma_1, method, flag_static,
                                  t_min, t_max,
                                  best_theta_1[rows(best_theta_1)],
                                  cache[best_idx_1].dY,
                                  cache[best_idx_1].times,
                                  cache[best_idx_1].uid),
                    *cache[best_idx_1].pZ, best_A, Omega,
                    cache[best_idx_1].n_rows))
        }
        else {
            xdpt2_kink_apply_joint(best_V,
                xdpt2_kink_joint_V(cache[best_idx_1].ZW,
                    xdpt2_jump_xg(units, best_gamma_1, best_theta_1, method,
                                  flag_static, t_min, t_max,
                                  cache[best_idx_1].dY, cache[best_idx_1].dW,
                                  cache[best_idx_1].times,
                                  cache[best_idx_1].uid, tpl_main,
                                  tpl_main_st),
                    *cache[best_idx_1].pZ, best_A, Omega,
                    cache[best_idx_1].n_rows))
        }
        stage2_converged = 0
        if (n_refine > 0) {
            ref_exhausted = 0
            ref_complete = 0
        }
        return
    }

    // Save the GLOBAL stage-2 solution before any local support search.
    // Every candidate below is evaluated under this same W2; W2 is never
    // reconstructed after the stage-1 residuals have selected it.
    stage2_global_gamma = best_gamma_2
    stage2_global_obj = best_obj_2
    if (W2_builds != 1) {
        errprintf("xtdpthresh: internal two-step error: W2 was not constructed exactly once\n")
        exit(498)
    }

    // Final jump-only local refinement.  The pool is fixed once from the
    // two adjacent W2-admitted GLOBAL anchors around the stage-2 argmin.
    // Up to 30 never-before-evaluated transformed-support points are added
    // per iteration.  Cache and objective values only expand; no previous
    // point is recomputed and the common W2 remains fixed.
    if (n_refine > 0) {
        real scalar ref_pos, ref_j, ref_first, ref_last
        real scalar ref_grid_lo0, ref_grid_hi0, ref_scale, ref_tol
        real scalar ref_final_pos, ref_nlo, ref_nhi
        real colvector ref_supp, ref_anchors, ref_pool, ref_cand
        real colvector ref_keep, ref_ix, ref_anchors_final

        ref_grid_lo0 = min(gamma_grid)
        ref_grid_hi0 = max(gamma_grid)
        ref_supp = q_split_supp
        if (rows(ref_supp) > 0) {
            ref_supp = select(ref_supp, ref_supp :< .)
            if (rows(ref_supp) > 0) ref_supp = uniqrows(sort(ref_supp, 1))
        }
        ref_anchors = gamma_adm2
        if (rows(ref_anchors) > 0) {
            ref_anchors = uniqrows(sort(ref_anchors, 1))
        }

        ref_pos = 0
        for (ref_j = 1; ref_j <= rows(ref_anchors); ref_j++) {
            if (ref_anchors[ref_j] == best_gamma_2) ref_pos = ref_j
        }

        if (rows(ref_supp) == 0 | rows(ref_anchors) < 2 | ref_pos == 0) {
            ref_remaining = .
            ref_exhausted = 0
            ref_complete = 0
        }
        else {
            ref_alo = (ref_pos > 1 ? ref_anchors[ref_pos - 1] : ref_grid_lo0)
            ref_ahi = (ref_pos < rows(ref_anchors) ?
                       ref_anchors[ref_pos + 1] : ref_grid_hi0)
            ref_pool = select(ref_supp,
                              (ref_supp :> ref_alo) :& (ref_supp :< ref_ahi))
            if (rows(ref_pool) > 0) {
                ref_pool = uniqrows(sort(ref_pool, 1))
                ref_keep = J(rows(ref_pool), 1, 1)
                for (ref_j = 1; ref_j <= rows(ref_pool); ref_j++) {
                    if (sum(gamma_grid :== ref_pool[ref_j]) > 0) {
                        ref_keep[ref_j] = 0
                    }
                }
                if (sum(ref_keep) > 0) ref_pool = select(ref_pool, ref_keep)
                else ref_pool = J(0, 1, .)
            }
            ref_pool_n = rows(ref_pool)

            while (ref_it < n_refine) {
                ref_cand = ref_pool
                if (rows(ref_cand) > 0) {
                    ref_keep = J(rows(ref_cand), 1, 1)
                    for (ref_j = 1; ref_j <= rows(ref_cand); ref_j++) {
                        if (sum(gamma_grid :== ref_cand[ref_j]) > 0) {
                            ref_keep[ref_j] = 0
                        }
                    }
                    if (sum(ref_keep) > 0) ref_cand = select(ref_cand, ref_keep)
                    else ref_cand = J(0, 1, .)
                }
                if (rows(ref_cand) == 0) break

                if (rows(ref_cand) > 30) {
                    // Endpoint-inclusive, deterministic coverage of the
                    // currently remaining pool; at most 30 additions/round.
                    ref_ix = floor((0::29) :*
                             ((rows(ref_cand) - 1) / 29)) :+ 1
                    ref_ix = uniqrows(ref_ix)
                    ref_cand = ref_cand[ref_ix]
                }

                xdpt2_append_candidates(units, ref_cand, method,
                                         flag_static, flag_kink,
                                         t_min, t_max, q_supp, min_user,
                                         gamma_grid, cache,
                                         ref_first, ref_last,
                                         tpl_main, tpl_main_st)
                if (ref_first == 0) break
                ref_added = ref_added + ref_last - ref_first + 1
                ref_it = ref_it + 1
                xdpt2_profile_update(gamma_grid, cache, ref_first, ref_last,
                                      1, W_n_2, best_idx_2, best_obj_2,
                                      best_gamma_2, best_theta_2, best_V_2,
                                      obj_hi2, n_prof2, gamma_adm2)
            }

            ref_remaining = 0
            for (ref_j = 1; ref_j <= rows(ref_pool); ref_j++) {
                if (sum(gamma_grid :== ref_pool[ref_j]) == 0) {
                    ref_remaining = ref_remaining + 1
                }
            }
            ref_exhausted = (ref_remaining == 0)
            ref_in_basin = (best_gamma_2 >= ref_alo & best_gamma_2 <= ref_ahi)

            // A stricter local diagnostic: count transformed-support values
            // still unevaluated between the nearest W2-admitted points around
            // the final argmin.  Failed but attempted candidates count as
            // evaluated, matching the historical refine() contract.
            ref_anchors_final = gamma_adm2
            if (rows(ref_anchors_final) > 0) {
                ref_anchors_final = uniqrows(sort(ref_anchors_final, 1))
            }
            ref_final_pos = 0
            for (ref_j = 1; ref_j <= rows(ref_anchors_final); ref_j++) {
                if (ref_anchors_final[ref_j] == best_gamma_2) {
                    ref_final_pos = ref_j
                }
            }
            if (rows(ref_anchors_final) < 2 | ref_final_pos == 0) {
                ref_neigh_rem = .
                ref_complete = 0
            }
            else {
                ref_nlo = (ref_final_pos > 1 ?
                           ref_anchors_final[ref_final_pos - 1] : ref_grid_lo0)
                ref_nhi = (ref_final_pos < rows(ref_anchors_final) ?
                           ref_anchors_final[ref_final_pos + 1] : ref_grid_hi0)
                ref_neigh_rem = 0
                for (ref_j = 1; ref_j <= rows(ref_supp); ref_j++) {
                    if (ref_supp[ref_j] <= ref_nlo |
                        ref_supp[ref_j] >= ref_nhi) continue
                    if (sum(gamma_grid :== ref_supp[ref_j]) == 0) {
                        ref_neigh_rem = ref_neigh_rem + 1
                    }
                }
                ref_complete = (ref_exhausted & ref_neigh_rem == 0)
            }

            ref_scale = max((abs(stage2_global_obj), abs(best_obj_2)))
            if (ref_scale == 0) ref_obj_gain = 0
            else ref_obj_gain = (stage2_global_obj - best_obj_2) / ref_scale
            if (ref_obj_gain < 0) ref_obj_gain = 0

            // Expanding a candidate set under one fixed objective cannot make
            // its minimum worse. Treat any violation as an internal error.
            ref_tol = xdpt2_objtol(best_obj_2, stage2_global_obj, 1e-10)
            if (best_obj_2 > stage2_global_obj + ref_tol) {
                errprintf("xtdpthresh: internal refinement error: fixed-W2 objective increased\n")
                exit(498)
            }
        }
    }

    // v0.9.34 (C3): kink -- the fixed-W2 profile is minimized between the
    // grid neighbours of the grid argmin (refine() is jump-only).
    if (flag_kink) {
        xdpt2_kink_refine(units, method, flag_static, t_min, t_max, q_supp,
                           min_user, gamma_adm2, 1, W_n_2, gamma_grid, cache,
                           best_idx_2, best_obj_2, best_gamma_2,
                           best_theta_2, best_V_2, tpl_main, tpl_main_st)
    }

    // Refresh the two-step admission diagnostics after local candidates.
    n_adm2 = n_prof2
    if (rows(gamma_adm2) > 0) {
        gamma_adm2_lo = min(gamma_adm2)
        gamma_adm2_hi = max(gamma_adm2)
    }

    // Final two-step fit and conditional cluster-robust VCE.
    best_gamma = best_gamma_2
    best_obj = best_obj_2
    best_theta = best_theta_2
    best_A = W_n_2
    best_twostep = 1

    real colvector r_2_final
    real matrix Omega_2, best_V_cr
    r_2_final = cache[best_idx_2].dY -
                cache[best_idx_2].dW * best_theta_2
    Omega_2 = xdpt2_build_cluster_omega(*cache[best_idx_2].pZ, r_2_final,
                                         cache[best_idx_2].uid)
    best_V_cr = xdpt2_gmm_sandwich(cache[best_idx_2].ZW, best_A, Omega_2,
                                    cache[best_idx_2].n_rows)
    if (rows(best_V_cr) == 0) {
        errprintf("xtdpthresh: cluster-robust variance failed at the final two-step estimate\n")
        exit(498)
    }
    best_V = best_V_cr
    best_V_influence = best_V

    // Reporting-only Windmeijer correction; W2 is not rebuilt.
    external real scalar xdpt_vce_wind, xdpt_wind_applied
    if (xdpt_vce_wind == 1) {
        real colvector gbar2
        real matrix V_wind
        gbar2 = cache[best_idx_2].ZY -
                cache[best_idx_2].ZW * best_theta_2
        external real scalar xdpt_expg
        if (xdpt_expg == 1) {
            external real matrix xdpt_w_ZW1, xdpt_w_X1, xdpt_w_Z, xdpt_w_Om1
            external real matrix xdpt_w_W1, xdpt_w_W2, xdpt_w_ZW2
            external real colvector xdpt_w_uid, xdpt_w_r1, xdpt_w_gbar2
            external real scalar xdpt_w_n
            xdpt_w_ZW1 = cache[best_idx_1].ZW
            xdpt_w_X1 = cache[best_idx_1].dW
            xdpt_w_Z = *cache[best_idx_1].pZ
            xdpt_w_Om1 = Omega
            xdpt_w_W1 = *cache[best_idx_1].pW1
            xdpt_w_W2 = W_n_2
            xdpt_w_ZW2 = cache[best_idx_2].ZW
            xdpt_w_uid = cache[best_idx_1].uid
            xdpt_w_r1 = r_1
            xdpt_w_gbar2 = gbar2
            xdpt_w_n = cache[best_idx_2].n_rows
        }
        V_wind = xdpt2_windmeijer(cache[best_idx_1].ZW,
                                   cache[best_idx_1].dW,
                                   *cache[best_idx_1].pZ,
                                   cache[best_idx_1].uid, r_1, Omega,
                                   *cache[best_idx_1].pW1, W_n_2,
                                   cache[best_idx_2].ZW, gbar2,
                                   cache[best_idx_2].n_rows)
        // v0.9.35 (R2): only a positive semidefinite corrected variance
        if (rows(V_wind) > 0) {
            if (xdpt2_psd_ok(V_wind)) {
                best_V = V_wind
                xdpt_wind_applied = 1
            }
        }
    }

    // v0.9.32: kink model, joint (theta, gamma) variance of the reported type
    // (robust sandwich, or Windmeijer-corrected when that correction was
    // applied above). v0.9.35: the jump model too (xdpt2_jump_xg).
    if (1) {
        real colvector xg_2, xg_1, gbar2_k
        real matrix V_joint
        external real matrix xdpt_V_joint_full
        if (flag_kink) {
            xg_2 = xdpt2_kink_xg(units, best_gamma_2, method, flag_static,
                                 t_min, t_max, best_theta_2[rows(best_theta_2)],
                                 cache[best_idx_2].dY, cache[best_idx_2].times,
                                 cache[best_idx_2].uid)
        }
        else {
            xg_2 = xdpt2_jump_xg(units, best_gamma_2, best_theta_2, method,
                                 flag_static, t_min, t_max,
                                 cache[best_idx_2].dY, cache[best_idx_2].dW,
                                 cache[best_idx_2].times,
                                 cache[best_idx_2].uid, tpl_main, tpl_main_st)
        }
        if (xdpt_wind_applied == 1) {
            if (flag_kink) {
                xg_1 = xdpt2_kink_xg(units, best_gamma_1, method, flag_static,
                                     t_min, t_max,
                                     best_theta_1[rows(best_theta_1)],
                                     cache[best_idx_1].dY,
                                     cache[best_idx_1].times,
                                     cache[best_idx_1].uid)
            }
            else {
                xg_1 = xdpt2_jump_xg(units, best_gamma_1, best_theta_1, method,
                                     flag_static, t_min, t_max,
                                     cache[best_idx_1].dY,
                                     cache[best_idx_1].dW,
                                     cache[best_idx_1].times,
                                     cache[best_idx_1].uid, tpl_main,
                                     tpl_main_st)
            }
            gbar2_k = cache[best_idx_2].ZY -
                      cache[best_idx_2].ZW * best_theta_2
            V_joint = xdpt2_kink_joint_wind(cache[best_idx_1].ZW,
                                            cache[best_idx_1].dW, xg_1,
                                            *cache[best_idx_1].pZ,
                                            cache[best_idx_1].uid, r_1, Omega,
                                            *cache[best_idx_1].pW1, W_n_2,
                                            cache[best_idx_2].ZW, xg_2,
                                            gbar2_k,
                                            cache[best_idx_2].n_rows)
            // v0.9.35 (R2): a corrected joint variance that is not positive
            // semidefinite is not reported: both corrections are dropped and
            // the cluster-robust variances are used (e(vce_applied) = 0), so
            // that e(V), e(V_cond), and the AR statistics stay consistent
            if (rows(V_joint) > 0) {
                if (!xdpt2_psd_ok(xdpt_V_joint_full)) {
                    best_V = best_V_cr
                    xdpt_wind_applied = 0
                    V_joint = xdpt2_kink_joint_V(cache[best_idx_2].ZW, xg_2,
                                                 *cache[best_idx_2].pZ, best_A,
                                                 Omega_2,
                                                 cache[best_idx_2].n_rows)
                }
            }
        }
        else {
            V_joint = xdpt2_kink_joint_V(cache[best_idx_2].ZW, xg_2,
                                         *cache[best_idx_2].pZ, best_A,
                                         Omega_2, cache[best_idx_2].n_rows)
        }
        xdpt2_kink_apply_joint(best_V, V_joint)
    }
}


// v0.8.0 (audit R5, finding #2 FIX): residual-bootstrap PERCENTILE CIs for
// the slope coefficients that account for THRESHOLD-SEARCH variability --
// the response to "analytic SEs are conditional on gamma-hat". Each draw
// rebuilds y* = W(gamma-hat)theta-hat + e-hat*eta_i (unit-level Mammen),
// RE-SEARCHES gamma* over the estimation grid (fast 1-step), and re-
// estimates theta* at that argmin; percentile bounds of theta* are
// returned (2 x k, rows lo/hi). This is a threshold-search-aware wild
// residual approximation, not the Gong-Seo coefficient bootstrap. Batched
// fast path only; returns J(0,0,.) if any admitted grid
// entry lacks it. Called after the other inference objects for reporting
// order; under rseed() it receives its own component-specific seed.
real matrix xdpt2_coef_bootstrap(struct xdpt2_gamma_cache rowvector cache,
                                  real colvector gamma_grid,
                                  real scalar best_gamma,
                                  real colvector best_theta,
                                  real scalar n_boot, real scalar alpha,
                                  real scalar best_2s,
                                  real scalar B_eff, real scalar B_2s,
                                  real scalar B_fb, real scalar n_skip_out,
                                  real scalar valid_out, real scalar n_g1,
                                  real scalar n_g2, real scalar B_att)
{
    real scalar idx, gl, b, n_u, n_rows, k, mb, gsel, j
    real colvector fastl, resid, fit, col_ok, okl, uid_draw
    real matrix ETA, Ymat, OBJ, ZYall, TH, out
    // v0.9.30: centering flag read once for the per-draw cluster Omega
    // (xdpt2_build_cluster_omega_c), instead of an -external- bind per draw.
    external real scalar xdpt_center
    real scalar cb_center
    cb_center = xdpt_center
    // v0.9.31: the reported stage 1 searched the initial grid only; refine()
    // points (appended after it) enter stage 2 alone, so the replay's stage-1
    // set is restricted the same way. Before, it also searched them.
    external real scalar xdpt_n_stage1
    real scalar cb_n1
    cb_n1 = (xdpt_n_stage1 < . ? xdpt_n_stage1 : cols(cache))

    B_eff = 0
    // v0.8.3 R12 (#1/#2): out-arg defaults -- valid_out=1 is set ONLY once
    // the CI matrix is guaranteed to be returned (B_eff alone cannot signal
    // that: 1 <= B_eff <= 9 bails with an empty matrix below).
    valid_out = 0
    n_g1 = 0
    n_g2 = 0
    // v0.8.4 R13 (#2): attempted stays 0 through every early bail (no idx,
    // bad cache entry, no fast-path points) -- nothing was attempted there.
    B_att = 0
    idx = 0
    for (gl = 1; gl <= rows(gamma_grid); gl++) {
        if (gamma_grid[gl] == best_gamma) {
            idx = gl
            gl = rows(gamma_grid) + 1
        }
    }
    external real scalar xdpt_verbose
    if (idx == 0) {
        if (xdpt_verbose) printf("  [coefboot] bail: idx==0\n")
        return(J(0,0,.))
    }
    if (!cache[idx].ok) {
        if (xdpt_verbose) printf("  [coefboot] bail: cache[idx] not ok\n")
        return(J(0,0,.))
    }
    n_rows = cache[idx].n_rows
    // v0.8.0: restrict the gamma* re-search to entries offering the batched
    // fast path with matched rows -- SKIP the others (boundary gammas with
    // ill-conditioned A typically lack it), mirroring the fast-only
    // convention of the CI's sample-side scan. Bail only if none qualify or
    // the best-gamma entry itself is unusable.
    real scalar n_skip
    n_skip = 0
    fastl = J(0, 1, 0)
    okl   = J(0, 1, 0)
    for (gl = 1; gl <= cols(cache); gl++) {
        if (!cache[gl].ok) continue
        if (cache[gl].n_rows != n_rows) {
            n_skip = n_skip + 1
            continue
        }
        // v0.8.2 R11 (#4): two search sets. Stage 1 needs the fast path
        // (its solver rejects cond(ZW'W1 ZW) > 1e12 -- exactly !fast_ok),
        // but the MAIN two-step search runs over all structurally-ok points
        // and tests conditioning under W2 itself; the replay must match.
        okl = okl \ gl
        if (cache[gl].fast_ok != 1) {
            n_skip = n_skip + 1
            continue
        }
        if (gl <= cb_n1) fastl = fastl \ gl
    }
    if (rows(fastl) == 0) {
        if (xdpt_verbose) printf("  [coefboot] bail: no fast-path grid entries\n")
        return(J(0,0,.))
    }
    // v0.8.3 R12 (#3): NO fast-path requirement at the reported gamma-hat.
    // Since R11 the main two-step search runs over all structurally-ok
    // points, so the selected gamma-hat can lack the one-step fast path;
    // the bootstrap DGP here needs only dY/dW/uid/pZ at idx -- never C_g
    // (the draw-level one-step fallback uses cache[gsel], a fast-path
    // point, and the stage-2 re-search runs on okl which includes idx).
    // A one-step-reported gamma-hat always has fast_ok by construction.
    if (xdpt_verbose & n_skip > 0) {
        printf("  [coefboot] gamma* search on %g fast entries (%g skipped)\n",
               rows(fastl), n_skip)
    }
    // v0.8.3 R12 (#2): replay search-space sizes returned via OUTPUT
    // ARGUMENTS. (R11 wrote them into r() from here, but xdpt2_run calls
    // st_rclear() before its own export block, so they never survived --
    // e(boot_grid_stage1/2) read back as 0 on every run.)
    n_g1 = rows(fastl)
    n_g2 = rows(okl)

    k     = rows(best_theta)
    fit   = cache[idx].dW * best_theta
    resid = cache[idx].dY - fit
    uid_draw = xdpt2_dense_uid(cache[idx].uid)
    n_u   = max(uid_draw)

    external real scalar xdpt_coefboot_2s
    real scalar do_2s
    // v0.8.1 R7 (#4.1): replay what was actually reported.
    do_2s = (xdpt_coefboot_2s == 1 & best_2s == 1)
    B_2s = 0
    B_fb = 0
    n_skip_out = n_skip
    B_att = 0
    TH = J(k, n_boot, .)
    real colvector r1_b, zy1_b, zy_keep
    real matrix Om_b, W2_b, A2_b, A2inv_b
    real scalar g2sel, mb2, obj2, gl2, ch, inv_ok_b
    real scalar nprof1_b, hi1_b, tol1_b, prof_scale_b
    real scalar nprof2_b, hi2_b, tol2_b
    // Draw exactly B replications. Failed two-step solves are not replaced
    // by one-step estimates and are not redrawn: replacement draws would
    // hide the original failure rate. If failures occur, quantiles over the
    // survivors remain conditional on solver success; the 90% gate and
    // exported fail rate make that limitation explicit.
    // B_fb counts failed draws and B_att equals the requested B once the
    // draw loop starts.
    ch = n_boot
        ETA = J(n_u, ch, 0)
        for (b = 1; b <= ch; b++) ETA[., b] = xdpt2_mammen_draw(n_u)
        Ymat = fit :+ resid :* ETA[uid_draw, .]
        // v0.9.31: objectives by the split (xdpt2_fast_obj_split_list)
        OBJ = xdpt2_fast_obj_split_list(fit, resid, uid_draw, ETA, cache, fastl)
        // Z is gamma-invariant (bitwise): one batched cross product serves
        // every grid point.
        ZYall = (*cache[idx].pZ)' * Ymat / n_rows
        for (b = 1; b <= ch; b++) {
            B_att = B_att + 1
            mb = .
            gsel = 0
            nprof1_b = 0
            hi1_b = .
            for (gl = 1; gl <= rows(fastl); gl++) {
                if (OBJ[gl, b] < .) {
                    nprof1_b = nprof1_b + 1
                    if (hi1_b >= . | OBJ[gl, b] > hi1_b) hi1_b = OBJ[gl, b]
                    if (mb >= .) {
                        mb = OBJ[gl, b]
                        gsel = fastl[gl]
                    }
                    else {
                        tol1_b = xdpt2_objtol(OBJ[gl, b], mb, 1e-12)
                        if (OBJ[gl, b] < mb - tol1_b |
                            (abs(OBJ[gl, b] - mb) <= tol1_b &
                             gamma_grid[fastl[gl]] < gamma_grid[gsel])) {
                            mb = OBJ[gl, b]
                            gsel = fastl[gl]
                        }
                    }
                }
            }
            if (gsel == 0 | nprof1_b < 2) {
                B_fb = B_fb + 1
                continue
            }
            prof_scale_b = max((abs(mb), abs(hi1_b)))
            if (abs(hi1_b - mb) <= 1e-12 * prof_scale_b) {
                B_fb = B_fb + 1
                continue
            }
            if (!do_2s) {
                B_eff = B_eff + 1
                TH[., B_eff] = cache[gsel].C_g * ZYall[., b]
                continue
            }
            // ---- stage 1 residuals at gamma1* ----
            r1_b = Ymat[., b] - cache[gsel].dW * (cache[gsel].C_g * ZYall[., b])
            Om_b = xdpt2_build_cluster_omega_c(*cache[gsel].pZ, r1_b,
                                               cache[gsel].uid, cb_center)
            xdpt2_syminv(Om_b, inv_ok_b, W2_b)
            if (!inv_ok_b) {
                B_fb = B_fb + 1
                continue
            }
            // ---- stage 2: grid pass with W2* fixed ----
            mb2 = .
            g2sel = 0
            nprof2_b = 0
            hi2_b = .
            zy_keep = J(k, 1, .)
            for (gl2 = 1; gl2 <= rows(okl); gl2++) {
                A2_b = cache[okl[gl2]].ZW' * W2_b * cache[okl[gl2]].ZW
                xdpt2_syminv(A2_b, inv_ok_b, A2inv_b)
                if (!inv_ok_b) continue
                zy1_b = A2inv_b * (cache[okl[gl2]].ZW' * (W2_b * ZYall[., b]))
                r1_b  = ZYall[., b] - cache[okl[gl2]].ZW * zy1_b
                obj2  = n_rows * (r1_b' * W2_b * r1_b)
                if (obj2 >= . | hasmissing(zy1_b)) continue
                nprof2_b = nprof2_b + 1
                if (hi2_b >= . | obj2 > hi2_b) hi2_b = obj2
                if (mb2 >= .) {
                    mb2 = obj2
                    g2sel = gl2
                    zy_keep = zy1_b
                }
                else {
                    tol2_b = xdpt2_objtol(obj2, mb2, 1e-12)
                    if (obj2 < mb2 - tol2_b |
                        (abs(obj2 - mb2) <= tol2_b &
                         gamma_grid[okl[gl2]] < gamma_grid[okl[g2sel]])) {
                        mb2 = obj2
                        g2sel = gl2
                        zy_keep = zy1_b
                    }
                }
            }
            if (g2sel == 0 | nprof2_b < 2) {
                B_fb = B_fb + 1
                continue
            }
            prof_scale_b = max((abs(mb2), abs(hi2_b)))
            if (abs(hi2_b - mb2) <= 1e-12 * prof_scale_b) {
                B_fb = B_fb + 1
                continue
            }
            B_eff = B_eff + 1
            B_2s = B_2s + 1
            TH[., B_eff] = zy_keep
        }
    // v0.9.3 R19 (#7): hard validity gate -- at least 90% of the request
    // and never fewer than 10 valid replications (quantiles from a handful
    // of draws are numerics, not inference).
    if (B_eff < 10 | B_eff < ceil(0.9 * n_boot)) {
        if (xdpt_verbose & n_boot > 0) {
            printf("  [coefboot] bail: valid=%g of %g requested (attempted %g)\n",
                   B_eff, n_boot, B_att)
        }
        return(J(0,0,.))
    }

    // v0.8.1 (audit R6, #5): SYMMETRIC percentile intervals by default --
    // theta_hat +/- c*, with c* the (1-alpha) quantile of |theta*-theta_hat|.
    // coefcitype(percentile) gives Efron's percentile interval
    // [q(alpha/2), q(1-alpha/2)] of theta*. The percentile interval studied
    // by Gong and Seo (2026, eq. 9) is the basic interval, and their
    // coverage results concern their own bootstrap, not this legacy scheme.
    external real scalar xdpt_coefci_sym
    out = J(2, k, .)
    for (j = 1; j <= k; j++) {
        col_ok = select(TH[j, .]', TH[j, .]' :< .)
        if (rows(col_ok) != B_eff) return(J(0,0,.))
        if (xdpt_coefci_sym == 1) {
            real scalar c_j
            real colvector dev_j
            dev_j = abs(col_ok :- best_theta[j])
            // A finite theta* and theta-hat can still overflow in their
            // difference. Do not let xdpt2_quantile() silently discard that
            // draw and manufacture a degenerate interval.
            if (hasmissing(dev_j)) return(J(0,0,.))
            c_j = xdpt2_quantile(dev_j, 1 - alpha)
            out[1, j] = best_theta[j] - c_j
            out[2, j] = best_theta[j] + c_j
        }
        else {
            out[1, j] = xdpt2_quantile(col_ok, alpha/2)
            out[2, j] = xdpt2_quantile(col_ok, 1 - alpha/2)
        }
    }
    if (hasmissing(out)) return(J(0,0,.))
    valid_out = 1
    return(out)
}

// v0.9.33: 1 if the normal matrix ZW' A ZW of the fixed-weight solve
// (xdpt2_solve_gmm_1step_pre, as run by the two-step search with A = W2)
// passes the gate of xdpt2_syminv, else 0.
real scalar xdpt2_w2_solvable(real matrix ZW, real matrix A)
{
    real scalar inv_ok
    real matrix M, Minv
    if (rows(A) == 0 | rows(A) != cols(A) | rows(A) != rows(ZW)) return(0)
    M = ZW' * A * ZW
    xdpt2_syminv(M, inv_ok, Minv)
    return(inv_ok)
}

// v0.9.34 (C1): fixed-W2 solve of every ok entry, with the expressions of
// xdpt2_solve_gmm_1step_pre (A = ZW' W ZW, theta = Ainv ZW' W ZY), so that
// C_g2 * ZY is bit for bit the stage-2 coefficient vector at that gamma.
void xdpt2_cache_w2(struct xdpt2_gamma_cache rowvector cache, real matrix W2)
{
    real scalar g, inv_ok
    real matrix A2, A2inv
    for (g = 1; g <= cols(cache); g++) {
        cache[g].fast2_ok = 0
        cache[g].C_g2 = J(0, 0, .)
        if (!cache[g].ok) continue
        if (rows(W2) != rows(cache[g].ZW) | cols(W2) != rows(cache[g].ZW)) continue
        A2 = cache[g].ZW' * W2 * cache[g].ZW
        xdpt2_syminv(A2, inv_ok, A2inv)
        if (!inv_ok) continue
        cache[g].C_g2 = A2inv * cache[g].ZW' * W2
        if (!hasmissing(cache[g].C_g2)) cache[g].fast2_ok = 1
    }
}

// v0.9.36: xdpt2_fast_gmm_boot with the fixed second-step weight when w2m = 1
// (the solve C_g2 of xdpt2_cache_w2 and the objective n g'W2 g), else the
// one-step pair (C_g, W1) exactly as xdpt2_fast_gmm_boot.
void xdpt2_fast_gmm_boot_w(real colvector Y_boot,
                            struct xdpt2_gamma_cache scalar gc,
                            real scalar w2m, real matrix W2,
                            real scalar ok, real colvector theta,
                            real scalar obj)
{
    real colvector ZY, r, g
    if (!w2m) {
        xdpt2_fast_gmm_boot(Y_boot, gc, ok, theta, obj)
        return
    }
    ok = 0
    if (gc.ok == 0 | gc.fast2_ok != 1) return
    if (rows(Y_boot) != gc.n_rows) return
    ZY = (*gc.pZ)' * Y_boot / gc.n_rows
    theta = gc.C_g2 * ZY
    if (hasmissing(theta)) return
    r = Y_boot - gc.dW * theta
    if (hasmissing(r)) return
    g = (*gc.pZ)' * r / gc.n_rows
    if (hasmissing(g)) return
    obj = gc.n_rows * (g' * W2 * g)
    if (obj >= .) return
    ok = 1
}

// Grid-bootstrap inversion. The default is an xthenreg-style cluster wild
// residual approximation; boottype(unit) is an experimental unit-resampling
// extension. Neither path is certified as the exact Gong-Seo Algorithm 1.
// Returns the convex-hull interval summary when inversion is complete.
void xdpt2_grid_bootstrap(struct xdpt2_unit rowvector units,
                           struct xdpt2_gamma_cache rowvector gamma_cache,
                           real colvector gamma_grid,
                           real colvector gamma_ci_grid,
                           real colvector q_supp, real scalar min_user,
                           real scalar best_obj, real scalar best_gamma,
                           string scalar method, real scalar flag_static,
                           real scalar flag_kink,
                           real scalar t_min, real scalar t_max,
                           real scalar n_boot, real scalar alpha,
                           real scalar gam_lo, real scalar gam_hi,
                           real scalar ci_empty, real scalar ci_nseg,
                           real scalar gci_adm, real scalar gci_lo,
                           real scalar gci_hi,
                           real scalar bt2s, real matrix bA,
                           real colvector rhat, real scalar gb_minB,
                           real matrix ci_tab, real matrix ci_seg,
                           real scalar ci_unres,
                           struct xdpt2_stack_tpl scalar tpl_main,
                           real scalar tpl_main_st,
                           | real colvector D_out)
{
    real scalar n_ci, l, b, ok_r, obj_r, D_sample, D_boot, crit, min_obj_b
    real scalar has_alt_b
    real scalar ok_b_r, obj_b_r, ok_b_u, obj_b_u, gb
    real matrix dY_r, dW_r, Z_r, V_r, V_dummy
    real colvector times_r, theta_r, resid_r, Y_boot, theta_b
    real colvector D_vec, accept, uid_r, uid_b, uid_draw
    real colvector eta_unit, eta
    real scalar n_u, n_draw, i, u
    real matrix dY_b, dW_b, Z_b
    real colvector times_b, theta_b_r, theta_b_u
    // Declarations for 1-step sample D (Gong-Seo Alg. 1 consistency fix)
    real scalar obj_r_1s, best_obj_1s, ok_1s, gb_s, obj_u_1s
    real colvector theta_s_dummy
    theta_s_dummy = J(0, 1, 0)

    n_ci = rows(gamma_ci_grid)
    accept = J(n_ci, 1, 0)
    // v0.9.3 R19 (#8): full inversion table -- the convex hull alone hides
    // disconnected acceptance regions.
    // v0.9.4 R20 (#3): the accepted column starts MISSING, not 0 --
    // "could not be evaluated" is not "rejected". Column 6 = status:
    //   1 valid inversion result        4 sample solve failed
    //   2 mechanical accept (D == 0)    5 insufficient valid draws
    //   3 not admissible (no solve      6 bootstrap quantile failed (a
    //                                     guard; unreachable when every draw
    //                                     is valid)
    //     under W1 nor, for a two-step fit, under W2)
    // Statuses 4-6 are UNRESOLVED: excluded from the reported set but
    // counted in ci_unres / e(ci_unresolved) and flagged loudly.
    ci_tab = J(n_ci, 6, .)
    if (n_ci > 0) ci_tab[., 1] = gamma_ci_grid
    ci_seg = J(0, 2, .)
    n_u = length(units)
    // v0.9.34: the unrestricted minimum runs over one fixed set on the sample
    // side and in every draw: the initial grid (entries 1..n1_ci), plus the
    // candidate itself. Points added around the sample's own minimum
    // (refine() support points, the kink refinement) would lower the sample
    // minimum only, since no draw gets the same search, and the test would
    // over-reject. best_obj is the minimum over that grid (the stage-2
    // minimum for a two-step fit); gamma-hat itself still has D = 0.
    external real scalar xdpt_n_stage1
    real scalar n1_ci
    n1_ci = (xdpt_n_stage1 < . ? min((xdpt_n_stage1, cols(gamma_cache))) :
                                 cols(gamma_cache))

    // === Per-gamma caches: the estimation-grid cache is passed in (built
    // once in xtdpthresh_run, v0.7.0 D1); only the CI-grid cache is built here ===
    struct xdpt2_gamma_cache rowvector gamma_ci_cache
    gamma_ci_cache = xdpt2_build_gamma_cache_t(units, gamma_ci_grid, method,
                                                flag_static, flag_kink, t_min,
                                                t_max, q_supp, min_user,
                                                tpl_main, tpl_main_st,
                                                gamma_cache)

    // v0.9.34 (C1): the wild inversion uses the criterion of the reported
    // estimator. With a two-step fit (bt2s = 1) that is the stage-2 criterion
    // with the second-step weight W2 (bA) held at its sample value, whose
    // minimum over the grid is attained at gamma-hat: D(gamma-hat) = 0 and
    // gamma-hat always belongs to its own set, as in Gong and Seo (2026),
    // whose statistic is built on the estimation criterion. Up to 0.9.33 the
    // wild inversion used the one-step criterion (W1), whose zero is the
    // one-step argmin, and the two-step gamma-hat could be rejected. After a
    // one-step fallback (bt2s = 0) the reported criterion is the one-step one.
    real scalar w2m
    w2m = (xdpt_boot_exact != 1 & bt2s == 1 & rows(bA) > 0)
    if (w2m) {
        xdpt2_cache_w2(gamma_cache, bA)
        xdpt2_cache_w2(gamma_ci_cache, bA)
    }
    // v0.9.34 (O1): common random numbers -- one set of Mammen weights for
    // all candidates (the clusters are the same at every gamma), drawn at
    // the first candidate that needs them.
    real matrix ETA_crn
    ETA_crn = J(0, 0, .)

    // Filled from the ACTUAL sample-side statuses after inversion. In
    // particular, boottype(unit) uses the fixed W2 solve and must not inherit
    // the wild path's original-sample fast_ok gate.
    real scalar gci_l
    gci_adm = 0
    gci_lo = .
    gci_hi = .

    // === v0.7.9 (A): hoist the sample-side unrestricted grid minimum ===
    // Inside the l-loop, D_sample needs the minimum over the gamma grid of
    // the 1-step objective evaluated on dY_r -- but dY is gamma-invariant,
    // so dY_r is the SAME vector at every CI point and the scan recomputed
    // the identical minimum n_ci times. It is computed ONCE here on the
    // first ok CI entry's dY; each l reuses it only after an EXACT bitwise
    // dY comparison (falling back to the original scan otherwise).
    // xdpt2_fast_gmm_boot is deterministic, so every objective the per-l
    // scan would produce is bitwise equal to the hoisted one, and min()
    // over bitwise-identical values is order-free -- D_sample is
    // bit-for-bit unchanged. No RNG is involved in the scan.
    real scalar smin_ready, smin_has, smin_val, smin_ref_l, l0s
    smin_ready = 0
    smin_has = 0
    smin_val = .
    smin_ref_l = 0
    for (l0s = 1; l0s <= n_ci; l0s++) {
        if (!gamma_ci_cache[l0s].ok) continue
        smin_ref_l = l0s
        break
    }
    if (smin_ref_l > 0 & !w2m) {
        for (gb_s = 1; gb_s <= n1_ci; gb_s++) {
            if (!gamma_cache[gb_s].ok) continue
            if (gamma_cache[gb_s].n_rows != gamma_ci_cache[smin_ref_l].n_rows) continue
            xdpt2_fast_gmm_boot(gamma_ci_cache[smin_ref_l].dY,
                                 gamma_cache[gb_s],
                                 ok_1s, theta_s_dummy, obj_u_1s)
            if (!ok_1s) continue
            if (!smin_has) {
                smin_val = obj_u_1s
                smin_has = 1
            }
            else if (obj_u_1s < smin_val) smin_val = obj_u_1s
        }
        smin_ready = 1
    }

    external real scalar xdpt_verbose
    external real scalar xdpt_boot_exact
    if (xdpt_verbose) {
        if (xdpt_boot_exact == 1) printf("  Grid bootstrap (B=%g, gridci=%g, EXPERIMENTAL unit resampling -- Alg. 1-oriented, not certified)...\n", n_boot, n_ci)
        else printf("  Grid bootstrap (B=%g, gridci=%g, unit-level Mammen; cached)...\n", n_boot, n_ci)
    }
    else {
        printf("  Grid bootstrap CI  (. per γ point, %g total)\n", n_ci)
        printf("  ")
        displayflush()
    }

    for (l = 1; l <= n_ci; l++) {
        if (!xdpt_verbose) {
            printf(".")
            displayflush()
        }
        // Sample: restricted at γ_ℓ — pull from cache
        if (!gamma_ci_cache[l].ok) {
            accept[l] = 0
            ci_tab[l, 6] = 3
            continue
        }
        dY_r    = gamma_ci_cache[l].dY
        dW_r    = gamma_ci_cache[l].dW
        Z_r     = *gamma_ci_cache[l].pZ
        times_r = gamma_ci_cache[l].times
        uid_r   = gamma_ci_cache[l].uid
        uid_draw = xdpt2_dense_uid(uid_r)
        n_draw = max(uid_draw)

        if (rows(dY_r) < 20) {
            accept[l] = 0
            ci_tab[l, 6] = 3
            continue
        }
        // v0.9.6 R22 (#2/#3): the ONE-STEP sample statistic, its D == 0
        // mechanical-accept shortcut, and its W_first solvability gate
        // belong to the WILD inversion only. Under boottype(unit) the
        // sample statistic is the fixed-W2 two-stage criterion computed
        // inside the unit branch -- running the one-step shortcut first
        // auto-accepted the one-step argmin (D_1step = 0 there by
        // construction, and that gamma is appended to the CI grid for the
        // wild inversion) even when the two-stage statistic is nonzero,
        // and the W_first gate could kill a point (status 4) whose
        // fixed-W2 solve is perfectly feasible.
        if (xdpt_boot_exact != 1) {
            // Sample D_sample. v0.9.34 (C1): with a two-step fit, the stage-2
            // restricted fit with the fixed W2 (the solve of the reported
            // search, so its objective is bit for bit the stage-2 profile) and
            // the stage-2 minimum best_obj; with a one-step fit, the one-step
            // criterion (W1) and its grid minimum as before.
            if (w2m) {
                xdpt2_solve_gmm_1step_pre(dY_r, dW_r, Z_r,
                                           gamma_ci_cache[l].ZW,
                                           gamma_ci_cache[l].ZY, bA,
                                           ok_r, theta_r, obj_r_1s, V_dummy)
                if (!ok_r) {
                    accept[l] = 0
                    // not admissible (3) when the reported estimator's normal
                    // matrix fails the gate there; any other failure is
                    // unresolved (4)
                    ci_tab[l, 6] = (gamma_ci_cache[l].fast2_ok == 1 ? 4 : 3)
                    continue
                }
                best_obj_1s = obj_r_1s
                if (best_obj < best_obj_1s) best_obj_1s = best_obj
            }
            else {
                xdpt2_fast_gmm_boot(dY_r, gamma_ci_cache[l],
                                     ok_r, theta_r, obj_r_1s)
                if (!ok_r) {
                    accept[l] = 0
                    // v0.9.31/0.9.33: not admissible (3) when the reported
                    // (one-step) estimator cannot use the point; any other
                    // failure of the sample solve is unresolved (4)
                    ci_tab[l, 6] = (gamma_ci_cache[l].fast_ok == 0 ? 3 : 4)
                    continue
                }
                // Sample unrestricted: min 1-step obj over γ_grid using dY_r
                // v0.7.9 (A): reuse the hoisted grid minimum under an exact
                // bitwise dY guard (see the smin_* block above); identical dY
                // implies an identical participating cache set and bitwise-
                // identical objectives, and min() is order-free.
                best_obj_1s = obj_r_1s
                if (smin_ready & dY_r == gamma_ci_cache[smin_ref_l].dY) {
                    if (smin_has) {
                        if (smin_val < best_obj_1s) best_obj_1s = smin_val
                    }
                }
                else {
                    for (gb_s = 1; gb_s <= n1_ci; gb_s++) {
                        if (!gamma_cache[gb_s].ok) continue
                        if (gamma_cache[gb_s].n_rows != rows(dY_r)) continue
                        xdpt2_fast_gmm_boot(dY_r, gamma_cache[gb_s],
                                             ok_1s, theta_s_dummy, obj_u_1s)
                        if (!ok_1s) continue
                        if (obj_u_1s < best_obj_1s) best_obj_1s = obj_u_1s
                    }
                }
            }
            D_sample = obj_r_1s - best_obj_1s
            // best_obj_1s is at most the restricted objective, so this
            // distance is nonnegative by construction. Only an EXACT zero
            // can be accepted mechanically; every positive statistic needs
            // its bootstrap critical value, regardless of outcome units.
            if (D_sample == 0) {
                accept[l] = 1
                ci_tab[l, 2] = D_sample
                ci_tab[l, 4] = 1
                ci_tab[l, 6] = 2
                continue
            }

            // Residuals under restricted DGP
            resid_r = dY_r - dW_r * theta_r
        }

        // Bootstrap loop  [FAST PATH: uses precomputed C_g, no cluster-Ω]
        // v0.7.6 SPEEDUP: when γ_ℓ uses the fast path, the whole
        // B-replication bootstrap is done in batched matrix form (0.9.31:
        // xdpt2_fast_obj_split_list). Mammen weights are drawn in the SAME
        // per-replication order as the scalar loop, so D_vec -- and the
        // resulting CI -- agree with the scalar loop to rounding error.
        // v0.9.33: an unrestricted γ without the one-step solve
        // (fast_ok = 0) is left out of the batch, as the scalar loop skips it
        // (xdpt2_fast_gmm_boot returns ok = 0 there). Before, one such γ
        // sent every CI point to the scalar loop, which is hundreds of times
        // slower (e.g., γ = 0 when q >= 0 is also a regressor, where
        // q*1(q > γ) equals q).
        D_vec = J(n_boot, 1, .)
        if (args() == 32) D_out = J(0, 1, .)
        real scalar n_rows_r_b, use_batch_b
        real colvector fast_gb_b
        real matrix ETA_b, OBJ_b
        real colvector F_b
        real rowvector objr_b, altmin_b, minobj_b, Dvec_b
        n_rows_r_b  = rows(dY_r)
        use_batch_b = (gamma_ci_cache[l].ok == 1 &
                       (w2m ? gamma_ci_cache[l].fast2_ok :
                              gamma_ci_cache[l].fast_ok) == 1 &
                       gamma_ci_cache[l].n_rows == n_rows_r_b)
        fast_gb_b = J(0, 1, 0)
        if (use_batch_b) {
            for (gb = 1; gb <= n1_ci; gb++) {
                if (!gamma_cache[gb].ok) continue
                if (gamma_cache[gb].n_rows != n_rows_r_b) continue
                // v0.9.34 (C1): the unrestricted set of the reported search,
                // the points solvable with the weight of the inversion
                if ((w2m ? gamma_cache[gb].fast2_ok : gamma_cache[gb].fast_ok) != 1) continue
                fast_gb_b = fast_gb_b \ gb
            }
        }

        // ============ v0.9.2 R18 (#1): boottype(unit) REWRITE ============
        // Reviewer-specified Alg. 1 structure (round 18):
        //  - SAMPLE statistic: two-stage criterion under the run's fixed
        //    second-step weight W_n_2 (bA), not the 1-step W_first;
        //  - DGP: Y* = W(gamma_l) alpha2(gamma_l) + eps-hat, where eps-hat
        //    are the UNRESTRICTED residuals at the reported
        //    (gamma-hat, theta-hat): restricted coefficients, unrestricted
        //    residuals;
        //  - RECENTERING: subtract the sample moment at theta-hat (S_hat),
        //    the same vector at every null point;
        //  - per draw: stage-1 argmin under W_first -> bootstrap residuals
        //    at that argmin -> recentered per-unit moments -> Omega* ->
        //    W2* -> stage-2 restricted and grid-min objectives (mirrors
        //    the sample's two-stage fixed-weight construction draw by
        //    draw). Labeled boottype(unit): Alg. 1-ORIENTED, NOT certified.
        if (xdpt_boot_exact == 1) {
            real colvector ex_cids, ex_w, ex_wrow, ex_pick, ex_Shat
            real colvector ex_Yb0, ex_ZYs, ex_m, ex_th, ex_rb, ex_th1b, ex_v
            real matrix ex_Zw, ex_Gt, ex_Om, ex_W2, ex_A, ex_Ainv, ex_W1ref
            real matrix ex_ZWall, ex_ZWg, ex_ZWl
            real rowvector ex_sb
            real scalar ex_nc, ex_nb, ex_g2, ex_j2, ex_obj, ex_objr, ex_min
            real scalar ex_o1, ex_b1, ex_ok2, ex_J2l, ex_k, ex_c0
            real scalar ex_att, ex_valid, ex_invok
            real scalar ex_n1, ex_hi1, ex_tol1, ex_scale
            real scalar ex_n2, ex_hi2, ex_tol2, ex_b2
            external real scalar xdpt_center
            // v0.9.31: stage 1 of the replay on the initial grid only, as in
            // the reported search (refine() points enter stage 2 alone).
            real scalar ex_n1max
            ex_n1max = n1_ci
            real colvector ex_okl
            ex_okl = J(0, 1, 0)
            for (gb = 1; gb <= n1_ci; gb++) {
                if (!gamma_cache[gb].ok) continue
                if (gamma_cache[gb].n_rows != n_rows_r_b) continue
                ex_okl = ex_okl \ gb
            }
            if (rows(ex_okl) == 0 | rows(rhat) != n_rows_r_b | bt2s != 1) {
                accept[l] = 0
                ci_tab[l, 6] = 4
                continue
            }
            ex_W1ref = *gamma_ci_cache[l].pW1
            ex_k = cols(dW_r)
            // --- two-stage SAMPLE statistic under the fixed W_n_2 ---
            xdpt2_solve_gmm_1step_pre(dY_r, dW_r, Z_r,
                                       gamma_ci_cache[l].ZW,
                                       gamma_ci_cache[l].ZY, bA,
                                       ex_ok2, ex_th, ex_J2l, V_dummy)
            if (!ex_ok2) {
                accept[l] = 0
                // v0.9.34: the rule of the wild path -- a point whose normal
                // matrix under W2 fails the gate (for example an exactly
                // collinear design, q*1(q > gamma) = q) is not admissible (3)
                // for the reported estimator; any other failure is
                // unresolved (4). Before, every failure was 4, so the unit
                // bootstrap never reported a set on such data.
                ci_tab[l, 6] = (xdpt2_w2_solvable(gamma_ci_cache[l].ZW, bA) ? 4 : 3)
                continue
            }
            // The unrestricted comparison must contain the null candidate
            // even when the CI point is not on the estimation grid.
            D_sample = ex_J2l - min((ex_J2l, best_obj))
            if (D_sample == 0) {
                accept[l] = 1
                ci_tab[l, 2] = D_sample
                ci_tab[l, 4] = 1
                ci_tab[l, 6] = 2
                continue
            }
            // --- DGP pieces: restricted 2-stage coefs + UNRESTRICTED resid
            ex_Yb0 = dW_r * ex_th + rhat
            ex_Shat = Z_r' * rhat
            ex_cids = uniqrows(uid_r)
            ex_nc = rows(ex_cids)
            // Use exactly the requested B random draws. Numerical failures
            // stay missing and count against the all-B unit validity rule;
            // drawing replacements would condition the bootstrap distribution
            // on solver success and hide the original failure rate.
            ex_valid = 0
            for (ex_att = 1; ex_att <= n_boot; ex_att++) {
                // floor(u*n)+1: runiform() can return exactly 0, which
                // ceil() would map to index 0 (user report, round 18)
                ex_pick = floor(runiform(ex_nc, 1) :* ex_nc) :+ 1
                ex_w = J(n_u, 1, 0)
                for (ex_j2 = 1; ex_j2 <= ex_nc; ex_j2++) {
                    ex_w[ex_cids[ex_pick[ex_j2]]] = ex_w[ex_cids[ex_pick[ex_j2]]] + 1
                }
                ex_wrow = ex_w[uid_r]
                ex_nb = sum(ex_wrow)
                if (ex_nb < 20) continue
                ex_Zw = Z_r :* ex_wrow
                ex_ZYs = ex_Zw' * ex_Yb0
                // weighted cross-products, cached once per draw for both stages
                ex_ZWall = J(cols(Z_r), ex_k * rows(ex_okl), 0)
                for (ex_g2 = 1; ex_g2 <= rows(ex_okl); ex_g2++) {
                    ex_c0 = (ex_g2 - 1) * ex_k
                    ex_ZWall[|1, ex_c0 + 1 \ cols(Z_r), ex_c0 + ex_k|] = ex_Zw' * gamma_cache[ex_okl[ex_g2]].dW
                }
                // stage 1: argmin over W_first-solvable candidates
                ex_o1 = .
                ex_b1 = 0
                ex_n1 = 0
                ex_hi1 = .
                for (ex_g2 = 1; ex_g2 <= rows(ex_okl); ex_g2++) {
                    if (ex_okl[ex_g2] > ex_n1max) continue
                    ex_c0 = (ex_g2 - 1) * ex_k
                    ex_ZWg = ex_ZWall[|1, ex_c0 + 1 \ cols(Z_r), ex_c0 + ex_k|]
                    ex_A = ex_ZWg' * (ex_W1ref * ex_ZWg)
                    xdpt2_syminv(ex_A, ex_invok, ex_Ainv)
                    if (!ex_invok) continue
                    ex_th1b = ex_Ainv * (ex_ZWg' * (ex_W1ref * (ex_ZYs - ex_Shat)))
                    ex_m = ex_ZYs - ex_ZWg * ex_th1b - ex_Shat
                    ex_v = ex_W1ref * ex_m
                    ex_obj = (ex_m' * ex_v) / ex_nb
                    if (ex_obj >= . | hasmissing(ex_th1b)) continue
                    ex_n1 = ex_n1 + 1
                    if (ex_hi1 >= . | ex_obj > ex_hi1) ex_hi1 = ex_obj
                    if (ex_o1 >= .) {
                        ex_o1 = ex_obj
                        ex_b1 = ex_g2
                        ex_th = ex_th1b
                    }
                    else {
                        ex_tol1 = xdpt2_objtol(ex_obj, ex_o1, 1e-12)
                        if (ex_obj < ex_o1 - ex_tol1 |
                            (abs(ex_obj - ex_o1) <= ex_tol1 &
                             gamma_grid[ex_okl[ex_g2]] <
                             gamma_grid[ex_okl[ex_b1]])) {
                            ex_o1 = ex_obj
                            ex_b1 = ex_g2
                            ex_th = ex_th1b
                        }
                    }
                }
                if (ex_b1 == 0 | ex_n1 < 2) continue
                ex_scale = max((abs(ex_o1), abs(ex_hi1)))
                if (abs(ex_hi1 - ex_o1) <= 1e-12 * ex_scale) continue
                // bootstrap residuals at the stage-1 argmin -> Omega* -> W2*
                ex_rb = ex_Yb0 - gamma_cache[ex_okl[ex_b1]].dW * ex_th
                ex_Gt = xdpt2_gsum_by_unit(Z_r :* ex_rb, uid_r, n_u)
                ex_Gt[ex_cids, .] = ex_Gt[ex_cids, .] :- (ex_Shat' / ex_nc)
                ex_Om = ((ex_Gt :* ex_w)' * ex_Gt) / ex_nb
                if (xdpt_center == 1) {
                    // Center over the ex_nc resampled clusters, while keeping
                    // the command's per-bootstrap-row normalization.
                    ex_sb = colsum(ex_Gt :* ex_w)
                    ex_Om = ex_Om - (ex_sb' * ex_sb) / (ex_nc * ex_nb)
                }
                ex_Om = (ex_Om + ex_Om') / 2
                xdpt2_syminv(ex_Om, ex_invok, ex_W2)
                if (!ex_invok) continue
                // stage 2 restricted at gamma_l ...
                ex_ZWl = ex_Zw' * dW_r
                ex_A = ex_ZWl' * (ex_W2 * ex_ZWl)
                xdpt2_syminv(ex_A, ex_invok, ex_Ainv)
                if (!ex_invok) continue
                ex_th = ex_Ainv * (ex_ZWl' * (ex_W2 * (ex_ZYs - ex_Shat)))
                ex_m = ex_ZYs - ex_ZWl * ex_th - ex_Shat
                ex_v = ex_W2 * ex_m
                ex_objr = (ex_m' * ex_v) / ex_nb
                if (ex_objr >= . | hasmissing(ex_th)) continue
                // ... and the stage-2 grid minimum over ALL ok candidates
                // Include the restricted candidate in the unrestricted set;
                // CI points need not belong to the estimation grid.
                ex_min = .
                ex_b2 = 0
                ex_n2 = 0
                ex_hi2 = .
                for (ex_g2 = 1; ex_g2 <= rows(ex_okl); ex_g2++) {
                    ex_c0 = (ex_g2 - 1) * ex_k
                    ex_ZWg = ex_ZWall[|1, ex_c0 + 1 \ cols(Z_r), ex_c0 + ex_k|]
                    ex_A = ex_ZWg' * (ex_W2 * ex_ZWg)
                    xdpt2_syminv(ex_A, ex_invok, ex_Ainv)
                    if (!ex_invok) continue
                    ex_th1b = ex_Ainv * (ex_ZWg' * (ex_W2 * (ex_ZYs - ex_Shat)))
                    ex_m = ex_ZYs - ex_ZWg * ex_th1b - ex_Shat
                    ex_v = ex_W2 * ex_m
                    ex_obj = (ex_m' * ex_v) / ex_nb
                    if (ex_obj >= . | hasmissing(ex_th1b)) continue
                    ex_n2 = ex_n2 + 1
                    if (ex_hi2 >= . | ex_obj > ex_hi2) ex_hi2 = ex_obj
                    if (ex_min >= .) {
                        ex_min = ex_obj
                        ex_b2 = ex_g2
                    }
                    else {
                        ex_tol2 = xdpt2_objtol(ex_obj, ex_min, 1e-12)
                        if (ex_obj < ex_min - ex_tol2 |
                            (abs(ex_obj - ex_min) <= ex_tol2 &
                             gamma_grid[ex_okl[ex_g2]] <
                             gamma_grid[ex_okl[ex_b2]])) {
                            ex_min = ex_obj
                            ex_b2 = ex_g2
                        }
                    }
                }
                // A restricted-only solve used to create D_boot=0 and count
                // the draw as valid even when every unrestricted grid solve
                // failed. Mirror the reported estimator's searchable,
                // non-flat stage-2 profile gate before adding the null point.
                if (ex_b2 == 0 | ex_n2 < 2) continue
                ex_scale = max((abs(ex_min), abs(ex_hi2)))
                if (abs(ex_hi2 - ex_min) <= 1e-12 * ex_scale) continue
                ex_min = min((ex_objr, ex_min))
                D_boot = ex_objr - ex_min
                ex_valid = ex_valid + 1
                D_vec[ex_valid] = D_boot
            }
        }
        else if (use_batch_b & rows(fast_gb_b) > 0) {
            // ---- batched path (agrees with the scalar loop to rounding) ----
            // v0.9.34 (O1): the common Mammen weights (contributing clusters
            // only), drawn once in the same per-replication order as before
            if (rows(ETA_crn) != n_draw | cols(ETA_crn) != n_boot) {
                ETA_crn = J(n_draw, n_boot, 0)
                for (b = 1; b <= n_boot; b++) {
                    ETA_crn[., b] = xdpt2_mammen_draw(n_draw)
                }
            }
            ETA_b = ETA_crn
            // v0.9.31: Y* = F + E (restricted fit plus reweighted residuals),
            // evaluated by the split (xdpt2_fast_obj_split_list); v0.9.33:
            // E = resid_r :* ETA_b[uid_draw, .] is passed by its factors;
            // v0.9.34: with the fixed W2 for a two-step fit.
            F_b    = dW_r * theta_r
            if (w2m) objr_b = xdpt2_fast_obj_split_list(F_b, resid_r, uid_draw, ETA_b,
                                                        gamma_ci_cache, l, bA)
            else     objr_b = xdpt2_fast_obj_split_list(F_b, resid_r, uid_draw, ETA_b,
                                                        gamma_ci_cache, l)
            // v0.7.9 (D): preallocated stack (the old append recopied the
            // accumulator per gamma); values and row order identical
            OBJ_b  = J(1 + rows(fast_gb_b), n_boot, .)
            OBJ_b[1, .] = objr_b                       // restricted = initial min
            if (w2m) OBJ_b[|2, 1 \ 1 + rows(fast_gb_b), n_boot|] =
                xdpt2_fast_obj_split_list(F_b, resid_r, uid_draw, ETA_b,
                                          gamma_cache, fast_gb_b, bA)
            else OBJ_b[|2, 1 \ 1 + rows(fast_gb_b), n_boot|] =
                xdpt2_fast_obj_split_list(F_b, resid_r, uid_draw, ETA_b,
                                          gamma_cache, fast_gb_b)
            minobj_b = colmin(OBJ_b)                   // per-sample min over γ
            Dvec_b   = objr_b - minobj_b
            D_vec    = Dvec_b'
            // Require at least one finite threshold-model alternative in
            // each draw. The restricted row alone must not manufacture
            // D*=0 after every unrestricted solve overflowed/failed.
            altmin_b = colmin(OBJ_b[|2, 1 \ 1 + rows(fast_gb_b), n_boot|])
            minobj_b = colmin(objr_b \ altmin_b)
            Dvec_b   = objr_b - minobj_b
            D_vec    = J(n_boot, 1, .)
            for (b = 1; b <= n_boot; b++) {
                if (objr_b[b] >= . | altmin_b[b] >= . | Dvec_b[b] >= .) continue
                D_vec[b] = Dvec_b[b]
            }
        }
        else if (!w2m) {
            // ---- original scalar loop (fallback; one-step criterion) ----
            // v0.9.34: reached only when no alternative is batched (every
            // draw is then invalid); with W2 the draws stay missing (status 5)
            if (rows(ETA_crn) != n_draw | cols(ETA_crn) != n_boot) {
                ETA_crn = J(n_draw, n_boot, 0)
                for (b = 1; b <= n_boot; b++) {
                    ETA_crn[., b] = xdpt2_mammen_draw(n_draw)
                }
            }
            for (b = 1; b <= n_boot; b++) {
                // UNIT-LEVEL Mammen weights: the common draws (O1)
                eta_unit = ETA_crn[., b]
                eta = eta_unit[uid_draw]
                Y_boot = dW_r * theta_r + resid_r :* eta

                // Bootstrap restricted at γ_ℓ — fast 1-step GMM
                xdpt2_fast_gmm_boot(Y_boot, gamma_ci_cache[l],
                                     ok_b_r, theta_b_r, obj_b_r)
                // v0.7.13 (audit): the old 2-step fallback here was
                // unreachable (this cache entry already solved for the
                // sample with the same rows and weight) and would have mixed
                // a cluster-Omega objective into a W_first criterion. Skip
                // the replication defensively instead.
                if (!ok_b_r) continue

                // Bootstrap unrestricted: grid search (fast 1-step per γ)
                min_obj_b = obj_b_r
                has_alt_b = 0
                for (gb = 1; gb <= n1_ci; gb++) {
                    if (!gamma_cache[gb].ok) continue
                    if (gamma_cache[gb].n_rows != rows(Y_boot)) continue
                    xdpt2_fast_gmm_boot(Y_boot, gamma_cache[gb],
                                         ok_b_u, theta_b_u, obj_b_u)
                    // v0.7.13 (audit): no 2-step fallback here. The sample
                    // statistic's unrestricted scan is fast-path-only, so
                    // the bootstrap must search the SAME γ feasibility set
                    // under the SAME W_first criterion. The old fallback
                    // admitted extra γ under a cluster-Omega objective,
                    // deflating min_obj_b and inflating crit (over-wide CI).
                    if (!ok_b_u) continue
                    has_alt_b = 1
                    if (obj_b_u < min_obj_b) min_obj_b = obj_b_u
                }

                if (!has_alt_b) continue
                D_boot = obj_b_r - min_obj_b
                D_vec[b] = D_boot
            }
        }

        // v0.9.2 R18 (user): skipped/singular draws leave missing entries;
        // track the smallest per-point valid count (exported) and refuse to
        // invert on an incomplete set of draws (rule below).
        real scalar n_valid_b, min_valid_b
        n_valid_b = sum(D_vec :< .)
        if (gb_minB >= . | n_valid_b < gb_minB) gb_minB = n_valid_b
        // v0.9.4 R20 (#2): survivors of numerical failure are a SELECTED
        // subsample. Points below the floor are UNRESOLVED (status 5), not
        // rejected.
        // v0.9.33: every draw must be valid, as under boottype(unit); the
        // 90% floor of the wild inversion is gone. A wild draw fails only
        // on a nonfinite value, so fits with all draws valid are unchanged.
        min_valid_b = n_boot
        if (n_valid_b < min_valid_b) {
            accept[l] = 0
            ci_tab[l, 5] = n_valid_b
            ci_tab[l, 6] = 5
            continue
        }
        crit = xdpt2_crit_orderstat(D_vec, 1 - alpha)
        if (crit == .) {
            accept[l] = 0
            ci_tab[l, 5] = n_valid_b
            ci_tab[l, 6] = 6
            continue
        }
        accept[l] = (D_sample <= crit)
        ci_tab[l, 2] = D_sample
        // v0.9.34: no finite critical value when boot() is too small for the
        // level (every candidate is accepted); stored as missing
        ci_tab[l, 3] = (crit >= maxdouble() ? . : crit)
        ci_tab[l, 4] = accept[l]
        ci_tab[l, 5] = n_valid_b
        ci_tab[l, 6] = 1
        // v0.9.36: the bootstrap statistics of the (last) point, for the
        // p-value of citest()
        if (args() == 32) D_out = D_vec

        if (xdpt_verbose) {
            printf("    γ_ℓ=%6.4f  D_n=%7.3f  crit=%7.3f  %s\n",
                   gamma_ci_grid[l], D_sample, crit,
                   (accept[l] ? "accept" : "reject"))
        }
    }
    if (!xdpt_verbose) {
        printf(" done\n")
        displayflush()
    }

    // Status 3 is not admissible; status 4 failed the relevant sample
    // solve. Statuses 1/2/5/6 all passed that solve, regardless of
    // whether the subsequent bootstrap was sufficiently complete.
    for (gci_l = 1; gci_l <= n_ci; gci_l++) {
        if (ci_tab[gci_l, 6] == 3 | ci_tab[gci_l, 6] == 4) continue
        gci_adm = gci_adm + 1
        if (gci_lo == .) gci_lo = gamma_ci_grid[gci_l]
        gci_hi = gamma_ci_grid[gci_l]
    }

    // Convex hull of accepted γ.
    // v0.7.0 (B3): an empty acceptance set is REPORTED (ci_empty = 1, bounds
    // left missing) instead of being silently collapsed to the degenerate
    // point CI {γ̂}, which looked like ultra-precise inference when the test
    // inversion in fact rejected everywhere. A disconnected acceptance set is
    // flagged via ci_nseg > 1 (the hull is still returned for continuity).
    gam_lo = .
    gam_hi = .
    ci_nseg = 0
    real scalar prev_acc
    prev_acc = 0
    for (l = 1; l <= n_ci; l++) {
        if (accept[l] == 1) {
            if (gam_lo == .) gam_lo = gamma_ci_grid[l]
            gam_hi = gamma_ci_grid[l]
            if (!prev_acc) ci_nseg = ci_nseg + 1
            prev_acc = 1
        }
        else prev_acc = 0
    }
    ci_empty = (gam_lo == .)
    // v0.9.3 R19 (#8): explicit accepted segments [lower, upper] -- the
    // hull stays as a summary, this is the real confidence set.
    real scalar seg_lo
    seg_lo = .
    for (l = 1; l <= n_ci; l++) {
        if (accept[l] == 1) {
            if (seg_lo == .) seg_lo = gamma_ci_grid[l]
        }
        else {
            if (seg_lo < .) {
                ci_seg = ci_seg \ (seg_lo, gamma_ci_grid[l - 1])
                seg_lo = .
            }
        }
    }
    if (seg_lo < .) ci_seg = ci_seg \ (seg_lo, gamma_ci_grid[n_ci])
    // v0.9.4 R20 (#3): unresolved = numerical/bootstrap failure, never a
    // statistical rejection.
    ci_unres = 0
    for (l = 1; l <= n_ci; l++) {
        if (ci_tab[l, 6] == 4 | ci_tab[l, 6] == 5 | ci_tab[l, 6] == 6) {
            ci_unres = ci_unres + 1
        }
    }
    // v0.9.5 R21 (blocker): an INCOMPLETE inversion must not masquerade as
    // a complete inversion result -- the R20 warning still shipped hull bounds,
    // ci_empty, segment counts, and boundary diagnostics built as if the
    // unresolved points had been REJECTED (a set like {0.2} u {0.4} with
    // 0.3 unresolved understates the truth, and "rejected ALL candidates"
    // with everything unresolved is simply false). With any unresolved
    // point the hull, the empty flag, and the segment count are WITHDRAWN
    // (missing). The acceptance runs over the points that WERE evaluated
    // stay in ci_seg -- the ado stores them as e(ci_segments_evaluated),
    // explicitly labeled evaluated-only -- and the full table is in e(ci_grid).
    if (ci_unres > 0) {
        gam_lo = .
        gam_hi = .
        ci_empty = .
        ci_nseg = .
    }
}

// Continuity test (Gong-Seo 2026, §4.3 and Theorem 7):
//   H0: model is continuous (kink) vs H1: discontinuous (jump)
//   Test stat T_n = n·(Q̂_kink(θ̃) - Q̂_jump(θ̂)) on sample
//   Bootstrap p-value under kink DGP.
real scalar xdpt2_continuity_test(struct xdpt2_unit rowvector units,
                                    real colvector gamma_grid,
                                    real colvector q_supp, real scalar min_user,
                                    struct xdpt2_gamma_cache rowvector cache_jump,
                                    real scalar best_obj_jump,
                                    string scalar method,
                                    real scalar flag_static,
                                    real scalar t_min, real scalar t_max,
                                    real scalar n_boot,
                                    real scalar valid_out,
                                    real scalar common_out,
                                    real scalar bt2s, real matrix bA)
{
    real scalar g_kink, obj_kink, T_sample, T_boot_b, count_exceed, valid_boot
    real scalar b, gl, ok, obj_cur, n_u, u, i, ok_b, obj_kink_b, obj_jump_b
    real scalar gl_j, min_obj_kink_b, min_obj_jump_b
    real scalar K, ci_C, tol_nest, obj_kcand, obj_jmatch, ok_jmatch
    real colvector theta_kink_sample, theta_kcand, r_kink, Y_boot, eta, eta_unit, uid_draw
    real matrix dY_k, dW_k, Z_k, V_dummy
    real colvector times_k, uid_k, theta_cur, times_cur, uid_cur
    real matrix dY_cur, dW_cur, Z_cur, V_cur
    real matrix W_first
    // 1-step T_sample (Gong-Seo Alg. 1 consistency)
    real scalar best_k_1s, best_j_1s, ok_1s, gl_1s
    real colvector theta_1s_dummy
    theta_1s_dummy = J(0, 1, 0)
    valid_out = .
    common_out = 0

    // Step 1: compute sample T_n = obj_kink - obj_jump
    // v0.7.0 (D1): the kink cache is built ONCE here and reused by the kink
    // grid search, the sample statistic, and the bootstrap. The jump cache is
    // the main estimation cache passed in by the caller (the continuity test
    // only runs when the estimated model is the jump model, flag_kink = 0).
    // v0.9.34: both models are searched over the initial grid on the sample
    // side and in every draw (see xdpt2_grid_bootstrap); refine() points of
    // the jump fit are left out.
    external real scalar xdpt_n_stage1
    real scalar n1_t
    n1_t = min((rows(gamma_grid), cols(cache_jump)))
    if (xdpt_n_stage1 < .) n1_t = min((n1_t, xdpt_n_stage1))
    if (n1_t < 2) return(.)
    struct xdpt2_gamma_cache rowvector cache_kink
    struct xdpt2_stack_tpl scalar tpl_k
    real scalar tpl_k_st
    tpl_k_st = 0
    cache_kink = xdpt2_build_gamma_cache_t(units, gamma_grid[|1 \ n1_t|],
                                            method, flag_static, 1, t_min,
                                            t_max, q_supp, min_user, tpl_k,
                                            tpl_k_st, cache_jump)

    // v0.9.36 (power): with a two-step fit both models are compared on the
    // criterion of the reported estimator -- the fixed second-step weight W2
    // of the jump fit (bA), common to the kink and jump solves so that the
    // nesting (T >= 0) is kept -- as the threshold CI does since 0.9.34. The
    // one-step weight used up to 0.9.35 (the MA(1) H-matrix under FD) is
    // efficient only for homoskedastic serially uncorrelated errors; the
    // distance statistic under it has a noisier null distribution and lower
    // power against a jump. After a one-step fallback W1 is kept.
    real scalar w2c
    w2c = (bt2s == 1 & rows(bA) > 0)
    if (w2c) {
        xdpt2_cache_w2(cache_kink, bA)
        xdpt2_cache_w2(cache_jump, bA)
    }

    // The computational comparison is nested only where BOTH specifications
    // solve on the same row sample with the same one-step criterion. The jump
    // design has more columns and can fail rank/conditioning at gamma values
    // where the kink design succeeds. Selecting the restricted minimum from
    // those kink-only points and then clamping a negative distance to zero
    // silently turned a nonnested numerical comparison into a p-value.
    // Gong-Seo (2026, eq. in sec. 2 and Theorem 4): the continuity statistic
    // is the GMM distance under the efficient second-step weight. If that
    // solve leaves fewer than two jointly feasible points the test is not
    // reported; it is never recomputed under W1, a different statistic that
    // their theory does not cover (v0.9.37 review).
    real colvector common_C
    common_C = J(0, 1, 0)
    for (gl = 1; gl <= min((cols(cache_kink), cols(cache_jump))); gl++) {
        if (!cache_kink[gl].ok |
            (w2c ? cache_kink[gl].fast2_ok : cache_kink[gl].fast_ok) != 1) continue
        if (!cache_jump[gl].ok |
            (w2c ? cache_jump[gl].fast2_ok : cache_jump[gl].fast_ok) != 1) continue
        if (cache_kink[gl].n_rows != cache_jump[gl].n_rows) continue
        if (any(cache_kink[gl].uid :!= cache_jump[gl].uid)) continue
        if (any(cache_kink[gl].times :!= cache_jump[gl].times)) continue
        if (any(cache_kink[gl].dY :!= cache_jump[gl].dY)) continue
        common_C = common_C \ gl
    }
    // Select the restricted DGP with the SAME one-step objective used by both
    // the sample statistic and bootstrap draws. fast_ok certifies the normal
    // matrix only; require the kink and its matching jump solve to be finite
    // on the observed sample before calling a gamma jointly feasible.
    real scalar idx_k
    real colvector common_eval_C
    idx_k = 0
    best_k_1s = .
    common_eval_C = J(0, 1, 0)
    for (ci_C = 1; ci_C <= rows(common_C); ci_C++) {
        gl_1s = common_C[ci_C]
        xdpt2_fast_gmm_boot_w(cache_kink[gl_1s].dY, cache_kink[gl_1s],
                               w2c, bA, ok_1s, theta_kcand, obj_kcand)
        if (!ok_1s) continue
        xdpt2_fast_gmm_boot_w(cache_kink[gl_1s].dY, cache_jump[gl_1s],
                               w2c, bA, ok_jmatch, theta_1s_dummy, obj_jmatch)
        if (!ok_jmatch) continue
        common_eval_C = common_eval_C \ gl_1s
        // v0.9.14 R33 (#5): deterministic tie-break -- the selected kink
        // model SEEDS the whole continuity bootstrap DGP, so near-flat
        // profiles must not resolve by grid insertion order.
        if (best_k_1s == .) {
            best_k_1s = obj_kcand
            idx_k = gl_1s
            theta_kink_sample = theta_kcand
        }
        else {
            real scalar tol_k
            tol_k = xdpt2_objtol(obj_kcand, best_k_1s, 1e-12)
            if (obj_kcand < best_k_1s - tol_k |
                (abs(obj_kcand - best_k_1s) <= tol_k &
                 gamma_grid[gl_1s] < gamma_grid[idx_k])) {
                best_k_1s = obj_kcand
                idx_k = gl_1s
                theta_kink_sample = theta_kcand
            }
        }
    }
    common_C = common_eval_C
    common_out = rows(common_C)
    if (common_out < 2 | idx_k == 0) return(.)
    dY_k    = cache_kink[idx_k].dY
    dW_k    = cache_kink[idx_k].dW
    Z_k     = *cache_kink[idx_k].pZ
    times_k = cache_kink[idx_k].times
    uid_k   = cache_kink[idx_k].uid
    r_kink = dY_k - dW_k * theta_kink_sample

    uid_draw = xdpt2_dense_uid(uid_k)
    n_u = max(uid_draw)
    count_exceed = 0
    valid_boot = 0

    // v0.9.37: jump entries that may enter the comparison are those on the
    // kink rows (same uid, times and dY), not merely the same row count; the
    // jump residuals of the bootstrap DGP come from one of them.
    real colvector jalign
    jalign = J(n1_t, 1, 0)
    for (gl_1s = 1; gl_1s <= n1_t; gl_1s++) {
        if (!cache_jump[gl_1s].ok) continue
        if (cache_jump[gl_1s].n_rows != rows(dY_k)) continue
        if (any(cache_jump[gl_1s].uid :!= uid_k)) continue
        if (any(cache_jump[gl_1s].times :!= times_k)) continue
        if (any(cache_jump[gl_1s].dY :!= dY_k)) continue
        jalign[gl_1s] = 1
    }

    best_j_1s = .
    real scalar idx_j
    real colvector theta_jump_sample
    idx_j = 0
    for (gl_1s = 1; gl_1s <= n1_t; gl_1s++) {
        if (!jalign[gl_1s]) continue
        xdpt2_fast_gmm_boot_w(dY_k, cache_jump[gl_1s], w2c, bA,
                               ok_1s, theta_1s_dummy, obj_cur)
        if (!ok_1s) continue
        if (best_j_1s == . | obj_cur < best_j_1s) {
            best_j_1s = obj_cur
            idx_j = gl_1s
            theta_jump_sample = theta_1s_dummy
        }
    }
    if (best_k_1s == . | best_j_1s == . | idx_j == 0) return(.)
    // v0.9.36 (power): the bootstrap DGP is the restricted (kink) fit plus
    // the UNRESTRICTED (jump) residuals, reweighted by cluster: restricted
    // coefficients, unrestricted residuals, as in Gong-Seo Alg. 1 and the
    // unit bootstrap of the CI. Under H0 both residual vectors estimate the
    // same errors; under H1 the kink residuals also carry the omitted jump,
    // delta*(1(q>gamma0) - kink fit), which inflated every bootstrap
    // statistic and with it the critical value -- the continuity test lost
    // power exactly where it should reject. dY and the row sample are the
    // same for both models (rows aligned by jalign above).
    r_kink = dY_k - cache_jump[idx_j].dW * theta_jump_sample
    if (hasmissing(r_kink)) return(.)
    T_sample = best_k_1s - best_j_1s
    tol_nest = xdpt2_objtol(best_k_1s, best_j_1s, 1e-10)
    if (T_sample < -tol_nest) return(.)
    if (T_sample < 0) T_sample = 0

    external real scalar xdpt_verbose
    if (xdpt_verbose) printf("  Continuity test (H0: kink, B=%g, unit-level Mammen; cached)...\n", n_boot)
    else {
        printf("  Continuity test  (. per bootstrap, %g total)\n  ", n_boot)
        displayflush()
    }

    // v0.7.7 SPEEDUP: batch all B replications when every participating γ in
    // BOTH the kink and jump caches uses the fast path. Same Mammen draw order
    // -> the scalar loop's p-value up to rounding; else original scalar loop.
    // v0.9.33: a jump γ without the one-step solve (fast_ok = 0) is left out
    // of the batch, as the scalar loop skips it; it no longer sends every
    // draw to the scalar loop. The common set holds fast entries only.
    real scalar use_batch_C, n_rows_k
    real colvector fast_k_C, fast_j_C
    real matrix ETA_C, OBJk_C, OBJj_C, OBJkc_C
    real colvector F_C
    real rowvector mink_C, minj_C, T_C
    real colvector jump_pos_C, pos_C
    n_rows_k = rows(dY_k)
    use_batch_C = 1
    // Restricted searches stay on the jointly feasible set. The jump
    // alternative may use its full feasible set; it necessarily contains
    // the matching jump model at every common gamma.
    fast_k_C = common_C
    fast_j_C = J(0, 1, 0)
    if (use_batch_C) {
        for (gl_j = 1; gl_j <= n1_t; gl_j++) {
            if (!jalign[gl_j]) continue
            if ((w2c ? cache_jump[gl_j].fast2_ok : cache_jump[gl_j].fast_ok) != 1) continue
            fast_j_C = fast_j_C \ gl_j
        }
    }
    jump_pos_C = J(rows(common_C), 1, 0)
    if (use_batch_C) {
        for (ci_C = 1; ci_C <= rows(common_C); ci_C++) {
            pos_C = selectindex(fast_j_C :== common_C[ci_C])
            if (rows(pos_C) != 1) {
                use_batch_C = 0
                break
            }
            jump_pos_C[ci_C] = pos_C[1]
        }
    }
    if (use_batch_C & rows(fast_k_C) > 0 & rows(fast_j_C) > 0) {
        ETA_C = J(n_u, n_boot, 0)
        for (b = 1; b <= n_boot; b++) ETA_C[., b] = xdpt2_mammen_draw(n_u)
        // v0.9.31: Y* = F + E, evaluated by the split; v0.9.33: E is
        // passed by its factors, r_kink and ETA_C[uid_draw, .]
        F_C = dW_k * theta_kink_sample
        // v0.7.9 (D): preallocated stacks; values and row order identical
        if (w2c) {
            OBJk_C = xdpt2_fast_obj_split_list(F_C, r_kink, uid_draw, ETA_C,
                                               cache_kink, fast_k_C, bA)
            OBJj_C = xdpt2_fast_obj_split_list(F_C, r_kink, uid_draw, ETA_C,
                                               cache_jump, fast_j_C, bA)
        }
        else {
            OBJk_C = xdpt2_fast_obj_split_list(F_C, r_kink, uid_draw, ETA_C,
                                               cache_kink, fast_k_C)
            OBJj_C = xdpt2_fast_obj_split_list(F_C, r_kink, uid_draw, ETA_C,
                                               cache_jump, fast_j_C)
        }
        // The restricted minimum may use a gamma only when the matching
        // unrestricted jump solve is finite in that draw. The unrestricted
        // minimum itself still uses its full feasible set.
        OBJkc_C = OBJk_C
        for (ci_C = 1; ci_C <= rows(common_C); ci_C++) {
            for (b = 1; b <= n_boot; b++) {
                if (OBJj_C[jump_pos_C[ci_C], b] >= .) OBJkc_C[ci_C, b] = .
            }
        }
        mink_C = colmin(OBJkc_C)
        minj_C = colmin(OBJj_C)
        T_C = mink_C - minj_C
        // A materially negative distance means numerical nesting failed for
        // that draw. Exclude it; clamp only roundoff-sized negatives.
        valid_boot = 0
        count_exceed = 0
        for (b = 1; b <= n_boot; b++) {
            if (mink_C[b] >= . | minj_C[b] >= . | T_C[b] >= .) continue
            tol_nest = xdpt2_objtol(mink_C[b], minj_C[b], 1e-10)
            if (T_C[b] < -tol_nest) continue
            T_boot_b = T_C[b]
            if (T_boot_b < 0) T_boot_b = 0
            valid_boot = valid_boot + 1
            if (T_boot_b >= T_sample) count_exceed = count_exceed + 1
        }
    }
    else {
        for (b = 1; b <= n_boot; b++) {
            if (!xdpt_verbose) {
                if (mod(b, 50) == 0) printf("+ %g\n  ", b)
                else                 printf(".")
                displayflush()
            }
            eta_unit = xdpt2_mammen_draw(n_u)
            eta = eta_unit[uid_draw]
            Y_boot = dW_k * theta_kink_sample + r_kink :* eta

            min_obj_kink_b = .
            for (ci_C = 1; ci_C <= rows(common_C); ci_C++) {
                gl = common_C[ci_C]
                xdpt2_fast_gmm_boot_w(Y_boot, cache_kink[gl], w2c, bA,
                                       ok_b, theta_cur, obj_kcand)
                if (!ok_b) continue
                xdpt2_fast_gmm_boot_w(Y_boot, cache_jump[gl], w2c, bA,
                                       ok_jmatch, theta_cur, obj_jmatch)
                if (!ok_jmatch) continue
                if (obj_kcand < min_obj_kink_b) min_obj_kink_b = obj_kcand
            }

            min_obj_jump_b = .
            for (gl_j = 1; gl_j <= n1_t; gl_j++) {
                if (!jalign[gl_j]) continue
                xdpt2_fast_gmm_boot_w(Y_boot, cache_jump[gl_j], w2c, bA,
                                       ok_b, theta_cur, obj_cur)
                if (!ok_b) continue
                if (obj_cur < min_obj_jump_b) min_obj_jump_b = obj_cur
            }

            if (min_obj_kink_b == . | min_obj_jump_b == .) continue
            T_boot_b = min_obj_kink_b - min_obj_jump_b
            tol_nest = xdpt2_objtol(min_obj_kink_b, min_obj_jump_b, 1e-10)
            if (T_boot_b < -tol_nest) continue
            if (T_boot_b < 0) T_boot_b = 0
            valid_boot = valid_boot + 1
            if (T_boot_b >= T_sample) count_exceed = count_exceed + 1
        }
    }
    if (!xdpt_verbose) {
        printf(" done\n")
        displayflush()
    }

    // v0.9.3 R19 (#7): a p-value from a handful of surviving draws is
    // noise. v0.9.33: every draw must be valid (was >= 90% and >= 10), so
    // the p-value is never computed from a subsample selected by failures.
    // v0.9.4 R20 (#6): the count is returned either way, so the display
    // can SAY why p is missing instead of printing a bare dot.
    valid_out = valid_boot
    if (valid_boot < n_boot) return(.)
    // v0.7.0 (B2): add-one correction (Davidson-MacKinnon 2000) — a valid
    // bootstrap p-value is never exactly zero.
    return((1 + count_exceed) / (1 + valid_boot))
}

// Linearity test (H0: no regime, δ=0). Wild bootstrap.
real scalar xdpt2_linearity_test(struct xdpt2_unit rowvector units,
                                  real colvector gamma_grid,
                                  struct xdpt2_gamma_cache rowvector gamma_cache,
                                  real scalar best_obj,
                                  string scalar method, real scalar flag_static,
                                  real scalar flag_kink,
                                  real scalar t_min, real scalar t_max,
                                  real scalar n_boot,
                                  real scalar valid_out)
{
    // Restricted: no regime — W has only K (β) cols; γ irrelevant
    real scalar K, ok, obj_r, n_rows, b, ok_b_r, obj_b_r, ok_b_u, obj_b_u
    real scalar gl_b, min_obj_b_u, supW_sample, supW_b_s, has_alt_b
    real scalar count_exceed, valid_boot
    real scalar n_u, u
    real matrix dY_s, dW_s, Z_s, W_beta_s, V_r, V_dummy
    real colvector times_s, theta_r, fit_r, resid_r, Y_boot, theta_b_r, theta_b_u
    real colvector eta_unit, eta, uid_s, uid_b, uid_draw
    real matrix dY_b, dW_b, Z_b
    real colvector times_b
    // 1-step sample obj (Gong-Seo Alg. 1 consistency)
    real matrix W_first_lin
    real scalar obj_r_1s, gl_s, obj_u_1s, ok_1s, min_obj_u_1s, has_alt_1s
    real colvector theta_1s_dummy
    theta_1s_dummy = J(0, 1, 0)
    // v0.9.34: the threshold model is searched over the initial grid on the
    // sample side and in every draw (see xdpt2_grid_bootstrap)
    external real scalar xdpt_n_stage1
    real scalar n1_l

    K = cols(units[1].X)
    n1_l = (xdpt_n_stage1 < . ? min((xdpt_n_stage1, cols(gamma_cache))) :
                                cols(gamma_cache))

    // v0.7.0 (D1): the per-γ cache is passed in (built once in
    // xtdpthresh_run). The restricted (no-regime) model uses only the
    // γ-invariant β columns, so any valid cache entry supplies the stacked
    // sample — no fresh stack_at_gamma / cache build needed here.
    real scalar idx_base
    idx_base = 0
    for (gl_s = 1; gl_s <= cols(gamma_cache); gl_s++) {
        if (gamma_cache[gl_s].ok) {
            idx_base = gl_s
            gl_s = cols(gamma_cache) + 1
        }
    }
    if (idx_base == 0) return(.)

    dY_s    = gamma_cache[idx_base].dY
    dW_s    = gamma_cache[idx_base].dW
    Z_s     = *gamma_cache[idx_base].pZ
    times_s = gamma_cache[idx_base].times
    uid_s   = gamma_cache[idx_base].uid
    uid_draw = xdpt2_dense_uid(uid_s)
    n_u = max(uid_draw)
    if (rows(dY_s) < 20) return(.)

    // Restricted regressors: the base β columns.
    W_beta_s = J(rows(dW_s), 0, 0)
    if (K > 0) W_beta_s = dW_s[., 1..K]

    // BUG 4b FIX: use 1-step theta_r with W_first (consistent with bootstrap).
    // Previously theta_r came from 2-step solve_gmm while supW uses 1-step obj,
    // creating a scale inconsistency between sample and bootstrap statistics.
    // Both sample and bootstrap now use 1-step with the same W_first weight.
    W_first_lin = *gamma_cache[idx_base].pW1
    if (cols(W_beta_s) == 0) {
        // Static with no RHS but valid external IVs has a zero-parameter
        // linear null; its GMM objective is still well defined.
        real colvector g_r0
        theta_r = J(0, 1, 0)
        V_r = J(0, 0, .)
        g_r0 = Z_s' * dY_s / rows(dY_s)
        obj_r_1s = rows(dY_s) * (g_r0' * W_first_lin * g_r0)
        ok = (obj_r_1s < .)
    }
    else {
        xdpt2_solve_gmm_1step(dY_s, W_beta_s, Z_s, W_first_lin,
                               ok, theta_r, obj_r_1s, V_r)
    }
    if (!ok) return(.)

    // Bootstrap under H0 (unit-level Mammen) — uses 1-step theta_r
    fit_r = (cols(W_beta_s) == 0 ? J(rows(dY_s), 1, 0) : W_beta_s * theta_r)
    resid_r = dY_s - fit_r
    count_exceed = 0
    valid_boot = 0

    // The unrestricted set contains the restricted model. Start from its
    // objective so the distance is nonnegative by construction, but still
    // require at least one numerically valid threshold-model solve.
    min_obj_u_1s = obj_r_1s
    has_alt_1s = 0
    for (gl_s = 1; gl_s <= n1_l; gl_s++) {
        if (!gamma_cache[gl_s].ok) continue
        if (gamma_cache[gl_s].n_rows != rows(dY_s)) continue
        xdpt2_fast_gmm_boot(dY_s, gamma_cache[gl_s],
                             ok_1s, theta_1s_dummy, obj_u_1s)
        if (!ok_1s) continue
        has_alt_1s = 1
        if (obj_u_1s < min_obj_u_1s) min_obj_u_1s = obj_u_1s
    }
    if (!has_alt_1s) return(.)
    supW_sample = obj_r_1s - min_obj_u_1s

    // Precompute C_beta for fast restricted bootstrap
    // (dY_s/Z_s/times_s/uid_s come from gamma_cache[idx_base]; Z and the row
    //  sample are γ-invariant, so any valid entry is equivalent)
    real matrix W_first_s, ZW_beta, A_beta, Ainv_beta, C_beta
    real scalar fast_r_ok, n_rows_s, inv_ok_beta
    n_rows_s = rows(dY_s)
    W_first_s = *gamma_cache[idx_base].pW1
    if (cols(W_beta_s) == 0) {
        ZW_beta = J(cols(Z_s), 0, 0)
        A_beta = J(0, 0, 0)
        C_beta = J(0, cols(Z_s), 0)
        fast_r_ok = 1
    }
    else {
        ZW_beta = Z_s' * W_beta_s / n_rows_s
        A_beta = ZW_beta' * W_first_s * ZW_beta
        xdpt2_syminv(A_beta, inv_ok_beta, Ainv_beta)
        fast_r_ok = inv_ok_beta
        if (fast_r_ok) {
            C_beta = Ainv_beta * ZW_beta' * W_first_s
            if (hasmissing(C_beta)) fast_r_ok = 0
        }
    }

    external real scalar xdpt_verbose
    if (xdpt_verbose) printf("  Linearity test (H0: δ=0, B=%g, unit-level Mammen; cached)...\n", n_boot)
    else {
        printf("  Linearity test   (. per bootstrap, %g total)\n  ", n_boot)
        displayflush()
    }

    real colvector ZY_b, r_b, g_b
    // v0.7.7 SPEEDUP: batch all B replications when the restricted (C_beta)
    // solve is fast. Mammen draws keep the same per-replication order -> the
    // scalar loop's p-value up to rounding. Otherwise fall back to the
    // original scalar loop.
    // v0.9.33: an unrestricted γ without the one-step solve (fast_ok = 0) is
    // left out of the batch, as the scalar loop skips it; it no longer sends
    // every draw to the scalar loop.
    real scalar use_batch_L
    real colvector fast_gl_L
    real matrix ETA_L, OBJu_L, G0_L
    real rowvector objr_L, altminu_L, minu_L, supW_L
    use_batch_L = fast_r_ok
    fast_gl_L = J(0, 1, 0)
    if (use_batch_L) {
        for (gl_b = 1; gl_b <= n1_l; gl_b++) {
            if (!gamma_cache[gl_b].ok) continue
            if (gamma_cache[gl_b].n_rows != n_rows_s) continue
            if (gamma_cache[gl_b].fast_ok != 1) continue
            fast_gl_L = fast_gl_L \ gl_b
        }
    }
    if (use_batch_L & rows(fast_gl_L) > 0) {
        ETA_L = J(n_u, n_boot, 0)
        for (b = 1; b <= n_boot; b++) ETA_L[., b] = xdpt2_mammen_draw(n_u)
        // v0.9.34: no n x B matrix. The restricted (linear) objectives go
        // through the split with a one-entry cache for the linear model;
        // with no regressor the fit is zero and g = Z'E/n = S'ETA/n.
        if (cols(W_beta_s) == 0) {
            G0_L = xdpt2_gsum_by_unit(Z_s :* resid_r, uid_draw, n_u)' * ETA_L / n_rows_s
            // v0.9.35 (R3): a draw with a missing term has a missing objective
            real matrix P_L
            real rowvector bad_L
            P_L = G0_L :* (W_first_s * G0_L)
            objr_L = n_rows_s :* colsum(P_L)
            if (hasmissing(P_L)) {
                bad_L = selectindex(colmissing(P_L) :> 0)
                objr_L[bad_L] = J(1, cols(bad_L), .)
            }
        }
        else {
            struct xdpt2_gamma_cache rowvector lin_c
            lin_c = xdpt2_gamma_cache(1, 1)
            lin_c[1].ok      = 1
            lin_c[1].fast_ok = 1
            lin_c[1].n_rows  = n_rows_s
            lin_c[1].dW      = W_beta_s
            lin_c[1].C_g     = C_beta
            lin_c[1].ZW      = ZW_beta
            lin_c[1].pZ      = gamma_cache[idx_base].pZ
            lin_c[1].pW1     = gamma_cache[idx_base].pW1
            objr_L = xdpt2_fast_obj_split_list(fit_r, resid_r, uid_draw, ETA_L,
                                               lin_c, 1)
        }
        // v0.7.9 (D): preallocated stack; values and row order identical
        OBJu_L = J(1 + rows(fast_gl_L), n_boot, .)
        OBJu_L[1, .] = objr_L                        // restricted = initial min
        // v0.9.31: objectives by the split (xdpt2_fast_obj_split_list)
        OBJu_L[|2, 1 \ 1 + rows(fast_gl_L), n_boot|] =
            xdpt2_fast_obj_split_list(fit_r, resid_r, uid_draw, ETA_L,
                                      gamma_cache, fast_gl_L)
        minu_L = colmin(OBJu_L)
        // Alternative-only minimum: the restricted row cannot stand in for
        // an unrestricted threshold solve when every alternative is nonfinite.
        altminu_L = colmin(OBJu_L[|2, 1 \ 1 + rows(fast_gl_L), n_boot|])
        minu_L = colmin(objr_L \ altminu_L)
        supW_L = objr_L - minu_L
        valid_boot = 0
        count_exceed = 0
        for (b = 1; b <= n_boot; b++) {
            if (objr_L[b] >= . | altminu_L[b] >= . |
                minu_L[b] >= . | supW_L[b] >= .) continue
            valid_boot = valid_boot + 1
            if (supW_L[b] >= supW_sample) count_exceed = count_exceed + 1
        }
    }
    else {
        for (b = 1; b <= n_boot; b++) {
            if (!xdpt_verbose) {
                if (mod(b, 50) == 0) printf("+ %g\n  ", b)
                else                 printf(".")
                displayflush()
            }
            eta_unit = xdpt2_mammen_draw(n_u)
            eta = eta_unit[uid_draw]
            Y_boot = fit_r + resid_r :* eta

            if (cols(W_beta_s) == 0) {
                g_b = Z_s' * Y_boot / n_rows_s
                obj_b_r = n_rows_s * (g_b' * W_first_s * g_b)
                ok_b_r = 1
            }
            else if (fast_r_ok) {
                ZY_b = Z_s' * Y_boot / n_rows_s
                theta_b_r = C_beta * ZY_b
                r_b = Y_boot - W_beta_s * theta_b_r
                g_b = Z_s' * r_b / n_rows_s
                obj_b_r = n_rows_s * (g_b' * W_first_s * g_b)
                ok_b_r = 1
            }
            else {
                // v0.7.13 (audit): unreachable — fast_r_ok tests the same
                // matrix the pre-loop sample solve already required, so if
                // that solve succeeded fast_r_ok==1. The old 2-step-solve
                // fallback would have mixed a cluster-Omega objective into
                // the W_first supW comparison; skip defensively instead.
                continue
            }
            if (obj_b_r >= .) continue

            min_obj_b_u = obj_b_r
            has_alt_b = 0
            for (gl_b = 1; gl_b <= n1_l; gl_b++) {
                if (!gamma_cache[gl_b].ok) continue
                if (gamma_cache[gl_b].n_rows != rows(Y_boot)) continue
                xdpt2_fast_gmm_boot(Y_boot, gamma_cache[gl_b],
                                     ok_b_u, theta_b_u, obj_b_u)
                if (!ok_b_u) continue
                has_alt_b = 1
                if (obj_b_u < min_obj_b_u) min_obj_b_u = obj_b_u
            }

            if (!has_alt_b) continue
            supW_b_s = obj_b_r - min_obj_b_u
            valid_boot = valid_boot + 1
            if (supW_b_s >= supW_sample) count_exceed = count_exceed + 1
        }
    }
    if (!xdpt_verbose) {
        printf(" done\n")
        displayflush()
    }

    // v0.9.3 R19 (#7): same validity rule as the linearity test (v0.9.33:
    // every draw must be valid).
    valid_out = valid_boot
    if (valid_boot < n_boot) return(.)
    // v0.7.0 (B2): add-one correction (Davidson-MacKinnon 2000).
    return((1 + count_exceed) / (1 + valid_boot))
}

// Hansen J over-identification test.
// Under 2-step GMM with efficient weight Ω̂⁻¹: J = n·g'Ω̂⁻¹g ~ χ²(df)
// df = # instruments - # parameters. Stored obj IS J under 2-step.
real rowvector xdpt2_hansen_j(real scalar obj, real scalar n_iv, real scalar k_W)
{
    real scalar df, pval
    real rowvector out_miss
    out_miss = (., ., .)
    // v0.8.0 (audit R5): gamma is an estimated parameter too -- the regular
    // jump/kink estimator has k_W + 1 parameters, so df = L - k_W - 1. Under
    // continuity the chi-square reference is itself only a diagnostic (the
    // Jacobian degenerates; see help).
    df = n_iv - k_W - 1
    if (df <= 0) return(out_miss)
    // Use the survival function directly. 1-chi2() cancels to exact zero
    // in the far right tail even while the representable p-value is positive.
    pval = chi2tail(df, obj)
    return((obj, df, pval))
}

// v0.7.12: removed dead helper xdpt2_recompute_cluster_j(). It was never
// called on any live path -- Hansen J is now taken from the paired best_obj
// when best_twostep (xdpt2_hansen_j at the run level). Kept out to avoid
// maintenance confusion (the old "BUG 1 FIX" comment described a path that no
// longer exists).

// Arellano-Bond AR(k) test for serial correlation in transformed residuals.
// m_k = Σ_i e_i / sqrt(Σ_i e_i^2), where e_i = Σ_t r_{it} · r_{i,t-k}
// Under H0 (no k-th order autocorrelation): m_k ~ N(0, 1)
// For FD: expect reject at k=1 (Δε has MA(1)); fail to reject at k=2 is GOOD.
// Full Arellano-Bond (1991, eq. 8) m_k statistic — v0.7.2 (B1 fix).
//   m_k = b0 / sqrt(T1 + T2 + T3)
//   b0 = ẽ_{-k}' ê          (lag-aligned residual cross product, raw sum)
//   T1 = Σ_i (ẽ_{i,-k}' ê_i)²
//   T2 = -2 g' (G A G')^{-1} G A s   — covariance with θ̂ estimation error
//   T3 = g' V̂(θ̂) g                   — direct θ̂-variance contribution
// where g = X_t' ẽ_{-k}  (∂ê_test/∂θ direction, zero-padded to dim(θ̂)),
//       G = X_s' Z_s  (raw sums, estimation equation),
//       s = Σ_i (Z_{s,i}' ê_{s,i}) · c_i  with c_i = ẽ_{i,-k}' ê_i,
//       A = the moment weight actually paired with θ̂ (exported by
//           xdpt2_grid_search), V̂ = the reported variance of θ̂.
// ẽ_{-k} is zero where the lag-k residual is unobserved — the AB trimming
// convention. T2/T3 use the ESTIMATION-equation pieces (which for
// method(fod) differs from the FD test residuals; the c_i scalars link
// the two within unit). If the full variance pieces are unavailable or
// nonpositive, the statistic and p-value are missing; a negative pair count
// records that diagnostic state. No simplified T1-only statistic is reported.
real rowvector xdpt2_ar_full(real scalar k,
                              real colvector e_t, real matrix X_t,
                              real colvector times_t, real colvector uid_t,
                              real colvector e_s, real matrix Z_s,
                              real colvector uid_s, real matrix X_s,
                              real matrix A, real matrix V)
{
    real scalar i, j, jj, t_jk, found, n_pairs, n_pairs_i, n_pair_units, n_units_t
    real scalar b0, T1, T2, T3, c_u, mk, pval, denom2, full_ok, k_par
    real scalar fast_t, fast_s, jlo, jhi, nt, ns
    real colvector rows_i, times_i, r_i, elk, uniq_t, g, g_pad, s, c_row
    real matrix G_s, GAG, B_map, info_t, info_s
    real rowvector out_miss
    transmorphic cmap

    // 7 elements: (mk, np, pval, b0, T1, TT, pair-cluster count).
    // callers read [4]..[6] unconditionally, so a 3-element failure return
    // crashed the whole command with Mata 3301 whenever the AR test could
    // not be computed (e.g. AR(2) with zero lag-2 pairs on T=5 panels).
    out_miss = (., ., ., ., ., ., .)
    if (rows(e_t) < 2) return((., ., ., ., ., ., 0))

    // --- Pass 1 (test equation): lag-aligned ẽ_{-k}, per-unit c_i, b0, T1 ---
    elk = J(rows(e_t), 1, 0)
    uniq_t = uniqrows(uid_t)
    n_units_t = rows(uniq_t)
    // v0.9.35: every stack is sorted by unit and, within unit, by time. The
    // rows of a unit are then one block (panelsetup, in the order of uniq_t),
    // and the lag-k partner of a row, if any, is among the k rows before it,
    // times being increasing integers. The rows, their order, and the
    // arithmetic are those of the scans, which remain for any other order.
    nt = rows(uid_t)
    fast_t = (rows(times_t) == nt & nt >= 2)
    if (fast_t) {
        fast_t = all((uid_t[|2 \ nt|] :> uid_t[|1 \ nt - 1|]) :|
                     ((uid_t[|2 \ nt|] :== uid_t[|1 \ nt - 1|]) :&
                      (times_t[|2 \ nt|] :> times_t[|1 \ nt - 1|])))
    }
    if (fast_t) {
        fast_t = all(times_t :== trunc(times_t))
    }
    if (fast_t) {
        info_t = panelsetup(uid_t, 1)
        fast_t = (rows(info_t) == n_units_t)
    }
    cmap = asarray_create("real", 1)
    asarray_notfound(cmap, 0)
    b0 = 0
    T1 = 0
    n_pairs = 0
    n_pair_units = 0
    for (i = 1; i <= n_units_t; i++) {
        if (fast_t) rows_i = (info_t[i, 1]::info_t[i, 2])
        else        rows_i = selectindex(uid_t :== uniq_t[i])
        if (rows(rows_i) < 1) continue
        times_i = times_t[rows_i]
        r_i = e_t[rows_i]
        c_u = 0
        n_pairs_i = 0
        for (j = 1; j <= rows(rows_i); j++) {
            t_jk = times_i[j] - k
            found = 0
            jlo = 1
            jhi = rows(rows_i)
            if (fast_t) {
                jlo = max((1, j - k))
                jhi = j - 1
            }
            for (jj = jlo; jj <= jhi; jj++) {
                if (times_i[jj] == t_jk) {
                    found = jj
                    break
                }
            }
            if (found > 0) {
                elk[rows_i[j]] = r_i[found]
                c_u = c_u + r_i[j] * r_i[found]
                n_pairs = n_pairs + 1
                n_pairs_i = n_pairs_i + 1
            }
        }
        if (n_pairs_i > 0) n_pair_units = n_pair_units + 1
        b0 = b0 + c_u
        T1 = T1 + c_u^2
        asarray(cmap, uniq_t[i], c_u)
    }
    // The N(0,1) reference is cluster-asymptotic. Five pairs contributed by
    // one long panel are not five independent pieces of information.
    if (T1 <= 0 | n_pairs < 5 | n_pair_units < 5) {
        return((., ., ., b0, T1, ., n_pair_units))
    }

    // --- Estimator-dependent variance pieces (T2, T3) ---
    // Requires: X_t aligned to e_t rows; estimation-equation residuals e_s,
    // instruments Z_s, regressors X_s; weight A; variance V.
    T2 = 0
    T3 = 0
    full_ok = 0
    k_par = rows(V)
    if (rows(X_t) == rows(e_t) & rows(e_s) > 0 &
        rows(Z_s) == rows(e_s) & rows(X_s) == rows(e_s) &
        rows(A) == cols(Z_s) & cols(A) == cols(Z_s) &
        k_par == cols(X_s) & cols(X_t) <= k_par) {

        // g = X_t' ẽ_{-k}, zero-padded to dim(θ̂) when needed.
        g = X_t' * elk
        if (rows(g) < k_par) g_pad = g \ J(k_par - rows(g), 1, 0)
        else                 g_pad = g

        // s = Σ_rows Z_s[r,.]' ê_s[r] · c(uid_s[r])  — units absent from the
        // test stack contribute c = 0 (asarray notfound default)
        s = J(cols(Z_s), 1, 0)
        c_row = J(rows(e_s), 1, 0)
        // v0.9.35: one lookup per unit block when uid_s is sorted
        ns = rows(uid_s)
        fast_s = (ns >= 2)
        if (fast_s) fast_s = all(uid_s[|2 \ ns|] :>= uid_s[|1 \ ns - 1|])
        if (fast_s) {
            info_s = panelsetup(uid_s, 1)
            for (i = 1; i <= rows(info_s); i++) {
                c_row[|info_s[i, 1] \ info_s[i, 2]|] =
                    J(info_s[i, 2] - info_s[i, 1] + 1, 1,
                      asarray(cmap, uid_s[info_s[i, 1]]))
            }
        }
        else {
            for (j = 1; j <= rows(e_s); j++) {
                c_row[j] = asarray(cmap, uid_s[j])
            }
        }
        s = Z_s' * (e_s :* c_row)

        G_s = X_s' * Z_s                       // k_par × k_iv (raw sums)
        GAG = G_s * A * G_s'
        real scalar inv_ok_ar
        real matrix GAGinv
        xdpt2_syminv(GAG, inv_ok_ar, GAGinv)
        if (inv_ok_ar) {
            B_map = GAGinv * G_s * A           // k_par × k_iv
            T2 = -2 * (g_pad' * B_map * s)
            T3 = g_pad' * V * g_pad
            if (T2 < . & T3 < .) full_ok = 1
        }
    }

    // Do not silently replace the Arellano-Bond variance with b0/sqrt(T1).
    // Dropping the parameter-estimation terms is a different statistic and
    // does not justify the displayed N(0,1) p-value.
    if (!full_ok) return((., -n_pairs, ., b0, T1, ., n_pair_units))
    denom2 = T1 + T2 + T3
    if (denom2 >= . | denom2 <= 0) {
        return((., -n_pairs, ., b0, T1, T2 + T3, n_pair_units))
    }
    mk = b0 / sqrt(denom2)
    // Direct lower-tail evaluation avoids catastrophic cancellation for
    // large |m| (1-normal(|m|) rounds to zero around eight standard errors).
    pval = 2 * normal(-abs(mk))
    // v0.7.2-d5: expose b0, T1, T2+T3 for external verification
    return((mk, n_pairs, pval, b0, T1, T2 + T3, n_pair_units))
}

// v0.9.29: for each variable in -vars-, 1 if its value is identical within
// every group of -grp- over the rows marked by -touse- (groups with a single
// row carry no variation and are ignored), 0 otherwise. With grp = panel
// variable this detects regressors removed exactly by FD/FOD; with grp = time
// variable it detects regressors common to all units in every period, which
// td removes. Exact comparisons: FOD of such a column would otherwise leave
// rounding noise (Mata's mean of k equal doubles is not always exact) that
// the estimator would fit.
real rowvector xdpt2_const_within(string scalar vars, string scalar grp,
                                  string scalar touse)
{
    real matrix X, info
    real colvector g, ord
    real rowvector out
    real scalar j, i, r0, r1
    X = st_data(., vars, touse)
    g = st_data(., grp, touse)
    out = J(1, cols(X), 1)
    if (rows(X) < 2) return(out)
    ord = order(g, 1)
    X = X[ord, .]
    g = g[ord]
    info = panelsetup(g, 1)
    for (j = 1; j <= cols(X); j++) {
        for (i = 1; i <= rows(info); i++) {
            r0 = info[i, 1]
            r1 = info[i, 2]
            if (r1 <= r0) continue
            if (max(X[|r0, j \ r1, j|]) != min(X[|r0, j \ r1, j|])) {
                out[j] = 0
                break
            }
        }
    }
    return(out)
}

// v0.9.33 (an identical copy lives in xtdpthresh_p.ado, whose Mata block
// cannot see this one; keep the two function texts identical -- the files
// differ only in their line endings): data
// signature for -predict- (e(p_dsig_type) "rowsig2") that ties every value to
// its (panel, time) key. The first two variables must be the panel and time
// variables. Rows are sorted by those keys and each column is summarized by
// its plain sum and two sums weighted by Park-Miller sequences in the row
// rank, so values moved between rows, relabelled keys, added or dropped
// observations, and changes to single values change the signature (unless
// smaller than the rounding of the column sums, about 1e-16 of a sum),
// while a change of storage type (compress) or of the sort order does not.
// Missing values enter as a fixed sentinel. rowsig1 (0.9.29-0.9.32) used the
// weights 1 + mod(r*a, m)/m, which are linear in r below r = m/a (about
// 30,000 rows), so its two weighted sums added only sum(r*x) to sum(x): a
// change such as +c, -2c, +c on three consecutive rows kept all three sums.
string scalar xdpt2_rowsig2(string scalar vars)
{
    real matrix X
    real colvector ord, w1, w2
    real scalar j, n, x1, x2
    string scalar s
    X = st_data(., vars)
    n = rows(X)
    if (n == 0) return("0")
    // v0.9.34: sort by every column, not only (panel, time): rows that
    // share a key (several rows of a panel with a missing time) are ties
    // that order() does not keep in a fixed order; identical rows are
    // interchangeable. For unique keys the order is unchanged.
    ord = order(X, (1..cols(X)))
    X = editmissing(X[ord, .], -9876543210.125)
    // Park-Miller sequences x(j) = a*x(j-1) mod (2^31 - 1), x(0) = 1, for
    // a = 48271 and a = 69621; a*x < 2^48, so every step is exact in double
    w1 = J(n, 1, .)
    w2 = J(n, 1, .)
    x1 = 1
    x2 = 1
    for (j = 1; j <= n; j++) {
        x1 = mod(48271 * x1, 2147483647)
        x2 = mod(69621 * x2, 2147483647)
        w1[j] = x1
        w2[j] = x2
    }
    w1 = 1 :+ w1 :/ 2147483647
    w2 = 1 :+ w2 :/ 2147483647
    s = "2:" + strofreal(n) + ":" + strofreal(cols(X))
    for (j = 1; j <= cols(X); j++) {
        s = s + ":" + strofreal(sum(X[., j]), "%21x") +
                "," + strofreal(sum(X[., j] :* w1), "%21x") +
                "," + strofreal(sum(X[., j] :* w2), "%21x")
    }
    return(s)
}

// v0.7.3 (D5): fill the -predict- target variable from the persisted
// AR-test rows. which: 1 = residuals (ê), 2 = xb (= Δy − ê). Writes only on
// rows whose (panelvar, timevar) key matches a stored estimation row AND that
// are marked by touse; everything else is left untouched (missing).
// serial/token/checksum identity guards against stale Mata state (mata clear,
// restart, or e() restored from a different run).
void xdpt2_p_fill(string scalar pvar, string scalar tvar,
                   string scalar outvar, string scalar touse,
                   real scalar source, real scalar which,
                   real scalar serial_expect, real scalar sig_expect,
                   string scalar token_expect)
{
    external transmorphic xdpt_p_store_r, xdpt_p_store_e, xdpt_p_store_sig
    external transmorphic xdpt_p_store_token
    real matrix D, S, S_r, S_e
    real scalar r, v, sig_have, sig_actual
    string scalar token_have
    transmorphic A

    // v0.9.3 R19 (#4): rows looked up BY SERIAL, so a model brought back
    // with -estimates restore- predicts from its own stored rows.
    if (serial_expect >= . | sig_expect >= . | token_expect == "" |
        eltype(xdpt_p_store_r) == "real" | eltype(xdpt_p_store_e) == "real" |
        eltype(xdpt_p_store_sig) == "real" |
        eltype(xdpt_p_store_token) == "real") {
        errprintf("xtdpthresh predict: stored estimation rows not found in Mata memory\n")
        errprintf("  (cleared by -mata: mata clear-, -discard-, or restarting Stata).\n")
        errprintf("  Re-run xtdpthresh, then predict.\n")
        exit(498)
    }
    if (!asarray_contains(xdpt_p_store_token, serial_expect)) {
        errprintf("xtdpthresh predict: this run's exact cache token is no longer in Mata memory\n")
        errprintf("  Re-run xtdpthresh, then predict.\n")
        exit(498)
    }
    token_have = asarray(xdpt_p_store_token, serial_expect)
    if (token_have != token_expect) {
        errprintf("xtdpthresh predict: cached rows belong to a different model/run\n")
        errprintf("  (the Mata serial was reused after state was cleared).\n")
        errprintf("  Re-run xtdpthresh, then predict.\n")
        exit(498)
    }
    if (!asarray_contains(xdpt_p_store_sig, serial_expect)) {
        errprintf("xtdpthresh predict: this run's cache guard is no longer in Mata memory\n")
        errprintf("  (evicted after 20 newer runs, or Mata was cleared).\n")
        errprintf("  Re-run xtdpthresh, then predict.\n")
        exit(498)
    }
    sig_have = asarray(xdpt_p_store_sig, serial_expect)
    if (sig_have != sig_expect) {
        errprintf("xtdpthresh predict: cached rows belong to a different model/run\n")
        errprintf("  (the Mata cache key was reused after state was cleared).\n")
        errprintf("  Re-run xtdpthresh, then predict.\n")
        exit(498)
    }
    // v0.9.19: validate the checksum against the ACTUAL row matrices, not
    // merely against its parallel stored copy. Otherwise an overwritten or
    // corrupted Mata row cache passes the old sig_have==sig_expect test.
    if (!asarray_contains(xdpt_p_store_r, serial_expect) |
        !asarray_contains(xdpt_p_store_e, serial_expect)) {
        errprintf("xtdpthresh predict: this run's stored rows are no longer in Mata\n")
        errprintf("  memory (evicted after 20 newer runs, or Mata was cleared).\n")
        errprintf("  Re-run xtdpthresh, then predict.\n")
        exit(498)
    }
    S_r = asarray(xdpt_p_store_r, serial_expect)
    S_e = asarray(xdpt_p_store_e, serial_expect)
    sig_actual = hash1(S_r, 2147483647) * 4194304 +
                 mod(hash1(S_e), 4194304)
    if (sig_actual != sig_have) {
        errprintf("xtdpthresh predict: cached estimation rows failed their checksum\n")
        errprintf("  (Mata cache contents were changed or corrupted).\n")
        errprintf("  Re-run xtdpthresh, then predict.\n")
        exit(498)
    }
    // source: 1 = AR-test (FD) series; 2 = estimation-equation series
    S = (source == 2 ? S_e : S_r)

    A = asarray_create("real", 2)
    asarray_notfound(A, .)
    for (r = 1; r <= rows(S); r++) {
        asarray(A, (S[r, 1], S[r, 2]),
                (which == 2 ? S[r, 3] - S[r, 4] : S[r, 4]))
    }
    st_view(D = ., ., (pvar, tvar, outvar), touse)
    for (r = 1; r <= rows(D); r++) {
        v = asarray(A, (D[r, 1], D[r, 2]))
        if (v < .) D[r, 3] = v
    }
}

// Main orchestrator
void xtdpthresh_run(string scalar depvar_name,
                      string scalar Ly_name,
                      string scalar exog_names,
                      string scalar endog_names,
                      string scalar predet_names,
                      string scalar inst_names,
                      string scalar q_name,
                      string scalar panelvar_name,
                      string scalar timevar_name,
                      string scalar method,
                      real scalar flag_static,
                      real scalar flag_kink,
                      real scalar flag_collapse,
                      real scalar maxlag_lo, real scalar maxlag_hi,
                      real scalar n_grid,
                      real scalar n_gridci,
                      real scalar trim_rate,
                      real scalar q_lo,
                      real scalar q_hi,
                      real scalar do_grid_ci,
                      real scalar n_boot,
                       real scalar alpha,
                       real scalar flag_iv_collapse,
                       real scalar flag_exportgmm,
                       real scalar flag_notest,
                       real scalar flag_cont_test)
{
    // Keep the numerical backend closed to unsupported transformations.
    if (method != "fd" & method != "fod") {
        errprintf("xtdpthresh: backend supports method(fd) and method(fod) only\n")
        exit(198)
    }
    external real scalar xdpt_collapse, xdpt_iv_collapse, xdpt_lag_lo, xdpt_lag_hi
    external real scalar xdpt_verbose, xdpt_trim_rate
    xdpt_collapse = flag_collapse
    xdpt_iv_collapse = flag_iv_collapse
    xdpt_lag_lo = maxlag_lo
    xdpt_lag_hi = maxlag_hi
    xdpt_trim_rate = trim_rate
    xdpt_verbose = strtoreal(st_local("flag_verbose"))
    // v0.7.13 (C1): Windmeijer flag read from the caller's frame; the
    // applied-flag is reset here and set inside the grid search only when
    // the correction actually replaced the reported V (two-step path with a
    // successful correction solve).
    external real scalar xdpt_vce_wind, xdpt_wind_applied
    xdpt_vce_wind = strtoreal(st_local("flag_vce_wind"))
    if (xdpt_vce_wind >= .) xdpt_vce_wind = 0
    xdpt_wind_applied = 0
    // v0.9.32: kink joint-variance flag and the conditional V kept for the
    // AR statistics; set inside the grid search for the kink model only.
    external real matrix xdpt_V_cond_ar
    external real scalar xdpt_kink_joint
    xdpt_V_cond_ar = J(0, 0, .)
    xdpt_kink_joint = .
    // v0.9.35: the full joint variance of (theta, gamma), for the AR tests
    external real matrix xdpt_V_joint_full
    xdpt_V_joint_full = J(0, 0, .)
    // v0.9.34 (C3): 1 if the kink refinement moved gamma-hat off the grid
    external real scalar xdpt_kref
    xdpt_kref = .
    // v0.9.35: residual diagnostics of the dropped instrument columns, set
    // by the final stack at gamma-hat (xdpt2_indep_cols)
    external real scalar xdpt_ivc_diag, xdpt_ivc_dep_res, xdpt_ivc_dep_near
    xdpt_ivc_diag = 0
    xdpt_ivc_dep_res = .
    xdpt_ivc_dep_near = 0
    // v0.9.35: kernel bandwidth of the jump model's joint variance
    external real scalar xdpt_bwscale, xdpt_gamma_bw
    xdpt_bwscale = strtoreal(st_local("bwscale"))
    if (xdpt_bwscale >= . | xdpt_bwscale <= 0) xdpt_bwscale = 1.5
    xdpt_gamma_bw = .
    // v0.9.9 R26: exportgmm flag visible to grid_search, plus a RESET of
    // the Windmeijer-input externals (populated inside grid_search when
    // the correction runs under exportgmm; used by the independent
    // numerical certification _cert_windmeijer.do -- finite-difference
    // dOmega/dtheta and component-level V2/DV2/V2D'/DV1rD' checks).
    external real scalar xdpt_expg
    xdpt_expg = strtoreal(st_local("flag_exportgmm"))
    if (xdpt_expg >= .) xdpt_expg = 0
    external real matrix xdpt_w_ZW1, xdpt_w_X1, xdpt_w_Z, xdpt_w_Om1
    external real matrix xdpt_w_W1, xdpt_w_W2, xdpt_w_ZW2
    external real colvector xdpt_w_uid, xdpt_w_r1, xdpt_w_gbar2
    external real scalar xdpt_w_n
    xdpt_w_ZW1 = J(0, 0, .)
    xdpt_w_X1 = J(0, 0, .)
    xdpt_w_Z = J(0, 0, .)
    xdpt_w_Om1 = J(0, 0, .)
    xdpt_w_W1 = J(0, 0, .)
    xdpt_w_W2 = J(0, 0, .)
    xdpt_w_ZW2 = J(0, 0, .)
    xdpt_w_uid = J(0, 1, .)
    xdpt_w_r1 = J(0, 1, .)
    xdpt_w_gbar2 = J(0, 1, .)
    xdpt_w_n = .
    // v0.9.29: first-step weight fallback flag (see xdpt2_build_W_ma1).
    external real scalar xdpt_w1_fallback
    xdpt_w1_fallback = 0
    // v0.9.30: instrument columns dropped by xdpt2_drop_cellconst (last stack).
    external real scalar xdpt_ivc_drop
    xdpt_ivc_drop = 0
    // v0.9.31: number of stage-1 grid points (the initial grid; refine()
    // appends its candidates after them), for the bootstrap replays.
    external real scalar xdpt_n_stage1
    xdpt_n_stage1 = .
    // v0.9.31: drop-mask cache of xdpt2_stack_at_gamma, reset per run.
    external real matrix xdpt_ivc_Zref
    external real rowvector xdpt_ivc_mask
    external real colvector xdpt_ivc_tref
    xdpt_ivc_Zref = J(0, 0, .)
    xdpt_ivc_mask = J(1, 0, .)
    xdpt_ivc_tref = J(0, 1, .)
    xdpt2_tpl_release()
    // v0.7.13 (C2): FWL time-dummy partialling flag for the stack builder.
    external real scalar xdpt_td_fwl
    xdpt_td_fwl = strtoreal(st_local("flag_td_fwl"))
    if (xdpt_td_fwl >= .) xdpt_td_fwl = 0
    // v0.7.13 (C3): quantile-spaced γ grids.
    external real scalar xdpt_grid_quant
    xdpt_grid_quant = strtoreal(st_local("flag_grid_quant"))
    if (xdpt_grid_quant >= .) xdpt_grid_quant = 0
    // v0.8.0 (R5): centered moment covariance flag.
    external real scalar xdpt_center
    xdpt_center = strtoreal(st_local("flag_center"))
    if (xdpt_center >= .) xdpt_center = 0
    // v0.8.1 (R6): coefficient-bootstrap controls.
    external real scalar xdpt_coefci_sym, xdpt_coefboot_2s
    xdpt_coefci_sym = strtoreal(st_local("flag_coefci_sym"))
    if (xdpt_coefci_sym >= .) xdpt_coefci_sym = 1
    // v0.9.0: exact-resampling grid bootstrap flag
    external real scalar xdpt_boot_exact
    xdpt_boot_exact = strtoreal(st_local("flag_boot_exact"))
    if (xdpt_boot_exact >= .) xdpt_boot_exact = 0
    // v0.9.3 R19 (#5): coefboot(none)
    external real scalar xdpt_coefboot_off
    xdpt_coefboot_off = strtoreal(st_local("flag_coefboot_off"))
    if (xdpt_coefboot_off >= .) xdpt_coefboot_off = 0
    xdpt_coefboot_2s = strtoreal(st_local("flag_coefboot_2s"))
    if (xdpt_coefboot_2s >= .) xdpt_coefboot_2s = 1
    real colvector y, q, pid, tid, gamma_grid
    real matrix Ly
    real matrix X_exog, X_endog, X_predet, X_inst
    struct xdpt2_unit rowvector units
    real scalar t_min, t_max, i, n_units

    y   = st_data(., depvar_name)
    q   = st_data(., q_name)
    pid = st_data(., panelvar_name)
    tid = st_data(., timevar_name)
    // v0.8.1 (R6 #1): equation-eligibility flag (1 = complete row that may
    // form a GMM equation; 0 = history-only instrument-source row).
    real colvector eqv
    eqv = st_data(., st_local("eqvar"))
    // v0.9.1 R17 (#3): global sorted vector of OBSERVED equation times.
    // Instrument blocks and td-fod dummy columns are indexed by RANK in
    // this vector, never by raw calendar offset -- a sparse delta-1 index
    // (daily or finer dates) can no longer force span-sized allocations
    // anywhere. On a gap-free index rank == offset, results bit-for-bit.
    external real colvector xdpt_teq
    xdpt_teq = select(tid, eqv :== 1)
    if (rows(xdpt_teq) == 0) xdpt_teq = tid
    xdpt_teq = uniqrows(xdpt_teq)
    Ly  = J(rows(y), 0, 0)
    if (!flag_static & Ly_name != "") Ly = st_data(., Ly_name)

    X_exog   = J(rows(y), 0, 0)
    X_endog  = J(rows(y), 0, 0)
    X_predet = J(rows(y), 0, 0)
    X_inst   = J(rows(y), 0, 0)
    if (exog_names   != "") X_exog   = st_data(., exog_names)
    if (endog_names  != "") X_endog  = st_data(., endog_names)
    if (predet_names != "") X_predet = st_data(., predet_names)
    if (inst_names   != "") X_inst   = st_data(., inst_names)

    t_min = min(tid)
    t_max = max(tid)

    // v0.9.3 R19 (#11): estimate the uncollapsed instrument width BEFORE
    // any allocation. The default maxlag() is unlimited; on long panels
    // the block-diagonal Z proliferates (weak Hansen, heavy bootstrap,
    // possible memory failure) -- say so while it is still cheap to stop.
    if (st_local("nowarn") == "") {
        real scalar est_hi, est_lo_y, est_lo_p, est_n_y, est_n_p
        real scalar est_lag, est_wid
        external real scalar xdpt_collapse
        est_hi = maxlag_hi
        if (est_hi > t_max - t_min) est_hi = t_max - t_min
        est_lo_y = (maxlag_lo > 2 ? maxlag_lo : 2)
        est_lo_p = (maxlag_lo > 1 ? maxlag_lo : 1)
        est_n_y = (est_hi >= est_lo_y ? est_hi - est_lo_y + 1 : 0)
        est_n_p = (est_hi >= est_lo_p ? est_hi - est_lo_p + 1 : 0)
        est_lag = max((est_n_y, est_n_p))
        est_wid = rows(xdpt_teq) * (1 + (flag_static ? 0 : est_n_y) +
                   cols(X_endog)*est_n_y + cols(X_predet)*est_n_p + cols(X_exog))
        if (cols(X_inst) > 0) {
            est_wid = est_wid + (xdpt_iv_collapse ? cols(X_inst) :
                      rows(xdpt_teq)*cols(X_inst))
        }
        if (xdpt_collapse != 1 & est_wid > 1500) {
            printf("{err}warning: ~%g instrument columns projected -- this will be slow and\n", est_wid)
            printf("{err}         memory-heavy. The unlimited maxlag() default is a RISKY\n")
            printf("{err}         choice on long panels, not a safe one: use maxlag(1 3) or\n")
            printf("{err}         -collapse- unless the full lag ladder is intentional.\n")
        }
        else if (xdpt_collapse != 1 & est_wid > 300) {
            printf("{txt}note: the uncollapsed instrument matrix may reach ~%g columns\n", est_wid)
            printf("{txt}      (maxlag depth %g over %g observed periods); consider maxlag(1 3)\n", est_lag, rows(xdpt_teq))
            printf("{txt}      or -collapse- to curb instrument proliferation.\n")
        }
    }
    units = xdpt2_build_units(y, Ly, X_exog, X_endog, X_predet, X_inst,
                               q, pid, tid, eqv, flag_static)
    n_units = length(units)
    if (n_units < 5) {
        errprintf("xtdpthresh: need >= 5 units (got %g)\n", n_units)
        exit(498)
    }

    // Rebuild the global time support from RETAINED units. Under
    // history(panel), panels wholly outside if/in were previously allowed to
    // stretch t_min/t_max and the IV lag width despite contributing no row.
    real scalar _nteq, _teqpos, _ui
    real colvector _tei
    _nteq = 0
    t_min = .
    t_max = .
    for (_ui = 1; _ui <= n_units; _ui++) {
        if (t_min >= . | units[_ui].t[1] < t_min) t_min = units[_ui].t[1]
        if (t_max >= . | units[_ui].t[rows(units[_ui].t)] > t_max) ///
            t_max = units[_ui].t[rows(units[_ui].t)]
        _nteq = _nteq + sum(units[_ui].eq)
    }
    if (_nteq == 0) {
        errprintf("xtdpthresh: retained units contain no equation-eligible rows\n")
        exit(498)
    }
    xdpt_teq = J(_nteq, 1, .)
    _teqpos = 1
    for (_ui = 1; _ui <= n_units; _ui++) {
        _tei = select(units[_ui].t, units[_ui].eq :== 1)
        if (rows(_tei) == 0) continue
        xdpt_teq[|_teqpos \ _teqpos + rows(_tei) - 1|] = _tei
        _teqpos = _teqpos + rows(_tei)
    }
    xdpt_teq = uniqrows(xdpt_teq)
    // v0.9.31: t_max only sets the lag width (xdpt2_transform_unit) and the
    // allocation gate below. Instruments are lags, so the widest lag any
    // equation can use is (last equation time) - t_min: under FD t - L >=
    // t_min, under FOD t + 1 - L >= t_min with t at most the last equation
    // time minus one. Trailing history rows (after the last equation, e.g.
    // from -if t <= 30- on a longer panel) only added all-zero lag columns,
    // which were dropped, but they inflated the gate and could reject a
    // valid large design.
    t_max = max(xdpt_teq)

    // v0.9.19: fail BEFORE allocating a nominally enormous instrument
    // matrix. Lag columns are pruned only after stacking; on a long or
    // sparse delta-1 calendar, an open/large maxlag() can
    // otherwise allocate a huge Z and then form an L x L weight matrix,
    // exhausting RAM before the existing proliferation warning can help.
    // These are conservative upper bounds before exact-zero columns drop.
    real scalar _iv_span, _iv_hi, _iv_lo_y, _iv_lo_p, _iv_ny, _iv_np
    real scalar _iv_tstart, _iv_toff, _iv_nt, _iv_per_t, _iv_tcols
    real scalar _iv_total, _iv_row_upper, _iv_maxcols, _iv_maxcells
    _iv_span = t_max - t_min
    _iv_hi = maxlag_hi
    if (_iv_hi > _iv_span) _iv_hi = _iv_span
    _iv_lo_y = (maxlag_lo > 2 ? maxlag_lo : 2)
    _iv_lo_p = (maxlag_lo > 1 ? maxlag_lo : 1)
    _iv_ny = (_iv_hi >= _iv_lo_y ? _iv_hi - _iv_lo_y + 1 : 0)
    _iv_np = (_iv_hi >= _iv_lo_p ? _iv_hi - _iv_lo_p + 1 : 0)
    _iv_tstart = (method == "fd" ? t_min + 1 : t_min)
    _iv_toff = xdpt2_tpos(xdpt_teq, _iv_tstart)
    _iv_nt = rows(xdpt_teq) - _iv_toff + 1
    if (_iv_nt < 1) _iv_nt = 1
    _iv_per_t = 1 + (flag_static ? 0 : _iv_ny) + cols(X_exog) +
                cols(X_endog) * _iv_ny + cols(X_predet) * _iv_np
    _iv_tcols = (flag_collapse ? _iv_per_t : _iv_nt * _iv_per_t) +
                (cols(X_inst) == 0 ? 0 :
                 (flag_iv_collapse ? cols(X_inst) : _iv_nt * cols(X_inst)))
    _iv_total = _iv_tcols
    _iv_row_upper = _nteq
    _iv_maxcols = 5000
    _iv_maxcells = 50000000
    if (_iv_total > _iv_maxcols |
        _iv_row_upper * _iv_total > _iv_maxcells) {
        errprintf("xtdpthresh: projected instrument allocation is unsafe before zero-column pruning\n")
        errprintf("  (up to %g columns and %g Z cells over retained span %g).\n",
                  _iv_total, _iv_row_upper * _iv_total, _iv_span)
        errprintf("  Use an explicit, tighter maxlag() and/or specify collapse.\n")
        errprintf("  This gate bounds the size of the instrument matrix; the per-gamma\n")
        errprintf("  caches need additional memory that it does not count.\n")
        exit(498)
    }

    // Recompute trim bounds on the rows that can actually enter the GMM
    // equations. The pre-Mata percentiles are only a safe initialization;
    // transformation/history/zero-IV filters can materially change q support.
    real matrix Y_eff, W_eff, Z_eff
    real colvector times_eff, uid_eff, q_eff
    real scalar n_units_eff
    xdpt2_stack_at_gamma(units, (q_lo + q_hi) / 2, method,
                          flag_static, flag_kink, t_min, t_max,
                          Y_eff, W_eff, Z_eff, times_eff, uid_eff)
    if (rows(Y_eff) < 20) {
        errprintf("xtdpthresh: fewer than 20 usable GMM rows after transformation\n")
        exit(498)
    }
    // v0.9.28: Z does not depend on gamma, so its conditioning is checked
    // once here. The grid loop applies the same Z'Z test at every candidate
    // and, when it fails, reported only that no gamma was admitted. The
    // decision is identical; only the message and the time to fail change.
    real matrix ZZ_eff, ZZ_eff_inv
    real scalar zz_eff_ok
    ZZ_eff = Z_eff' * Z_eff / rows(Z_eff)
    xdpt2_syminv(ZZ_eff, zz_eff_ok, ZZ_eff_inv)
    if (!zz_eff_ok) {
        // v0.9.34: dependent columns are dropped in the stack
        // (xdpt2_indep_cols), so what reaches this gate is near dependence.
        // v0.9.35: columns within a relative residual of 3.2e-7 are dropped
        // too, so this is near dependence beyond that.
        errprintf("xtdpthresh: the instrument matrix is ill-conditioned: some instrument\n")
        errprintf("  columns are nearly, though not exactly, linear combinations of others\n")
        errprintf("  (columns closer than a relative residual of 3.2e-7 are dropped). Common\n")
        errprintf("  causes: an iv() variable close to a combination of instruments already\n")
        errprintf("  present (for example one built from lags of y), or lags of a variable\n")
        errprintf("  common to all units only up to rounding. A smaller instrument set\n")
        errprintf("  (collapse, maxlag()) can help, as can rescaling nearly equal\n")
        errprintf("  instruments (for example, (z2 - z1)/c in place of z2).\n")
        // v0.9.30: variables common to all units only up to rounding are not
        // recognized as common; name them, since they are a likely cause.
        string scalar _ncv
        _ncv = xdpt2_near_common(st_local("_cc_vars"), st_local("_cc_labs"),
                                 st_local("timevar"), st_local("touse"))
        if (_ncv != "") {
            errprintf("  Nearly common to all units in every period, but not exactly: %s.\n", _ncv)
            errprintf("  Make such a variable exactly equal across units (for example, merge it\n")
            errprintf("  by period) so that xtdpthresh can treat it as a common variable.\n")
        }
        exit(498)
    }
    // v0.9.31: a regressor removed by the unit and time effects together
    // (a_i + g_t) under td -- see xdpt2_stack_at_gamma. FOD used to fit the
    // rounding residue (coefficients near 1e25, rc 0); FD failed generically.
    external real rowvector xdpt_td_gone
    if (cols(xdpt_td_gone) > 0) {
        string rowvector _bl
        string scalar _gone
        real scalar _gi
        _bl = J(1, 0, "")
        if (!flag_static) _bl = _bl, ("L." + st_local("depvar"))
        _bl = _bl, tokens(st_local("all_exog_lab")), tokens(st_local("endog_lab")),
              tokens(st_local("predet_lab"))
        _gone = ""
        for (_gi = 1; _gi <= cols(xdpt_td_gone); _gi++) {
            _gone = _gone + (_gone == "" ? "" : " ") +
                (xdpt_td_gone[_gi] <= cols(_bl) ? _bl[xdpt_td_gone[_gi]] :
                 "regressor " + strofreal(xdpt_td_gone[_gi]))
        }
        errprintf("xtdpthresh: %s is removed by the unit and time effects together: it\n", _gone)
        errprintf("  has the form a_i + g_t (for example, firm age = year - founding year),\n")
        errprintf("  which the %s transformation and td remove exactly, so its coefficient\n", strupper(method))
        errprintf("  is not identified. Remove it, or drop td.\n")
        exit(498)
    }
    // v0.9.31: after the columns that repeat the constant instruments are
    // dropped, the threshold needs at least one instrument more than the
    // coefficients; otherwise every candidate would be rejected with a
    // generic message.
    external real scalar xdpt_ivc_drop
    if (xdpt_ivc_drop > 0 & cols(Z_eff) < cols(W_eff) + 1) {
        errprintf("xtdpthresh: %g instrument column(s) that repeat the constant instruments\n", xdpt_ivc_drop)
        errprintf("  (variables common to all units in each period) were dropped, leaving %g\n", cols(Z_eff))
        errprintf("  instruments for %g coefficients; at least %g are needed. Add\n", cols(W_eff), cols(W_eff) + 1)
        errprintf("  instruments (a wider maxlag(), or iv()) or remove such variables.\n")
        exit(498)
    }

    // v0.8.1 R7 (audit): the grid/trim support is the set of q values whose
    // indicators enter the RETAINED transformed rows (q_t and q_{t-1} under
    // FD; q_t and future equation-row q under FOD) -- not only current-row
    // q (v0.8.0 and earlier: missed breakpoints) and not the full history
    // (v0.8.1 first cut: history-only rows and dropped units could stretch
    // the trim range with values that never touch the criterion).
    if (st_local("gridsample") == "observed") {
        // v0.8.2 R10 (#5): xthenreg-style support -- current-row q of the
        // retained equation rows only.
        // Deduplicate by (unit,time) before mapping back to q.
        real matrix obs_keys
        obs_keys = uniqrows((uid_eff, times_eff))
        q_eff = xdpt2_q_at_rows(units, obs_keys[., 2], obs_keys[., 1])
        q_eff = select(q_eff, q_eff :< .)
    }
    else {
        q_eff = xdpt2_q_support(units, times_eff, uid_eff, method)
    }
    // v0.8.2 R9 (#7): the deduplicated effective support is passed to the
    // cache builders as an ARGUMENT (grid construction and grid admission
    // share one definition) -- no global Mata state to go stale.
    real scalar minreg_user
    minreg_user = strtoreal(st_local("minregime"))
    if (minreg_user >= .) minreg_user = 0
    // v0.8.2 R10 (#3): fail fast with a SPECIFIC message when minregime()
    // makes every threshold split inadmissible, instead of the generic
    // no-gamma-admitted error downstream.
    if (minreg_user > 0 & 2*minreg_user > rows(q_eff)) {
        errprintf("xtdpthresh: minregime(%g) leaves no admissible threshold split\n", minreg_user)
        errprintf("  (effective support has only %g observations; each regime needs\n", rows(q_eff))
        errprintf("  at least minregime() of them at every candidate gamma)\n")
        exit(498)
    }
    // v0.8.2 R11 (#7): reproducibility metadata -- the default trimming
    // floor and the floor actually applied (max of the two); mirrors
    // xdpt2_build_gamma_cache exactly.
    real scalar minreg_def, minreg_applied
    if (rows(q_eff) < 20) {
        errprintf("xtdpthresh: threshold variable has too few usable values\n")
        exit(498)
    }
    q_lo = xdpt2_quantile(q_eff, trim_rate / 2)
    q_hi = xdpt2_quantile(q_eff, 1 - trim_rate / 2)
    if (q_lo >= q_hi) {
        errprintf("xtdpthresh: qx() has insufficient variation on the effective GMM sample\n")
        exit(498)
    }
    // v0.9.34: the default floor is the number of observations in the
    // smaller trimmed tail, min(#{q <= q_lo}, #{q > q_hi}): the counts at the
    // two ends of the grid, so both ends are admitted. It was ceil(trim*n/2);
    // when q_hi was a data value only ceil(trim*n/2) - 1 observations lay
    // above it, and the top grid point was admitted or not by rounding.
    external real scalar xdpt_minreg_def
    minreg_def = min((sum(q_eff :<= q_lo), sum(q_eff :> q_hi)))
    if (minreg_def < 2) minreg_def = 2
    minreg_applied = (minreg_user > minreg_def ? minreg_user : minreg_def)
    xdpt_minreg_def = minreg_def
    n_units_eff = rows(uniqrows(uid_eff))
    if (n_units_eff < 5) {
        errprintf("xtdpthresh: need >= 5 contributing units (got %g)\n", n_units_eff)
        exit(498)
    }
    printf("  %g units built (n_obs=%g, T range [%g, %g])\n",
           n_units_eff, rows(y), t_min, t_max)

    // Grid search over γ
    // v0.7.13 (C3): estimation grid — uniform values (default) or empirical
    // quantiles of q over the effective sample (gridtype(quantile)).
    if (xdpt_grid_quant == 1) {
        gamma_grid = xdpt2_quantile_grid(q_eff, trim_rate/2, 1 - trim_rate/2,
                                          n_grid)
    }
    else gamma_grid = xdpt2_rangen(q_lo, q_hi, n_grid)

    real scalar best_obj, best_gamma
    real colvector best_theta
    real matrix best_V, best_V_influence
    // v0.8.7 hotfix: out-args passed into TYPED (real scalar) parameters
    // must be initialized scalars -- an auto-created 0x0 local fails the
    // callee's argument check (Mata 3204) even though the callee assigns
    // them first thing.
    real scalar n_adm2, grid2_adm_lo, grid2_adm_hi
    n_adm2 = .
    grid2_adm_lo = .
    grid2_adm_hi = .

    printf("  Grid search over %g γ points in [%8.4f, %8.4f]...\n",
           n_grid, q_lo, q_hi)
    // v0.7.0 (D1): build the estimation-grid cache ONCE; it is shared by the
    // grid search, the grid-bootstrap CI, the linearity test, and (as the
    // jump cache) the continuity test below.
    struct xdpt2_gamma_cache rowvector cache_main
    // v0.9.34 (SPEEDUP): the template of this build serves the later builds
    // of the same model (refinements, the confidence-set grid)
    struct xdpt2_stack_tpl scalar tpl_main
    real scalar tpl_main_st
    tpl_main_st = 0
    cache_main = xdpt2_build_gamma_cache_t(units, gamma_grid, method,
                                            flag_static, flag_kink, t_min,
                                            t_max, q_eff, minreg_user,
                                            tpl_main, tpl_main_st)
    // v0.8.2 R10 (#2/#4): the search space users should reason about is the
    // ADMITTED grid, not the nominal trim bounds -- minregime,
    // ties, rank/condition failures and solver availability all prune
    // points. Track requested vs distinct vs admitted and the admitted span.
    real colvector gamma_admitted
    real scalar grid_adm_lo, grid_adm_hi, gg_a, n_struct
    gamma_admitted = J(0, 1, .)
    n_struct = 0
    // Structural admission is cache ok. Search admission is exported by the
    // full fixed-W1 solve below: fast_ok alone checks C_g but not finite theta,
    // residuals, moments, or objective.
    for (gg_a = 1; gg_a <= rows(gamma_grid); gg_a++) {
        if (!cache_main[gg_a].ok) continue
        n_struct = n_struct + 1
    }
    real matrix best_A
    real scalar best_twostep
    // v0.9.25 adaptive-search bookkeeping.  q_split_supp always uses the
    // transformed design support, even when gridsample(observed) defines the
    // initial grid, because stopping is about the realized regime split.
    real colvector q_split_supp
    q_split_supp = xdpt2_q_support(units, times_eff, uid_eff, method)
    if (rows(q_split_supp) > 0) {
        q_split_supp = uniqrows(sort(q_split_supp, 1))
    }
    real scalar search_l1_n, search_l2_n, search_l3_n
    real scalar search_s1_level, search_s2_level
    real scalar search_s1_n, search_s2_n
    real scalar search_s1_same, search_s2_same
    real scalar search_s1_gain, search_s2_gain
    real scalar search_s1_conv, search_s2_conv
    real scalar search_hit_max, search_W2_builds
    real scalar search_g1, search_o1, search_g2_global, search_o2_global
    real scalar n_refine, ref_it, ref_added, ref_obj_gain
    real scalar ref_pool_n, ref_remaining, ref_exhausted
    real scalar ref_alo, ref_ahi, ref_in_basin, ref_neigh_rem, ref_complete
    string scalar search_mode
    real scalar search_tol, search_max_level
    search_mode = st_local("searchmode")
    search_tol = strtoreal(st_local("searchtol"))
    search_max_level = strtoreal(st_local("search_max_level"))
    n_refine = strtoreal(st_local("refine"))
    if (n_refine >= .) n_refine = 0
    search_l1_n = .
    search_l2_n = .
    search_l3_n = .
    search_s1_level = .
    search_s2_level = .
    search_s1_n = .
    search_s2_n = .
    search_s1_same = .
    search_s2_same = .
    search_s1_gain = .
    search_s2_gain = .
    search_s1_conv = .
    search_s2_conv = .
    search_hit_max = .
    search_W2_builds = .
    search_g1 = .
    search_o1 = .
    search_g2_global = .
    search_o2_global = .
    ref_it = .
    ref_added = .
    ref_obj_gain = .
    ref_pool_n = .
    ref_remaining = .
    ref_exhausted = .
    ref_alo = .
    ref_ahi = .
    ref_in_basin = .
    ref_neigh_rem = .
    ref_complete = .

    xdpt2_grid_search(units, gamma_grid, cache_main, method, flag_static,
                       flag_kink, t_min, t_max, q_eff, q_split_supp,
                       minreg_user, search_mode, search_tol,
                       search_max_level, n_refine,
                       best_gamma, best_obj, best_theta, best_V,
                       best_V_influence, best_A,
                       best_twostep, n_adm2, grid2_adm_lo, grid2_adm_hi,
                       gamma_admitted, search_l1_n, search_l2_n,
                       search_l3_n, search_s1_level, search_s2_level,
                       search_s1_n, search_s2_n, search_s1_same,
                       search_s2_same, search_s1_gain, search_s2_gain,
                       search_s1_conv, search_s2_conv, search_hit_max,
                       search_W2_builds, search_g1, search_o1,
                       search_g2_global, search_o2_global,
                       ref_it, ref_added, ref_pool_n, ref_remaining,
                       ref_exhausted, ref_alo, ref_ahi, ref_in_basin,
                       ref_neigh_rem, ref_complete, ref_obj_gain,
                       tpl_main, tpl_main_st)
    // v0.9.31: stage 1 of the search runs on the initial grid only
    // (entries 1..search_s1_n of cache_main); the replays must match.
    xdpt_n_stage1 = search_s1_n
    // The engine may have appended global level-2/3 and local support entries.
    n_struct = 0
    for (gg_a = 1; gg_a <= rows(gamma_grid); gg_a++) {
        if (cache_main[gg_a].ok) n_struct = n_struct + 1
    }
    if (rows(gamma_admitted) > 0) {
        grid_adm_lo = min(gamma_admitted)
        grid_adm_hi = max(gamma_admitted)
    }
    else {
        grid_adm_lo = .
        grid_adm_hi = .
    }
    printf("  Grid: initial %g, cached %g, structural %g, stage-1 admitted %g; span [%8.4f, %8.4f]\n",
           n_grid, rows(gamma_grid), n_struct, rows(gamma_admitted), grid_adm_lo, grid_adm_hi)
    if (search_mode == "adaptive") {
        printf("  Nested search: stage 1 level %g (%g points, converged=%g); stage 2 level %g (%g points, converged=%g); W2 builds=%g\n",
               search_s1_level, search_s1_n, search_s1_conv,
               search_s2_level, search_s2_n, search_s2_conv,
               search_W2_builds)
    }
    displayflush()

    if (best_gamma == . & rows(gamma_admitted) >= 2) {
        // v0.9.34: stage 1 admitted candidates, but the profile is flat
        errprintf("xtdpthresh: the GMM criterion is flat over the admitted grid: every\n")
        errprintf("  candidate threshold gives the same fit (for example, q takes few values,\n")
        errprintf("  so all admitted candidates give the same regime split). The threshold is\n")
        errprintf("  not identified on this grid.\n")
        exit(498)
    }
    if (best_gamma == .) {
        errprintf("xtdpthresh: point estimation failed (no γ admitted a valid GMM solve)\n")
        errprintf("  Hard requirements include >=20 usable rows, at least K+1 moments,\n")
        errprintf("  nonsingular moment matrices, and a non-flat profiled criterion.\n")
        errprintf("  Possible causes: too few units/periods, too many instruments relative\n")
        errprintf("  to N, weak threshold identification, or no q variation in the grid.\n")
        errprintf("  Remedies: collapse, tighter maxlag(), larger trim(), or method(fod).\n")
        exit(498)
    }
    if (xdpt_verbose) printf("  γ̂ = %8.4f, obj = %8.4f\n", best_gamma, best_obj)

    if (n_refine > 0) {
        printf("  Final fixed-W2 refine: %g iteration(s), %g support point(s) added; complete=%g\n",
               ref_it, ref_added, ref_complete)
        if (best_twostep != 1 & st_local("nowarn") == "") {
            printf("{err}warning: refine() was not applied because the estimator fell back to one step.\n")
        }
        else if (ref_complete != 1 & st_local("nowarn") == "") {
            printf("{err}warning: final fixed-W2 refinement is incomplete; see e(refine_remaining)\n")
            printf("{err}         and e(refine_neigh_unevaluated).\n")
        }
        displayflush()
    }

    // The unit-resampling scheme remains available as an experimental
    // extension on unbalanced FD stacks, but the common-T panel theory does
    // not certify that extension. Diagnose it before any seed is applied or
    // bootstrap draw is taken.
    if (do_grid_ci & xdpt_boot_exact == 1) {
        real scalar _ubg, _ubb
        _ubb = .
        for (_ubg = 1; _ubg <= cols(cache_main); _ubg++) {
            if (!cache_main[_ubg].ok) continue
            if (cache_main[_ubg].gamma == best_gamma) {
                _ubb = xdpt2_is_strongly_balanced(cache_main[_ubg].uid,
                                                   cache_main[_ubg].times)
                break
            }
        }
        if (_ubb == 0 & st_local("nowarn") == "") {
            printf("{err}warning: boottype(unit) is running on an unbalanced effective FD stack;\n")
            printf("{err}         this is an experimental extension beyond the common-T theory.\n")
        }
    }

    // Grid bootstrap CI + linearity test + continuity test
    real scalar gam_lo, gam_hi, pval_lin, pval_cont, ci_empty, ci_nseg
    real scalar gci_adm, gci_lo, gci_hi, gci_eff_n, gci_eval, gb_minB
    real colvector gamma_ci_grid
    // v0.8.0 (#2): coefficient percentile bootstrap outputs
    real matrix bci
    real scalar bci_B, bci_2s, bci_fb, bci_skip
    real scalar bci_valid, bci_g1, bci_g2, bci_att
    bci = J(0, 0, .)
    bci_B = 0
    bci_2s = 0
    bci_fb = 0
    bci_skip = 0
    bci_valid = 0
    bci_g1 = 0
    bci_g2 = 0
    bci_att = 0
    // Initialize CI bounds to missing; overwritten below only when CI computed.
    // This allows downstream users to detect "no CI" via missing(e(gamma_lo)).
    gam_lo = .
    // v0.8.2 R11 (#3): CI-grid bookkeeping defaults (stay missing
    // under -noboot-)
    gci_adm = .
    gci_lo = .
    gci_hi = .
    gci_eff_n = .
    gci_eval = .
    gb_minB = .
    real matrix ci_tab_r, ci_seg_r
    real scalar ci_unres_r
    ci_tab_r = J(0, 0, .)
    ci_seg_r = J(0, 0, .)
    ci_unres_r = .
    // v0.9.36: citest(#) result (gamma, D, crit, accept, draws, status, p)
    real matrix cit_row
    real scalar seed_citest
    cit_row = J(0, 7, .)
    seed_citest = .
    gam_hi = .
    pval_lin = .
    pval_cont = .
    real scalar lin_valid, cont_valid, cont_common
    real scalar seed_threshold, seed_linearity, seed_continuity, seed_coefficient
    lin_valid = .
    cont_valid = .
    cont_common = .
    seed_threshold = .
    seed_linearity = .
    seed_continuity = .
    seed_coefficient = .
    ci_empty = .
    ci_nseg = .

    if (do_grid_ci) {
        if (xdpt_grid_quant == 1) {
            gamma_ci_grid = xdpt2_quantile_grid(q_eff, trim_rate/2,
                                                 1 - trim_rate/2, n_gridci)
        }
        else gamma_ci_grid = xdpt2_rangen(q_lo, q_hi, n_gridci)
        // v0.7.13 (audit R4): the inverted confidence set must contain the
        // point at which the test statistic is exactly zero, so that it can
        // never be empty by discretization alone (Gong-Seo property).
        // v0.9.34 (C1): with a two-step fit the wild inversion uses the
        // stage-2 criterion, whose zero is the reported gamma-hat (appended
        // below); only after a one-step fallback is the one-step argmin over
        // the estimation grid appended as well.
        real scalar _g1s_obj, _g1s_gamma, _gok, _gobj, _gl, _w2ci
        real colvector _gtheta
        real matrix _gV
        _w2ci = (xdpt_boot_exact != 1 & best_twostep == 1 & rows(best_A) > 0)
        _g1s_obj = .
        _g1s_gamma = .
        for (_gl = 1; _gl <= (_w2ci ? 0 :
             (xdpt_n_stage1 < . ? min((xdpt_n_stage1, cols(cache_main))) :
                                   cols(cache_main))); _gl++) {
            if (!cache_main[_gl].ok) continue
            if (rows(cache_main[_gl].dY) < 20) continue
            xdpt2_solve_gmm_1step_pre(cache_main[_gl].dY, cache_main[_gl].dW,
                                       *cache_main[_gl].pZ, cache_main[_gl].ZW,
                                       cache_main[_gl].ZY, *cache_main[_gl].pW1,
                                       _gok, _gtheta, _gobj, _gV)
            if (!_gok) continue
            if (_g1s_obj == . | _gobj < _g1s_obj) {
                _g1s_obj = _gobj
                _g1s_gamma = cache_main[_gl].gamma
            }
        }
        // v0.9.6 R22 (#3): the one-step argmin is appended for the WILD
        // (one-step) inversion, where D(gamma_1step) = 0 guarantees a
        // non-empty acceptance set. The unit inversion is two-stage: its
        // zero point is best_gamma (already appended below).
        if (xdpt_boot_exact != 1 & !_w2ci & _g1s_gamma < .) gamma_ci_grid = gamma_ci_grid \ _g1s_gamma
        gamma_ci_grid = sort(uniqrows(gamma_ci_grid \ best_gamma), 1)
        gci_eff_n = rows(gamma_ci_grid)
        // v0.9.2 R18 (#1): the unit bootstrap needs the reported two-step
        // pieces -- W_n_2 (best_A) for the two-stage sample statistic and
        // the UNRESTRICTED residuals at (gamma-hat, theta-hat) for the
        // bootstrap DGP (restricted coefficients + unrestricted residuals).
        real colvector resid_hat_v
        real scalar _rh_g
        resid_hat_v = J(0, 1, .)
        if (xdpt_boot_exact == 1) {
            if (best_twostep != 1) {
                errprintf("xtdpthresh: boottype(unit) requires the two-step estimator\n")
                errprintf("  (this run fell back to one-step)\n")
                exit(498)
            }
            for (_rh_g = 1; _rh_g <= cols(cache_main); _rh_g++) {
                if (cache_main[_rh_g].ok == 0) continue
                if (cache_main[_rh_g].gamma == best_gamma) {
                    resid_hat_v = cache_main[_rh_g].dY - cache_main[_rh_g].dW * best_theta
                    break
                }
            }
            if (rows(resid_hat_v) == 0) {
                errprintf("xtdpthresh: boottype(unit): gamma-hat not found in the cache\n")
                exit(498)
            }
        }
        // Keep the threshold-CI stream pinned to the historical rseed() draw
        // sequence; later inference objects receive component-specific seeds.
        seed_threshold = xdpt2_component_seed(0)
        // v0.9.34: the sample side of the inversion takes its minimum over
        // the initial grid, as every draw does (see xdpt2_grid_bootstrap)
        real scalar ci_gmin
        ci_gmin = ((best_twostep == 1 & search_o2_global < .) ?
                   search_o2_global : best_obj)
        xdpt2_grid_bootstrap(units, cache_main, gamma_grid, gamma_ci_grid,
                              q_eff, minreg_user,
                              ci_gmin, best_gamma,
                              method, flag_static, flag_kink,
                              t_min, t_max, n_boot, alpha,
                              gam_lo, gam_hi, ci_empty, ci_nseg,
                              gci_adm, gci_lo, gci_hi,
                              best_twostep, best_A, resid_hat_v, gb_minB,
                              ci_tab_r, ci_seg_r, ci_unres_r,
                              tpl_main, tpl_main_st)
        if (rows(ci_tab_r) > 0) {
            gci_eval = sum((ci_tab_r[., 6] :== 1) :| (ci_tab_r[., 6] :== 2))
        }
        if (xdpt_verbose) printf("  Grid CI = [%8.4f, %8.4f]\n", gam_lo, gam_hi)

        // v0.7.7: -notest- skips the linearity + continuity bootstraps
        // entirely (they are independent of the CI). For a coverage study that
        // only needs e(gamma_lo)/e(gamma_hi) this removes the dominant cost
        // after the CI itself is batched; pval_lin/pval_cont stay missing.
        if (!flag_notest) {
            seed_linearity = xdpt2_component_seed(104729)
            pval_lin = xdpt2_linearity_test(units, gamma_grid, cache_main,
                                             best_obj,
                                             method, flag_static, flag_kink,
                                             t_min, t_max, n_boot, lin_valid)
            if (xdpt_verbose) printf("  Linearity p-value = %6.4f\n", pval_lin)

            // Continuity test only when unrestricted (jump) model is estimated;
            // cache_main is then exactly the jump cache it needs.
            // v0.9.12 R31: note that refine() refines the JUMP search only.
            // The restricted KINK cache below is built on the same grid,
            // but the kink criterion varies continuously in gamma, so its
            // minimum remains a finite-grid approximation (documented in
            // the help; a dense numerical kink grid is future work).
            if (flag_cont_test) {
                seed_continuity = xdpt2_component_seed(224737)
                pval_cont = xdpt2_continuity_test(units, gamma_grid,
                                                    q_eff, minreg_user, cache_main,
                                                    best_obj,
                                                    method, flag_static,
                                                    t_min, t_max, n_boot,
                                                    cont_valid, cont_common,
                                                    best_twostep, best_A)
                if (xdpt_verbose) printf("  Continuity p-value = %6.4f\n", pval_cont)
            }
        }

        // Threshold-search-aware coefficient CIs. A component-specific seed
        // makes this replay invariant to whether earlier tests were requested.
        // v0.8.4 R13 (#2): attempted is stamped INSIDE the helper, only
        // once the bootstrap samples exist -- an early bail (no idx, bad
        // cache, no fast-path entries) attempted nothing, and stamping
        // n_boot out here would overstate it.
        // v0.9.3 R19 (#5): coefboot(none) -> zero requested draws (the
        // helper then attempts nothing and reports attempted = 0).
        external real scalar xdpt_coefboot_off
        if (xdpt_coefboot_off != 1) {
            seed_coefficient = xdpt2_component_seed(350377)
            bci = xdpt2_coef_bootstrap(cache_main, gamma_grid, best_gamma,
                                        best_theta, n_boot, alpha, best_twostep,
                                        bci_B, bci_2s, bci_fb, bci_skip,
                                        bci_valid, bci_g1, bci_g2, bci_att)
        }
        if (xdpt_verbose & bci_B > 0) {
            printf("  Coef bootstrap: B_eff = %g draws\n", bci_B)
        }

        // v0.9.36: citest(#), the grid-bootstrap test of H0: gamma = # with
        // the statistic, bootstrap and weight of the confidence set. Run
        // last, with its own component seed, so every other result is
        // unchanged by it. p = (1 + #{D* >= D}) / (1 + B), the add-one rule
        // whose "p > alpha" is exactly the set's "D <= crit".
        real scalar cit_g, cit_un, cit_mb
        real scalar cit_lo, cit_hi, cit_emp, cit_ns, cit_ad, cit_gl, cit_gh
        real matrix cit_tab, cit_seg
        real colvector cit_D
        cit_g = strtoreal(st_local("citest"))
        if (cit_g < .) {
            seed_citest = xdpt2_component_seed(477377)
            cit_mb = .
            cit_D = J(0, 1, .)
            xdpt2_grid_bootstrap(units, cache_main, gamma_grid, (cit_g),
                                  q_eff, minreg_user,
                                  ci_gmin, best_gamma,
                                  method, flag_static, flag_kink,
                                  t_min, t_max, n_boot, alpha,
                                  cit_lo, cit_hi, cit_emp, cit_ns,
                                  cit_ad, cit_gl, cit_gh,
                                  best_twostep, best_A, resid_hat_v, cit_mb,
                                  cit_tab, cit_seg, cit_un,
                                  tpl_main, tpl_main_st, cit_D)
            if (rows(cit_tab) == 1) {
                cit_row = (cit_tab, .)
                if (cit_tab[1, 6] == 2) cit_row[1, 7] = 1
                else if (cit_tab[1, 6] == 1 & rows(cit_D) > 0) {
                    cit_row[1, 7] = (1 + sum((cit_D :< .) :& (cit_D :>= cit_tab[1, 2]))) /
                                    (1 + sum(cit_D :< .))
                }
            }
        }
    }

    // === Count sample sizes and instruments at best γ ===
    real matrix dY_f, dW_f, Z_f
    real colvector times_f, uid_f
    // v0.9.35: with the joint variance, a template of this stack (and under
    // FOD of the FD stack of the AR test) gives the derivative column of the
    // AR statistics (xdpt2_ar_xg)
    external real scalar xdpt_tpl_rec
    struct xdpt2_stack_tpl scalar tp_f, tp_fd
    real scalar ar_try
    ar_try = (xdpt_kink_joint == 1)
    xdpt_ivc_diag = 1
    xdpt_tpl_rec = ar_try
    xdpt2_stack_at_gamma(units, best_gamma, method, flag_static, flag_kink,
                          t_min, t_max,
                          dY_f, dW_f, Z_f, times_f, uid_f)
    xdpt_tpl_rec = 0
    xdpt_ivc_diag = 0
    if (ar_try) {
        tp_f = xdpt2_tpl_build(units, method, flag_kink, dY_f, Z_f, times_f,
                               uid_f)
        if (!tp_f.ok) ar_try = 0
    }
    // v0.9.34: units whose regime changes within their equations at gamma-hat
    real scalar n_switch
    n_switch = xdpt2_n_switch(units, best_gamma, method, uid_f, times_f)
    // v0.9.30: instrument columns dropped as multiples of the per-period
    // constants in the estimation Z (read before the FD restack for the AR
    // test overwrites the counter).
    external real scalar xdpt_ivc_drop
    real scalar ivc_final
    ivc_final = xdpt_ivc_drop
    // v0.9.34: numerically dependent instrument columns dropped; v0.9.35:
    // those not linear combinations of the kept ones, and the largest
    // relative residual
    external real scalar xdpt_ivc_dep
    real scalar ivdep_final, ivdep_res_final, ivdep_near_final
    ivdep_final = xdpt_ivc_dep
    ivdep_res_final = xdpt_ivc_dep_res
    ivdep_near_final = xdpt_ivc_dep_near

    real scalar n_raw, n_trans, n_usable, n_iv, balanced_eff
    // v0.8.1: complete-case count = equation-eligible rows (rows(y) now
    // counts the full history sample).
    n_raw    = sum(eqv)
    n_usable = rows(dY_f)
    n_iv     = cols(Z_f)
    balanced_eff = xdpt2_is_strongly_balanced(uid_f, times_f)

    n_trans = n_usable

    // === Hansen J over-identification test ===
    // Hansen J uses the same second-step weight and criterion minimized by the
    // reported estimator. A first-step fallback has no efficient robust weight,
    // so its Hansen statistic is deliberately reported missing.
    real rowvector hj
    real scalar k_W_final, hansen_stat, hansen_df, hansen_p
    k_W_final = cols(dW_f)
    if (best_twostep) hj = xdpt2_hansen_j(best_obj, n_iv, k_W_final)
    else              hj = (., ., .)
    hansen_stat = hj[1]
    hansen_df = hj[2]
    hansen_p = hj[3]

    // === Arellano-Bond AR(1) and AR(2) tests on FD residuals ===
    // xtabond2 (Roodman 2009) convention: AR test is ALWAYS computed on
    // first-difference residuals, regardless of the transformation used for
    // estimation. This makes AR(1)/AR(2) interpretation consistent across
    // FD and FOD.
    // v0.7.2 (B1 FIX): full Arellano-Bond (1991, eq. 8) statistic including
    // the estimated-parameter variance terms (see xdpt2_ar_full). The
    // estimation-equation pieces (residuals/instruments/regressors at γ̂, the
    // weight A actually paired with θ̂, and V̂) feed Terms 2-3; for
    // method(fod) the test residuals are the FD restack while the
    // estimator pieces remain those of the FOD stack, linked within
    // unit via the c_i scalars.
    // v0.9.29: Term 3 uses the REPORTED variance of theta-hat (best_V), i.e.
    // the Windmeijer-corrected V under vce(windmeijer), as xtabond2 does with
    // twostep robust; the default vce(robust) V is unchanged. Certified
    // against xtabond2 at gamma fixed at gamma-hat (identical Z, W2, and
    // residuals): with the same V the AR(1)/AR(2) statistics coincide under
    // both FD and FOD (_dev_0928/tests/ar_diag.do, _dev_0929/tests).
    real matrix dY_fd_ar, dW_fd_ar, Z_fd_ar, X_ar
    real colvector times_fd_ar, uid_fd_ar, resid_trans, times_trans, uid_trans
    real colvector resid_est, dy_test
    resid_est = dY_f - dW_f * best_theta
    if (method == "fd") {
        resid_trans = resid_est
        times_trans = times_f
        uid_trans = uid_f
        X_ar = dW_f
        dy_test = dY_f
    }
    else {
        // FOD: re-stack with FD to get proper AR-test residuals
        xdpt_tpl_rec = ar_try
        xdpt2_stack_at_gamma(units, best_gamma, "fd", flag_static, flag_kink,
                              t_min, t_max,
                              dY_fd_ar, dW_fd_ar, Z_fd_ar, times_fd_ar, uid_fd_ar)
        xdpt_tpl_rec = 0
        if (ar_try) {
            tp_fd = xdpt2_tpl_build(units, "fd", flag_kink, dY_fd_ar, Z_fd_ar,
                                    times_fd_ar, uid_fd_ar)
            if (!tp_fd.ok) ar_try = 0
        }
        resid_trans = dY_fd_ar - dW_fd_ar * best_theta
        times_trans = times_fd_ar
        uid_trans = uid_fd_ar
        X_ar = dW_fd_ar
        dy_test = dY_fd_ar
    }

    real rowvector ar1, ar2, ar1_c, ar2_c, ar1_j, ar2_j
    real scalar ar1_stat, ar1_p, ar2_stat, ar2_p, ar_joint, K_base
    real colvector xg_est, xg_test
    // v0.9.32: under kink the reported V includes the estimation error of
    // gamma-hat; the statistics that treat gamma-hat as known (Terms 2-3
    // without a gamma column) use the conditional V. v0.9.35: in the jump
    // model too. They equal xtabond2's with gamma fixed at gamma-hat and are
    // kept in e(ar1_cond), e(ar2_cond).
    real matrix V_ar
    V_ar = best_V
    if (xdpt_kink_joint == 1 & rows(xdpt_V_cond_ar) == rows(best_V) &
        cols(xdpt_V_cond_ar) == cols(best_V)) V_ar = xdpt_V_cond_ar
    ar1_c = xdpt2_ar_full(1, resid_trans, X_ar, times_trans, uid_trans,
                          resid_est, Z_f, uid_f, dW_f, best_A, V_ar)
    ar2_c = xdpt2_ar_full(2, resid_trans, X_ar, times_trans, uid_trans,
                          resid_est, Z_f, uid_f, dW_f, best_A, V_ar)
    ar1 = ar1_c
    ar2 = ar2_c
    // v0.9.35: with the joint variance, the reported statistics include the
    // estimation of gamma-hat. The derivative column x_g of the residuals in
    // gamma (xdpt2_ar_xg: on the estimation rows, and under FOD on the FD
    // rows of the test) is appended to the test and estimation designs, and V
    // is the joint variance of (theta, gamma) of the reported type, so Terms
    // 2-3 treat gamma as a parameter. If a joint statistic cannot be formed
    // where the conditional one can, both conditional statistics are kept.
    ar_joint = 0
    xg_est = J(0, 1, .)
    xg_test = J(0, 1, .)
    if (ar_try & rows(xdpt_V_joint_full) == cols(dW_f) + 1 &
        cols(xdpt_V_joint_full) == cols(dW_f) + 1) {
        K_base = cols(units[1].X)
        xg_est = xdpt2_ar_xg(tp_f, best_gamma, best_theta, K_base, flag_kink,
                             dW_f, xdpt_gamma_bw)
        if (method == "fd") xg_test = xg_est
        else xg_test = xdpt2_ar_xg(tp_fd, best_gamma, best_theta, K_base,
                                   flag_kink, dW_fd_ar, xdpt_gamma_bw)
        if (rows(xg_est) == rows(dW_f) & rows(xg_test) == rows(X_ar) &
            rows(xg_est) > 0 & rows(xg_test) > 0) {
            ar1_j = xdpt2_ar_full(1, resid_trans, (X_ar, xg_test), times_trans,
                                  uid_trans, resid_est, Z_f, uid_f,
                                  (dW_f, xg_est), best_A, xdpt_V_joint_full)
            ar2_j = xdpt2_ar_full(2, resid_trans, (X_ar, xg_test), times_trans,
                                  uid_trans, resid_est, Z_f, uid_f,
                                  (dW_f, xg_est), best_A, xdpt_V_joint_full)
            if (!((ar1_j[1] >= . & ar1_c[1] < .) |
                  (ar2_j[1] >= . & ar2_c[1] < .))) {
                ar1 = ar1_j
                ar2 = ar2_j
                ar_joint = 1
            }
        }
    }
    // v0.9.35: jump model with the joint variance: distinct values of q
    // within two bandwidths of gamma-hat (a warning below 10)
    real scalar q_nvals_bw
    q_nvals_bw = .
    if (!flag_kink & xdpt_kink_joint == 1) {
        q_nvals_bw = xdpt2_q_nvals(units, best_gamma, xdpt_gamma_bw)
    }

    // v0.7.3 (D5 rework): persist the EXACT AR-test row set for -predict-.
    // Mata globals survive -restore-, so xtdpthresh_p can merge residuals
    // back by actual (panelvar, timevar) keys instead of re-deriving the
    // transformation/trim/zero-Z filters in ado (which cannot be done
    // exactly: the t-2 membership and B4 zero-instrument filters depend on
    // the estimation-time unit structs). uid in the stacked vectors is the
    // unit INDEX, so map through units[.].id to recover panelvar values.
    external real matrix xdpt_p_serial_m
    // v0.9.3 R19 (#4): per-run storage keyed by serial in asarrays so that
    // -estimates store/restore- workflows can predict from an OLDER run
    // in the same session (the single-global design errored on any
    // restored model). Entries beyond the 20 most recent runs are evicted;
    // a Mata clear or a restart still requires re-running (documented).
    external transmorphic xdpt_p_store_r, xdpt_p_store_e, xdpt_p_store_sig
    external transmorphic xdpt_p_store_token
    real colvector p_id_map, p_id_est
    real matrix p_store_r_m, p_store_e_m
    real scalar pr, pe, p_cache_sig
    // AR-test (FD) series: keyed (id, time) -> (Δy_FD, ê_FD)
    p_id_map = J(rows(uid_trans), 1, .)
    for (pr = 1; pr <= rows(uid_trans); pr++) {
        p_id_map[pr] = units[uid_trans[pr]].id
    }
    // Estimation-equation series: y/resid on the ACTUAL estimation rows
    // (FOD rows under method(fod); identical to the FD series under fd).
    p_id_est = J(rows(uid_f), 1, .)
    for (pe = 1; pe <= rows(uid_f); pe++) {
        p_id_est[pe] = units[uid_f[pe]].id
    }
    p_store_r_m = (p_id_map, times_trans, dy_test, resid_trans)
    p_store_e_m = (p_id_est,  times_f,     dY_f,     resid_est)
    // A 53-bit checksum over the cached outputs supplements the exact per-fit
    // token. The token supplies identity; the checksum is a corruption guard
    // only (Jenkins hash1 is not a cryptographic/digital signature). Neither
    // mechanism consumes the statistical RNG.
    p_cache_sig = hash1(p_store_r_m, 2147483647) * 4194304 +
                  mod(hash1(p_store_e_m), 4194304)
    if (rows(xdpt_p_serial_m) == 0) xdpt_p_serial_m = 0
    xdpt_p_serial_m = xdpt_p_serial_m[1, 1] + 1
    if (eltype(xdpt_p_store_r) == "real") xdpt_p_store_r = asarray_create("real", 1)
    if (eltype(xdpt_p_store_e) == "real") xdpt_p_store_e = asarray_create("real", 1)
    if (eltype(xdpt_p_store_sig) == "real") xdpt_p_store_sig = asarray_create("real", 1)
    if (eltype(xdpt_p_store_token) == "real") xdpt_p_store_token = asarray_create("real", 1)
    asarray(xdpt_p_store_r, xdpt_p_serial_m[1, 1], p_store_r_m)
    asarray(xdpt_p_store_e, xdpt_p_serial_m[1, 1], p_store_e_m)
    asarray(xdpt_p_store_sig, xdpt_p_serial_m[1, 1], p_cache_sig)
    asarray(xdpt_p_store_token, xdpt_p_serial_m[1, 1], st_local("_p_cache_token"))
    if (asarray_elements(xdpt_p_store_r) > 20) {
        real colvector _pks
        _pks = sort(asarray_keys(xdpt_p_store_r), 1)
        asarray_remove(xdpt_p_store_r, _pks[1])
        asarray_remove(xdpt_p_store_e, _pks[1])
        asarray_remove(xdpt_p_store_sig, _pks[1])
        asarray_remove(xdpt_p_store_token, _pks[1])
    }
    st_numscalar("r(xdpt2_p_serial)", xdpt_p_serial_m[1, 1])
    st_numscalar("r(xdpt2_p_sig)", p_cache_sig)

    // v0.7.5: exportgmm — expose the GMM pieces of the reported estimate for
    // external verification (e.g. re-computing the full AB AR statistic
    // outside the package). The externals are created/RESET on EVERY run:
    // 0x0 when the option is off (so stale matrices from a previous
    // exportgmm run can never leak), populated when on.
    //   xdpt_best_A    — moment weight paired with theta_hat/V_hat (k_iv x k_iv)
    //   xdpt_best_Z_f  — instruments at gamma_hat, estimation stack (n x k_iv)
    //   xdpt_best_X_f  — regressors at gamma_hat, estimation stack (n x k_par)
    //   xdpt_best_Xar  — FD AR-test-equation regressors (rows match
    //                    xdpt_p_resid); equals X_f under method(fd), the FD
    //                    restack under FOD — needed (with xdpt_p_resid
    //                    and xdpt_p_est) to reproduce e(ar*) externally for
    //                    those methods.
    //   v0.9.35, when the AR statistics include gamma-hat (e(ar_joint) = 1):
    //   xdpt_best_xg_f, xdpt_best_xg_ar -- their derivative columns on the
    //   estimation and test rows; xdpt_best_Vj -- the joint variance of
    //   (theta, gamma) they use.
    external real matrix xdpt_best_A, xdpt_best_Z_f, xdpt_best_X_f, xdpt_best_Xar
    external real matrix xdpt_best_Vj, xdpt_best_xg_f, xdpt_best_xg_ar
    xdpt_best_A   = J(0, 0, .)
    xdpt_best_Z_f = J(0, 0, .)
    xdpt_best_X_f = J(0, 0, .)
    xdpt_best_Xar = J(0, 0, .)
    xdpt_best_Vj    = J(0, 0, .)
    xdpt_best_xg_f  = J(0, 0, .)
    xdpt_best_xg_ar = J(0, 0, .)
    if (flag_exportgmm) {
        xdpt_best_A   = best_A
        xdpt_best_Z_f = Z_f
        xdpt_best_X_f = dW_f
        xdpt_best_Xar = X_ar
        if (ar_joint) {
            xdpt_best_Vj    = xdpt_V_joint_full
            xdpt_best_xg_f  = xg_est
            xdpt_best_xg_ar = xg_test
        }
    }
    ar1_stat = ar1[1]
    ar1_p = ar1[3]
    ar2_stat = ar2[1]
    ar2_p = ar2[3]
    real scalar ar1_b0, ar1_T1, ar1_TT, ar2_b0, ar2_T1, ar2_TT
    ar1_b0 = ar1[4]; ar1_T1 = ar1[5]; ar1_TT = ar1[6]
    ar2_b0 = ar2[4]; ar2_T1 = ar2[5]; ar2_TT = ar2[6]

    if (xdpt_verbose) {
        printf("  Hansen J=%6.3f (df=%g) p=%6.4f\n", hansen_stat, hansen_df, hansen_p)
        printf("  AR(1): m=%6.3f p=%6.4f   AR(2): m=%6.3f p=%6.4f\n",
               ar1_stat, ar1_p, ar2_stat, ar2_p)
    }

    // v0.9.31: release the drop-mask cache (a copy of Z)
    xdpt_ivc_Zref = J(0, 0, .)
    xdpt_ivc_mask = J(1, 0, .)
    xdpt_ivc_tref = J(0, 1, .)
    // v0.9.34: and the template records of the per-gamma caches
    xdpt2_tpl_release()

    // Return results
    st_rclear()
    // v0.7.5 hotfix (restored): persist serial AFTER st_rclear so e(p_serial)
    // survives for -predict-. The 0.7.6/0.7.7 speedup refactor (branched from
    // 0.7.3) dropped this re-set; without it st_rclear() wipes the scalar set
    // above and e(p_serial) returns missing, breaking predict (r301).
    st_numscalar("r(xdpt2_p_serial)", xdpt_p_serial_m[1, 1])
    st_numscalar("r(xdpt2_p_sig)", p_cache_sig)
    st_matrix("r(xdpt2_theta)", best_theta')
    st_matrix("r(xdpt2_V)",     best_V)
    // v0.9.32: kink joint variance flag and the conditional V it replaced
    st_numscalar("r(xdpt2_kink_joint)", xdpt_kink_joint)
    st_numscalar("r(xdpt2_kink_refined)", xdpt_kref)
    st_numscalar("r(xdpt2_gamma_bw)", (flag_kink ? . : xdpt_gamma_bw))
    if (xdpt_kink_joint == 1) st_matrix("r(xdpt2_V_cond)", xdpt_V_cond_ar)
    st_numscalar("r(xdpt2_gamma)", best_gamma)
    st_numscalar("r(xdpt2_obj)",   best_obj)
    st_numscalar("r(xdpt2_nused)",     n_usable)
    st_numscalar("r(xdpt2_n_raw)",     n_raw)
    st_numscalar("r(xdpt2_n_trans)",   n_trans)
    st_numscalar("r(xdpt2_n_iv)",      n_iv)
    st_numscalar("r(xdpt2_w1_fallback)", xdpt_w1_fallback)
    st_numscalar("r(xdpt2_iv_common)", ivc_final)
    st_numscalar("r(xdpt2_iv_dep)", ivdep_final)
    st_numscalar("r(xdpt2_iv_dep_res)", ivdep_res_final)
    st_numscalar("r(xdpt2_iv_dep_near)", ivdep_near_final)
    st_numscalar("r(xdpt2_n_units)",   rows(uniqrows(uid_f)))
    st_numscalar("r(xdpt2_n_switch)",  n_switch)
    st_numscalar("r(xdpt2_balanced_eff)", balanced_eff)
    st_numscalar("r(xdpt2_wind_applied)", xdpt_wind_applied)
    st_numscalar("r(xdpt2_bci_B)", bci_B)
    st_numscalar("r(xdpt2_bci_2s)", bci_2s)
    st_numscalar("r(xdpt2_bci_fb)", bci_fb)
    st_numscalar("r(xdpt2_bci_skip)", bci_skip)
    st_numscalar("r(xdpt2_bci_valid)", bci_valid)
    st_numscalar("r(xdpt2_bci_att)",   bci_att)
    st_numscalar("r(xdpt2_bci_g1)",    bci_g1)
    st_numscalar("r(xdpt2_bci_g2)",    bci_g2)
    st_numscalar("r(xdpt2_twostep)", best_twostep)
    st_numscalar("r(xdpt2_grid_req)", n_grid)
    st_numscalar("r(xdpt2_grid_eff)", rows(gamma_grid))
    st_numscalar("r(xdpt2_grid_adm)", rows(gamma_admitted))
    st_numscalar("r(xdpt2_grid_lo)",  grid_adm_lo)
    st_numscalar("r(xdpt2_grid_hi)",  grid_adm_hi)
    st_numscalar("r(xdpt2_minreg)",   (minreg_user > 0 ? minreg_user : 0))
    st_numscalar("r(xdpt2_grid_struct)", n_struct)
    st_numscalar("r(xdpt2_search_l1_n)", search_l1_n)
    st_numscalar("r(xdpt2_search_l2_n)", search_l2_n)
    st_numscalar("r(xdpt2_search_l3_n)", search_l3_n)
    st_numscalar("r(xdpt2_search_s1_level)", search_s1_level)
    st_numscalar("r(xdpt2_search_s2_level)", search_s2_level)
    st_numscalar("r(xdpt2_search_s1_n)", search_s1_n)
    st_numscalar("r(xdpt2_search_s2_n)", search_s2_n)
    st_numscalar("r(xdpt2_search_s1_same)", search_s1_same)
    st_numscalar("r(xdpt2_search_s2_same)", search_s2_same)
    st_numscalar("r(xdpt2_search_s1_gain)", search_s1_gain)
    st_numscalar("r(xdpt2_search_s2_gain)", search_s2_gain)
    st_numscalar("r(xdpt2_search_s1_conv)", search_s1_conv)
    st_numscalar("r(xdpt2_search_s2_conv)", search_s2_conv)
    st_numscalar("r(xdpt2_search_hit_max)", search_hit_max)
    st_numscalar("r(xdpt2_search_W2_builds)", search_W2_builds)
    st_numscalar("r(xdpt2_search_g1)", search_g1)
    // v0.9.29: the Windmeijer correction is derived for a second-step weight
    // built on the same regressor design as the final estimate; flag whether
    // the stage-1 threshold (where W2 was built) gives the reported design.
    // v0.9.34: in the jump model the design depends on gamma only through
    // the regime split, so the flag compares splits over the effective
    // support (refine() moves gamma-hat to a support point with the grid
    // point's split, and the flag used to drop to 0 with an identical fit);
    // the kink design changes with gamma itself. Under kink with the joint
    // variance the correction differentiates W2 in the stage-1 threshold as
    // well, so it does not require the two thresholds to agree (missing).
    // v0.9.35: the same holds for the jump model's joint variance.
    real scalar same_thr
    same_thr = .
    if (search_g1 < . & best_gamma < . & xdpt_kink_joint != 1) {
        if (flag_kink) same_thr = (search_g1 == best_gamma)
        else same_thr = (sum((q_eff :> min((search_g1, best_gamma))) :&
                             (q_eff :<= max((search_g1, best_gamma)))) == 0)
    }
    st_numscalar("r(xdpt2_same_threshold)", same_thr)
    st_numscalar("r(xdpt2_search_o1)", search_o1)
    st_numscalar("r(xdpt2_search_g2_global)", search_g2_global)
    st_numscalar("r(xdpt2_search_o2_global)", search_o2_global)
    st_numscalar("r(xdpt2_ref_it)",  ref_it)
    st_numscalar("r(xdpt2_ref_add)", ref_added)
    st_numscalar("r(xdpt2_ref_pool)", ref_pool_n)
    st_numscalar("r(xdpt2_ref_rem)",  ref_remaining)
    st_numscalar("r(xdpt2_ref_exh)",  ref_exhausted)
    st_numscalar("r(xdpt2_ref_lo)",   ref_alo)
    st_numscalar("r(xdpt2_ref_hi)",   ref_ahi)
    st_numscalar("r(xdpt2_ref_inb)",  ref_in_basin)
    st_numscalar("r(xdpt2_ref_nrem)", ref_neigh_rem)
    st_numscalar("r(xdpt2_ref_comp)", ref_complete)
    st_numscalar("r(xdpt2_ref_obj_gain)", ref_obj_gain)
    st_numscalar("r(xdpt2_grid_adm2)",   n_adm2)
    st_numscalar("r(xdpt2_grid2_lo)",    grid2_adm_lo)
    st_numscalar("r(xdpt2_grid2_hi)",    grid2_adm_hi)
    st_numscalar("r(xdpt2_minreg_def)",  minreg_def)
    st_numscalar("r(xdpt2_minreg_app)",  minreg_applied)
    if (rows(bci) > 0) st_matrix("r(xdpt2_bci)", bci)
    st_numscalar("r(xdpt2_q_lo)",      q_lo)
    st_numscalar("r(xdpt2_q_hi)",      q_hi)
    st_numscalar("r(xdpt2_gam_lo)", gam_lo)
    st_numscalar("r(xdpt2_gam_hi)", gam_hi)
    st_numscalar("r(xdpt2_gci_eff)", gci_eff_n)
    st_numscalar("r(xdpt2_gci_adm)", gci_adm)
    st_numscalar("r(xdpt2_gci_eval)", gci_eval)
    st_numscalar("r(xdpt2_gci_lo)",  gci_lo)
    st_numscalar("r(xdpt2_gci_hi)",  gci_hi)
    st_numscalar("r(xdpt2_ci_minB)", gb_minB)
    if (rows(ci_tab_r) > 0) st_matrix("r(xdpt2_ci_grid)", ci_tab_r)
    if (rows(ci_seg_r) > 0) st_matrix("r(xdpt2_ci_segments)", ci_seg_r)
    st_numscalar("r(xdpt2_ci_unres)", ci_unres_r)
    if (rows(cit_row) == 1) st_matrix("r(xdpt2_citest)", cit_row)
    st_numscalar("r(xdpt2_seed_citest)", seed_citest)
    st_numscalar("r(xdpt2_ci_empty)", ci_empty)
    st_numscalar("r(xdpt2_ci_nseg)",  ci_nseg)
    st_numscalar("r(xdpt2_pval_lin)", pval_lin)
    st_numscalar("r(xdpt2_pval_cont)", pval_cont)
    st_numscalar("r(xdpt2_lin_valid)",  lin_valid)
    st_numscalar("r(xdpt2_cont_valid)", cont_valid)
    st_numscalar("r(xdpt2_cont_common)", cont_common)
    st_numscalar("r(xdpt2_seed_threshold)", seed_threshold)
    st_numscalar("r(xdpt2_seed_linearity)", seed_linearity)
    st_numscalar("r(xdpt2_seed_continuity)", seed_continuity)
    st_numscalar("r(xdpt2_seed_coefficient)", seed_coefficient)
    st_numscalar("r(xdpt2_hansen)",    hansen_stat)
    st_numscalar("r(xdpt2_hansen_df)", hansen_df)
    st_numscalar("r(xdpt2_hansen_p)",  hansen_p)
    st_numscalar("r(xdpt2_ar1)",       ar1_stat)
    st_numscalar("r(xdpt2_ar1_p)",     ar1_p)
    st_numscalar("r(xdpt2_ar2)",       ar2_stat)
    st_numscalar("r(xdpt2_ar2_p)",     ar2_p)
    st_numscalar("r(xdpt2_ar1_b0)",    ar1_b0)
    st_numscalar("r(xdpt2_ar1_T1)",    ar1_T1)
    st_numscalar("r(xdpt2_ar1_TT)",    ar1_TT)
    st_numscalar("r(xdpt2_ar2_b0)",    ar2_b0)
    st_numscalar("r(xdpt2_ar2_T1)",    ar2_T1)
    st_numscalar("r(xdpt2_ar2_TT)",    ar2_TT)
    st_numscalar("r(xdpt2_ar1_np)",    ar1[2])
    st_numscalar("r(xdpt2_ar2_np)",    ar2[2])
    st_numscalar("r(xdpt2_ar1_nclust)", ar1[7])
    st_numscalar("r(xdpt2_ar2_nclust)", ar2[7])
    // v0.9.35
    st_numscalar("r(xdpt2_ar_joint)", ar_joint)
    st_numscalar("r(xdpt2_ar1_cond)", ar1_c[1])
    st_numscalar("r(xdpt2_ar2_cond)", ar2_c[1])
    st_numscalar("r(xdpt2_q_nvals_bw)", q_nvals_bw)
}

end
// End of xtdpthresh.ado


* ============================================================================
* CHANGELOG  (plain comments -- not shown by -which-; full detail in SSC history)
* ----------------------------------------------------------------------------
* 0.8.2   10jul2026 ROUND-9 audit release. GUARD (#2, result-changing for
*                   high-K/small-support designs): the per-gamma regime
*                   guard is now a pure TRIMMING rule on the effective
*                   support -- max(ceil(trim*n/2), 2), user-overridable via
*                   new minregime(#) -- the old cols(dW)+1 per-side floor
*                   demanded 2K+2 observations each side and could exclude
*                   valid thresholds before the objective was evaluated
*                   (identification is rank-based; conditioning is enforced
*                   by cond/fast_ok/stage guards). MEMORY (#4): q_support
*                   uses per-unit markers (O(NT)) instead of the R8 key
*                   multiset (~O(NT^2) preallocation and a per-row deep copy
*                   of the unit struct); the support SET is identical.
*                   HYGIENE (#7): q_supp/minregime passed as ARGUMENTS to
*                   the cache builders (no external Mata state); dead
*                   q_cur/q_ref work removed from the cache path; stale
*                   "one-step" comments fixed. SCOPE (#3, documented): if/in
*                   bounds the history sample too -- excluded rows are not
*                   instrument sources; history(panel) is future work.
*                   e(cmdversion)=0.8.2.
*                   R10 REFINEMENTS: ADMITTED-GRID tracking -- e(gamma_grid_
*                   lo/hi) (span of cache-admitted points), e(grid_requested/
*                   effective/admitted), run-time "requested/distinct/
*                   admitted" line, and the boundary-pin warning now
*                   references the ADMITTED span (minregime/ties/rank pruning
*                   can make an interior-looking endpoint the true search
*                   edge). minregime() is a FLOOR (max with the default trim
*                   rule, option semantics match its name) with a fail-fast
*                   validation when 2*minregime exceeds the support. NEW
*                   gridsample(effective|observed): observed = current-row q
*                   of retained equation rows (closer to the xthenreg
*                   convention for replication); docs no longer imply the
*                   default reproduces xthenreg's grid. Coefficient-bootstrap
*                   mixture POLICY: <=1% fallback note, 1-5% warning, >5%
*                   e(b_bootci) SUPPRESSED (a one-step/two-step mixture is no
*                   single estimator's distribution); e(boot_coef_attempted/
*                   fallback_rate/suppressed) stored. Comment/doc exactness:
*                   td-FOD is exact on ANY panel (stale approximation claims
*                   removed in ado+help); if/in-history comment rewritten;
*                   xdpt2_q_at_rows live again via gridsample(observed).
* 0.8.1   10jul2026 ROUND-6 audit release (result-changing; all anchors
*                   re-pinned). SAMPLE SPLIT (#1, adjudicated against the
*                   xthenreg source, which reshapes the FULL N x T panel and
*                   builds L.y internally so y_i1 IS an instrument): the
*                   command now keeps every in-scope row as a HISTORY row
*                   (instrument source, value-guarded) while GMM equations
*                   form only on strict complete-case EQUATION rows (u.eq).
*                   Restores y_i1-type instruments and the earliest FD
*                   equations that keep-if-touse used to delete. FOD forward
*                   means, level equations, min-row gates, the td-fod dummy
*                   operator, and e(N_raw) are all defined on equation rows,
*                   so estimators are IDENTICAL when no history-only rows
*                   exist. GRID SUPPORT (#2): trim bounds and quantile grids
*                   now use the FULL q history (FD steps at q_t AND q_{t-1};
*                   FOD at future q; xthenreg's grid_con convention), not
*                   only current-row q. COEF BOOTSTRAP (#3): each draw now
*                   replays the REPORTED two-step estimator (stage-1 argmin
*                   -> Omega* -> W2* -> stage-2 grid pass) by default;
*                   coefboot(onestep) keeps the fast replay. (#5) SYMMETRIC
*                   percentile intervals are the default (Gong-Seo report
*                   raw percentile under-coverage); coefcitype(percentile)
*                   available. (#6) CENTERED moment covariance (Seo-Shin
*                   eq. 11 / xthenreg) is now the DEFAULT; -nocenter- gives
*                   the xtabond2 convention for diagnostic cross-checks.
*                   (#7) centered-covariance denominator counts CONTRIBUTING
*                   clusters (uniqrows), not max(unit_id). iv_avail now also
*                   requires a non-missing instrument VALUE.
*                   R7 HOTFIXES (same-day audit of the R6 changes): the
*                   unit prefilter now requires only >=2 equation rows (a
*                   T=4 dynamic panel legitimately contributes FD equations
*                   at t=3,4 with y_1,y_2 instruments; length-based
*                   prefiltering also selected units systematically);
*                   grid/trim support = q values entering the RETAINED
*                   transformed rows (q_t and q_{t-1} under FD; + future
*                   equation-row q under FOD) -- neither current-row-only
*                   (pre-R6, missed breakpoints) nor full-history (R6 first
*                   cut, stretched by rows that never touch the criterion);
*                   level-equation iv_here moved INSIDE the value guards
*                   (partially-missing history rows could otherwise keep
*                   constant-only level rows); coefficient bootstrap honours
*                   best_twostep (one-step fallback estimates replay
*                   one-step), reports its composition via
*                   e(boot_coef_twostep)/e(boot_coef_fallback)/
*                   e(boot_grid_skipped) and warns when >5% of draws mixed
*                   estimators; e(cmdversion)=0.8.1; b_bootci wording no
*                   longer claims continuity robustness.
*                   R8 HOTFIXES (same-day audit of R7): xdpt2_q_support()
*                   moved AFTER the struct definition (it was declared
*                   before struct xdpt2_unit and would fail to compile);
*                   support deduplicated by (unit,time) KEY -- the R7
*                   multiset overweighted late-period observations in FOD
*                   forward means and biased quantile trims under trending
*                   q; the cache regime guard now counts the SAME support
*                   (current-row q alone can mis-declare an empty regime
*                   when only lagged/future indicators move the design);
*                   stack_at_gamma tags rows (1=transformed, 2=level) so
*                   system level rows contribute current-q only to the
*                   support; e(coefboot) reports the replay actually used
*                   (request preserved in e(coefboot_requested);
*                   e(estimator_twostep) stored).
* 0.8.0   10jul2026 ROUND-5 audit release (breaking version bump: the 0.7.13
*                   line had accumulated result-changing fixes under one
*                   version string -- a reproducibility hazard by itself).
*                   INFERENCE LABELLING: e(vcetype)="Conditional on estimated
*                   threshold"; run-time note that analytic slope SEs treat
*                   gamma-hat as fixed and are not continuity-robust; the
*                   Windmeijer correction explicitly labelled conditional on
*                   the selected threshold; e(vce) = VCE actually delivered,
*                   e(vce_requested) preserves the request. Full-Jacobian
*                   (G_gamma kernel) joint VCE and the Gong-Seo coefficient
*                   bootstrap are the documented roadmap.
*                   BOOTSTRAP METADATA: e(bootstrap_method)/e(resampling_
*                   unit)/e(moment_recentering) stored; package description
*                   corrected (wild scheme = computational approximation of
*                   Alg. 1, validity via seeded MC, not the exact theorem).
*                   HANSEN J: df = L - k_W - 1 (gamma is estimated too);
*                   chi-square reference documented as diagnostic under
*                   continuity. NEW OPTION center: centered clustered moment
*                   covariance (Seo-Shin eq. 11 / xthenreg convention);
*                   default stays uncentered (AB/xtabond2). td: run-time
*                   approximation note for method(fod) on unbalanced panels
*                   (exact under fd; exact under fod when balanced);
*                   singleton time cells are DROPPED after FWL demeaning
*                   (they carried zero information but inflated n_rows and
*                   the residual pools). TIER NOTE printed for fod/system
*                   covering analytic VCE, Hansen J, AR and the CI. DISPLAY:
*                   sample line shows obs-used (=e(N)) vs complete-case vs
*                   stacked. e(cmdversion) stored. API CLEANUP (pre-release):
*                   citype() REMOVED (duplicated noboot exactly); tdpurge
*                   REMOVED (the legacy pre-demeaning construction the review
*                   identified as not-time-dummies; td/FWL is the only
*                   treatment; system+td stays blocked).
* 0.7.13  03jul2026 audit fixes (5-agent code audit). VALIDATION: reject
*                   duplicate variables WITHIN one regressor group / iv()
*                   (they survive -syntax- via ts aliases like L.x l1.x and
*                   silently produce collinear columns); reject the depvar as
*                   its own regressor/instrument/threshold (endogenous(y),
*                   qx(y)); panel-only xtset now gets the right error; boot()
*                   no longer validated under noboot; auto-L.y duplicate check
*                   is case-sensitive on the variable name. PREDICT: -r- now
*                   abbreviates residuals (was: silently matched regime);
*                   regime requires -reg-. INFERENCE (result-changing on edge
*                   cases): CI/linearity bootstraps are fast-1-step-only on
*                   BOTH sample and bootstrap sides (removed 2-step fallbacks
*                   that searched a larger gamma set under an incomparable
*                   cluster-Omega objective, inflating crit); transformed-eq
*                   zero-IV row filter now actually fires (block constant
*                   written AFTER the filter, mirroring the level equation);
*                   history gate anchors on any lag row in [hist_req, lag_hi]
*                   instead of exactly t-hist_req (gapless panels unchanged;
*                   gapped panels keep rows with valid deeper-lag IVs and
*                   drop instrument-less rows). DISPLAY/RETURN: notest no
*                   longer prints "p = ."; e(level) stored; boundary-pin
*                   guard also checks missing q_lo/q_hi; xdpt2_solve_gmm now
*                   uncalled (kept as reference).
*                   FOLLOW-UP audit of the 0.7.13 fixes: (perf) the B6 history
*                   gate is bounded by the unit's earliest observation, not by
*                   xdpt_lag_hi (=9999 default) -- the unbounded scan was a
*                   large-panel regression; (correctness) the B5 zero-IV row
*                   filter, and the pre-existing level-equation BUG 4a filter,
*                   now drop on a structural-availability MASK (iv_avail /
*                   iv_here) instead of rowsum(abs(Z)) -- a valid instrument
*                   equal to 0 or on a tiny scale is no longer treated as
*                   absent; (fix) panel-only xtset now matched (r(timevar) is
*                   "." not ""); gridci()/boot()/rseed() gated on whether a
*                   bootstrap will actually run (citype(grid) & !noboot), so
*                   citype(none) boot(0) is accepted and rseed() no longer
*                   perturbs the RNG on point-estimate-only calls.
*                   ROUND-3 audit: (validation) a static model with no
*                   regressors is rejected with a clear message -- the
*                   availability mask would otherwise drop the whole sample;
*                   (correctness) the all-zero instrument-COLUMN drop now uses
*                   an exact-zero test (colsum > 0) instead of a 1e-12 floor,
*                   so a valid tiny-scaled instrument column is not deleted;
*                   (doc) the xtdpthresh_p header syntax line reflects the
*                   Residuals/REGime abbreviations.
*                   ROUND-4 audit: (CI) the test-inversion grid now always
*                   contains BOTH the reported gamma-hat and the 1-step
*                   argmin over the estimation grid (where the inversion
*                   statistic is exactly 0), so the confidence set can no
*                   longer be empty purely by discretization (Gong-Seo
*                   property); (memory) per-variable instrument lag windows
*                   are capped at the maxlag() upper bound instead of the
*                   full global time span -- results identical (the excess
*                   columns were all-zero and deleted post hoc), allocation
*                   down from O(T^2 K) toward O(T*lag_hi*K) uncollapsed;
*                   (semantics) e(N) = raw panel-time obs in e(sample) on all
*                   methods; stacked equation rows move to new e(N_stack);
*                   (honesty) td documented and announced as a PRE-purging
*                   convention, NOT equivalent to time dummies in the
*                   threshold model (interactions are formed from purged
*                   inputs; regime intercept not demeaned) -- an FWL-correct
*                   per-gamma demeaning is future work; bootstrap comments
*                   relabelled as the xthenreg-style fast cluster wild
*                   residual bootstrap, not Gong-Seo exact Algorithm 1.
*                   NEW OPTION vce(uncorrected|windmeijer): default reports
*                   the uncorrected asymptotic two-step cluster-robust
*                   sandwich (unchanged behaviour); vce(windmeijer) applies
*                   the Windmeijer (2005) finite-sample correction (robust
*                   variant V_c = V2 + DV2 + V2D' + DV1rD', Omega-derivative
*                   at the stage-1 residuals that built W2), computed once at
*                   the final estimate -- no effect on grid search/bootstrap.
*                   e(vce) records the request; e(vce_applied) flags whether
*                   the correction replaced e(V) (two-step path only).
*                   td REDESIGNED (FWL-correct): td now partials common-
*                   across-regime time dummies out of the STACKED system --
*                   dY, every dW(gamma) column (incl. 1(q>gamma) and the
*                   interactions), and every Z column are demeaned within
*                   each time cell after stacking (xdpt2_demean_bytime);
*                   algebraically = dummies in both W and Z + FWL. Exact
*                   under FD; exact under FOD on balanced panels. Blocked
*                   for method(system) (level constant collinear). The old
*                   pre-purging behaviour survives as tdpurge (legacy, with
*                   note); e(td_mode) = "fwl"|"purge". Instrument columns
*                   constant within every time cell are snapped to zero
*                   (relative scale test) and dropped.
*                   NEW OPTION gridtype(uniform|quantile): quantile places
*                   grid points on empirical quantiles of q over the
*                   effective sample (both estimation and CI grids; ties
*                   collapsed); default uniform preserves the xthenreg
*                   convention. xdpt2_solve_gmm() (orphaned) deleted.
*                   MEMORY (C5): the gamma-invariant Z (n x L) and W_first
*                   (L x L) are now HEAP OBJECTS shared across cache entries
*                   via pointers (pZ/pW1) under the existing bitwise reuse
*                   guards -- Mata struct assignment real-copies matrices
*                   (verified: 50 assignments of an 80MB matrix peaked at
*                   ~4GB), so the by-value cache duplicated Z across every
*                   grid/CI/kink entry. Cache memory drops from O(G x n x L)
*                   to O(n x L) for the dominant objects; values bit-for-bit
*                   unchanged (pointer deref of the same object).
* 0.7.12  03jul2026 usability: under kink, warn when qx() is absent from the
*                   regressor list / endogenous() / predetermined() -- without
*                   the baseline q slope the level term is a one-sided hinge,
*                   not a two-sided kink. Warning only (respects nowarn); no
*                   change to estimation. Renamed -nosearch- to -noboot-: the
*                   old name wrongly implied it skipped the gamma grid search,
*                   which always runs; it only turns off the bootstrap CI and
*                   the linearity/continuity tests (point estimate only). Hard
*                   rename (pre-release, no alias). Help/paper clarify qx()-in-
*                   RHS, the endogenous(q) semantics, and noboot. method() and
*                   citype() are now case-insensitive (method(FOD) accepted).
*                   Post-estimation notes (respect nowarn): #IV > #units, #units
*                   < 30, and boot() < 999 for publication. Kink note prints γ
*                   directly (not raw SMCL {&gamma}); dead helper
*                   xdpt2_recompute_cluster_j() removed.
* 0.7.11  02jul2026 speedup/cleanup, no result change: gamma-invariant q
*                   vector reused under the exact times/uid guard instead of
*                   an interpreted per-row rebuild per gamma; dead qrow struct
*                   member and dead pos_tm1 lookup dropped; new e(q_lo)/
*                   e(q_hi) expose the effective-sample trim bounds.
* 0.7.10  02jul2026 audit fixes: immutable q under td; effective-sample trim
*                   and regime guard; correct static/FOD history and exogenous
*                   instruments; exact e(sample)/N_units; reject delta!=1;
*                   fixed-weight sandwich V and paired AR weight; one-step
*                   continuity DGP; TS-overlap/auto-L.y/name guards.
* 0.7.9.1 02jul2026 bugfix: xdpt2_ar_full failure return was 3 elements but
*                   callers read [4]..[6] -- panels where the AR test cannot
*                   be computed (e.g. AR(2) with zero lag-2 pairs on T=5)
*                   crashed the whole command with Mata 3301 instead of
*                   reporting AR as missing. Pre-existing since 0.7.2.
* 0.7.9  02jul2026  speedup: exact-guarded reuse of gamma-invariant work
*                   (CI-min hoist, W_first & Z'Y reuse, cached cross-products)
* 0.7.8  02jul2026  speedup: large-N cache build (preallocated stacking,
*                   run-based cluster-moment aggregation)
* 0.7.7  10jun2026  speedup: batched linearity/continuity bootstraps; -notest-
* 0.7.6  10jun2026  speedup: batched grid-bootstrap CI (typically 10-40x)
* 0.7.5  10jun2026  -exportgmm-; e(ar*_np) full-formula path flags
* 0.7.4  10jun2026  -predict- returns estimation-eq residuals; -arresiduals-
* 0.7.3  10jun2026  -predict- merge architecture (exact AR-test rows)
* 0.7.2  10jun2026  full Arellano-Bond (1991) AR m-statistic
* 0.7.1  10jun2026  hotfix: method(system) conformability
* 0.7.0  10jun2026  audit: e(sample), system level constant, unified grid
*                   search, add-one bootstrap p-values, rseed(), e(cmdline)
* 0.6.1  26apr2026  doc fixes
* 0.6.0  25apr2026  initial SSC release
* ============================================================================

* ---------------------------------------------------------------------------
* v0.8.2 R11 (audit round 11):
*   #1  BLOCKER: bci_rate/bci_suppressed were used by the post-estimation
*       display block before being defined (an undefined local expands to
*       nothing, so -if `bci_suppressed'- became -if {- whenever the coef
*       bootstrap succeeded). Both are now computed at retrieval time, and
*       the display decides suppression BEFORE advertising e(b_bootci).
*   #2  "admitted" now means ok & fast_ok (the stage-1 searchable set; the
*       solver rejects !fast_ok points with the same conditioning test).
*       New: e(grid_structural) (ok-only) and e(grid_twostep_admitted)
*       (solvable under W_n_2; missing on one-step paths). Display line
*       reports requested/distinct/structural/admitted.
*   #3  boundary-pin warning now compares the CI endpoints against the CI
*       grid's OWN admitted span (new e(gamma_ci_grid_lo/hi),
*       e(gridci_requested/effective/admitted)), not the estimation grid's:
*       the two grids can admit different ranges.
*   #4  coefboot stage-2 replay searches ALL structurally-ok grid points
*       (conditioning tested per point under W2_b), matching the main
*       two-step search, instead of only the stage-1 fast list. New
*       e(boot_grid_stage1/stage2); draws whose stage 2 finds no point
*       remain counted in e(boot_coef_fallback).
*   #5  Windmeijer derivative gains the centering term
*       +(h_j' s + s' h_j)/(n_contrib*n_rows) when center (default) is on,
*       matching xdpt2_build_cluster_omega's Omega. vce(windmeijer) SEs
*       change slightly under the default; identical under -nocenter-.
*   #6  gridsample(observed) deduplicates (unit,time) before extracting q:
*       method(system) stacks each observation in a transformed AND a level
*       row, double-counting q in bounds/grid/minregime. No-op under fd/fod.
*   #7  reproducibility metadata: e(minregime_default/applied), e(trim),
*       e(gridtype), e(gridsample), e(gridci_*) -- grid config no longer
*       recoverable only from e(cmdline).
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.8.3 R12 (audit round 12):
*   #1  BLOCKER: with 1-9 successful coef-bootstrap draws, bci_B > 0 but
*       r(xdpt2_bci) does not exist -> -matrix ... = r(xdpt2_bci)- died.
*       New valid flag (e(boot_coef_valid)) gates every matrix consumer.
*   #2  BLOCKER (R11 regression): coefboot's in-helper r() exports were
*       wiped by xdpt2_run's st_rclear() before the export block --
*       e(boot_grid_stage1/2) were 0 on every run. Now threaded as output
*       arguments and exported after st_rclear().
*   #3  removed the fast_ok bail at the reported gamma-hat: since R11 the
*       two-step search selects over ALL ok points, so gamma-hat-2 can
*       legitimately lack the one-step fast path; the DGP needs only
*       dY/dW/uid/pZ there. One-step gamma-hat has fast_ok by construction.
*   #4  stage-labeled grid spans: e(gamma_grid1_lo/hi) (one-step search
*       space, ex gamma_grid_lo/hi) and e(gamma_grid2_lo/hi) (W_n_2-solvable
*       span from the stage-2 loop) -- a two-step gamma-hat can lie outside
*       the stage-1 span without being an error.
*   #5  attempt accounting: e(boot_coef_requested/attempted/success/valid);
*       attempted no longer reads 0 when all draws failed.
*   #6  version bump 0.8.2 -> 0.8.3 (R11+R12 change results/metadata).
*   #7  defensive markout of panelvar/timevar on the master sample before
*       the history copy (novarlist left the join keys unmarked).
*   #8  fixed stale comment: under system GMM e(N_stack) counts stacked
*       rows; e(N) remains the raw panel-time union.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.8.4 R13 (audit round 13):
*   #2  e(boot_coef_attempted) is stamped inside xdpt2_coef_bootstrap only
*       once ETA/Ymat exist and the draw loop is committed -- early bails
*       (no idx, bad cache entry, no fast-path grid points) now correctly
*       report attempted = 0 instead of n_boot.
*   #3  e(boot_coef_suppressed) requires bci_valid: with 1-9 successful
*       draws there is no CI to suppress. The mixture-policy locals are
*       computed after bci_valid retrieval (they previously preceded it).
*   AR  replaced the stale pre-release TODO with the parity record: FD path
*       certified vs -abar- to 1e-6 (suite M.6); fod/system documented as an
*       extension with no external reference implementation.
*   version bump 0.8.3 -> 0.8.4.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.8.5 R14 (audit round 14):
*   #1  method(system): the two-equation-row floor in xdpt2_build_units
*       wrongly excluded units that contribute only a LEVEL equation (one
*       eligible row + instrument history). Floor is now method-dependent
*       (system: 1; fd/fod: 2). Can change system-GMM samples and gamma-hat
*       on short/irregular panels.
*   #2  method(system) with zero usable level equations now hard-errors
*       (498) instead of silently reporting an FOD-only fit as system.
*   #3  fod+td: rows saturated by the time-dummy projection (leverage ~= 1
*       on unbalanced panels) are dropped via a gamma-invariant leverage
*       mask -- they carried no moments but inflated e(N_stack), the Omega
*       normalization, residual pools, and AR/bootstrap diagnostics.
*   #4  e(boot_coef_requested) = 0 under -noboot- (the boot() default no
*       longer masquerades as a request).
*   version bump 0.8.4 -> 0.8.5.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.8.6 R15 (audit round 15):
*   #1  method(system) can no longer degenerate to LEVEL-ONLY either (the
*       R14 min_eq=1 floor made that reachable when no unit has two
*       complete rows): system now hard-requires both equation blocks.
*   #2  the check fires EARLY in Mata on the initial effective stack
*       (eqtype_eff; row availability is gamma-invariant) -- before the
*       grid build, estimation, grid bootstrap, and coefficient bootstrap
*       -- instead of erroring after all of them ran. The ado-side check
*       remains as an internal-regression backstop on the final stack.
*   version bump 0.8.5 -> 0.8.6.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.8.7 R16 (audit round 16):
*   #4  fod+td: the time-dummy matrix is allocated by OBSERVED time
*       columns instead of the raw numeric span (delta-1 daily-date
*       indexes could demand a rows x 10^4 matrix for 30 observed
*       periods). Two-pass build; compact column set equals the old
*       post-_dcols set, so the projection is bit-for-bit unchanged.
*   #5  per-block unit participation: e(N_units_trans/level/both), plus a
*       level-dominated warning under method(system) when fewer than
*       max(5, 10% of level units) contribute transformed equations.
*   version bump 0.8.6 -> 0.8.7.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.0 (boottype(exact)):
*   New option boottype(wild|exact). wild (default) keeps the fast cluster
*   wild residual bootstrap. exact implements Gong-Seo (2026) Alg. 1-style
*   resampling for the threshold-CI test inversion: contributing units
*   drawn iid with replacement; per-unit moments recentered at the
*   restricted (gamma_l, theta_r) fit so H0 holds in the bootstrap world;
*   weight matrix recomputed per draw from the recentered resampled
*   moments; candidate set and 1-step D functional identical to D_sample.
*   Applies to the threshold CI only -- linearity/continuity tests and the
*   coefficient bootstrap keep the wild scheme. e(boottype) records the
*   choice. Roughly 30-60x the wild runtime (documented; note printed).
*   RNG: draws honor rseed(); wild-path results are bit-for-bit unchanged.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.1 R17 (audit round 17):
*   #2  method(system) cluster-participation gates on the initial stack:
*       >= 5 transformed clusters, >= 5 level clusters, and >= 1 unit in
*       BOTH blocks (hard errors); plus symmetric warnings (level-dominated
*       / thin-level / low-overlap) on the final-stack counts.
*   #3  span-free time indexing everywhere: global sorted observed-eq-time
*       vector xdpt_teq + binary-search rank (xdpt2_tpos). Instrument
*       blocks in transform_unit/level_unit and td-fod dummy columns are
*       rank-indexed -- no allocation anywhere scales with the raw
*       calendar span (R16's fix still used a span-length marker vector,
*       and the IV blocks were untouched and bigger). Gap-free index:
*       rank == offset, bit-for-bit unchanged.
*   #4  history(panel|sample): panel (default) = if/in restricts the
*       equation sample only, matching the pre-restriction materialization
*       of L.y / ts-operator regressors; sample = hard boundary, with the
*       auto L.y nulled where its source is out-of-history (user ts terms
*       documented as non-retro-restrictable). e(history).
*   #5  (rebutted) the 0.9.0 build already carried *! 0.9.0 and
*       e(cmdversion)=0.9.0; the report cited a stale copy.
*   #6  Difference-in-Hansen for the level block under method(system):
*       e(diffhansen_level/_df/_p), J_sys - J_fod at gamma-hat, each with
*       its own efficient weight; df = L_sys - L_fod - 1; displayed under
*       the Hansen line; missing on one-step paths.
*   version bump 0.9.0 -> 0.9.1.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.2 R18 (audit round 18):
*   #1  boottype(exact) -> boottype(unit) and REWRITTEN per the reviewer's
*       Alg. 1 specification: two-stage sample statistic under the run's
*       W_n_2; DGP = restricted 2-stage coefficients + UNRESTRICTED
*       residuals at (gamma-hat, theta-hat); recentering subtracts the
*       sample moment at theta-hat; per draw stage-1 argmin -> bootstrap
*       residuals -> recentered per-unit Omega* -> W2* -> stage-2
*       restricted/grid-min. fd-without-kink only; requires the two-step
*       path; labeled experimental/not-certified everywhere.
*   #2  Difference-in-Hansen: the FOD side is now FULLY re-estimated with
*       its own grid search (J_fod = min_gamma), fixing the downward bias
*       of pinning at the system gamma-hat. e(hansen_fod/_df/_p),
*       e(gamma_fod) expose the reduced fit. #2.1: negative C no longer
*       clamped -- p missing + e(diffhansen_negative)=1 + warning.
*   #3  bootstrap metadata split per inference object:
*       e(threshold_bootstrap/_resampling/_recentering),
*       e(coefficient_bootstrap), e(linearity_bootstrap),
*       e(continuity_bootstrap).
*   #4  history(sample) now REJECTS user ts-operator terms (they are
*       materialized on the full panel and cannot be retro-restricted);
*       the auto L.y null-enforcement stays.
*   user reports: per-point valid-draw floor (>=10) and tracking
*       (e(gridboot_min_draws) + note); floor(u*n)+1 index draw (ceil
*       could map u==0 to index 0); duplicate dh_df line NOT present in
*       this build (stale copy); dev file renamed xtdpthresh_dev.ado.
*   version bump 0.9.1 -> 0.9.2.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.3 R19 (audit round 19):
*   #4  predict storage keyed by serial in asarrays (20 most recent runs):
*       -estimates store/restore- + predict now works within a session.
*   #5  coefboot(none) added; coefboot(gs) rejected with an honest message
*       (the wild scheme is NOT the Gong-Seo coefficient bootstrap);
*       e(coefficient_bootstrap) says so explicitly.
*   #6  coefficient bootstrap NEVER mixes estimators: failed two-step draws
*       are discarded and redrawn (chunked, up to 5x the request); the old
*       one-step fallback substitution and the >5% suppression policy are
*       gone. e(boot_coef_failed)/e(boot_coef_fail_rate) replace
*       fallback/fallback_rate/suppressed. Zero-failure runs reproduce the
*       old draws bit-for-bit (same RNG order in the first chunk).
*   #7  validity floors: coefficient CI requires >= 90% of the request and
*       >= 10 valid draws (warning under 100); linearity and continuity
*       p-values are missing below the same 90%/>=10 floor (they previously
*       accepted a single surviving draw); the grid CI already refuses
*       points below 10 valid draws (R18) and reports e(gridboot_min_draws).
*   #8  e(ci_segments) (accepted [lower,upper] runs) and e(ci_grid)
*       (gamma, D_stat, crit, accepted, B_valid) exported -- the convex
*       hull is now labeled a summary, not the confidence set.
*   #9  explicit -version 15.0- before the Mata block; matastrict off
*       retained deliberately (documented).
*   #10 misleading "EXACT unit resampling" banner replaced with
*       "EXPERIMENTAL ... Alg. 1-oriented, not certified".
*   #11 pre-build advisory when the uncollapsed instrument matrix is
*       projected to exceed ~300 columns under the unlimited maxlag()
*       default (recommend maxlag(1 3) or collapse).
*   version bump 0.9.2 -> 0.9.3.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.4 R20 (audit round 20):
*   #1  e(coefboot)="none" under coefboot(none) (was "onestep",
*       contradicting e(coefficient_bootstrap)).
*   #2  grid-CI per-point validity floor raised to max(10, 90% of the
*       request), matching every other bootstrap object.
*   #3  numerical failure is no longer coded as rejection: e(ci_grid) adds
*       a status column (1 valid / 2 mechanical accept / 3 structural /
*       4 solve failed / 5 too few draws / 6 quantile failed); accepted
*       starts missing; statuses 4-6 are UNRESOLVED -- excluded from the
*       set, counted in e(ci_unresolved), flagged via e(ci_incomplete)=1
*       and a loud warning that the set may be understated.
*   #4  mechanical accepts (D ~ 0) carry status 2, so a missing B_valid
*       there is a legitimate shortcut, not a bootstrap failure.
*   #5  the coefboot failure warning names the replay mode that actually
*       ran (two-step vs one-step).
*   #6  linearity/continuity: valid-draw counts returned and stored
*       (e(boot_linearity_requested/valid), e(boot_continuity_*)); a
*       missing p is now explained ("only X of B draws were valid")
*       instead of printing a bare dot.
*   #8  instrument-proliferation advisory escalates to a strong warning
*       above ~1500 projected columns; help now states the unlimited
*       maxlag() default is risky on long panels.
*   version bump 0.9.3 -> 0.9.4.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.5 R21 (audit round 21, blocker):
*   With ANY unresolved gamma point (status 4-6) the grid-bootstrap
*   inversion is INCOMPLETE and no formal confidence set is reported:
*   e(gamma_lo)/e(gamma_hi)/e(ci_empty)/e(ci_nseg) are missing, boundary
*   diagnostics stay silent, and the acceptance runs over the evaluated
*   points are stored as e(ci_segments_evaluated) (explicitly non-formal)
*   with the full table in e(ci_grid). R20 had flagged unresolved points
*   but still shipped hull bounds / segment counts / "rejected ALL
*   candidates" built as if they were rejections. Minor: a missing
*   linearity/continuity p from an early return (sample statistic not
*   computable) is now explained separately from the too-few-draws case.
*   version bump 0.9.4 -> 0.9.5.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.6 R22 (audit round 22):
*   #2/#3 boottype(unit): the one-step sample statistic, its D<1e-6
*       mechanical-accept shortcut, and the W_first solvability gate now
*       run ONLY for the wild inversion -- under unit they auto-accepted
*       the one-step argmin (appended to the CI grid, D_1step = 0 there)
*       even when the two-stage statistic is nonzero, and could status-4 a
*       point whose fixed-W2 solve is feasible. The one-step argmin is no
*       longer appended to the CI grid under unit (best_gamma is the
*       two-stage zero point). Wild results are unchanged.
*   #5  predict data-integrity: the source columns of every cached series
*       are signed at estimation (_datasignature over keys, depvar, qx,
*       and the base variables of all regressors/instruments; stored in
*       e(p_dsig)/e(p_dsig_vars)); predict refuses (rc 459) when the data
*       changed or a coincident-key dataset is loaded.
*   #6  predict now REQUIRES an explicit statistic (no silent residuals
*       default -- xtabond2-family commands default to xb, so a silent
*       default either way misleads); xb documented as the cached
*       transformed-equation fit, not a current-data linear prediction.
*   #7  predictor: header/notes updated (store/restore supported for the
*       20 most recent runs), version synced to the package, top-level
*       -version 15.0- before its Mata block, legacy-globals comment
*       corrected.
*   version bump 0.9.5 -> 0.9.6.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.7 R23 (proactive fix from external probe):
*   the auto lag L.y was materialized as an untyped (FLOAT) variable --
*   precision loss for large-magnitude depvars fed every FD/FOD
*   difference. Now -gen double-. iv() nested-syntax parsing audited
*   against the same probe: outer varlist/if/in/maxlag/collapse are saved
*   and restored correctly (no fix needed).
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.8 R24 (full code audit):
*   Front end: standard replay + sort preservation; strict missing-value,
*       seed, method and bootstrap validation; all tsrevar temporaries use
*       double precision and preserve c(type); empty marked samples and
*       percentile failures stop cleanly; overlapping invalid external IVs
*       and nonnested automatic continuity comparisons are rejected/omitted.
*   Sample/moments: external IVs can identify early/static equations without
*       blanket history-row deletion; level lag allocation is span-capped;
*       retained panels define the time support; panel construction uses
*       panelsetup(); the auto FD td pre-demeaning reset was removed.
*   Identification/search: each gamma needs at least K+1 moments; flat
*       one- and two-stage profiles no longer select an arbitrary endpoint.
*       Fixed-grid and joint-rank limitations are exposed in e().
*   Inference: AR tests require the full AB variance and use the uncorrected
*       influence sandwich paired with A; Difference-in-Hansen is suppressed
*       when system and reduced FOD cluster universes differ.
*   Bootstrap: contributing clusters receive dense stable draw IDs; unit and
*       coefficient paths use fixed-B draws (no hidden replacement draws);
*       unit covariance honors centering and all-B validity; coefboot(none)
*       truly skips work. Actual/requested metadata are gated consistently.
*   Reporting: effective-stack balance is separate from raw xtset balance;
*       evaluated vs sample-admitted CI-grid counts are distinct; heuristic
*       bootstrap intervals/tests and regular-rank limits are labeled plainly.
*   Performance/data integrity: O(N) panel grouping/balance scans, robust base-
*       variable data signatures, delayed e(cmd), and honest inference metadata.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.9 R25 (audit round 25):
*   vce(uncorrected) renamed vce(robust): the default was always the
*   cluster-robust two-step sandwich, only without the Windmeijer
*   small-sample correction -- "uncorrected" wrongly suggested a
*   model-based/nonrobust VCE. NO alias kept: vce() only ever existed in
*   unreleased builds (v0.7.13+), so the old name is simply gone.
*   e(vce)/e(vce_requested) now report "robust".
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.9 R26 (audit round 26): Windmeijer numerical certification support.
*   Under exportgmm the exact correction inputs (stage-1 ZW/X/Z/uid/r1/
*   Omega1/W1 and stage-2 W2/ZW2/gbar2/n) are exposed as xdpt_w_*
*   externals; _cert_windmeijer.do recomputes the correction from scratch
*   (fresh analytic derivative AND central finite differences of the
*   centered clustered Omega) and compares dOmega/dtheta_j, D, and each
*   component V2 / DV2 / V2D' / DV1rD' plus the total against e(V).
*   xtabond2 parity is infeasible by design (no threshold model there);
*   the FD axis independently certifies sign, orientation, scaling, and
*   the centering term.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.10 R27 (audit round 27):
*   refine(#): opt-in local grid refinement. Each iteration collects the
*   OBSERVED q support values strictly between gamma-hat's two grid
*   neighbors (<= 30 quantile-spaced per iteration), appends them to the
*   estimation grid AND its cache (so the bootstraps, the CI-grid union,
*   and the admission bookkeeping see the refined points natively), and
*   re-runs the two-stage fixed-weight search; stops when no new support
*   values remain or the cap is hit. gamma-hat then effectively lands on
*   observed split points near the optimum. Default 0 (off) keeps results
*   grid()-comparable and every fixed-grid anchor unchanged.
*   e(refine_requested/iterations/added).
*   R28 (same build): the Hansen J and Diff-Hansen display lines are
*   labeled "Diagnostic ... (conditional on gamma-hat)" -- the chi-square
*   reference conditions on a grid-selected, possibly irregular gamma-hat
*   and is not a fully standard specification test. Scalar names unchanged.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.11 R30 (audit round 30):
*   BLOCKER 1: xdpt_wind_applied latched across grid_search calls -- a
*     failed correction on a refine() re-search inherited the coarse
*     pass's applied=1 (e(vce) misreported), and the diff-Hansen reduced
*     FOD re-search could both latch the flag off the REDUCED model and
*     overwrite the xdpt_w_* certification exports with reduced-model
*     inputs. Fix: reset before every refine re-search; hold/zero/restore
*     xdpt_vce_wind, xdpt_expg, xdpt_wind_applied around the reduced
*     search.
*   BLOCKER 2: refine() rejected (198) under kink -- (q-gamma)*1(q>gamma)
*     varies continuously in gamma, so observed-support refinement cannot
*     bracket its optimum. Jump-only until a numerical variant exists.
*   Stage-2 admitted span tracked order-free (grid unsorted after refine
*     appends). refine() brackets now come from the ORIGINAL coarse grid
*     (invalid refined points can no longer stop refinement early) and
*     candidates from the TRANSFORM support even under
*     gridsample(observed).
*   (Reviewer correction acknowledged: <30-unit and IV>units warnings and
*     the projected-width advisory were already present.)
*   version bump 0.9.10 -> 0.9.11.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.12 R31 (audit round 31):
*   refine() brackets now come from the ADMITTED coarse points (ok &
*   fast_ok) plus gamma-hat -- an invalid coarse point is not a criterion
*   evaluation and no longer walls off a possibly better basin behind it.
*   Windmeijer certification exports (xdpt_w_*) are cleared before every
*   refine re-search (a one-step-fallback final pass can no longer leave
*   stale coarse-pass matrices). Deterministic smaller-gamma tie-break in
*   both grid-search stages (post-refine iteration order is arbitrary and
*   a uniform point can tie a support point with an identical design).
*   Continuity diagnostic documented as finite-grid for its restricted
*   kink search even when the jump search is refined.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.13 R32 (audit round 32):
*   refine(): the candidate pool is fixed ONCE from the initial coarse
*   optimum's stage-1-searchable bracket; iterations consume <=30
*   unevaluated pool points each until the pool is exhausted -- the
*   30-point batch cap plus moving brackets previously made the outcome
*   path-dependent (a provisional best in the first batch could orphan
*   never-evaluated support on the far side). Comment corrected: the
*   admitted base is the stage-1 searchable set; excluded points may be
*   stage-2 solvable, and exclusion only widens brackets.
*   e(threshold_search) reflects refine(); e(continuity_kink_search)
*   documents the finite-grid restricted-kink comparison, plus an
*   on-screen note when refine() > 0. Tie-break tolerances declared and
*   precomputed (matastrict-on ready).
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.14 R33 (audit round 33):
*   refine(): endpoint-inclusive 30-point batches (ceil(k*n/30) skipped the
*   smallest candidate -- right-end bias). Bracket base matched to the
*   reported estimator's stage: under two-step, rebuilt from the stage-2
*   admission test with the run's fixed W_n_2 (a stage-1-searchable but
*   stage-2-singular coarse point no longer walls off the basin); one-step
*   fallback keeps the stage-1 set. Pool accounting exported:
*   e(refine_pool/remaining/exhausted/lo/hi), and e(threshold_search) says
*   when the pool was NOT exhausted. Windmeijer correction + xdpt_w_*
*   certification exports now computed ONCE at the final estimate under
*   refine() (disabled during coarse/intermediate passes; one final
*   enabled search). Deterministic smaller-gamma tie-break for the
*   restricted kink selection that seeds the continuity bootstrap DGP.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.15 R34 (audit round 34):
*   refine(): EXPAND-ONLY pool -- when the re-search (with its rebuilt W2)
*   moves gamma-hat outside pooled coverage, the new coarse basin's
*   support is unioned in (old basins never dropped), so refinement
*   follows the estimator instead of stopping at the first basin's edge.
*   Completeness split into two exported statements:
*   e(refine_final_in_basin) and e(refine_neigh_unevaluated), with
*   e(refine_complete) = pool consumed AND final neighbourhood fully
*   evaluated; on-screen warnings when the pool stops unconsumed or the
*   reported gamma-hat's neighbourhood holds unevaluated support. Stage-2
*   admission reconstruction now mirrors the solver's rows>=20 gate. The
*   final enabled re-search runs only under vce(windmeijer) (exportgmm
*   alone never populated the certification inputs).
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.16 R35 (audit round 35):
*   refine(): the coarse-anchor base is REBUILT for the current estimator
*   state after every re-search (helper xdpt2_ref_base_current over the
*   first n_coarse anchors; two-step uses the CURRENT W2 admission,
*   one-step the stage-1 set) -- a base frozen at the initial W2
*   mis-bracketed migrated optima, could reopen the WRONG cell, and
*   produced FALSE completeness against stale anchors. The current basin
*   is now ALWAYS unioned into the pool (the old hull test missed interior
*   never-pooled basins once coverage went non-contiguous).
*   e(refine_hull_lo/hi) renamed to say what they are (convex hull, not
*   refined coverage); e(refine_final_in_initial_basin) renamed for
*   precision; e(threshold_search) is three-state (complete / pool
*   unconsumed / exhausted-but-neighbourhood-unrefined); the unexhausted
*   warning is cap-aware at refine(20).
*   R36 hardening (same build): with fewer than two valid coarse anchors
*   under the final criterion the neighbourhood bracket is degenerate --
*   e(refine_complete)=0, e(refine_neigh_unevaluated)=., and a warning,
*   instead of certifying against a zero-width cell.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.17 R37 (user decision, 11jul2026): default grid(30) -> grid(100),
*   matching the xthenreg convention. Motivated by the R29-F grid-
*   convergence certificate on the package's own demo data: grid(30)
*   selected the wrong basin outright (gamma .0623/obj 231.8 vs
*   gamma .1339/obj 166.6 at grid(100); grid(1000) converges to
*   gamma .1327/obj 158.5), and refine() -- being local -- polishes the
*   coarse argmin's basin, it cannot recover a basin the coarse grid never
*   saw. All default-grid anchors repinned (S.1, V.9, pin2, smoke Z0).
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.18 R38 (full production audit, 11jul2026):
*   Predict cache identity now combines the serial with an exact Stata-
*   generated per-fit token and a 53-bit cache checksum; serial reuse after
*   mata clear fails closed, and a failed predict leaves no output variable.
*   Continuity testing is restricted to the jointly feasible nested grid and
*   rejects materially negative numerical distances instead of clamping them.
*   Static zero-RHS/external-IV models receive the linearity bootstrap.
*   rseed() uses deterministic component-specific seeds, decoupling threshold,
*   linearity, continuity, and coefficient inference objects.
*   refine() covers support beyond an edge anchor and counts distinct q values.
*   AR p-values require at least five pair-contributing panel clusters.
*   Sparse time lookup is memory-bounded with binary fallback; transformed-IV
*   blocks allocate only the effective maxlag(lo hi) interval. Regression
*   tests certify bit-for-bit point estimates/VCEs on FD, FOD+TD, and system.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.19 R39 (second independent production audit, 11jul2026):
*   Predict now recomputes the combined 53-bit checksum from BOTH actual
*   cached row matrices; changing a matrix while leaving its parallel stored
*   checksum intact fails closed (498). Hansen/AR p-values use chi2tail() and
*   normal(-abs()), preserving representable extreme tails. The omitted
*   maxlag() upper bound is truly open (the old 9999 sentinel silently capped
*   long calendars), paired with a pre-allocation gate at 5,000 nominal IV
*   columns / 50 million Z cells to prevent OOM. A maxlag() interval ending
*   below lag 2 in a dynamic model is no longer a hard error: a note is
*   printed and identification must come from exogenous/predetermined
*   moments or external iv(), with the downstream rank/conditioning gates
*   failing closed otherwise (suite J.7 repinned to this contract).
*   Pure L/F operator chains are
*   reduced to their net shift before checking against auto L.y. Help timing,
*   centering, history, coefficient-count, bootstrap, and predict guidance
*   synchronized with the runtime.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.20 R41 (scale-invariance audit, 12jul2026):
*   GMM objective tie, flat-profile, and continuity-nesting tolerances are
*   now purely relative to the compared objectives; the previous unit floor
*   made search and test decisions depend on outcome units below scale one.
*   Symmetric GMM normal, instrument, and moment matrices are admitted after
*   Jacobi equilibration and, when raw conditioning is unit-driven, inverted
*   on that balanced scale then mapped back. Rescaling y (including dynamic
*   L.y regressors and lag-y instruments) or another design column no longer
*   creates a false rank/conditioning failure solely from its measurement
*   units; genuinely ill-conditioned balanced systems still fail closed.
*   Grid-CI status 2 is now reserved for an exact zero distance. Positive
*   statistics, however small in absolute units, receive their bootstrap
*   critical value. Restricted objectives are included explicitly in the
*   unit-CI and linearity unrestricted sets, making those distances
*   nonnegative by construction instead of relying on an absolute clamp.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.21 R42 (bootstrap replay audit, 12jul2026):
*   Coefficient-bootstrap stage 1 and stage 2 now use the estimator's
*   relative objective tie rule, smaller-gamma tie break, two-point search
*   gate, and non-flat-profile gate in every draw. Unit-bootstrap draws use
*   the same identification gates; a restricted-only stage-2 solve can no
*   longer manufacture D_boot=0 and count as valid when the unrestricted
*   grid failed. Bootstrap inference now fails fast for td with
*   boottype(unit), and for td with coefboot(twostep): those paths reused the
*   original-sample FWL projection although resampling changes the required
*   draw-specific projection. Wild threshold/test inference and
*   coefboot(onestep|none) remain available with td; noboot is unaffected.
*   GMM solvers and raw bootstrap passes also reject nonfinite theta,
*   residual, moment, objective, or variance results instead of marking an
*   overflowed numerical solve as valid.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.22 R43 (fail-closed inference audit, 13jul2026):
*   Wild threshold-CI and linearity draws now require at least one finite
*   unrestricted threshold-model objective. A restricted-only solve can no
*   longer manufacture D*=0 after every alternative overflowed or failed.
*   Flat profiles remain admissible in these tests because gamma is a
*   nuisance parameter under their nulls.
*   Cluster-sandwich overflow returns the helper's documented empty failure
*   result; AR variance-component overflow preserves the signed negative
*   pair-count failure contract. Refinement coarse anchors now pass the full
*   fixed-weight solve, not merely the normal-matrix/fast-path precheck.
*   boottype(unit) metadata now states its actual weighting contract: fixed
*   sample W1 with per-draw recentered Omega and W2*.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.23 R44 (extended fail-closed audit, 13jul2026):
*   Symmetric covariance/weight/normal matrices must be positive definite;
*   indefinite matrices can no longer pass a singular-value cond() gate.
*   Cluster-sandwich failure hard-errors instead of attaching Ainv/n under
*   vce(robust). Stage-1 admission metadata comes from the full fixed-W1
*   solve. Refinement certifies completeness from the number of valid coarse
*   anchors, excluding an appended refined optimum from that count.
*   Coefficient intervals fail closed on deviation/bound overflow. Continuity
*   bootstrap minima use a draw-specific jointly finite nested set. noboot
*   bypasses inactive boottype(unit) compatibility gates, while levmaxlag()
*   is rejected outside method(system) instead of being silently ignored.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.24 R45 (CI-grid default update, 16jul2026):
*   Raise the default threshold-CI inversion grid from gridci(25) to
*   gridci(100), matching the seeded MC diagnostic evidence that the coarser
*   grid materially under-covered while the finer grid moved coverage close
*   to nominal. A non-fatal note is printed for gridci()<100 unless nowarn is
*   specified.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.25 R46 (search-architecture correction, 15aug2026):
*   Point estimation uses one incremental nested global cache
*   (g -> 2g-1 -> 4g-3). Stage 1 holds W1 fixed; W2 is constructed exactly
*   once from its selected residuals; stage 2 profiles the retained nested
*   levels under that same W2. For the jump model, refine() is now a final
*   support-point search inside the selected stage-2 basin under the fixed
*   W2: it appends at most 30 new candidates per round and never restarts the
*   estimator or rebuilds W2. Stage-specific gamma/objective, convergence,
*   grid, pool, and W2-build diagnostics are posted. searchmax(199|397)
*   exposes an explicit cap for the default 100-point adaptive schedule;
*   capped nonconvergence remains visible in e(search_incomplete). Omitting
*   searchmax() preserves g -> 2g-1 -> 4g-3 for custom initial grids.
*   e(search_cap_exhausted) identifies a cap reached without stage stability;
*   e(search_hit_max) is its compatibility alias and no longer fires merely
*   because a converged level equals the configured cap. The experimental
*   System GMM API, level-equation builder, moments, diagnostics, returns, and predict
*   branches were removed; the implementation now contains FD and FOD only.
*   The transformed-equation timing classes remain intact: exogenous(),
*   predetermined(), endogenous(), iv(), maxlag(), and collapse.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.26 R47 (public-scope freeze, 16aug2026):
*   The reported numerical estimator is now defined on the finite grid chosen
*   explicitly by grid(), followed by optional final fixed-W2 refine() for the
*   jump model. Adaptive escalation and its user controls searchmode(),
*   searchmax(), and searchtol() are retired with informative compatibility
*   errors; legacy e(search_*) fields remain temporarily for MC-harness
*   compatibility and are non-applicable on the fixed path. The coefficient
*   bootstrap is no longer automatic: coefboot(none) is the default and the
*   onestep/twostep implementations are retained only as opt-in legacy paths.
*   boottype(unit) remains callable solely as a verification-only,
*   Algorithm-1-oriented threshold-CI path; it is not certified as exact and
*   e(boottype_status) records that status. The default wild threshold
*   bootstrap, the FD/FOD estimators, instrument controls, diagnostics, and
*   fixed-W2 refinement are otherwise unchanged. The boundary warning now
*   describes Gong-Seo's p10-p90 application grid as trim(.20) with quantile
*   spacing rather than incorrectly calling trim(.15) their convention.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.27 (25sep2026): method(fod) now dates internal lag instruments from t+1,
*   the xtabond2 convention (Roodman 2009), so maxlag() admits the FOD-valid
*   lags: L.y and endogenous() from t-1, predetermined() from t. Earlier FOD
*   reused the FD dates (valid, but one usable lag short per variable, and an
*   endogenous() variable could drop the first FOD equation). method(fd) is
*   bit-for-bit unchanged.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.28 (25sep2026): exact duplicate instrument columns are dropped (first
*   kept) after stacking. Declaring a variable together with its own lag
*   (x and L.x, or L2.depvar next to the automatic L.depvar) produced
*   identical (variable, date) columns, made Z'Z singular, and rejected every
*   candidate threshold; such models are now estimable. A lag of the
*   dependent variable in indepvars or exogenous() is rejected (error 198),
*   because it is not strictly exogenous; it belongs in predetermined(). The
*   conditioning of Z, which does not depend on gamma, is checked once before
*   the grid search, with a specific message. Fits that succeeded under
*   0.9.27 had no duplicate column and are bit-for-bit unchanged.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.29 (25sep2026): the AR(1)/AR(2) statistics now use the reported
*   variance of theta-hat in their third term, so under vce(windmeijer) they
*   use the Windmeijer-corrected variance, as xtabond2 does with twostep
*   robust; under the default vce(robust) nothing changes. iv() variables no
*   longer define the equation sample: a missing value contributes a zero
*   instrument in that row, as for the internal lag instruments (previously it
*   removed the equation, and under FD also the next period's equation). Fits
*   with vce(robust) and without missing iv() values are bit-for-bit
*   unchanged.
*   Second round (code read-through): regressors or q with no variation within
*   units, or common to all units under td, are rejected (498) instead of being
*   fitted on rounding noise; any operator on depvar in indepvars/exogenous(),
*   and a lag of an endogenous/predetermined variable placed in indepvars, are
*   rejected (198); predict checks a key-tied data signature (459 on swapped
*   rows); e(wind_same_threshold), e(W1_fallback), fallback warnings, and a
*   corrected rank-failure message; display and help fixes.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.30 (25sep2026): instrument columns that are exact multiples of the
*   constant instrument columns are dropped (xdpt2_drop_cellconst). A
*   regressor or iv() variable that takes one value for every unit in a period
*   (a macro variable, a trend) produced such columns in every period block,
*   made Z'Z singular, and stopped the command at the rank check; it is now
*   instrumented by the per-period constants, as xtabond2 and xthenreg do
*   implicitly through a generalized inverse, and a note names it
*   (e(N_iv_common), e(iv_common)). Under td such a regressor is still not
*   identified (498). The rank-failure message names variables that are common
*   only up to rounding. Fits that succeeded under 0.9.29 had no such column
*   and are bit-for-bit unchanged.
*   Speed (bit-for-bit): xdpt2_transform_unit, called once per unit per grid
*   point, no longer declares -external- variables. Mata binds externals at
*   every call at a cost that grows with the number of live objects, which
*   took about 90% of a point estimate and made it O(N^2); the settings are
*   now read once per stack (xdpt2_unit_cfg). Point estimates with the Monte
*   Carlo settings: N = 400 51 s -> 8 s, N = 800 ~3 min -> 14 s, N = 1600
*   ~15-20 min -> 34 s (same PC). The bootstraps and tests compute Z'Y once
*   per shared Z instead of once per grid point (xdpt2_fast_obj_batch_list),
*   and the coefficient bootstrap no longer binds an external per draw.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.31 (25sep2026): the bootstrap one-step objectives (threshold
*   confidence set, linearity and continuity tests, coefficient bootstrap) are
*   evaluated by xdpt2_fast_obj_split_list. Each bootstrap sample is F + E,
*   the fit that generates every draw plus the reweighted residuals; the
*   moment vector is linear, so g(F + E) = g(F) + g(E). g(F) is formed from
*   residuals as before (once per candidate threshold); g(E) from cross
*   products, Z'E/n - (Z'dW/n)(C Z'E/n), with Z'E/n formed once per shared Z;
*   no n-row residual matrix is formed per draw and candidate threshold. The
*   sample statistics and the scalar fallback paths keep the 0.9.30
*   arithmetic; point estimates, SEs, Hansen J and AR tests are unchanged, and
*   bootstrap critical values agree with 0.9.30 to rounding (against a 40-digit
*   reference the split is as accurate as the residual form).
*   Audit fixes: a regressor of the form a_i + g_t (firm age) under td is
*   rejected (498) -- FOD used to return coefficients near 1e25; a CI point with
*   a singular one-step normal matrix is status 3 (not admissible) instead of
*   4, which withheld the whole set; the stage-1 replay of the coefficient
*   bootstrap and of boottype(unit) searches the initial grid only, as the
*   reported stage 1 does; operator temporaries are recomputed in double
*   (tsrevar keeps a float source's type); the common-column rule compares
*   across periods to a relative 1e-10 (a trend's differences), reuses its
*   mask for an unchanged Z, and a specific error replaces the generic failure
*   when it leaves too few instruments; real-valued search diagnostics
*   (e(gamma_stage1), e(obj_stage1), e(gamma_stage2_global), ...) are carried
*   in scalars instead of locals, which rounded them; t_max is the last
*   equation time, so trailing history rows no longer inflate the allocation
*   gate; iv() keeps Stata's own error codes; the dropped-column note is worded
*   exactly; dead code (xdpt2_fast_obj_batch) removed. Non-fatal "Warning:"
*   and "Note:" labels are coloured with {err} inside -as text- output: shown
*   -as err-, they escaped -quietly-, so -qui xtdpthresh- printed fragments.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.32 (25sep2026): two changes from the full audit of 0.9.31.
*   (1) The critical value of the grid-bootstrap test is the empirical
*   (1 - alpha) quantile of the valid bootstrap statistics, the
*   ceil((1 - alpha)*B_v)-th order statistic (Gong and Seo 2026, eq. 7),
*   instead of the type-7 interpolated quantile, which fell up to one rank
*   lower and made the test liberal by about 0.9/(B_v + 1). It now agrees
*   with the add-one rule of the linearity and continuity tests.
*   (2) Kink model: e(V) is the slope block of the joint variance of the
*   slopes and gamma-hat. The kink term is continuous in gamma, with
*   derivative -delta_k*1(q > gamma), so the design is augmented by
*   x_g = -delta_k*T(1(q > gamma-hat)) and the existing sandwich and
*   Windmeijer code applies. The conditional variance is kept in e(V_cond)
*   and used by the AR tests, which stay conditional on gamma-hat;
*   e(kink_joint_vce) = 1 (0 if the joint variance could not be computed,
*   missing for the jump model). The jump model is unchanged.
*   Text only: the disconnected-set note says "rejected or inadmissible";
*   the trim warning no longer attributes trim(0.40) to Seo-Shin or
*   trim(0.10) to xthreg2; e(threshold_bootstrap) no longer says
*   "xthenreg-style"; two code comments corrected.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.33 (26sep2026): speed only; no change in method.
*   (1) The batched bootstraps (threshold confidence set, linearity and
*   continuity tests) leave out a candidate threshold without the one-step
*   solve (fast_ok = 0), as the scalar loops always did. Before, one such
*   point sent every draw to the scalar loop: with qx(debt) and debt >= 0 also
*   a regressor, gamma = 0 makes debt*1(debt > 0) equal debt, and the
*   default confidence set on the Hansen investment data would have taken
*   several hours instead of seconds.
*   (2) xdpt2_fast_obj_split_list receives the reweighted residuals by their
*   factors (residuals, cluster index, cluster weights) and forms
*   Z'E/n = S'ETA/n from the within-cluster sums S of z_i*r_i, without the
*   n x B matrix E; g(F) is formed for all candidate thresholds that share Z
*   by one product Z'[r_1, ..., r_m] instead of one matrix-vector product
*   each. The arithmetic of each term is unchanged; only the order of
*   summation differs, so bootstrap objectives equal those of 0.9.32 to
*   rounding. Point estimates, SEs, J and AR are bit-for-bit unchanged.
*   (3) The per-gamma caches no longer restack the data at every candidate
*   threshold. Only the regime columns of dW depend on gamma; the rows, dY,
*   Z, the base columns of dW and the time-effect partialling do not. The
*   first stack of each cache build records them (xdpt2_tpl_build), and dW
*   is then rebuilt with the arithmetic of xdpt2_transform_unit
*   (xdpt2_tpl_dW): W(t) - W(t-1) under FD; under FOD, the mean of the later
*   equation rows as a quad-precision sum divided by their number, which is
*   how mean() forms it. The second gamma is also stacked in full and
*   compared bit for bit; on any difference the build restacks at every
*   gamma. Results are bit-for-bit those of (1)-(2).
*   Example 4 of the help at the defaults (gridci(100), boot(299)): about
*   7 hours under 0.9.32 (profiled), 23 seconds now.
*   (4) predict: the data signature (e(p_dsig_type) "rowsig2") weights
*   each column by Park-Miller sequences in the row rank. The weights of
*   rowsig1, 1 + mod(r*a, m)/m, are linear in r below about 30,000 rows,
*   so rowsig1 checked only sum(x) and sum(r*x): a change of +c, -2c, +c
*   on three consecutive rows kept it unchanged, and predict then served
*   results of the old data. predict still checks rowsig1 results.
*   (5) Threshold confidence set: a candidate threshold without the one-step
*   solve (fast_ok = 0) is status 3 (not admissible) only if the reported
*   estimator cannot use it either. The two-step search runs over every
*   admissible point with the second-step weight W2, so a point whose
*   one-step normal matrix fails the numerical gate while Z_W'W2 Z_W passes
*   it can be the reported threshold; the one-step inversion cannot
*   evaluate it, so it is now unresolved (status 4) and the set is not
*   reported. Exactly collinear points (q*1(q > gamma) = q at gamma = 0 when
*   q >= 0) fail under both weights and stay status 3. The inversion and the
*   linearity and continuity tests use a candidate or report a p-value only
*   if every bootstrap draw is valid (before: at least 90% and 10), so no
*   result rests on draws selected by numerical failures; a wild draw can
*   fail only on a nonfinite value, so fits with all draws valid are
*   unchanged.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.34 (26sep2026): fixes from the independent audit of 0.9.33, and speed.
*   Bugs.
*   (1) A model whose only regressor is L.depvar stopped with a Mata error
*   (r(3598)): the list of variables checked for constancy within units was
*   blank but not empty. Such models now run.
*   (2) iv(): a lag of depvar or of an endogenous() variable must be of order
*   2 or more under FD and 1 or more under FOD, and a lag of a
*   predetermined() variable of order 1 or more under FD (error 198
*   otherwise). iv(L.y) under FD was accepted and biased every estimate,
*   and the Hansen test did not detect it.
*   (3) predict: the data signature sorts the rows by every column. Sorted by
*   (panel, time) only, rows that share a key (two or more rows of a panel
*   with a missing time) had no fixed order, and predict could refuse
*   unchanged data (459).
*   (4) boottype(unit): a candidate threshold whose normal matrix fails the
*   numerical gate under W2 (for example q*1(q > gamma) = q) is not
*   admissible (status 3), as in the wild inversion since 0.9.33. It was
*   unresolved (status 4), so the unit bootstrap reported no set on the
*   main model of the help.
*   (5) predict refuses results whose data signature is not rowsig2 (fits of
*   0.9.32 or earlier still in memory after an update) with error 498.
*   (6) The template of the per-gamma caches forms the FOD later-row sums in
*   forward order, as mean() does. 0.9.33 summed them in reverse order,
*   which differs by one unit in the last place for values whose exponents
*   span 2^35 or more with exact cancellation, so the rebuild equalled a
*   full restack only for data that passed the check at the second gamma;
*   it now does for all data.
*   (7) The template records are released at the start and end of each run.
*   Method.
*   (8) The threshold confidence set inverts the test on the criterion of the
*   reported estimator: for two-step estimates, the second-step criterion
*   with W2 held at its sample value (e(ci_criterion) "twostep"); after a
*   one-step fallback, the one-step criterion as before ("onestep"). The
*   statistic is then zero at gamma-hat, which always belongs to its set, as
*   in Gong and Seo (2026). Up to 0.9.33 the zero was at the one-step
*   argmin, and gamma-hat could be rejected (help example 4). All candidate
*   thresholds use the same Mammen weights (common random numbers), which
*   also removes most of the fragmentation of the set on flat profiles.
*   (9) The critical value is the ceil(p*(B + 1))-th order statistic, the
*   rule of the add-one p-values of the tests; with B < p/(1 - p) every
*   candidate is accepted. ceil(p*B) (0.9.32) is liberal when
*   (1 - p)*(B + 1) is not an integer (6% instead of 5% at B = 100).
*   Correction to the 0.9.32 entry: the two rules agree only when
*   (1 - p)*(B + 1) is an integer.
*   (10) Kink model: the criterion is continuous in gamma, and the joint
*   variance of the slopes and gamma-hat presumes its exact minimizer. After
*   the grid search, the second-step criterion (the one-step criterion after
*   a fallback) is minimized between the admitted grid points next to the
*   grid minimum: 40 equally spaced points, then rounds of 10 around the
*   best point so far, down to a spacing below 1e-6 of that interval. The
*   best point is added to the estimation grid; e(kink_refined) = 1 if it
*   lies between grid points. In a well-identified design (N = 3000) the
*   continuous minimizer lay 0.6-1.8 standard errors of gamma-hat from the
*   grid point, and the slope on q moved by up to 0.7 of its SE.
*   (11) The confidence set and the linearity and continuity tests take the
*   unrestricted minimum over the initial grid, on the sample side and in
*   every draw. Points added around the sample's own minimum (refine(), the
*   kink refinement) lowered only the sample minimum, since no draw gets the
*   same search. gamma-hat keeps a zero statistic.
*   (12) Instrument columns that are linear combinations of the others (a
*   period with fewer units than block columns; FOD with td when units share
*   their last period) are dropped, as the generalized inverse of xtabond2
*   does (e(N_iv_dep), shown in the output), instead of stopping the command
*   with 498. The rank check that follows concerns ill-conditioning only.
*   (13) The default regime floor is the number of observations in the
*   smaller trimmed tail, min(#{q <= q_lo}, #{q > q_hi}), and the last point
*   of a uniform grid is q_hi exactly. The floor was ceil(trim*n/2), and the
*   top grid point was admitted or not by rounding when q_hi is a data
*   value.
*   (14) e(N_switch): the number of units whose regime changes within their
*   equations at gamma-hat; a note is shown when it is below 10.
*   (15) e(wind_same_threshold) compares, in the jump model, the regime
*   splits of the stage-1 threshold and gamma-hat. Under kink with the joint
*   variance, whose correction differentiates W2 in the stage-1 threshold as
*   well, it is missing.
*   (16) The linearity test forms no n x B matrix (the split with a
*   one-entry cache for the linear model).
*   Messages: the note on invalid draws, the one-step fallback warning (it
*   named only a singular W2), the note on the approximate threshold set (it
*   referred to "intervals below"), the error for a flat criterion (it was
*   reported as "no gamma admitted"), and the grid count after refine().
*   Correction to the 0.9.33 entry, item (5): a continuity-test draw is also
*   discarded when the kink criterion falls below the jump criterion beyond
*   rounding, so a nonfinite value is not the only way a draw fails.
*   Speed, bit for bit: xdpt2_build_units allocates the unit vector once
*   (appending copied it at every unit, O(N^2): 16 of 22 seconds of a fit
*   at N = 3000); later cache builds of the same model reuse the template
*   and the first-step weight of the main build and compare the template
*   pieces with their reference once; the FOD later-row sums of all units
*   with the same number of equations are formed by one quadcross. Study A
*   settings at N = 1600: 7.8 to 2.6 seconds per fit.
* ---------------------------------------------------------------------------

* ---------------------------------------------------------------------------
* v0.9.37 (28sep2026): continuity-test hardening (code review of 0.9.36).
*   (a) The jump entries of the continuity comparison (its minimum, the
*   batched and scalar bootstrap alternatives, and the source of the jump
*   residuals) must lie on the kink rows -- same uid, times and dY -- not
*   merely have the same row count. The rows do not depend on gamma by
*   construction, so results are unchanged; the guard makes it explicit.
*   (b) No fallback to the first-step weight when the fixed-W2 solve leaves
*   fewer than two jointly feasible points: the test is then not reported.
*   Gong-Seo's statistic (sec. 3.2) uses the efficient weight W_n and the
*   limit of their Theorem 4 is built on Omega^{-1}; a W1 distance is a
*   different statistic their theory does not cover (unchanged from 0.9.36).
*   (c) The draws returned for citest() are reset at each point, so they can
*   only belong to the point just evaluated.
*   (d) The citest() line reads "rejected/not rejected at the 5% level"
*   (was "reject at 95%").
*   Estimates, confidence sets and p-values equal 0.9.36.
* ---------------------------------------------------------------------------
* v0.9.36 (28sep2026): continuity-test power; threshold-test diagnostic.
*   (a) Continuity-test power. The bootstrap DGP used the restricted (kink)
*   residuals; under a jump they carry the omitted discontinuity, which
*   inflated every bootstrap statistic and the critical value. It now uses
*   the kink fit plus the unrestricted (jump) residuals, as Gong-Seo Alg. 1
*   and boottype(unit). With a two-step fit both models are compared under
*   the fixed second-step weight W2 of the jump fit (the criterion of the
*   reported estimator, as the threshold CI since 0.9.34) instead of the
*   one-step weight; after a one-step fallback W1 is kept.
*   (b) New diagnostic citest(#): the grid-bootstrap test of H0: gamma = #
*   alone, with the statistic, bootstrap and weight of the confidence set,
*   run last under its own component seed (other results unchanged). Returns
*   e(citest_D), e(citest_crit), e(citest_accept), e(citest_p) (add-one
*   bootstrap p-value; p > alpha iff accepted), e(citest_status) (status
*   code of e(ci_grid)) and e(citest_draws). At the true threshold of a
*   simulation its rejection rate is the size of the test that the set
*   inverts; its acceptance rate is the coverage proved by Gong-Seo.
*   (c) Coverage convention. Gong and Seo (2026, eq. 7 and Theorem 5) define
*   the set as {gamma in the grid : accepted} and prove
*   P(gamma0 in set) -> 1 - tau, i.e. acceptance of the test AT gamma0; their
*   Monte Carlo (Table 1) scores exactly that, and the set may be convexified
*   (the hull e(gamma_lo), e(gamma_hi)). Scoring coverage by the union of
*   e(ci_segments) is not a coverage of the inference: a gamma0 between an
*   accepted and a rejected grid point was never tested and counts as a miss.
*   An independent prototype of the 0.9.35 wild procedure on the Gong-Seo
*   benchmark (FD, T=6, 24 lag instruments, 46-point grid, R=300-400)
*   accepts gamma0 in 94-96% of samples (kappa = 0, 1; N = 400, 800), with
*   hull coverage 97-100% and union-of-segments coverage 85-94%. The wild
*   threshold bootstrap is therefore kept and citest(#) scores the proved
*   property directly. A boundary refinement of the set (extra test points
*   where acceptance changes) was prototyped and dropped: union coverage
*   91.7 -> 93.0% at three times the CI cost. The same prototype gives the
*   continuity change of (a): size 1-2% (0.9.35: 3%) and power
*   45 -> 66% at kappa = 2 and 74 -> 94% at kappa = 3 (N = 400); at kappa = 1
*   no variant exceeds 15%, since a kink at gamma - kappa/delta3 reproduces
*   the jump regime above gamma.
* ---------------------------------------------------------------------------
* v0.9.35 (27sep2026): joint variance of the slopes and gamma-hat in the jump
* model.
*   Under the GMM asymptotics of Seo and Shin (2016), gamma-hat and the slopes
*   of the jump model are jointly sqrt(n)-normal, so a variance that treats
*   gamma-hat as known (0.9.34 and earlier, the Hansen convention) omits a
*   first-order term; vce(windmeijer) did not add it either. The sample
*   moments are step functions of gamma, but their expectation is smooth, with
*   derivative -E[z (1, x')delta f(gamma | .)]. As in xthenreg (Seo, Kim, and
*   Kim 2019), it is estimated with a Gaussian kernel, with a bandwidth of the
*   form of xthenreg's, h = 1.06 s_q n^(-1/5) (s_q: standard deviation of q
*   over the panel rows, n: units), times the new option bwscale(), whose
*   default 1.5 is xthenreg's multiplier (see item k). The derivative column
*   x_g = -T(phi_h(q - gamma-hat) (1, x')delta-hat) is transformed by a
*   template recorded at gamma-hat (FD or FOD, rows, td partialling), which
*   must reproduce the regime columns of the design there, and enters the
*   sandwich and Windmeijer code of the kink model's joint variance (0.9.32).
*   e(V) is the slope block; e(V_cond) keeps the conditional variance, as
*   under kink; e(joint_vce) (both models; e(kink_joint_vce) is kept for
*   kink), e(gamma_bw), e(bwscale). e(wind_same_threshold) is missing when the
*   joint variance is reported, since the correction then differentiates W2
*   in the stage-1 threshold as well. Point estimates, gamma-hat, the
*   confidence set, the tests and J are unchanged; the AR statistics change by
*   item (e) below.
*   Checks: against xthenreg's jump covariance (its own estimator replayed on
*   a common grid, 6 configurations incl. static, discrete q and gamma-hat
*   differing between the steps) the bandwidth is identical and the efficient
*   form built from xtdpthresh's design and derivative column equals
*   xthenreg's full covariance, gamma included, to 3e-14 - 1.3e-12; under FOD,
*   td and gaps the template column equals the stacker's own transformation of
*   the level kernel column to 4e-16.
*   Also in 0.9.35 (external review):
*   (a) Dropped instrument columns. The rule of 0.9.34 (Cholesky of the
*   unit-diagonal Gram matrix, pivot below 1e-13, i.e. relative residual below
*   3.2e-7) is unchanged, but rounding in the Gram matrix is of order 1e-8 in
*   relative residual, so it drops nearly dependent columns as well and cannot
*   tell them from exact combinations. The final stack at gamma-hat now
*   measures the relative residual of each dropped column on Z itself against
*   the kept columns (Householder QR): e(iv_dep_res) is the largest, and
*   e(N_iv_dep_near) counts those above 1e-10. The output calls the columns
*   linear combinations only when none is above 1e-10; otherwise it warns that
*   the directions they add are lost and that the estimates depend on how the
*   instruments are written (iv(z1 z2) with z2 = z1 + 1e-8*q^3 loses the q^3
*   direction that iv(z1 zc), zc = (z2 - z1)/1e-8, keeps). Measured residuals:
*   exact combinations (a thin FD period, FOD + td with a shared last period,
*   an iv() equal to L2.y + L3.y) 0 to 9e-13; nearly dependent columns (that
*   example, a variable common to all units up to 1e-9) 7e-10 to 3e-7.
*   Estimates are unchanged.
*   (b) bwscale() with a missing value is rejected (198); it was accepted and
*   treated as 1 while e(bwscale) stored the missing value.
*   (c) The note on e(N_switch) says that the change in the intercept rests on
*   the switching units; the slope changes also use units that stay in the
*   upper regime.
*   (d) Help: the chi-square reference of Hansen J presumes regular
*   identification, and a large p-value does not show that the instruments
*   are valid; the AR(2) and Hansen tests are diagnostics for the assumptions.
*   (e) The AR statistics include the estimation of gamma-hat when the joint
*   variance is reported (e(ar_joint) = 1). The derivative column of the
*   residuals in gamma, as in the joint variance (kernel under jump, exact
*   under kink; under FOD also on the FD rows of the test, from a template of
*   that stack), is appended to the test and estimation designs of the
*   Arellano-Bond (1991, eq. 8) statistic, and V is the joint variance of
*   (theta, gamma), of the reported type. The statistics with gamma-hat treated
*   as known, which equal xtabond2's with gamma fixed at gamma-hat, are kept
*   in e(ar1_cond), e(ar2_cond). On 10 fits (FD, jump and kink) AR(2) moved by
*   at most 0.024.
*   (f) e(ar_vcetype) follows e(ar_joint). (A Hansen J p-value with gamma not
*   counted in df was added during development and withdrawn: the argument
*   J(gamma-hat) <= J(gamma_0) needs the true threshold on the grid and a
*   consistent second-step weight, which is built at the first-step estimate;
*   Stock and Wright's result concerns the continuously updated criterion.)
*   (g) xdpt2_ar_full runs in O(n) on sorted stacks (every stack is sorted by
*   unit and time): the rows of a unit are one block (panelsetup) and the
*   lag-k partner of a row is among the k rows before it; c_i is looked up
*   once per unit. It scanned all rows for every unit (O(n x units)). Same
*   rows, order and arithmetic; the scans remain for any other order. The
*   four AR statistics of a fit now cost less than the two of 0.9.34.
*   (h) A discrete threshold variable. The joint variance of the jump model and
*   the AR statistics with gamma-hat assume that q has a continuous density,
*   positive at the threshold (Seo and Shin 2016, Assumption 2); with few
*   values of q, every threshold between two adjacent values gives the same
*   split and the kernel derivative has no density to estimate, yet the
*   variance is computed. e(q_nvals_bw) counts the distinct values of q within
*   two bandwidths of gamma-hat, and the output warns below 10; the help
*   states the limitation. The bwscale() help now says that the AR statistics
*   with gamma-hat depend on the bandwidth as well.
*   (i) Help: the unbalanced-panel assumptions are stated for the unit's
*   observation pattern and the equations actually retained (FD needs t-1,
*   FOD a later observation, the instruments depend on the periods observed,
*   and rows without a non-constant instrument are dropped); (U2) is the
*   identification and rank condition on the retained equations, which
*   selection independent of the regime does not replace. The title no longer
*   says "continuity-robust inference" (the default bootstrap is an
*   approximation that Gong and Seo's results do not cover), as in the package
*   description.
*   (j) Help: vce(windmeijer) with the joint variance is described as a
*   Windmeijer-type correction applied to the linearized joint moments (the
*   2005 result is derived for linear moment conditions); the standard-error
*   note limits the unreliability near continuity to the unrestricted jump
*   model (under kink, with the restriction true, a nonzero change in the
*   slope of q and full rank, the estimator is asymptotically normal); and the
*   bandwidth is said to have the form of xthenreg's rule.
*   (k) The default of bwscale() is 1.5, the default multiplier of xthenreg
*   (it was 1 during development), so the default bandwidth, and the jump
*   model's default joint SEs, can be compared with xthenreg directly. The
*   help cites Windmeijer (2005): the correction improves the variance
*   estimate for moment conditions linear in the parameters (sec. 2.1); for
*   nonlinear moments its term is of the same order as other remainders
*   (sec. 2), and the joint moments here are not linear in gamma.
*   (l) Help: collapse -- fewer instruments can reduce the problems caused by
*   too many instruments (it said that they make the Hansen test more
*   reliable).
*   (m) External review (faults), two computational errors:
*   R1. xdpt2_syminv checked the conditioning of the equilibrated matrix but,
*   when the raw matrix also passed, returned invsym() of the raw matrix,
*   which drops pivots below an absolute tolerance: a matrix of tiny scale
*   (condition number 3) came back as a zero "inverse" with ok = 1. With y
*   multiplied by 1e6 the command selected another threshold (-1.38 -> 0.32
*   in the review's fixture). The inverse is now always that of the
*   equilibrated matrix, mapped back, and checked (no dropped pivot; residual
*   of C * inv(C) within 10 k eps cond(C)). Estimates change at the rounding
*   level; they no longer depend on the units of the data.
*   R2. The Windmeijer correction V2 + D V2 + V2 D' + D V1 D' assumes equal
*   stage-1 and stage-2 Jacobians; they differ when gamma-hat_1 != gamma-hat_2
*   and, in the joint variance, through the derivative columns, and the
*   formula could give negative variances (a variance of -231 in the review's
*   FOD fixture; an indefinite matrix with positive diagonal in a correctly
*   specified design). It is now the variance of the linearized two-step
*   estimator, (L2 + D L1) Omega1 (L2 + D L1)' / n, which equals the former
*   when the Jacobians are equal and is positive semidefinite. A corrected
*   variance is used only if it is positive semidefinite (xdpt2_psd_ok);
*   otherwise the cluster-robust variances are reported (e(vce_applied) = 0),
*   for e(V), e(V_cond), and the AR statistics alike.
*   (n) External review (R4): the stacking template refused K = 0, a static
*   model without base regressors whose only regime term is the intercept
*   (the hinge under kink). The jump model's joint variance and the AR
*   statistics with gamma-hat rely on the template, so such models reported
*   the conditional variance (e(joint_vce) = 0, correctly labeled) and, under
*   kink, conditional AR statistics, although the joint Jacobian has full
*   rank. The template now accepts K = 0 (two slicings guarded). Its cache
*   path is checked against a full stack at the second grid point, as for any
*   K, so the estimates are unchanged.
*   (o) External review (R3): colsum() skips missing values. In the batched
*   bootstrap objectives (xdpt2_fast_obj_split_list, and the linearity test
*   without regressors) a term that overflowed, or a missing weight, left a
*   partial sum -- zero when the other terms were zero -- and the draw passed
*   as valid, whereas the scalar path reports it as failed. Such draws now have
*   a missing objective and go through the invalid-draw handling; the
*   per-unit moment sums (xdpt2_gsum_by_unit) return missing for a missing
*   input. No natural fit reproduced a wrong confidence set.
*   (p) Help and output shortened, with every point of the review kept.
* ---------------------------------------------------------------------------
