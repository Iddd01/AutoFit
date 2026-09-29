* preflight.do -- run_tonight.ps1 checks that the Seo-Shin command used by the
* XTH cells is installed before anything starts.
capture which xthenreg
if _rc {
    di as err "PREFLIGHT_FAIL xthenreg missing: ssc install xthenreg"
    exit 111
}
capture mata: mm_quantile((0\1), 1, .5)
if _rc {
    di as err "PREFLIGHT_FAIL moremata missing: ssc install moremata"
    exit 111
}
which xthenreg
di "PREFLIGHT_PASS"
