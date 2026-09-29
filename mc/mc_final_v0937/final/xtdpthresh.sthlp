{smcl}
{* *! version 0.9.37  28sep2026}{...}
{vieweralsosee "[XT] xtset" "help xtset"}{...}
{vieweralsosee "[XT] xtabond" "help xtabond"}{...}
{vieweralsosee "xthenreg (if installed)" "help xthenreg"}{...}
{vieweralsosee "xtabond2 (if installed)" "help xtabond2"}{...}
{viewerjumpto "Syntax" "xtdpthresh##syntax"}{...}
{viewerjumpto "Description" "xtdpthresh##description"}{...}
{viewerjumpto "Options" "xtdpthresh##options"}{...}
{viewerjumpto "Remarks" "xtdpthresh##remarks"}{...}
{viewerjumpto "  Moment conditions and assumptions" "xtdpthresh##model"}{...}
{viewerjumpto "  Bootstrap algorithm" "xtdpthresh##bootstrap"}{...}
{viewerjumpto "  Tests and diagnostics" "xtdpthresh##tests"}{...}
{viewerjumpto "  Notes and limitations" "xtdpthresh##notes"}{...}
{viewerjumpto "  Differences from xthenreg" "xtdpthresh##xthenreg"}{...}
{viewerjumpto "Examples" "xtdpthresh##examples"}{...}
{viewerjumpto "Postestimation" "xtdpthresh##postest"}{...}
{viewerjumpto "Stored results" "xtdpthresh##results"}{...}
{viewerjumpto "References" "xtdpthresh##references"}{...}
{viewerjumpto "Authors" "xtdpthresh##authors"}{...}
{viewerjumpto "Version" "xtdpthresh##versionhistory"}{...}
{cmd:help xtdpthresh}
{hline}

{title:Title}

{p2colset 5 19 21 2}{...}
{p2col:{cmd:xtdpthresh} {hline 2}}Dynamic panel threshold regression for
unbalanced panels, with endogenous regressors and grid-bootstrap confidence
sets for the threshold{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:xtdpthresh} {depvar} [{indepvars}] {ifin}{cmd:,}
{opt qx(varname)} [{it:options}]

{synoptset 32 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{p2coldent:* {opt qx(varname)}}threshold variable{p_end}
{synopt:{opt endo:genous(varlist)}}endogenous regressors; instrumented by levels from t-2 (t-1 under FOD){p_end}
{synopt:{opt pred:etermined(varlist)}}predetermined regressors; instrumented by levels from t-1 (t under FOD){p_end}
{synopt:{opt exo:genous(varlist)}}additional strictly exogenous regressors{p_end}
{synopt:{opt kink}}impose continuity of the regression function at the threshold{p_end}
{synopt:{opt td}}partial common time effects out of the estimating equations{p_end}
{synopt:{opt hist:ory(panel|sample)}}whether {cmd:if} and {cmd:in} also bound the lag history; default is {cmd:panel}{p_end}

{syntab:Instruments}
{synopt:{opt maxlag(# [#])}}lag range of the internal instruments; default is all available lags{p_end}
{synopt:{opt collapse}}collapse the instruments across periods{p_end}
{synopt:{opt iv(varlist [, collapse])}}additional external instruments{p_end}

{syntab:Transformation and threshold search}
{synopt:{opt method(fd|fod)}}first differences or forward orthogonal deviations; default is {cmd:fd}{p_end}
{synopt:{opt grid(#)}}number of grid points for the threshold search; default is {cmd:grid(100)}{p_end}
{synopt:{opt gridt:ype(uniform|quantile)}}placement of the grid points; default is {cmd:uniform}{p_end}
{synopt:{opt trim(#)}}total trimming rate of the threshold variable; default is {cmd:trim(0.10)}{p_end}
{synopt:{opt grids:ample(effective|observed)}}observations that define the trimming bounds; default is {cmd:effective}{p_end}
{synopt:{opt minreg:ime(#)}}raise the minimum number of observations in each regime{p_end}
{synopt:{opt ref:ine(#)}}rounds of local refinement of the threshold estimate, 0-20; default is {cmd:refine(0)}{p_end}

{syntab:Inference}
{synopt:{opt gridci(#)}}number of candidate thresholds in the confidence-set grid; default is {cmd:gridci(100)}{p_end}
{synopt:{opt boot(#)}}number of bootstrap replications; default is {cmd:boot(299)}{p_end}
{synopt:{opt rseed(#)}}seed for reproducible bootstrap draws{p_end}
{synopt:{opt noboot}}skip all bootstrap inference{p_end}
{synopt:{opt notest}}skip the bootstrap linearity test{p_end}
{synopt:{opt citest(#)}}run the pointwise bootstrap test H0: gamma = #{p_end}
{synopt:{opt vce(robust|windmeijer)}}variance estimator for the slope coefficients; default is {cmd:robust}{p_end}
{synopt:{opt bw:scale(#)}}scale of the kernel bandwidth of the jump model's joint variance; default is {cmd:bwscale(1.5)}{p_end}
{synopt:{opt nocenter}}use the uncentered moment covariance{p_end}
{synopt:{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}

{syntab:Reporting}
{synopt:{opt nowarn}}suppress nonfatal warnings and most notes{p_end}
{synopt:{opt verbose}}display intermediate results, including the test of each candidate threshold{p_end}
{synopt:{opt exportgmm}}export the GMM matrices to Mata for external verification{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
* {opt qx(varname)} is required.{p_end}
{p 4 6 2}
You must {cmd:xtset} your data before using {cmd:xtdpthresh}; see
{manhelp xtset XT}. The time variable must have delta 1.{p_end}
{p 4 6 2}
{it:indepvars} and the variable lists in {opt endogenous()},
{opt predetermined()}, {opt exogenous()}, and {opt iv()} may contain
time-series operators; {it:depvar} may not. Lags of {it:depvar} used as
regressors belong in {opt predetermined()}. Differences,
seasonal differences, and leads of {it:depvar} (D., S., F.) contain its
current or a future value and are rejected in every list; {opt iv()} may
contain only lags of {it:depvar}, subject to the rule given under
{opt iv()}. A lag of a variable declared in
{opt endogenous()} or {opt predetermined()} cannot be declared strictly
exogenous ({it:indepvars} or {opt exogenous()}); declaring a variable
strictly exogenous and a further lag of it predetermined is allowed. These
combinations are rejected with error 198. Each variable may appear in only one
of {it:indepvars}, {opt exogenous()}, {opt endogenous()}, and
{opt predetermined()}.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:xtdpthresh} fits the single-threshold dynamic panel model of Seo and
Shin (2016) by two-step GMM. The individual effects are removed either by
first differences (FD; Arellano and Bond 1991) or by forward orthogonal
deviations (FOD; Arellano and Bover 1995). When periods are missing in the
interior of a unit's time series, FOD discards fewer observations than FD
(Roodman 2009). Regressors may be strictly exogenous, predetermined, or
endogenous. The command extends {cmd:xthenreg} (Seo, Kim, and Kim 2019),
which implements the FD estimator for balanced panels.

{pstd}
The model is

{p 8 8 2}
y_it = x_it'β + (1, x_it')δ·1(q_it > γ) + η_i + ε_it,
i = 1, ..., n,  t = 1, ..., T_i,

{pstd}
where n is the number of panel units, T_i is the number of periods observed
for unit i, x_it is the p x 1 vector of regressors, q_it is the threshold
variable, γ is the threshold, η_i is a unit effect, and ε_it is an
idiosyncratic error. In the default dynamic model, x_it contains
y_{i,t-1}, which is added automatically, and the variables in
{it:indepvars}, {opt endogenous()}, {opt predetermined()}, and
{opt exogenous()}. The vector δ = (δ_1, δ_2')' has p + 1 elements:
δ_1 is the shift in the intercept and δ_2 the shift in the slopes when
q_it > γ. The regressors therefore have slopes β in the lower regime
(q_it ≤ γ) and β + δ_2 in the upper regime.

{pstd}
The slope coefficients are reported with cluster-robust standard errors
that include the estimation error of γ̂. No standard error is reported for
γ̂. Instead, a confidence set for γ is obtained by inverting a bootstrap test
over a grid of candidate thresholds, following the grid bootstrap of Gong and
Seo (2026), which they show to be valid whether or not the regression
function is continuous at γ. The default implementation is a faster
approximation to their algorithm that their results do not cover: it draws
one wild weight per panel unit and holds the weight matrix of the estimator
fixed; see {help xtdpthresh##bootstrap:Bootstrap algorithm}. The command also reports a
bootstrap linearity test, the Hansen J statistic, and the Arellano-Bond AR(1)
and AR(2) tests.

{pstd}
The methods, a Monte Carlo study, and worked examples are described in
{browse "https://ssrn.com/abstract=6619058":Nguyen and Lai (2026)}.


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opt qx(varname)} specifies the threshold variable q_it and is required. To
use q also as a regressor, whose slope can change at γ, list it in
{it:indepvars}, {opt endogenous()}, or {opt predetermined()}, as in
{cmd:qx(debt) endogenous(debt)}.

{phang}
{opt endogenous(varlist)} specifies regressors correlated with the current
error but not with later errors. They are instrumented by their levels from
t-2 under {cmd:method(fd)} and from t-1 under {cmd:method(fod)}.

{phang}
{opt predetermined(varlist)} specifies regressors uncorrelated with current
and later errors that may respond to earlier errors (Arellano and Bond 1991;
Bond 2002). They are instrumented by their levels from t-1 under
{cmd:method(fd)} and from t under {cmd:method(fod)}. Lags of {it:depvar}
used as regressors are declared here.

{phang}
{opt exogenous(varlist)} specifies additional strictly exogenous regressors,
treated like {it:indepvars}.

{phang}
{opt kink} imposes continuity of the regression function at γ (Seo and Shin
2016; Gong and Seo 2026). The model becomes
y_it = x_it'β + δ_q·(q_it - γ)·1(q_it > γ) + η_i + ε_it, where q should also
be a regressor so that it has a slope below γ; otherwise a note is
displayed. The criterion is then continuous in γ, and γ̂ minimizes it over
the interval between the grid points next to the grid minimum
({cmd:e(kink_refined)} = 1 if γ̂ lies between grid points). The standard
errors include the estimation error of γ̂; see {opt vce()}.

{phang}
{opt td} removes common time effects by partialling time dummies out of the
transformed dependent variable, regressors (including the regime
interactions), and instruments. Under {cmd:method(fod)}, and under
{cmd:method(fd)} when every unit has an equation in every period, the
estimates equal those obtained with the dummies among the regressors and
instruments; otherwise the two are asymptotically equivalent. The Windmeijer standard errors and the
Arellano-Bond statistics can differ slightly from that specification.

{phang}
{opt history(panel|sample)} determines how {cmd:if} and {cmd:in} are
applied. Under {cmd:history(panel)}, the default, they restrict the
estimating equations, and observations outside the restriction remain
available as instruments. An equation uses only rows inside the restriction:
under {cmd:method(fd)}, periods t and t-1 must both satisfy it. Under
{cmd:history(sample)}, the restriction also bounds the lag history, and
time-series operators are not allowed; create such lags as variables, as in
{cmd:generate Lx = L.x if L.(}{it:condition}{cmd:)}.

{dlgtab:Instruments}

{pstd}
{it:Instrument set.} The instruments for the equation of unit i in period t
form the vector z_it, and Z is the matrix with rows z_it'. By default, z_it
contains

{phang2}
(i) a constant for each period;{p_end}
{phang2}
(ii) in a dynamic model, the levels of y from t-2 under {cmd:method(fd)} and
from t-1 under {cmd:method(fod)};{p_end}
{phang2}
(iii) the strictly exogenous regressors, transformed as the equation;{p_end}
{phang2}
(iv) the lagged levels of the variables in {opt endogenous()} and
{opt predetermined()}, from the periods given above;{p_end}
{phang2}
(v) the period-t values of the variables in {opt iv()}.{p_end}

{pstd}
Each instrument has one column per period unless {opt collapse} is
specified. A row with no instrument other than the constant is dropped. A
threshold variable that appears only in {opt qx()} supplies no instruments;
its lags can be added with {opt iv()}, as in {cmd:iv(L2.q L3.q)}.

{pstd}
{it:Variables common to all units.} A regressor or {opt iv()} variable that
takes one value for every unit in a period, such as a macro variable or a
trend, is identified through the constants (i). Its own instrument columns
repeat the constants and are dropped; the output reports how many
({cmd:e(N_iv_common)}). The variable must be exactly equal across units, for
example merged by period, and is best declared among the regressors. Under
{opt td}, it is removed by the time effects and rejected.

{pstd}
{it:Dependent instrument columns.} Other instrument columns whose residual
on the preceding columns is below 3.2·10^-7 of their length are dropped
({cmd:e(N_iv_dep)}): linear combinations, as in a period with fewer units
than instrument columns, and also nearly dependent columns. The generalized
inverse of {cmd:xtabond2} also drops linear combinations. Nearly dependent
columns add directions that are lost when dropped; the output then warns
({cmd:e(N_iv_dep_near)}), since the estimates depend on how the instruments
are written. Entering (z2 - z1)/c instead of z2 keeps such a direction.

{phang}
{opt maxlag(# [#])} sets the lags used as internal instruments, like the
{cmd:lag()} suboption of {cmd:xtabond2}. {cmd:maxlag(}{it:L}{cmd:)} uses
lags 1 to {it:L} for predetermined variables and 2 to {it:L} for endogenous
variables and y; {cmd:maxlag(}{it:a b}{cmd:)} uses lags {it:a} to {it:b}
with the same minimums. By default, all available lags are used, and the
number of instruments grows quadratically with T. Under {cmd:method(fod)},
lags are dated as in {cmd:xtabond2} (Roodman 2009), so the same
{opt maxlag()} gives the same number of lags under FD and FOD. With
{it:a} ≥ 2, {opt predetermined()} and {opt endogenous()} receive the same
instruments. A variable may be declared together with its own lag, as in
{cmd:predetermined(x L.x)}; duplicate instrument columns are kept once.

{phang}
{opt collapse} collapses the instruments across periods (Roodman 2009): one
column per lag instead of one per lag and period. Fewer instruments can
reduce the problems caused by too many instruments but can weaken the
identification of γ.

{phang}
{opt iv(varlist [, collapse])} specifies external instruments, entered with
their period-t values; the {cmd:collapse} suboption collapses only these
columns. The period-t value must be uncorrelated with the transformed error:
with ε_it and ε_{i,t-1} under FD, and with ε_it and later errors under FOD.
Accordingly, a lag of {it:depvar} or of a variable in {opt endogenous()} must
be of order 2 or more under FD and 1 or more under FOD, and a lag of a
variable in {opt predetermined()} of order 1 or more under FD; other lags are
rejected with error 198. A missing value contributes a zero instrument in
that row.

{dlgtab:Transformation and threshold search}

{phang}
{opt method(fd|fod)} selects the transformation that removes η_i.
{cmd:fd}, the default, uses first differences, as in Seo and Shin (2016) and
{cmd:xthenreg}, with the first-step weight matrix of Arellano and Bond
(1991). {cmd:fod} uses forward orthogonal deviations (Arellano and Bover
1995) with the first-step weight (Z'Z)^(-1). An interior gap removes three
FD equations and two FOD equations.

{phang}
{opt grid(#)} sets the number of grid points for the threshold search; the
default is 100 and the minimum 10. γ̂ minimizes the GMM criterion over this
grid (under {opt kink}, also between its points). Uniform grids of 100, 199,
and 397 points are nested, so a finer grid can be used to check the
sensitivity of γ̂.

{phang}
{opt gridtype(uniform|quantile)} places the grid points equally spaced in q
({cmd:uniform}, the default) or at empirical quantiles of q
({cmd:quantile}), for both the estimation and the confidence-set grids.

{phang}
{opt trim(#)} sets the total trimming rate of q, split between the two
tails: {cmd:trim(0.10)}, the default, restricts the grid to the 5th to 95th
percentiles of q. The rate must be in [0.01, 0.45].

{phang}
{opt gridsample(effective|observed)} sets which values of q define the
trimming bounds and the regime counts of {opt minregime()}: those that enter
the transformed equations ({cmd:effective}, the default) or those in the
retained rows only ({cmd:observed}).

{phang}
{opt minregime(#)} raises the minimum number of observations in each regime
at every candidate threshold. The default, and the least allowed, is the
number of observations in the smaller trimmed tail, so that both ends of the
grid are admitted.

{phang}
{opt refine(#)} requests # rounds (0 to 20) of local refinement after the
grid search: the observed values of q next to γ̂ are added as candidates and
the second-step criterion is minimized again. The default, {cmd:refine(0)},
keeps γ̂ on the grid. The confidence set and the tests search the initial
grid only. {opt refine()} is not allowed with {opt kink}, whose γ̂ is refined
continuously.

{dlgtab:Inference}

{phang}
{opt gridci(#)} sets the number of candidate thresholds for the confidence
set; the default is 100 and the minimum 10.

{phang}
{opt boot(#)} sets the number of bootstrap replications for the confidence
set and each bootstrap test; the default is 299 and the minimum 10.

{phang}
{opt rseed(#)} sets the seed, from which separate seeds are derived for the
confidence set, the linearity test, and {opt citest()}, so that
changing one of them leaves the draws of the others unchanged. The seed is
not restored after the command. Under {cmd:set rng mt64s}, results also
depend on {cmd:set rngstream}.

{phang}
{opt noboot} skips the confidence set and the linearity test.

{phang}
{opt notest} skips the linearity test.

{phang}
{opt citest(#)} runs the confidence-set bootstrap test only at gamma = #.
It is evaluated after all other bootstrap components with a separate derived
seed, so requesting it does not change the estimate, confidence set, or other
tests. In Monte Carlo work, {cmd:e(citest_accept)} at the true threshold is
the direct pointwise coverage indicator. Statuses 3--6 are not rejections;
they are unavailable or unresolved tests.

{phang}
{opt vce(robust|windmeijer)} selects the variance of the two-step estimates:
the cluster-robust sandwich ({cmd:robust}, the default) or a Windmeijer-type
correction applied to the joint moments linearized in γ. Windmeijer (2005,
sec. 2.1) derives the correction for moment conditions linear in the
parameters; it need not improve the variance estimate for nonlinear moments.
The reported variance includes the estimation error of γ̂ (Seo and Shin
2016). The moment derivative in γ is computed analytically under {opt kink}
and estimated with a kernel in the jump model (see {opt bwscale()}). The
variance conditional on γ̂ is kept in {cmd:e(V_cond)}. See
{help xtdpthresh##notes:Notes and limitations}.

{phang}
{opt bwscale(#)} sets the kernel bandwidth of the jump model's joint
variance, h = #·1.06·s_q·n^(-1/5), with s_q the standard deviation of q and n
the number of units ({cmd:e(gamma_bw)}); the default, 1.5, is that of
{cmd:xthenreg}. Only {cmd:e(V)} and the AR statistics that include γ̂ depend
on it.

{phang}
{opt nocenter} uses the uncentered moment covariance, as {cmd:xtabond2}
does. The default is the centered covariance of Seo and Shin (2016, eq. 11)
and {cmd:xthenreg}.

{phang}
{opt level(#)} sets the confidence level; the default is {cmd:level(95)} or
as set by {helpb set level}.

{dlgtab:Reporting}

{phang}
{opt nowarn} suppresses nonfatal warnings and most notes. Notes on how the
reported results were formed, such as dropped instrument columns, a
confidence set with several intervals, or a missing Hansen J, are always
shown. The flag for a confidence-set bound at the edge of the grid is stored
in {cmd:e(boundary_warn)} regardless.

{phang}
{opt verbose} displays intermediate results, including the statistic,
critical value, and decision for each candidate threshold of the confidence
set.

{phang}
{opt exportgmm} stores the GMM matrices in Mata ({bf:xdpt_best_A},
{bf:xdpt_best_Z_f}, {bf:xdpt_best_X_f}, {bf:xdpt_best_Xar}) for external
verification.


{marker remarks}{...}
{title:Remarks}

{pstd}
Remarks are presented under the following headings:

{phang2}{help xtdpthresh##model:Moment conditions and assumptions}{p_end}
{phang2}{help xtdpthresh##bootstrap:Bootstrap algorithm for the threshold confidence set}{p_end}
{phang2}{help xtdpthresh##tests:Tests and diagnostics}{p_end}
{phang2}{help xtdpthresh##notes:Notes and limitations}{p_end}
{phang2}{help xtdpthresh##xthenreg:Differences from xthenreg}{p_end}
{phang2}{help xtdpthresh##labels:Coefficient labels}{p_end}


{marker model}{...}
{title:Moment conditions and assumptions}

{pstd}
Write the model as y_it = W_it(γ)'θ + η_i + ε_it, where
W_it(γ) = (x_it', (1, x_it')·1(q_it > γ))' and θ = (β', δ')'. With T(·) the
transformation selected by {opt method()} and z_it the instruments listed
under {help xtdpthresh##options:Options}, the estimator uses the moment
conditions E[z_it·T(ε)_it] = 0. Let
g_i(θ,γ) = Σ_t z_it·[T(y)_it - T(W(γ))_it'θ] and
ḡ_n(θ,γ) = n^(-1) Σ_i g_i(θ,γ). For a weight matrix A, the criterion is
Q(θ,γ) = ḡ_n' A ḡ_n; θ̂(γ) has a closed form, and γ̂ minimizes
Q(θ̂(γ),γ) over the grid. The first step uses the weight W_1 given under
{opt method()} and yields (θ̂_1, γ̂_1); the second step uses the inverse W_2
of the clustered covariance of g_i(θ̂_1, γ̂_1) and repeats the grid search.
Under {opt kink}, the second-step criterion is then minimized between the
grid points next to its grid minimum.

{pstd}
{bf:Assumptions.} The moment conditions hold if (a) ε_it has mean zero
conditional on η_i and on y, x, and q up to period t-1, so that ε_it is
serially uncorrelated; and (b) each regressor is correctly classified as
strictly exogenous, predetermined, or endogenous. The AR(2) and Hansen
tests are diagnostics for these assumptions. Seo and Shin (2016) establish the
asymptotic theory of the FD estimator when the regression function is
discontinuous at γ: γ̂ and the slopes are √n-consistent and jointly normal.
When it is continuous, the unrestricted estimator of γ converges at the rate
n^(1/4) to a nonnormal limit (Gong and Seo 2026, Theorem 2), whereas the
estimator with {opt kink} is asymptotically normal. No such results exist
for FOD; Nguyen and Lai (2026) report Monte Carlo evidence for both
transformations.

{pstd}
{bf:Unbalanced panels.} Let s_i = (s_i1, ..., s_iT) mark the periods observed
for unit i. The retained equations and their instruments depend on s_i as a
whole: under FD the equation for t needs period t-1, and an equation without
an instrument other than the constant is dropped. Estimation further assumes
(U1) E[z_it·T(ε)_it | s_i] = 0 for every retained equation, as when s_i is
independent of the initial value, the unit effect, the errors, and the
regressors; (U2) the identification and rank conditions of Seo and Shin
(2016) for the retained moments, with observed values of q_it near the
threshold on both sides; and (U3) under FOD, a later observation for every
retained equation, which the command imposes. (U1) fails, for example, when
units leave after a large shock, and (U2) when attrition leaves few
observations on one side of the threshold. Either failure can bias γ̂, and
neither can be checked from the data used by the command.

{pstd}
{bf:One threshold.} The model has one threshold and two regimes; the theory
cited above covers only this case.


{marker bootstrap}{...}
{title:Bootstrap algorithm for the threshold confidence set}

{pstd}
The confidence set collects the candidate thresholds γ_ℓ that a bootstrap
test does not reject (Hansen 1999b; Gong and Seo 2026, sec. 4.1). The
candidates are the {opt gridci()} grid and γ̂. The test uses the criterion Q
of the reported estimator with its weight matrix held at the sample value:
W_2 for the two-step estimates, and W_1 after a one-step fallback, when the
minimizer of Q over the grid is also a candidate ({cmd:e(ci_criterion)}).
Let Γ be the initial estimation grid, without the points added by
{opt refine()} or under {opt kink}. For each γ_ℓ:

{phang2}
1. {it:Restricted fit.} Estimate θ̂(γ_ℓ) = argmin_θ Q(θ, γ_ℓ) and the
residuals e_it = T(y)_it - T(W(γ_ℓ))_it'θ̂(γ_ℓ). If this fit fails the
numerical check of the command, γ_ℓ is not admissible when the reported
estimator cannot use it either (as when q·1(q > γ_ℓ) = q) and unresolved
otherwise.{p_end}

{phang2}
2. {it:Statistic.} D_n(γ_ℓ) = Q(θ̂(γ_ℓ), γ_ℓ) - m_n, where m_n is the
minimum of Q(θ̂(γ), γ) over Γ and γ_ℓ. D_n(γ̂) = 0, so γ̂ belongs to the
set.{p_end}

{phang2}
3. {it:Bootstrap samples under the null.} For b = 1, ..., B, draw one weight
η_i^b per panel unit from the two-point distribution of Mammen (1993) and
set T(y)*_it = T(W(γ_ℓ))_it'θ̂(γ_ℓ) + η_i^b·e_it for every period of unit i;
regressors and instruments keep their sample values. The same weights serve
every γ_ℓ.{p_end}

{phang2}
4. {it:Bootstrap statistic.} Compute D*_n,b(γ_ℓ) as in step 2, with the same
weight matrix and grid.{p_end}

{phang2}
5. {it:Decision.} Accept γ_ℓ if D_n(γ_ℓ) does not exceed the k-th smallest
of the D*_n,b(γ_ℓ), with k = ceil(p·(B + 1)) and p = {opt level()}/100, that
is, if the bootstrap p-value (1 + R)/(1 + B), defined as for the linearity test below,
exceeds 1 - p. If B < p/(1 - p), every candidate is accepted. If a draw is
not valid (a statistic is not finite), γ_ℓ is unresolved. With any
unresolved γ_ℓ, no confidence set is reported ({cmd:e(ci_incomplete)} =
1).{p_end}

{pstd}
A single weight multiplies all residuals of a unit, which preserves their
correlation within the unit and heteroskedasticity across units. Because the
data are generated from the restricted fit and E*[η_i] = 0, the null
γ = γ_ℓ holds in the bootstrap population without recentering.
{cmd:e(gamma_lo)} and {cmd:e(gamma_hi)} are the smallest and largest
accepted values; if the set has several intervals ({cmd:e(ci_nseg)} > 1),
the range between them also contains rejected or inadmissible values. {cmd:e(ci_grid)} holds
the full inversion table.

{pstd}
{it:Relation to Gong and Seo (2026).} Their Algorithm 1 resamples whole
units, generates the bootstrap data from the restricted coefficients and the
unrestricted residuals with explicit recentering, re-estimates the weight
matrix in every draw, and is proved valid (their Theorem 5; uniformly, for a
simplified model, their Theorem I.1). The default scheme above, with wild
weights, the restricted residuals, and the weight matrix held fixed, is a
faster approximation that these results do not cover; its coverage is
examined by simulation in Nguyen and Lai (2026).


{marker tests}{...}
{title:Tests and diagnostics}

{phang}
{bf:Linearity test} (H0: δ = 0). The statistic is the one-step criterion of
the linear model minus the minimum one-step criterion of the threshold model
over the initial grid; bootstrap samples are drawn from the linear fit, as in
the algorithm above. This GMM distance test differs from the sup-Wald test of
Seo and Shin (2016) reported by {cmd:xthenreg}.

{pstd}
The p-value is (1 + R)/(1 + B), where R is the number of
draws whose statistic is at least as large as the sample statistic; if a draw
is not valid, no p-value is reported.

{phang}
{bf:Hansen J statistic}, with N_iv - k - 1 degrees of freedom, where k is the
number of slope coefficients and γ counts as one estimated parameter. The
chi-square reference presumes regular identification of the slopes and γ,
which fails near a continuous threshold and under weak identification. Treat
J as a diagnostic: it loses power with many instruments (Roodman 2009), and
a large p-value does not show that the instruments are valid.

{phang}
{bf:Arellano-Bond AR(1) and AR(2) tests} on the first-difference residuals,
for both transformations, using the full statistic of Arellano and Bond
(1991, eq. 8). With the joint variance, the statistics include the
estimation of γ̂ ({cmd:e(ar_joint)} = 1); those that treat γ̂ as known, equal
to the statistics of {cmd:xtabond2} with γ fixed at γ̂, are {cmd:e(ar1_cond)}
and {cmd:e(ar2_cond)}. Near a continuous threshold in the jump model the tests
are diagnostics. AR(1) is expected to reject; a
rejection of AR(2) indicates serial correlation in ε_it, which invalidates
the lagged instruments.


{marker notes}{...}
{title:Notes and limitations}

{phang}
{bf:Sample size.} The theory is for many units and a fixed number of
periods. With few units, γ̂ and its confidence set can be unstable.

{phang}
{bf:Standard errors.} In the unrestricted (jump) model, the joint variance
follows the asymptotics of Seo and Shin (2016), in which the jump is not
close to zero; near a continuous threshold neither it nor the conditional
variance is reliable, and inference on γ rests on the confidence set. Under
{opt kink}, with the restriction true, a nonzero kink, and full rank, the
estimator is asymptotically normal (Gong and Seo 2026). If the joint variance
cannot be computed, the conditional one is reported ({cmd:e(joint_vce)} = 0);
its Windmeijer correction is then approximate when
{cmd:e(wind_same_threshold)} = 0.

{phang}
{bf:Discrete threshold variable.} In the jump model, the joint variance and
the AR statistics with γ̂ assume that q has a continuous density, positive at
the threshold (Seo and Shin 2016, Assumption 2). The output warns when q
takes fewer than 10 distinct values within two bandwidths of γ̂
({cmd:e(q_nvals_bw)}); {cmd:e(V_cond)}, {cmd:e(ar1_cond)}, and
{cmd:e(ar2_cond)} do not use that assumption.

{phang}
{bf:Regressors removed by the transformation.} A regressor that does not
vary over time within units, or, under {opt td}, one that is common to all
units in each period or has the form a_i + g_t (such as firm age), is
removed by the transformation and the time effects and is rejected with
error 498.

{phang}
{bf:Fallbacks.} If the second step cannot be completed (a singular
second-step weight matrix, no candidate threshold that can be solved with
it, or a flat second-step criterion), the one-step estimates are reported
({cmd:e(estimator_twostep)} = 0) without Hansen J. If the first-step weight
is singular, (Z'Z)^(-1) or the identity replaces it ({cmd:e(W1_fallback)}).
If the Windmeijer-corrected variance cannot be computed or is not positive
semidefinite, the cluster-robust one is reported ({cmd:e(vce_applied)} = 0).
These cases are noted in the output unless {opt nowarn} is specified.


{marker xthenreg}{...}
{title:Differences from xthenreg}

{pstd}
{cmd:xtdpthresh} follows the conventions of {cmd:xthenreg} (Seo, Kim, and
Kim 2019) where possible. Under {cmd:method(fd)} it uses the same moment
conditions and weight matrices and, on the same grid, gives the same
estimates. The differences are as follows:

{phang}
{cmd:*} the threshold variable is given in {opt qx()};{p_end}
{phang}
{cmd:*} {cmd:method(fod)}, {opt td}, and the instrument options
{opt maxlag()}, {opt collapse}, and {opt iv()} are added, and regressors can
be declared predetermined or endogenous; strictly exogenous regressors are
instrumented by their transformed values, whereas the {cmd:exogenous()}
option of {cmd:xthenreg} uses their levels;{p_end}
{phang}
{cmd:*} the default trimming rate is 0.10 instead of 0.4;{p_end}
{phang}
{cmd:*} under {opt kink}, γ̂ minimizes the criterion between grid points, not
only on the grid;{p_end}
{phang}
{cmd:*} no standard error is reported for γ̂; the slope standard errors are
those of the joint variance, in the jump model with the kernel and default
bandwidth of {cmd:xthenreg} ({opt bwscale()});{p_end}
{phang}
{cmd:*} inference on γ uses the grid-bootstrap confidence set; the
linearity test uses a GMM distance statistic with Mammen weights instead of
the sup-Wald statistic with normal weights.{p_end}


{marker labels}{...}
{title:Coefficient labels}

{pstd}
The coefficients in {cmd:e(b)} are labeled as in {cmd:xthenreg}.

{pstd}
{bf:Unrestricted (jump) model}, 2p + 1 coefficients:

{phang2}{cmd:Lag_y_b}{space 7}slope on L.{it:depvar} in the lower regime (dynamic model){p_end}
{phang2}{it:var}{cmd:_b}{space 9}slope on {it:var} in the lower regime{p_end}
{phang2}{cmd:cons_d}{space 8}intercept shift δ_1{p_end}
{phang2}{cmd:Lag_y_d}{space 7}shift in the slope on L.{it:depvar} (dynamic model){p_end}
{phang2}{it:var}{cmd:_d}{space 9}shift in the slope on {it:var}{p_end}

{pstd}
The slope in the upper regime is the sum of the {cmd:_b} and {cmd:_d}
coefficients; for example, {cmd:Lag_y_b} + {cmd:Lag_y_d}.

{pstd}
{bf:Kink model}, p + 1 coefficients:

{phang2}{cmd:Lag_y_b}{space 7}slope on L.{it:depvar} (dynamic model){p_end}
{phang2}{it:var}{cmd:_b}{space 9}slope on {it:var}{p_end}
{phang2}{cmd:kink_slope}{space 4}change in the slope on q at γ, δ_q{p_end}


{marker examples}{...}
{title:Examples}

{pstd}
The examples use the firm investment panel of Hansen (1999a), 565 firms
observed from 1973 to 1987, which is installed as an ancillary file by
{cmd:ssc install xtdpthresh, all}. The variables are {cmd:invest}
(investment over capital), {cmd:tobin_q}, {cmd:cashflow}, and {cmd:debt}
(leverage).

{pstd}Setup{p_end}
{phang2}{cmd:. use invest, clear}{p_end}
{phang2}{cmd:. xtset firm year}{p_end}

{pstd}
Point estimates and diagnostics only, FD, instruments restricted to lags 2 to
4{p_end}
{phang2}{cmd:. xtdpthresh invest tobin_q cashflow, qx(debt) maxlag(2 4) noboot}{p_end}

{pstd}
Full inference: threshold confidence set and linearity test, with a seed
for reproducibility{p_end}
{phang2}{cmd:. xtdpthresh invest tobin_q cashflow debt, qx(debt) maxlag(2 4) rseed(12345)}{p_end}


{pstd}
Cash flow treated as endogenous{p_end}
{phang2}{cmd:. xtdpthresh invest tobin_q, qx(debt) endogenous(cashflow) maxlag(2 4) rseed(12345)}{p_end}

{pstd}
Kink model; debt is also a regressor, so that its slope can change at the
threshold{p_end}
{phang2}{cmd:. xtdpthresh invest tobin_q cashflow debt, qx(debt) maxlag(2 4) kink rseed(12345)}{p_end}

{pstd}
Unbalanced panel: drop 15% of the observations at random, keep firms with
at least four remaining observations, and compare FD with FOD{p_end}
{phang2}{cmd:. preserve}{p_end}
{phang2}{cmd:. set seed 20260418}{p_end}
{phang2}{cmd:. drop if runiform() < 0.15}{p_end}
{phang2}{cmd:. bysort firm: drop if _N < 4}{p_end}
{phang2}{cmd:. xtdpthresh invest tobin_q cashflow debt, qx(debt) method(fd) maxlag(2 4) rseed(12345)}{p_end}
{phang2}{cmd:. xtdpthresh invest tobin_q cashflow debt, qx(debt) method(fod) maxlag(2 4) rseed(12345)}{p_end}
{phang2}{cmd:. restore}{p_end}

{pstd}
External instrument: the yearly cross-sectional mean of Tobin's q
instruments {cmd:tobin_q}, treated as endogenous{p_end}
{phang2}{cmd:. bysort year: egen double peer_q = mean(tobin_q)}{p_end}
{phang2}{cmd:. xtdpthresh invest cashflow, qx(debt) endogenous(tobin_q) iv(peer_q) maxlag(2 4) collapse noboot}{p_end}


{marker postest}{...}
{title:Postestimation: predict}

{p 8 17 2}
{cmd:predict} {dtype} {newvar} {ifin}{cmd:,} {it:statistic}

{synoptset 22 tabbed}{...}
{synopthdr:statistic}
{synoptline}
{synopt:{opt r:esiduals}}residuals of the estimated equation: FD residuals under {cmd:method(fd)}, FOD residuals under {cmd:method(fod)}{p_end}
{synopt:{opt arresiduals}}the first-difference residuals used by the AR(1) and AR(2) tests; equal to {cmd:residuals} under {cmd:method(fd)}{p_end}
{synopt:{opt xb}}fitted values of the transformed equation, so that {cmd:residuals} + {cmd:xb} equals the transformed {it:depvar}{p_end}
{synopt:{opt reg:ime}}regime indicator 1(q_it > γ̂) for all observations in {cmd:if} and {cmd:in}{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
A statistic must be specified. {cmd:residuals} and {cmd:xb} are defined for
the estimation rows ({cmd:e(sample)}), {cmd:arresiduals} for the rows of the
first-difference equations. They are the series computed during estimation,
not recomputed, so use {cmd:arresiduals} to reproduce the AR statistics under
{cmd:method(fod)}. Under {opt td}, they are those of the equations after the
time effects are partialled out. The 20 most recent fits are kept, so {cmd:estimates restore}
works. {cmd:predict} exits with error 459 if the data have changed since the
fit, and with error 498 if the fit is no longer in memory (for example after
{cmd:mata clear}).{p_end}

{pstd}Example{p_end}
{phang2}{cmd:. xtdpthresh invest tobin_q cashflow debt, qx(debt) maxlag(2 4) noboot}{p_end}
{phang2}{cmd:. predict double ehat, residuals}{p_end}
{phang2}{cmd:. predict byte upper, regime}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:xtdpthresh} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of transformed equations{p_end}
{synopt:{cmd:e(N_units)}}number of panel units{p_end}
{synopt:{cmd:e(N_clust)}}number of clusters (panel units) of the VCE and the bootstrap{p_end}
{synopt:{cmd:e(N_iv)}}number of instruments{p_end}
{synopt:{cmd:e(N_iv_common)}}number of instrument columns dropped as combinations of the constant instruments{p_end}
{synopt:{cmd:e(N_iv_dep)}}number of other instrument columns dropped as numerically dependent{p_end}
{synopt:{cmd:e(N_iv_dep_near)}}number of those that are not linear combinations of the kept columns (relative residual above 10^-10){p_end}
{synopt:{cmd:e(iv_dep_res)}}largest relative residual of a dropped column on the kept columns{p_end}
{synopt:{cmd:e(N_switch)}}number of units whose regime changes within their equations at γ̂{p_end}
{synopt:{cmd:e(gamma)}}threshold estimate γ̂{p_end}
{synopt:{cmd:e(obj)}}GMM criterion at the estimates, N·ḡ'Aḡ with ḡ the mean moment over the {cmd:e(N)} equations{p_end}
{synopt:{cmd:e(gamma_lo)}}smallest accepted threshold{p_end}
{synopt:{cmd:e(gamma_hi)}}largest accepted threshold{p_end}
{synopt:{cmd:e(ci_nseg)}}number of intervals in the confidence set{p_end}
{synopt:{cmd:e(ci_empty)}}1 if no candidate is accepted{p_end}
{synopt:{cmd:e(ci_incomplete)}}1 if some candidates are unresolved; no confidence set is then reported{p_end}
{synopt:{cmd:e(boundary_warn)}}bound at the edge of the confidence-set grid: 0 none, 1 lower, 2 upper, 3 both{p_end}
{synopt:{cmd:e(level)}}confidence level{p_end}
{synopt:{cmd:e(pval_lin)}}p-value of the linearity test{p_end}
{synopt:{cmd:e(citest_gamma)}}null threshold supplied in {opt citest()}{p_end}
{synopt:{cmd:e(citest_D)}}pointwise threshold-test statistic{p_end}
{synopt:{cmd:e(citest_crit)}}bootstrap critical value; missing for mechanical acceptance{p_end}
{synopt:{cmd:e(citest_accept)}}1 if the pointwise null is accepted, 0 if rejected{p_end}
{synopt:{cmd:e(citest_p)}}add-one bootstrap p-value{p_end}
{synopt:{cmd:e(citest_status)}}1 valid, 2 mechanical accept, 3 structural, 4 solve failed, 5 incomplete draws, 6 quantile failed{p_end}
{synopt:{cmd:e(citest_draws)}}number of valid draws; missing for mechanical acceptance{p_end}
{synopt:{cmd:e(seed_citest)}}component seed used by {opt citest()}{p_end}
{synopt:{cmd:e(hansen)}}Hansen J statistic{p_end}
{synopt:{cmd:e(hansen_df)}}degrees of freedom of J{p_end}
{synopt:{cmd:e(hansen_p)}}p-value of J{p_end}
{synopt:{cmd:e(ar1)}, {cmd:e(ar2)}}Arellano-Bond AR(1) and AR(2) statistics{p_end}
{synopt:{cmd:e(ar1_p)}, {cmd:e(ar2_p)}}p-values of the AR tests{p_end}
{synopt:{cmd:e(ar_joint)}}1 if the AR statistics include the estimation of γ̂{p_end}
{synopt:{cmd:e(ar1_cond)}, {cmd:e(ar2_cond)}}AR(1) and AR(2) statistics with γ̂ treated as known{p_end}
{synopt:{cmd:e(q_lo)}, {cmd:e(q_hi)}}trimming bounds of the threshold grid{p_end}
{synopt:{cmd:e(grid_admitted)}}number of admitted grid points{p_end}
{synopt:{cmd:e(estimator_twostep)}}1 if the estimates are two-step, 0 for the one-step fallback{p_end}
{synopt:{cmd:e(W1_fallback)}}1 or 2 if (Z'Z)^(-1) or the identity replaced a singular first-step weight{p_end}
{synopt:{cmd:e(joint_vce)}}1 if {cmd:e(V)} includes the estimation error of γ̂, 0 if it is conditional on γ̂{p_end}
{synopt:{cmd:e(gamma_bw)}}jump model: kernel bandwidth of the joint variance{p_end}
{synopt:{cmd:e(q_nvals_bw)}}jump model: number of distinct values of q within two bandwidths of γ̂{p_end}
{synopt:{cmd:e(kink_refined)}}under {opt kink}, 1 if γ̂ lies between grid points{p_end}
{synopt:{cmd:e(rseed)}}seed specified in {opt rseed()}; missing under {opt noboot}{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:xtdpthresh}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(cmdversion)}}version of {cmd:xtdpthresh}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(q_var)}}name of threshold variable{p_end}
{synopt:{cmd:e(indepvars)}}regressors in {it:indepvars}{p_end}
{synopt:{cmd:e(exog_extra)}}regressors in {opt exogenous()}{p_end}
{synopt:{cmd:e(endog)}}endogenous regressors{p_end}
{synopt:{cmd:e(predet)}}predetermined regressors{p_end}
{synopt:{cmd:e(inst)}}external instruments{p_end}
{synopt:{cmd:e(iv_common)}}variables common to all units in each period{p_end}
{synopt:{cmd:e(method)}}{cmd:fd} or {cmd:fod}{p_end}
{synopt:{cmd:e(panelvar)}}name of panel variable{p_end}
{synopt:{cmd:e(timevar)}}name of time variable{p_end}
{synopt:{cmd:e(vce)}}{cmd:robust} or {cmd:windmeijer}{p_end}
{synopt:{cmd:e(clustvar)}}panel variable (the cluster variable){p_end}
{synopt:{cmd:e(vcetype)}}title used to label Std. err.{p_end}
{synopt:{cmd:e(boottype)}}{cmd:wild}{p_end}
{synopt:{cmd:e(ci_criterion)}}criterion of the confidence-set test: {cmd:twostep} or {cmd:onestep}{p_end}
{synopt:{cmd:e(predict)}}program used to implement {cmd:predict}{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators; see {opt vce()}{p_end}
{synopt:{cmd:e(V_cond)}}variance conditional on γ̂, when {cmd:e(joint_vce)} = 1{p_end}
{synopt:{cmd:e(ci_segments)}}intervals of the confidence set{p_end}
{synopt:{cmd:e(ci_grid)}}inversion table: candidate γ, statistic (on the scale of {cmd:e(obj)}), critical value, accepted flag, valid draws, and status (1 resolved; 2 accepted with a zero statistic; 3 not admissible; 4-6 unresolved){p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks the period of each transformed equation{p_end}
{p2colreset}{...}

{pstd}
Further results record the grid, the refinement, the bootstrap draws, the
seeds, and the data signature used by
{cmd:predict}; type {cmd:ereturn list} to see them.


{marker references}{...}
{title:References}

{phang}
Arellano, M., and S. Bond. 1991. Some tests of specification for panel
data: Monte Carlo evidence and an application to employment equations.
{it:Review of Economic Studies} 58: 277-297.

{phang}
Arellano, M., and O. Bover. 1995. Another look at the instrumental variable
estimation of error-components models. {it:Journal of Econometrics} 68:
29-51.

{phang}
Bond, S. R. 2002. Dynamic panel data models: A guide to micro data methods
and practice. {it:Portuguese Economic Journal} 1: 141-162.

{phang}
Gong, W., and M. H. Seo. 2026. Bootstraps for dynamic panel threshold
models. {it:Journal of Econometrics} 253: 106153.

{phang}
Hansen, B. E. 1999a. Threshold effects in non-dynamic panels: Estimation,
testing, and inference. {it:Journal of Econometrics} 93: 345-368.

{phang}
------. 1999b. The grid bootstrap and the autoregressive model.
{it:Review of Economics and Statistics} 81: 594-607.

{phang}
Mammen, E. 1993. Bootstrap and wild bootstrap for high dimensional linear
models. {it:Annals of Statistics} 21: 255-285.

{phang}
Nguyen, D. C., and N. D. Lai. 2026. xtdpthresh: Dynamic panel threshold
regression for unbalanced panels, with endogenous regressors and
continuity-robust inference. SSRN Working Paper 6619058.
{browse "https://ssrn.com/abstract=6619058"}.

{phang}
Roodman, D. 2009. How to do xtabond2: An introduction to difference and
system GMM in Stata. {it:Stata Journal} 9: 86-136.

{phang}
Seo, M. H., S. Kim, and Y.-J. Kim. 2019. Estimation of dynamic panel
threshold model using Stata. {it:Stata Journal} 19: 685-697.

{phang}
Seo, M. H., and Y. Shin. 2016. Dynamic panels with threshold effect and
endogeneity. {it:Journal of Econometrics} 195: 169-186.

{phang}
Windmeijer, F. 2005. A finite sample correction for the variance of linear
efficient two-step GMM estimators. {it:Journal of Econometrics} 126: 25-51.


{marker authors}{...}
{title:Authors}

{pstd}
Duy Chinh Nguyen{break}
School of Business, International University, Ho Chi Minh City, Vietnam{break}
Vietnam National University Ho Chi Minh City, Vietnam{break}
{browse "mailto:ndchinh@hcmiu.edu.vn":ndchinh@hcmiu.edu.vn}{break}
ORCID: {browse "https://orcid.org/0000-0002-9157-9358":0000-0002-9157-9358}

{pstd}
Nhat Duy Lai (corresponding author){break}
Faculty of Finance and Accounting, Saigon University{break}
273 An Duong Vuong St, Cho Quan Ward, Ho Chi Minh City, Vietnam{break}
{browse "mailto:lnduy@sgu.edu.vn":lnduy@sgu.edu.vn}{break}
ORCID: {browse "https://orcid.org/0009-0008-5365-2893":0009-0008-5365-2893}


{title:Also see}

{psee}
Manual:  {manlink XT xtabond}

{psee}
Online:  {manhelp xtset XT}, {help xthenreg} (Seo, Kim, and Kim 2019, if
installed), {help xtabond2} (Roodman 2009, if installed){p_end}


{marker versionhistory}{...}
{title:Version}

{pstd}
This help file documents {cmd:xtdpthresh} 0.9.37; type {cmd:which xtdpthresh}
to see the installed version. Earlier releases, including the one posted on
SSC, contain errors that have since been corrected and should not be used.
The list of changes is kept at the end of the file {cmd:xtdpthresh.ado}.
{p_end}
