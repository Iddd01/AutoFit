# Python prototypes for the xtdpthresh 0.9.36 review

Independent numpy re-implementations on the Gong-Seo (2026) benchmark DGP
(FD, T=6, rho_y=.6, beta_q=1, delta3=2, gamma0=.25, sigma=.5, rho_q=.7,
rho_eu=.5; jump size kappa = delta1 + delta3*gamma0). They do not call
xtdpthresh; they reproduce its algorithms to separate properties of the
method from properties of the Stata code.

| Script | Question |
|---|---|
| `gs_sim.py` | Is gamma weakly identified in this DGP under FD-GMM (full vs short instrument ladder)? |
| `gs_boot.py` | Size of the grid-bootstrap test at gamma0 under four schemes: M0 wild fixed-W2 (= 0.9.35), M0u wild unrestricted residuals, M1 wild re-estimated W2, M2 unit resampling (Alg. 1-oriented) |
| `gs_set_cont.py` | `set`: pointwise vs hull vs union-of-segments coverage of the same inverted set; `cont`: continuity-test size/power for W1/W2 x kink/jump residuals |
| `gs_refine.py` | Does `cirefine()`-style boundary refinement recover union-of-segments coverage? |

Usage (arguments: N kappa R B):

    python3 gs_boot.py 400 1 400 199
    python3 gs_set_cont.py set 400 1 300 199
    python3 gs_set_cont.py cont 400 2 300 199
    python3 gs_refine.py 400 1 300 199

Results obtained in this review (R = 300-400, B = 199):

- Test at gamma0, B6 design: M0 94-96%, M2 95-99.8% (Gong-Seo Table 1:
  95.5-99.2%). The 0.9.35 wild procedure has correct size at gamma0.
- Same samples, 0.9.35 wild: pointwise 95.3-96.0%, hull 96.7-99.7%,
  union of segments 85.3-93.7%.
- Continuity (N=400): size 3.3% (0.9.35) vs 1.7% (0.9.36); power at
  kappa=1 9.7% vs 9.3%, kappa=2 45% vs 67%, kappa=3 74% vs 94%.
