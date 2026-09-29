"""Generate the registries of the final Monte Carlo (xtdpthresh 0.9.37).
Main axis: the Gong-Seo (2026) design and geometry (all available lags,
46-point p5-p95 quantile grid, trim .10, observed support, no refinement,
B = 500); every other cell changes one factor: FD -> FOD, balanced -> MCAR
gaps / monotone attrition, N, kappa, the DGP, the instrument set, T, or the
model (kink). Writes final_point_cells.csv and final_inf_cells.csv."""
import csv

HDR = ["block", "cell_id", "pair_id", "dgp", "spec", "method", "vce", "iv", "N", "T",
       "miss", "missp", "kappa", "c", "mode", "R", "B", "grid", "gridci",
       "gridtype", "gridsample", "trim", "refine"]
MISS = {"bal": ("balanced", 0), "gap30": ("mcar", .3), "attr15": ("attrition", .15)}
R = 500


def cell(rows, block, dgp, spec, vce, iv, N, T, miss, kappa, c, mode, B, gridci):
    m, p = MISS[miss]
    k = "%03d" % round(100 * kappa)
    cc = "%03d" % round(100 * c)
    pair = f"{block}_{dgp}_{spec}_{vce}_{iv}_n{N}_t{T}_{miss}_k{k}_c{cc}"
    for method in ("fd", "fod"):
        rows.append(dict(block=block, cell_id=0, pair_id=pair, dgp=dgp, spec=spec,
                         method=method, vce=vce, iv=iv, N=N, T=T, miss=m, missp=p,
                         kappa=kappa, c=c, mode=mode, R=R, B=B, grid=46, gridci=gridci,
                         gridtype="quantile", gridsample="observed", trim=.10, refine=0))


def point():
    r = []
    for N in (400, 800, 1600):                       # P1 main response surface
        for miss in MISS:
            for k in (0, .1, .2, .5, 1):
                cell(r, "P1", "base", "jump", "robust", "all", N, 6, miss, k, 0, "POINT", 0, 0)
    for miss in ("bal", "gap30"):
        for k in (0, 1):
            cell(r, "P2", "base", "jump", "robust", "all", 200, 6, miss, k, 0, "POINT", 0, 0)
            for dgp in ("endog", "persist", "heavy"):  # P3 robustness
                cell(r, "P3", dgp, "jump", "robust", "all", 400, 6, miss, k, 0, "POINT", 0, 0)
            for iv in ("L3", "collapse"):             # P4 instrument set, T = 6
                cell(r, "P4", "base", "jump", "robust", iv, 400, 6, miss, k, 0, "POINT", 0, 0)
            for N in (200, 400):                     # P5 T = 10, three instrument sets
                for iv in ("all", "L3", "collapse"):
                    cell(r, "P5", "base", "jump", "robust", iv, N, 10, miss, k, 0, "POINT", 0, 0)
        # P6 covariance: Windmeijer for the jump; robust and Windmeijer for the kink
        cell(r, "P6", "base", "jump", "windmeijer", "all", 400, 6, miss, 1, 0, "POINT", 0, 0)
        for vce in ("robust", "windmeijer"):
            cell(r, "P6", "base", "kink", vce, "all", 400, 6, miss, 0, 0, "POINT", 0, 0)
    return r


def inf():
    r = []
    for N in (400, 800):                             # I1 coverage + linearity power
        for miss in MISS:
            for k in (0, .1, .2, .5, 1):
                cell(r, "I1", "base", "jump", "robust", "all", N, 6, miss, k, 0, "FULL", 500, 46)
    for miss in ("bal", "gap30"):
        for k in (0, 1):
            cell(r, "I2", "base", "jump", "robust", "all", 200, 6, miss, k, 0, "FULL", 500, 46)
            cell(r, "I2", "endog", "jump", "robust", "all", 400, 6, miss, k, 0, "FULL", 500, 46)
    for N in (400, 800):                             # I3 linearity size
        for miss in MISS:
            cell(r, "I3", "linear", "jump", "robust", "all", N, 6, miss, 0, 0, "LIN", 500, 10)
        for miss in ("bal", "gap30"):                # I4 kink model
            cell(r, "I4", "base", "kink", "robust", "all", N, 6, miss, 0, 0, "CI", 500, 46)
    for miss in MISS:                                # I5 power, N = 400
        for k in (0, .5, 1):
            for c in (.1, .25, .5):
                cell(r, "I5", "base", "jump", "robust", "all", 400, 6, miss, k, c, "CI", 500, 10)
    for k in (0, .5, 1):                             # I5 power, N = 800 balanced
        for c in (.1, .25, .5):
            cell(r, "I5", "base", "jump", "robust", "all", 800, 6, "bal", k, c, "CI", 500, 10)
    for miss in ("bal", "gap30"):                    # I6 restricted lags
        for k in (0, 1):
            cell(r, "I6", "base", "jump", "robust", "L3", 400, 6, miss, k, 0, "FULL", 500, 46)
            for c in (.25, .5):
                cell(r, "I6", "base", "jump", "robust", "L3", 400, 6, miss, k, c, "CI", 500, 10)
    return r


for name, rows in (("final_point_cells.csv", point()), ("final_inf_cells.csv", inf())):
    for i, row in enumerate(rows, 1):
        row["cell_id"] = i
    assert len({(x["pair_id"], x["method"]) for x in rows}) == len(rows)
    with open(name, "w", newline="\n") as f:
        w = csv.DictWriter(f, HDR, lineterminator="\n")
        w.writeheader()
        for row in rows:
            row = {k: (("%g" % v) if isinstance(v, float) else v) for k, v in row.items()}
            w.writerow(row)
    print(name, len(rows), "cells", len(rows) * R, "fits")
