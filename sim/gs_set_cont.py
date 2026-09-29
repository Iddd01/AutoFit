"""(1) Set coverage (union of accepted-run segments on the CI grid, as Study B
scores it) vs pointwise acceptance at gamma0, M0 = xtdpthresh 0.9.35 wild.
(2) Continuity test variants: size (kappa=0) / power (kappa>0)."""
import numpy as np, sys, time
from gs_boot import dgp, Data, solve_all, omega_inv, mammen, estimate

def m0_accept(d, est_W2, i, rng, B, ZX, ZY, o2):
    """wild fixed-W2 test at candidate index i of d.gams (grid entries 1..)."""
    N = d.N
    th = solve_all(ZX[i:i+1], ZY, est_W2)[1][0][:,0]
    F = np.einsum('nik,k->ni', d.Xs[i], th); e = d.dY - F
    D = o2[i] - min(o2[i], o2[1:].min())
    if D == 0: return True
    ZF = np.einsum('nil,ni->nl', d.Z, F).sum(0); Se = np.einsum('nil,ni->nl', d.Z, e)
    eta = mammen(rng, (N,B)); ZYb = (ZF[:,None] + Se.T@eta)/N
    o, _ = solve_all(ZX, ZYb, est_W2)
    Dst = N*(o[i] - np.minimum(o[i], o[1:].min(0)))
    k = int(np.ceil(.95*(B+1)-1e-9))
    return D <= np.sort(Dst)[k-1]

def set_experiment(N, kappa, R, B, rng):
    pw = cu = hull = cl = 0; disc = 0
    for r in range(R):
        y, q = dgp(N, kappa, rng); d = Data(y, q)
        d.prep(np.r_[.25, d.grid]); est = estimate(d)
        cg = np.unique(np.r_[d.grid, d.gams[est['ih']]])     # CI grid + gamma-hat
        d.prep(np.r_[.25, cg])
        ZX = d.ZXu.mean(1); ZY = d.ZYu.mean(0)[:,None]
        o2 = d.N*solve_all(ZX, ZY, est['W2'])[0][:,0]
        # the estimator's grid minimum (initial grid only) -- same as sample side
        acc0 = m0_accept(d, est['W2'], 0, rng, B, ZX, ZY, o2)
        acc = np.array([m0_accept(d, est['W2'], i, rng, B, ZX, ZY, o2) for i in range(1, len(cg)+1)])
        pw += acc0
        runs = []; j = 0
        while j < len(cg):
            if acc[j]:
                k = j
                while k+1 < len(cg) and acc[k+1]: k += 1
                runs.append((cg[j], cg[k])); j = k+1
            else: j += 1
        disc += len(runs) > 1
        cu += any(lo <= .25 <= hi for lo, hi in runs)
        if runs: hull += runs[0][0] <= .25 <= runs[-1][1]
        # jump-cell closure of each run (support = effective q values)
        qe = np.unique(np.concatenate([q[:,t] for t in d.ts]+[q[:,t-1] for t in d.ts]))
        clo = []
        for lo, hi in runs:
            a = qe[qe <= lo]; b = qe[qe > hi]
            clo.append((a.max() if len(a) else lo, b.min() if len(b) else hi))
        cl += any(lo <= .25 <= hi for lo, hi in clo)
    print(f'SET N={N} kappa={kappa} R={R}: pointwise@g0={pw/R:.3f} union-of-segments={cu/R:.3f} '
          f'hull={hull/R:.3f} union+cell-closure={cl/R:.3f} disconnected={disc/R:.3f}', flush=True)

# ---------------- continuity test ----------------
def Xk(d, gam):
    y, q = d.y, d.q
    def lev(t):
        r = (q[:,t]>gam).astype(float)
        return np.column_stack([y[:,t-1], q[:,t], (q[:,t]-gam)*r])
    return np.stack([lev(t)-lev(t-1) for t in d.ts], axis=1)

def cont_experiment(N, kappa, R, B, rng):
    names = ['W1 kinkres (0.9.35)', 'W1 jumpres', 'W2 kinkres', 'W2 jumpres (0.9.36)']
    rej = np.zeros(4)
    for r in range(R):
        y, q = dgp(N, kappa, rng); d = Data(y, q); d.prep(d.grid)
        est_ZX = d.ZXu.mean(1); ZY = d.ZYu.mean(0)[:,None]
        XK = np.stack([Xk(d, g) for g in d.grid])
        ZXk = np.einsum('nil,gnik->glk', d.Z, XK)/N
        # W2 of the jump fit (stage 1 W1 argmin)
        o1, th1 = solve_all(est_ZX, ZY, d.W1); i1 = int(np.argmin(o1[:,0]))
        e1 = d.dY - np.einsum('nik,k->ni', d.Xs[i1], th1[i1][:,0])
        W2 = omega_inv(np.einsum('nil,ni->nl', d.Z, e1))
        eta = mammen(rng, (N,B)); k95 = int(np.ceil(.95*(B+1)-1e-9))
        for v, (W, res) in enumerate([(d.W1,'k'), (d.W1,'j'), (W2,'k'), (W2,'j')]):
            ok, thk = solve_all(ZXk, ZY, W); oj, thj = solve_all(est_ZX, ZY, W)
            ik, ij = int(np.argmin(ok[:,0])), int(np.argmin(oj[:,0]))
            Tn = N*(ok[ik,0] - oj[ij,0])
            F = np.einsum('nik,k->ni', XK[ik], thk[ik][:,0])
            e = d.dY - (F if res == 'k' else np.einsum('nik,k->ni', d.Xs[ij], thj[ij][:,0]))
            ZF = np.einsum('nil,ni->nl', d.Z, F).sum(0); Se = np.einsum('nil,ni->nl', d.Z, e)
            ZYb = (ZF[:,None] + Se.T@eta)/N
            Tb = N*(solve_all(ZXk, ZYb, W)[0].min(0) - solve_all(est_ZX, ZYb, W)[0].min(0))
            rej[v] += Tn > np.sort(Tb)[k95-1]
    print(f'CONT N={N} kappa={kappa} R={R}: ' + '  '.join(f'{n}={x/R:.3f}' for n, x in zip(names, rej)), flush=True)

if __name__ == '__main__':
    what, N, kappa, R, B = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    rng = np.random.default_rng(777 + N + int(10*kappa))
    (set_experiment if what == 'set' else cont_experiment)(N, kappa, R, B, rng)
