"""Kink-model threshold test at gamma0 (Study B supplement geometry, kappa=0):
the model imposes continuity, x = [dy_{-1}, dq, d((q-g)1(q>g))].
Scores: acceptance at gamma0 (citest), union of accepted grid runs, hull.
0.9.35/0.9.36 wild scheme: restricted fit + restricted residuals * unit Mammen
weights, fixed regressors/instruments and fixed sample W2."""
import numpy as np, sys
from gs_boot import dgp, solve_all, omega_inv, mammen

def build(y, q, ylag, qlag):
    N, T = y.shape; ts = list(range(2, T)); m = len(ts)
    blocks = []
    for t in ts:
        yl = [y[:, s] for s in range(0, t-1)]; ql = [q[:, s] for s in range(0, t)]
        if ylag: yl = yl[-ylag:]
        if qlag: ql = ql[-qlag:]
        blocks.append(np.column_stack(yl + ql + [np.ones(N)]))
    L = sum(b.shape[1] for b in blocks)
    Z = np.zeros((N, m, L)); c = 0
    for j, b in enumerate(blocks): Z[:, j, c:c+b.shape[1]] = b; c += b.shape[1]
    dY = np.column_stack([y[:, t]-y[:, t-1] for t in ts])
    H = 2*np.eye(m)-np.eye(m, k=1)-np.eye(m, k=-1)
    W1 = np.linalg.pinv(np.einsum('nil,ij,njk->lk', Z, H, Z)/N)
    return Z, dY, W1, ts

def Xk(y, q, g, ts):
    def lev(t):
        r = (q[:, t] > g).astype(float)
        return np.column_stack([y[:, t-1], q[:, t], (q[:, t]-g)*r])
    return np.stack([lev(t)-lev(t-1) for t in ts], axis=1)

def one(N, rng, B, ylag, qlag, trim, ngrid, nci):
    y, q = dgp(N, 0.0, rng)
    Z, dY, W1, ts = build(y, q, ylag, qlag)
    qe = np.concatenate([q[:, t] for t in ts] + [q[:, t-1] for t in ts])
    lo, hi = np.quantile(qe, [trim/2, 1-trim/2])
    grid = np.linspace(lo, hi, ngrid); cig = np.linspace(lo, hi, nci)
    gams = np.r_[.25, grid, cig]
    Xs = np.stack([Xk(y, q, g, ts) for g in gams])
    ZX = np.einsum('nil,gnik->glk', Z, Xs)/N; ZY = (np.einsum('nil,ni->l', Z, dY)/N)[:, None]
    G = slice(1, 1+ngrid)
    o1, th1 = solve_all(ZX[G], ZY, W1); i1 = 1+int(np.argmin(o1[:, 0]))
    e1 = dY - np.einsum('nik,k->ni', Xs[i1], th1[i1-1][:, 0])
    W2 = omega_inv(np.einsum('nil,ni->nl', Z, e1))
    o2, th2 = solve_all(ZX, ZY, W2); o2 = N*o2[:, 0]
    gmin = o2[G].min(); ih = 1+int(np.argmin(o2[G]))
    k = int(np.ceil(.95*(B+1)-1e-9))
    def accept(i):
        D = o2[i] - min(o2[i], gmin)
        if D == 0: return True
        F = np.einsum('nik,k->ni', Xs[i], th2[i][:, 0]); e = dY - F
        ZF = np.einsum('nil,ni->nl', Z, F).sum(0); Se = np.einsum('nil,ni->nl', Z, e)
        ZYb = (ZF[:, None] + Se.T@mammen(rng, (N, B)))/N
        ob = N*solve_all(ZX[i:i+1], ZYb, W2)[0][0]
        og = N*solve_all(ZX[G], ZYb, W2)[0].min(0)
        return D <= np.sort(ob - np.minimum(ob, og))[k-1]
    a0 = accept(0)
    idx = list(range(1+ngrid, 1+ngrid+nci)); cg = gams[idx]
    acc = np.array([accept(i) for i in idx])
    # gamma-hat is appended to the CI grid (D = 0 there)
    cg = np.r_[cg, gams[ih]]; acc = np.r_[acc, True]; o = np.argsort(cg); cg, acc = cg[o], acc[o]
    runs = []; j = 0
    while j < len(cg):
        if acc[j]:
            kk = j
            while kk+1 < len(cg) and acc[kk+1]: kk += 1
            runs.append((cg[j], cg[kk])); j = kk+1
        else: j += 1
    union = any(a <= .25 <= b for a, b in runs)
    hull = bool(runs) and runs[0][0] <= .25 <= runs[-1][1]
    return a0, union, hull, len(runs) > 1, (runs[-1][1]-runs[0][0]) if runs else np.nan, gams[ih]

if __name__ == '__main__':
    N, R, B, lab = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    cfg = {'supp': (2, 3, .15, 199, 100), 'gs': (None, None, .10, 46, 46),
           'gs100': (None, None, .10, 46, 100), 'gs200': (None, None, .10, 46, 200)}[lab]
    rng = np.random.default_rng(99 + N + len(lab))
    res = np.array([one(N, rng, B, *cfg) for _ in range(R)], dtype=float)
    pw, un, hu, dc, ln, gh = res.mean(0)[0], res[:, 1].mean(), res[:, 2].mean(), res[:, 3].mean(), np.nanmean(res[:, 4]), res[:, 5]
    print(f'KINK {lab} N={N} R={R} B={B}: at-gamma0={pw:.3f} union={un:.3f} hull={hu:.3f} '
          f'disconnected={dc:.3f} hull_len={ln:.2f} rmse_ghat={np.sqrt(np.mean((gh-.25)**2)):.3f}', flush=True)
