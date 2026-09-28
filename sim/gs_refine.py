"""Does cirefine() (xtdpthresh 0.9.36 boundary refinement) close the gap
between union-of-segments coverage and pointwise coverage? M0 wild test."""
import numpy as np, sys
from gs_boot import dgp, Data, solve_all, mammen, estimate

class Tester:
    def __init__(s, d, W2, gridZX, rng, B):
        s.d, s.W2, s.G, s.rng, s.B = d, W2, gridZX, rng, B
        s.ZY = d.ZYu.mean(0)[:,None]
        s.gmin = d.N*solve_all(gridZX, s.ZY, W2)[0][:,0].min()
        s.k = int(np.ceil(.95*(B+1)-1e-9))
    def accept(s, gam):
        d = s.d; N = d.N; X = d.X(gam)
        ZXc = np.einsum('nil,nik->lk', d.Z, X)[None]/N
        oc, th = solve_all(ZXc, s.ZY, s.W2); oc = N*oc[0,0]
        D = oc - min(oc, s.gmin)
        if D == 0: return True
        F = np.einsum('nik,k->ni', X, th[0][:,0]); e = d.dY - F
        ZF = np.einsum('nil,ni->nl', d.Z, F).sum(0); Se = np.einsum('nil,ni->nl', d.Z, e)
        ZYb = (ZF[:,None] + Se.T@mammen(s.rng, (N,s.B)))/N
        ob = N*solve_all(ZXc, ZYb, s.W2)[0][0]
        og = N*solve_all(s.G, ZYb, s.W2)[0].min(0)
        Dst = ob - np.minimum(ob, og)
        return D <= np.sort(Dst)[s.k-1]

def runs_cover(g, a, x):
    o = np.argsort(g); g, a = g[o], a[o]; j = 0
    while j < len(g):
        if a[j]:
            k = j
            while k+1 < len(g) and a[k+1]: k += 1
            if g[j] <= x <= g[k]: return True
            j = k+1
        else: j += 1
    return False

def refine_pts(g, a, supp, m=10):
    o = np.argsort(g); g, a = g[o], a[o]; out = []
    for j in range(len(g)-1):
        if a[j] == a[j+1]: continue
        lbv = supp[supp <= g[j+1]]
        if not len(lbv): continue
        c = supp[(supp > g[j]) & (supp < lbv.max())]
        if len(c) > m: c = c[np.unique([max(0, int(np.ceil(k*len(c)/(m+1)))-1) for k in range(1, m+1)])]
        out += list(c)
    return np.setdiff1d(np.unique(out), g)

if __name__ == '__main__':
    N, kappa, R, B = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    rng = np.random.default_rng(4242 + N + int(10*kappa))
    pw = 0; cov = np.zeros(4); npts = np.zeros(4)
    for r in range(R):
        y, q = dgp(N, kappa, rng); d = Data(y, q); d.prep(np.r_[.25, d.grid]); est = estimate(d)
        t = Tester(d, est['W2'], d.ZXu.mean(1)[1:], rng, B)
        supp = np.unique(np.concatenate([q[:,s] for s in d.ts]+[q[:,s-1] for s in d.ts]))
        pw += t.accept(.25)
        g = np.unique(np.r_[d.grid, d.gams[est['ih']]]); a = np.array([t.accept(x) for x in g])
        for rnd in range(4):
            cov[rnd] += runs_cover(g, a, .25); npts[rnd] += len(g)
            if rnd == 3: break
            new = refine_pts(g, a, supp)
            if len(new):
                g = np.r_[g, new]; a = np.r_[a, [t.accept(x) for x in new]]
    print(f'REFINE N={N} kappa={kappa} R={R}: pointwise={pw/R:.3f}  union coverage by cirefine(0..3)=' +
          ' '.join(f'{c/R:.3f}' for c in cov) + '  mean points=' + ' '.join(f'{p/R:.0f}' for p in npts), flush=True)
