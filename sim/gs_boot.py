"""Prototype: grid-bootstrap test of H0: gamma = gamma0 on the Gong-Seo DGP
(B6 geometry: FD, T=6, 24 lag instruments + 4 block constants, 46-point
p5-p95 quantile grid). Compares bootstrap schemes for the threshold CI:
  M0  wild, restricted residuals, W2 fixed           (= xtdpthresh 0.9.35)
  M0u wild, unrestricted residuals, W2 fixed
  M1  wild, restricted residuals, W2 re-estimated per draw (stage 1 -> Omega* -> W2*)
  M2  unit resampling, recentred moments, W2 re-estimated (Alg.1-oriented,
      = xtdpthresh boottype(unit))
Coverage of the inverted set at gamma0 = 1 - rejection rate here."""
import numpy as np, sys, time
MAM_P = (np.sqrt(5)+1)/(2*np.sqrt(5)); PHI = (np.sqrt(5)-1)/2
def mammen(rng, shape):
    d = rng.random(shape) < MAM_P
    return np.where(d, -PHI, 1/PHI)

def dgp(N, kappa, rng, T=6, burn=20, gam=.25):
    TT = T+burn
    es = rng.standard_normal((N,TT)); w = rng.standard_normal((N,TT))
    u = np.empty((N,TT)); u[:,0] = w[:,0]; u[:,1:] = .5*es[:,:-1]+np.sqrt(.75)*w[:,1:]
    q = np.empty((N,TT)); q[:,0] = rng.standard_normal(N)
    for t in range(1,TT): q[:,t] = .7*q[:,t-1]+u[:,t]
    g = (kappa-.5+2*q)*(q>gam)
    y = np.empty((N,TT)); y[:,0] = rng.standard_normal(N)
    for t in range(1,TT): y[:,t] = .6*y[:,t-1]+q[:,t]+g[:,t]+.5*es[:,t]
    return y[:,burn:], q[:,burn:]

class Data:
    def __init__(s, y, q, grid_n=46):
        N,T = y.shape; ts = list(range(2,T)); m = len(ts)
        blocks = [np.column_stack([y[:,:t-1], q[:,:t], np.ones(N)]) for t in ts]
        L = sum(b.shape[1] for b in blocks)
        Z = np.zeros((N,m,L)); c = 0
        for j,b in enumerate(blocks): Z[:,j,c:c+b.shape[1]] = b; c += b.shape[1]
        s.N, s.m, s.L, s.ts, s.y, s.q, s.Z = N, m, L, ts, y, q, Z
        s.dY = np.column_stack([y[:,t]-y[:,t-1] for t in ts])
        H = 2*np.eye(m)-np.eye(m,k=1)-np.eye(m,k=-1)
        s.W1 = np.linalg.inv(np.einsum('nil,ij,njk->lk',Z,H,Z)/N)
        qe = np.concatenate([q[:,t] for t in ts]+[q[:,t-1] for t in ts])
        s.grid = np.unique(np.quantile(qe, np.linspace(.05,.95,grid_n)))
        s.ZYu = np.einsum('nil,ni->nl', Z, s.dY)            # per-unit Z'dY
    def X(s, gam):
        y,q = s.y, s.q
        def lev(t):
            r = (q[:,t]>gam).astype(float)
            return np.column_stack([y[:,t-1],q[:,t],r,q[:,t]*r,y[:,t-1]*r])
        return np.stack([lev(t)-lev(t-1) for t in s.ts], axis=1)
    def prep(s, gams):
        s.gams = np.asarray(gams)
        s.Xs = np.stack([s.X(g) for g in s.gams])            # G,N,m,k
        s.ZXu = np.einsum('nil,gnik->gnlk', s.Z, s.Xs)       # G,N,L,k per unit

def solve_all(ZX, ZY, W):
    """ZX: G,L,k (means); ZY: L,B (means); W: L,L. Returns obj (G,B), theta (G,k,B)."""
    A = np.einsum('glk,lj,gjh->gkh', ZX, W, ZX)
    rhs = np.einsum('glk,lj,jb->gkb', ZX, W, ZY)
    th = np.linalg.solve(A, rhs)
    g = ZY[None] - np.einsum('glk,gkb->glb', ZX, th)
    obj = np.einsum('glb,lj,gjb->gb', g, W, g)
    return obj, th

def omega_inv(Gu):
    """Gu: N,L per-unit moments -> centred cluster Omega^{-1}."""
    N = Gu.shape[0]; Gc = Gu - Gu.mean(0)
    return np.linalg.pinv(Gc.T@Gc/N)

def estimate(d):
    """Two-step fixed-grid estimator. Index 0 of d.gams is gamma0, 1.. is grid."""
    N = d.N; ZX = d.ZXu.mean(1); ZY = d.ZYu.mean(0)[:,None]
    o1, th1 = solve_all(ZX[1:], ZY, d.W1)
    i1 = 1+int(np.argmin(o1[:,0]))
    e1 = d.dY - np.einsum('nik,k->ni', d.Xs[i1], th1[i1-1][:,0])
    W2 = omega_inv(np.einsum('nil,ni->nl', d.Z, e1))
    o2, th2 = solve_all(ZX, ZY, W2)
    o2 = N*o2[:,0]
    D = o2[0] - min(o2[0], o2[1:].min())
    ih = 1+int(np.argmin(o2[1:]))
    return dict(W2=W2, D=D, th0=th2[0][:,0], thh=th2[ih][:,0], ih=ih, o2=o2)

def boot_wild(d, est, rng, B, resid='restricted', reW2=False):
    N = d.N; ZX = d.ZXu.mean(1)
    F = np.einsum('nik,k->ni', d.Xs[0], est['th0'])
    th_e = est['th0'] if resid=='restricted' else est['thh']
    ie = 0 if resid=='restricted' else est['ih']
    e = d.dY - np.einsum('nik,k->ni', d.Xs[ie], th_e)
    ZF = np.einsum('nil,ni->nl', d.Z, F)                  # per unit
    Se = np.einsum('nil,ni->nl', d.Z, e)
    eta = mammen(rng, (N,B))
    ZYb = (ZF.sum(0)[:,None] + Se.T@eta)/N               # L,B
    if not reW2:
        o, _ = solve_all(ZX, ZYb, est['W2'])
        return N*(o[0] - np.minimum(o[0], o[1:].min(0)))
    o1, th1 = solve_all(ZX[1:], ZYb, d.W1)
    i1 = 1+np.argmin(o1, axis=0)
    Dst = np.empty(B)
    for b in range(B):
        Gu = ZF + Se*eta[:,b:b+1] - d.ZXu[i1[b]]@th1[i1[b]-1][:,b]
        W2b = omega_inv(Gu)
        o, _ = solve_all(ZX, ZYb[:,b:b+1], W2b)
        Dst[b] = N*(o[0,0] - min(o[0,0], o[1:,0].min()))
    return Dst

def boot_unit(d, est, rng, B):
    """Alg.1-oriented unit resampling (mirrors xtdpthresh boottype(unit))."""
    N = d.N
    ehat = d.dY - np.einsum('nik,k->ni', d.Xs[est['ih']], est['thh'])
    Y0 = np.einsum('nik,k->ni', d.Xs[0], est['th0']) + ehat
    ZY0u = np.einsum('nil,ni->nl', d.Z, Y0)
    Shat = np.einsum('nil,ni->nl', d.Z, ehat).sum(0)
    Dst = np.empty(B)
    for b in range(B):
        w = np.bincount(rng.integers(0,N,N), minlength=N).astype(float)
        nb = w.sum()
        ZX = np.einsum('n,gnlk->glk', w, d.ZXu)
        ZY = (w@ZY0u - Shat)[:,None]
        o1, th1 = solve_all(ZX[1:], ZY, d.W1)
        i1 = 1+int(np.argmin(o1[:,0]))
        Gu = ZY0u - d.ZXu[i1]@th1[i1-1][:,0] - Shat/N
        Gc = Gu - (w@Gu)/nb
        Om = (Gc*w[:,None]).T@Gc/nb
        W2b = np.linalg.pinv(Om)
        o, _ = solve_all(ZX, ZY, W2b)
        o = o[:,0]/nb
        Dst[b] = o[0] - min(o[0], o[1:].min())
    return Dst

def reject(D, Dst, alpha=.05):
    k = int(np.ceil((1-alpha)*(len(Dst)+1)-1e-9))
    return D > np.sort(Dst)[k-1]

METHODS = {'M0 wild fixW2 (0.9.35)': lambda d,e,r,B: boot_wild(d,e,r,B),
           'M0u wild unres.resid':    lambda d,e,r,B: boot_wild(d,e,r,B,'unrestricted'),
           'M1 wild re-est W2':       lambda d,e,r,B: boot_wild(d,e,r,B,reW2=True),
           'M2 unit Alg1-oriented':   lambda d,e,r,B: boot_unit(d,e,r,B)}

if __name__ == '__main__':
    N, kappa, R, B = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    only = sys.argv[5].split(',') if len(sys.argv) > 5 else None
    rng = np.random.default_rng(20260928 + N + int(10*kappa))
    rej = {k: [] for k in METHODS}; tm = {k: 0. for k in METHODS}
    for r in range(R):
        y, q = dgp(N, kappa, rng); d = Data(y, q); d.prep(np.r_[.25, d.grid])
        est = estimate(d)
        for k, f in METHODS.items():
            if only and k.split()[0] not in only: continue
            t0 = time.time(); Dst = f(d, est, rng, B); tm[k] += time.time()-t0
            rej[k].append(reject(est['D'], Dst))
    print(f'N={N} kappa={kappa} R={R} B={B}')
    for k in METHODS:
        if not rej[k]: continue
        p = np.mean(rej[k]); se = np.sqrt(p*(1-p)/len(rej[k]))
        print(f'  {k:26s} reject={p:.3f} (mcse {se:.3f})  coverage={1-p:.3f}  {tm[k]/len(rej[k]):.2f}s/rep')
