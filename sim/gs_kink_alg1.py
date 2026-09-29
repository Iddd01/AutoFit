"""Kink-model threshold test at gamma0 = .25 (kappa = 0, continuity imposed):
does Gong-Seo Algorithm 1 fix the mild undercoverage of the command's wild
scheme seen in the formal kink supplement (citest coverage 91-94%)?
Same samples for every scheme (supplement geometry: FD, T=6, y lags 2-3,
q lags 1-3, block constants, 199-point uniform grid, trim .15):
  M0  wild, restricted residuals, sample W2 fixed        (= xtdpthresh 0.9.36)
  M1  wild, restricted residuals, W2 re-estimated per draw
  M2  unit resampling, recentred moments, W2 re-estimated (Gong-Seo Alg. 1)
Usage: python3 gs_kink_alg1.py N R B SEED"""
import numpy as np, sys, time
from gs_boot import dgp, Data, estimate, boot_wild, boot_unit, reject


class KinkData(Data):
    def __init__(s, y, q, ylag=2, qlag=3, trim=.15, ngrid=199):
        N, T = y.shape; ts = list(range(2, T)); m = len(ts)
        blocks = []
        for t in ts:
            yl = [y[:, u] for u in range(0, t-1)][-ylag:]
            ql = [q[:, u] for u in range(0, t)][-qlag:]
            blocks.append(np.column_stack(yl + ql + [np.ones(N)]))
        L = sum(b.shape[1] for b in blocks)
        Z = np.zeros((N, m, L)); c = 0
        for j, b in enumerate(blocks):
            Z[:, j, c:c+b.shape[1]] = b; c += b.shape[1]
        s.N, s.m, s.L, s.ts, s.y, s.q, s.Z = N, m, L, ts, y, q, Z
        s.dY = np.column_stack([y[:, t]-y[:, t-1] for t in ts])
        H = 2*np.eye(m)-np.eye(m, k=1)-np.eye(m, k=-1)
        s.W1 = np.linalg.pinv(np.einsum('nil,ij,njk->lk', Z, H, Z)/N)
        qe = np.concatenate([q[:, t] for t in ts]+[q[:, t-1] for t in ts])
        lo, hi = np.quantile(qe, [trim/2, 1-trim/2])
        s.grid = np.linspace(lo, hi, ngrid)
        s.ZYu = np.einsum('nil,ni->nl', Z, s.dY)

    def X(s, gam):
        y, q = s.y, s.q
        def lev(t):
            r = (q[:, t] > gam).astype(float)
            return np.column_stack([y[:, t-1], q[:, t], (q[:, t]-gam)*r])
        return np.stack([lev(t)-lev(t-1) for t in s.ts], axis=1)


METHODS = {'M0 wild, W2 fixed (command)': lambda d, e, r, B: boot_wild(d, e, r, B),
           'M1 wild, W2 re-estimated':    lambda d, e, r, B: boot_wild(d, e, r, B, reW2=True),
           'M2 unit, Alg. 1':             lambda d, e, r, B: boot_unit(d, e, r, B)}

if __name__ == '__main__':
    N, R, B, seed = (int(a) for a in sys.argv[1:5])
    rng = np.random.default_rng(seed)
    rej = {k: [] for k in METHODS}; tm = {k: 0. for k in METHODS}
    for r in range(R):
        y, q = dgp(N, 0.0, rng); d = KinkData(y, q); d.prep(np.r_[.25, d.grid])
        est = estimate(d)
        for k, f in METHODS.items():
            t0 = time.time(); Dst = f(d, est, rng, B); tm[k] += time.time()-t0
            rej[k].append(bool(reject(est['D'], Dst)))
        if (r+1) % 25 == 0:
            print(f'rep {r+1}: ' + '  '.join(f'{k.split()[0]}={1-np.mean(v):.3f}' for k, v in rej.items()), flush=True)
    print(f'KINK-ALG1 N={N} R={R} B={B} seed={seed}')
    for k, v in rej.items():
        p = np.mean(v)
        print(f'  {k:30s} coverage={1-p:.3f} (mcse {np.sqrt(p*(1-p)/len(v)):.3f})  '
              f'rejections={int(np.sum(v))}/{len(v)}  {tm[k]/len(v):.2f}s/rep', flush=True)
