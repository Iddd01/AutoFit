import numpy as np, sys
def dgp(N, kappa, rng, T=6, burn=20, gam=.25):
    TT=T+burn
    es=rng.standard_normal((N,TT)); w=rng.standard_normal((N,TT))
    u=np.empty((N,TT)); u[:,0]=w[:,0]; u[:,1:]=.5*es[:,:-1]+np.sqrt(.75)*w[:,1:]
    q=np.empty((N,TT)); q[:,0]=rng.standard_normal(N)
    for t in range(1,TT): q[:,t]=.7*q[:,t-1]+u[:,t]
    g=(kappa-.5+2*q)*(q>gam)
    y=np.empty((N,TT)); y[:,0]=rng.standard_normal(N)
    for t in range(1,TT): y[:,t]=.6*y[:,t-1]+q[:,t]+g[:,t]+.5*es[:,t]
    return y[:,burn:], q[:,burn:]
def build(y,q,ylags,qlags):
    # FD equations t=3..6 (0-based 2..5); instrument blocks per t
    N,T=y.shape; ts=range(2,T)
    blocks=[]
    for t in ts:
        yl=[y[:,s] for s in range(0,t-1)][-ylags:] if ylags else [y[:,s] for s in range(0,t-1)]
        ql=[q[:,s] for s in range(0,t)][-qlags:] if qlags else [q[:,s] for s in range(0,t)]
        blocks.append(np.column_stack(yl+ql))
    L=sum(b.shape[1] for b in blocks); m=len(blocks)
    Z=np.zeros((N,m,L)); c=0
    for j,b in enumerate(blocks):
        Z[:,j,c:c+b.shape[1]]=b; c+=b.shape[1]
    dY=np.column_stack([y[:,t]-y[:,t-1] for t in ts])
    H=2*np.eye(m)-np.eye(m,k=1)-np.eye(m,k=-1)
    ZHZ=np.einsum('nil,ij,njk->lk',Z,H,Z)/N
    W1=np.linalg.pinv(ZHZ)
    return Z,dY,W1,ts
def design(y,q,gam,ts,kink):
    N,T=y.shape
    def lev(t):
        r=(q[:,t]>gam).astype(float)
        if kink: return np.column_stack([y[:,t-1],q[:,t],(q[:,t]-gam)*r])
        return np.column_stack([y[:,t-1],q[:,t],r,q[:,t]*r,y[:,t-1]*r])
    return np.stack([lev(t)-lev(t-1) for t in ts],axis=1)  # N,m,k
def fit(Z,dY,X,W):
    N=Z.shape[0]
    ZX=np.einsum('nil,nik->lk',Z,X)/N; ZY=np.einsum('nil,ni->l',Z,dY)/N
    A=ZX.T@W@ZX; b=np.linalg.solve(A,ZX.T@W@ZY)
    g=ZY-ZX@b; return N*g@W@g, b
def run(N,kappa,R,ylags,qlags,seed=1):
    rng=np.random.default_rng(seed); gh=[]; T_c=[]
    for r in range(R):
        y,q=dgp(N,kappa,rng)
        Z,dY,W1,ts=build(y,q,ylags,qlags)
        qe=q[:,1:].ravel(); grid=np.quantile(qe,np.linspace(.05,.95,46))
        objs=[fit(Z,dY,design(y,q,g,ts,False),W1)[0] for g in grid]
        i1=int(np.argmin(objs)); X1=design(y,q,grid[i1],ts,False)
        _,b1=fit(Z,dY,X1,W1)
        e=dY-np.einsum('nik,k->ni',X1,b1)
        m=np.einsum('nil,ni->nl',Z,e); m=m-m.mean(0)
        W2=np.linalg.pinv(m.T@m/N)
        o2=[fit(Z,dY,design(y,q,g,ts,False),W2)[0] for g in grid]
        gh.append(grid[int(np.argmin(o2))])
        ok=[fit(Z,dY,design(y,q,g,ts,True),W2)[0] for g in grid]
        T_c.append(min(ok)-min(o2))
    gh=np.array(gh)
    return np.sqrt(np.mean((gh-.25)**2)), gh.std(), np.mean(abs(gh-.25)<.05), np.median(T_c)
if __name__=='__main__':
    R=int(sys.argv[1])
    for yl,ql,lab in [(None,None,'GS full ladder (24 IV)'),(2,3,'short ladder')]:
        for kappa in (0,1):
            for N in (400,1600):
                rm,sd,w05,tc=run(N,kappa,R,yl,ql)
                print(f'{lab:24s} kappa={kappa} N={N:5d} rmse={rm:.3f} sd={sd:.3f} P(|err|<.05)={w05:.2f} median T_cont(W2)={tc:.1f}',flush=True)
