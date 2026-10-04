import json,math
from pathlib import Path
import numpy as np

T=math.pi/2
def weights(length,n):
    ts=np.arange(n)*length/n
    W=np.zeros((n,n))
    for i in range(n):
        for j in range(i):
            W[i,j]=math.cos(ts[i]-(j+1)*length/n)-math.cos(ts[i]-j*length/n)
    end=np.array([math.cos(length-(j+1)*length/n)-math.cos(length-j*length/n) for j in range(n)])
    return ts,W,end

rows=[]
for a in [.03,.01,.003]:
  h=a**(19/30); ep=a**.9
  for nq,ns in [(2,1),(8,3),(32,11),(128,43)]:
    mq,ms=3,2
    size=2+mq*nq+ms*ns+1
    A=np.zeros((size,size)); C=np.zeros((size,3)); d=np.ones(size)
    # A literal admitted toy old-state block with nonlinear state row O(a).
    A[1,0]=.2*a; C[0,0]=.7; C[1,0]=.3
    oldrow=np.array([.2*a,.3*a])
    tq,Wq,wq=weights(T,nq); ts,Ws,ws=weights(h,ns)
    qgroups=[]; sgroups=[]
    for k in range(mq):
        ids=np.arange(2+k*nq,2+(k+1)*nq); qgroups.append(ids); d[ids]=1/math.sqrt(nq)
        C[ids,0]=np.cos(tq); C[ids,1]=np.sin(tq)
        A[np.ix_(ids,np.arange(2))]=np.cos(tq)[:,None]*oldrow
        if k: A[np.ix_(ids,qgroups[k-1])]=-a*Wq
    off=2+mq*nq
    for k in range(ms):
        ids=np.arange(off+k*ns,off+(k+1)*ns); sgroups.append(ids); d[ids]=1/math.sqrt(ns)
        C[ids,1]=np.cos(ts); C[ids,2]=np.sin(ts)
        A[np.ix_(ids,qgroups[-1])]=-a*np.cos(ts)[:,None]*wq[None,:]
        if k: A[np.ix_(ids,sgroups[k-1])]=-a*Ws
    A[-1,qgroups[-1]]=-a*math.cos(h)*wq
    A[-1,sgroups[-1]]=-a*ws
    C[-1]=[0,math.cos(h),math.sin(h)]
    Aw=d[:,None]*A/d[None,:]; Cw=d[:,None]*C
    M=np.abs(Aw)
    row_schur=np.max((M@d)/d); col_schur=np.max((M.T@d)/d)
    bound=math.sqrt(row_schur*col_schur)
    c_norm=math.sqrt(float(np.linalg.eigvalsh(Cw.T@Cw)[-1]))
    l_norm=float(np.linalg.norm(Cw[:-1,2]))
    terminal=float(np.linalg.norm(Aw[-1]))
    assert bound/a<7
    assert c_norm<=math.sqrt(2+mq+ms+1)
    assert l_norm<=math.sqrt(ms)*math.sin(h)+1e-14
    assert terminal/a<2
    # Weighted source realizes exactly the unweighted original-gradient DAG.
    record=np.array([.3,-.7,.4]); p=np.zeros(size); pw=np.zeros(size)
    for i in range(size):
        z=C[i]@record+A[i]@p
        p[i]=.6*z+.2*math.sin(math.sqrt(a)*z)/math.sqrt(a)
        zw=(Cw[i]@record+Aw[i]@pw)/d[i]
        pw[i]=d[i]*(.6*zw+.2*math.sin(math.sqrt(a)*zw)/math.sqrt(a))
    assert np.max(np.abs(pw/d-p))<2e-13
    # Project every nonterminal L row; the terminal stays intact.
    Cp=Cw.copy(); Cp[:-1,2]=0
    Ctw=np.zeros((2*size,4)); Ctw[:size,:3]=Cp; Ctw[size:,:3]=Cp; Ctw[size+size-1,3]=ep
    gamma=math.sin(h)**2; beta=(1/gamma+1/(ep*ep))**-.5
    B=np.array([0.,0.,beta/math.sin(h),-beta/ep])
    target=np.zeros(2*size); target[size-1]=beta
    assert abs(B@B-1)<1e-13
    assert np.max(np.abs(B@Ctw.T-target))<1e-13
    rows.append({'a':a,'quarter_nodes':nq,'short_nodes':ns,'total_original_nodes':size,
                 'weighted_C_norm':c_norm,'unweighted_C_norm':math.sqrt(float(np.linalg.eigvalsh(C.T@C)[-1])),
                 'absolute_edge_schur_bound_over_a':bound/a,
                 'weighted_nonterminal_L_norm_over_h':l_norm/h,
                 'terminal_row_norm_over_a':terminal/a,
                 'same_value_max_error':float(np.max(np.abs(pw/d-p))),
                 'coisometry_error':abs(B@B-1)})

out={'status':'PASS: weighted rows, edge/terminal bounds, identical VALUE map, and projected twin row',
     'scope':'Finite scalar block fixtures; tensor-I lifting preserves bounds. General old-block admission uses its supplied operator/readout contract.',
     'fixtures':rows}
Path(__file__).with_name('weighted_source_graph_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
