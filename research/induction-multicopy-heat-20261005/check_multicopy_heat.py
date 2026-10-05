#!/usr/bin/env python3
"""Independent scalar/covariance diagnostics, not a replacement for the proofs."""
import itertools, json, math
from functools import lru_cache
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss

OUT=Path(__file__).resolve().parent
rng=np.random.default_rng(20261005)

def prufer_tree(seq,n):
    deg=[1]*n
    for x in seq: deg[x]+=1
    out=[]
    for x in seq:
        v=next(i for i in range(n) if deg[i]==1)
        out.append((v,x)); deg[v]-=1;deg[x]-=1
    rem=[i for i,d in enumerate(deg) if d==1]
    out.append(tuple(rem)); return out

def covariance(edges,t,n):
    adj=[[] for _ in range(n)]
    for (i,j),v in zip(edges,t):adj[i].append((j,v));adj[j].append((i,v))
    K=np.eye(n)
    for src in range(n):
        stack=[(src,-1,1.)]
        while stack:
            i,p,m=stack.pop(); K[src,i]=m
            for j,v in adj[i]:
                if j!=p:stack.append((j,i,min(m,v)))
    delta=np.ones(n)
    for (i,j),v in zip(edges,t):delta[i]=min(delta[i],1-v);delta[j]=min(delta[j],1-v)
    c=min(t);H=K-c*np.ones((n,n))-np.diag(delta)
    return K,c,delta,H

def root_psd(M):
    vals,U=np.linalg.eigh((M+M.T)/2)
    assert vals.min()>-2e-11
    return (U*np.sqrt(np.maximum(0,vals)))@U.T

def clock_history():
    x=np.exp(rng.uniform(np.log(1e-8),np.log(.49),3))
    zq,zs=np.exp(rng.uniform(np.log(1e-8),np.log(.99),2))
    r=np.sqrt(1-2*x);q=np.sqrt(1-zq);s=np.sqrt(1-zs)
    A=np.array([[r[0]*np.sqrt(zq),0,np.sqrt(x[0]),0,0],
                [r[1]*s*np.sqrt(zq),r[1]*np.sqrt(zs),0,np.sqrt(x[1]),0],
                [r[2]*s*np.sqrt(zq),r[2]*np.sqrt(zs),0,0,np.sqrt(x[2])]])
    d=np.array([r[0]*q,r[1]*s*q,r[2]*s*q])
    m=[min(x[1],x[2],zq),min(x[0],x[2],zs),min(x[0],x[1],zs)]
    return A,x,d,np.array(m)

stats={}
minH=1.;max_reconstruction=0.
for _ in range(12000):
    n=int(rng.integers(2,13));edges=prufer_tree(rng.integers(0,n,n-2),n)
    t=1-np.exp(rng.uniform(np.log(1e-10),0,n-1))
    K,c,delta,H=covariance(edges,t,n)
    minH=min(minH,np.linalg.eigvalsh(H).min())
    R=root_psd(H)
    max_reconstruction=max(max_reconstruction,np.max(abs(K-(c+R@R.T+np.diag(delta)))))
stats['n_copy_geometry']={'cases':12000,'n_max':12,'minimum_H_eigenvalue':minH,'maximum_covariance_error':max_reconstruction}
assert minH>-2e-11 and max_reconstruction<2e-11

min_ratio=1e100;row_err=0.;triple_err=0.;mixed_err=0.
for _ in range(6000):
    histories=[clock_history() for _ in range(4)]
    for A,x,d,m in histories:
        G=A@A.T
        for j in range(3):min_ratio=min(min_ratio,np.linalg.eigvalsh(G)[0]/m[j])
        row_err=max(row_err,np.max(abs(d*d+np.diag(G)-(1+(1-2*x))/2)))
    edges=[(0,1),(0,2)];t=1-np.exp(rng.uniform(np.log(1e-9),0,2))
    K,c,delta,H=covariance(edges,t,3)
    # Full triple with different histories on the SAME five-root original bank.
    AA=[v[0] for v in histories[:3]];DD=[];private=[]
    for A,x,d,m in histories[:3]:
        DD.append(np.eye(3)*max(m)/32);private.append(x)
    true=np.block([[K[i,j]*(AA[i]@AA[j].T) for j in range(3)] for i in range(3)])
    reconstructed=np.block([[(c+H[i,j])*(AA[i]@AA[j].T) for j in range(3)] for i in range(3)])
    for i in range(3):
        sl=slice(3*i,3*i+3);G=AA[i]@AA[i].T
        root_psd(G-DD[i])
        true[sl,sl]+=np.diag(private[i]);reconstructed[sl,sl]+=delta[i]*(G-DD[i])+np.diag(private[i])+delta[i]*DD[i]
    triple_err=max(triple_err,np.max(abs(true-reconstructed)))
    # Actual six-packet bank: five common roots, plus two 3-query residual roots.
    s_old=float(np.exp(rng.uniform(np.log(1e-9),0)))
    A6=np.zeros((6,11));oldprivate=[];D6=np.zeros((6,6))
    for ell,(A,x,d,m) in enumerate(histories[2:]):
        G=A@A.T;mold=m[1];sl=slice(3*ell,3*ell+3)
        A6[sl,:5]=np.sqrt(1-s_old)*A
        A6[sl,5+3*ell:8+3*ell]=np.sqrt(s_old)*root_psd(G-np.eye(3)*mold/32)
        oldprivate.extend(x+s_old*mold/32)
        D6[sl,sl]=np.eye(3)*s_old*max(m)/64
    root_psd(A6@A6.T-D6)
    Amix=[np.pad(histories[i][0],((0,0),(0,6))) for i in range(2)]+[A6]
    Dmix=DD[:2]+[D6];pv=[histories[i][1] for i in range(2)]+[np.array(oldprivate)]
    sizes=[3,3,6];starts=np.cumsum([0]+sizes)
    true=np.block([[K[i,j]*(Amix[i]@Amix[j].T) for j in range(3)] for i in range(3)])
    new=np.block([[(c+H[i,j])*(Amix[i]@Amix[j].T) for j in range(3)] for i in range(3)])
    for i in range(3):
        sl=slice(starts[i],starts[i+1]);G=Amix[i]@Amix[i].T
        true[sl,sl]+=np.diag(pv[i]);new[sl,sl]+=delta[i]*(G-Dmix[i])+np.diag(pv[i])+delta[i]*Dmix[i]
    mixed_err=max(mixed_err,np.max(abs(true-new)))
assert min_ratio>=1/16-1e-6 and row_err<1e-12 and triple_err<2e-11 and mixed_err<2e-11
stats['cubic_and_mixed_geometry']={'histories':24000,'minimum_lambda_min_over_m':min_ratio,'full_caller_row_identity_error':row_err,'triple_joint_covariance_error':triple_err,'mixed_joint_covariance_error':mixed_err}

# Exact finite exponent enumeration; fractions are halves, so integer tests suffice.
local_cases=0;diagonal_cases=0
for d in range(4):
    for hits in itertools.product(range(3),repeat=d):
        p=[0,1,0]
        for h in hits:p[h]+=1
        if max(p)>2:
            j=int(np.argmax(p));two_a=p[j]-2
            assert two_a<=2
            assert all(p[v]+two_a<=2 for v in range(3) if v!=j)
        local_cases+=1
for old_a,old_b in itertools.product(range(3),repeat=2):
    for d in range(4):
        for hits in itertools.product(range(6),repeat=d):
            powers=[0,1,0,0,1,0];powers[old_a]+=1;powers[3+old_b]+=1
            for h in hits:powers[h]+=1
            p=[powers[v]+powers[v+3] for v in range(3)]
            assert sum(p)==4+d
            if p[1]>4:
                two_a=p[1]-4
                assert two_a<=3
                assert p[0]+two_a<=3 and p[2]+two_a<=3
                assert two_a<=2 or d==3
            else:
                assert max(p[0],p[2])<=5
                assert sum(max(0,p[v]-4) for v in [0,2])<=1
                if d<=2:assert max(p[0],p[2])<=4
            diagonal_cases+=1
stats['power_census']={'local_d_at_most_3':local_cases,'diagonal_old_hit_extra_hit_cases':diagonal_cases,'max_diagonal_inverse_eta_power':.5}

# Weighted Prüfer counts.
census={}
for Ns in [(3,3,3),(3,3,6),(2,3,4,5),(1,2,2,3,3)]:
    n=len(Ns);literal=0
    for seq in itertools.product(range(n),repeat=n-2):
        edges=prufer_tree(seq,n);deg=[0]*n
        for i,j in edges:deg[i]+=1;deg[j]+=1
        literal+=math.prod(Ns[i]**deg[i] for i in range(n))
    formula=math.prod(Ns)*sum(Ns)**(n-2);assert literal==formula
    census[str(Ns)]={'literal':literal,'formula':formula}
stats['weighted_prufer']=census

# Positive dyadic bridge integrals and literal node maxima.
def dyadic_rule(eta,m=5):
    panels=[];lo=eta
    panels.append((0.,eta))
    while lo<1:
        hi=min(1.,2*lo);panels.append((lo,hi));lo=hi
    nodes=[];weights=[]
    for a,b in panels:
        step=(b-a)/m
        nodes.extend(a+(np.arange(m)+.5)*step);weights.extend([step]*m)
    return np.array(nodes),np.array(weights)
bridge=[]
for J in [4,8,16,24,32]:
    eta=2.**(-J);s,w=dyadic_rule(eta)
    x=s[:,None];y=s[None,:];u=w[:,None]*w[None,:];delta=np.minimum(x,y)
    main=np.sum(u/(delta+eta)**.5)
    central=np.sum(u/(delta+eta))
    leaf=np.sum(u/np.sqrt((delta+eta)*(x+eta)))
    mx_main=np.max(u/(delta+eta)**.5);mx_caller=max(np.max(u/(delta+eta)),np.max(u/np.sqrt((delta+eta)*(x+eta))))
    assert main<8/3+1e-9 and central<=2*np.log((1+eta)/eta)+1e-8 and leaf<=np.log((1+eta)/eta)+2+1e-8
    assert mx_main<=1 and mx_caller<=1
    bridge.append({'J':J,'nodes_per_edge':len(s),'main':main,'central_caller':central,'leaf_caller':leaf,'node_main_max':mx_main,'node_caller_max':mx_caller})
stats['positive_bridge_rule']=bridge

# Polynomial scalar verification of exact three-replica forest/cumulant formula.
def gaussian_moment(counts,K):
    @lru_cache(None)
    def rec(c):
        if sum(c)==0:return 1.
        if sum(c)%2:return 0.
        i=next(k for k,v in enumerate(c) if v)
        rem=list(c);rem[i]-=1;out=0.
        for j,v in enumerate(rem):
            if v:
                rr=rem.copy();rr[j]-=1;out+=v*K[i,j]*rec(tuple(rr))
        return out
    return rec(tuple(counts))
def mom1(p):return 0 if p%2 else math.prod(range(1,p,2))
xx,ww=leggauss(7);xx=(xx+1)/2;ww=ww/2
Ts=[[(0,1),(0,2)],[(0,1),(1,2)],[(0,2),(1,2)]]
max_poly_error=0.;poly_cases=0
for powers in itertools.product(range(1,5),repeat=3):
    target=mom1(sum(powers))-sum(mom1(powers[i]+powers[j])*mom1(powers[k]) for i,j,k in [(0,1,2),(0,2,1),(1,2,0)])+2*math.prod(mom1(p) for p in powers)
    val=0.
    for edges in Ts:
        deg=[0]*3
        for i,j in edges:deg[i]+=1;deg[j]+=1
        if any(p<d for p,d in zip(powers,deg)):continue
        fac=math.prod(math.factorial(p)//math.factorial(p-d) for p,d in zip(powers,deg))
        rem=[p-d for p,d in zip(powers,deg)]
        for order in [(0,1),(1,0)]:
            for a,wa in zip(xx,ww):
                for b,wb in zip(xx,ww):
                    tt=[0.,0.];tt[order[0]]=a;tt[order[1]]=a+(1-a)*b
                    K,_,_,_=covariance(edges,tt,3)
                    val+=wa*wb*(1-a)*fac*gaussian_moment(rem,K)
    max_poly_error=max(max_poly_error,abs(val-target));poly_cases+=1
assert max_poly_error<2e-8
stats['three_replica_polynomial_identity']={'cases':poly_cases,'max_absolute_error':max_poly_error}
# Dyadic majorant comparison and min-path Lipschitz for the sharper rule.
max_comparison_ratio=0.;max_min_path_ratio=0.;max_min_path_absolute_excess=0.
for _ in range(10000):
    n=int(rng.integers(2,13));edges=prufer_tree(rng.integers(0,n,n-2),n)
    deg=np.zeros(n,dtype=int)
    for i,j in edges:deg[i]+=1;deg[j]+=1
    for j in rng.integers(0,n,2):deg[j]+=1
    a=(deg-1)/2;assert abs(a.sum()-n/2)<1e-12
    J=int(rng.integers(3,25));eta=2.**(-J)
    levels=rng.integers(0,J+1,n-1)
    lo=np.where(levels==J,0.,2.**(-levels-1))
    hi=np.where(levels==J,eta,2.**(-levels))
    sA=lo+(hi-lo)*rng.random(n-1);sB=lo+(hi-lo)*rng.random(n-1)
    KA,_,dA,_=covariance(edges,1-sA,n);KB,_,dB,_=covariance(edges,1-sB,n)
    logratio=abs(float(np.dot(a,np.log((dA+eta)/(dB+eta)))))
    ratio=math.exp(logratio)/(2**(n/2))
    max_comparison_ratio=max(max_comparison_ratio,ratio)
    numerator=float(np.max(abs(KA-KB)));denominator=float(np.max(abs(sA-sB)))
    lip=numerator/max(1e-300,denominator)
    max_min_path_ratio=max(max_min_path_ratio,lip)
    max_min_path_absolute_excess=max(max_min_path_absolute_excess,numerator-denominator)
assert max_comparison_ratio<=1+1e-10 and max_min_path_absolute_excess<=1e-14
stats['integrated_majorant_cubature']={'cases':10000,'n_max':12,'maximum_ratio_to_proved_2pow_n_over_2':max_comparison_ratio,'maximum_min_path_Lipschitz_ratio':max_min_path_ratio,'maximum_absolute_Lipschitz_excess':max_min_path_absolute_excess,'absolute_tolerance':1e-14,'note':'Relative excess near 1e-8 gaps is floating subtraction roundoff; the absolute excess is the checked quantity.'}

(OUT/'checks.json').write_text(json.dumps(stats,indent=2)+'\n')
print(json.dumps(stats,indent=2))
