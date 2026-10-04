#!/usr/bin/env python3
"""Finite diagnostics; the companion report contains the proof and scope."""
import hashlib, json, math
from fractions import Fraction as F
from pathlib import Path
import numpy as np

OUT=Path(__file__).parent
ROOT=OUT.parents[2]

def coeff(A, lam, T, N, M):
    t=np.arange(N+1)*T/N
    W=np.zeros((N+1,N+1))
    for i in range(1,N+1):
        q=np.arange(i)
        W[i,q]=np.cos(t[i]-t[q+1])-np.cos(t[i]-t[q])
    b=np.column_stack((np.cos(t),np.sin(t)))
    u=b.copy()
    for _ in range(M):
        u=b-A*lam*(W@u)
    return u[-1],float(np.max(np.sum(W,axis=1)))

def packet(A,lam,j):
    R=j+.5
    Nq=max(1,math.ceil(A**(-(j-1))))
    Ns=max(1,math.ceil(A**(-max(0,R-3.4))))
    Mq=max(1,j-1)
    Ms=max(1,math.ceil((R-.5)/(34/15))-1)
    q,_=coeff(A,lam,math.pi/2,Nq,Mq)
    h=A**(19/30)
    s,_=coeff(A,lam,h,Ns,Ms)
    return float(s[0]*q[0]),float(s[0]*q[1]),float(s[1]),Nq,Ns,Mq,Ms

def level_vector(A,lam,j,seedscale):
    # Complete source coefficients on (seed,Z1,L1,...,Zj,Lj), divided by sqrt(A).
    if j==0:return np.array([seedscale])
    c,b,d,*_=packet(A,lam,j)
    v=np.zeros(1+2*j);v[0]=c**j*seedscale
    for i in range(1,j+1):v[2*i-1:2*i+1]=c**(j-i)*np.array([b,d])
    return v

fixtures=[]
for A in [.3,.2,.1,.05]:
    for lam in [.2,.6,.9]:
        for j in [1,2,3]:
            sig=1/math.sqrt(1+A*lam)
            f=level_vector(A,lam,j,sig)
            c=level_vector(A,lam,j-1,sig)
            cp=np.zeros_like(f);cp[0]=c[0]
            if j>1:cp[3:]=c[1:]
            gap=float(np.linalg.norm(f-cp))
            law=math.sqrt(A)*abs(float(np.linalg.norm(f))-sig)
            h=A**(19/30);Kh=1-math.cos(h)
            contraction=A/((1-A)*(1-A*Kh))
            c_j,b_j,d_j,*bill=packet(A,lam,j)
            assert abs(c_j)<=contraction+1e-12
            # Mean/zero telescope on a seed with known nonzero displacement.
            mu=.37
            vals=[(1 if k==0 else packet(A,lam,k)[0]**k)*mu for k in range(j+1)]
            inc=[vals[k]-vals[k-1] for k in range(1,j+1)]
            assert abs(vals[0]+sum(inc)-vals[-1])<1e-13
            fixtures.append(dict(A=A,lam=lam,j=j,strong_ratio=gap/A**(j-1),law_ratio=law/A**(j+.5),map_contraction=abs(c_j),contraction_bound=contraction,counts=bill))
assert max(x['strong_ratio'] for x in fixtures)<4
assert max(x['law_ratio'] for x in fixtures)<3

# Coarse numerical floor cannot be omitted. Identical reset maps perturbed by
# eta_j at the LAST common step leave exactly |eta_j-eta_(j-1)|.
coarse_floor=[]
for A in [.2,.1,.05]:
    for j in [2,3,4]:
        fine=A**(j+.5);coarse=A**(j-.5)
        gap=abs(coarse-fine)
        coarse_floor.append(dict(A=A,j=j,gap=gap,ratio_to_coarse=gap/coarse))
        assert abs(gap/coarse-(1-A))<1e-12

# Exact continuous constrained resource allocation, followed by legal rounding.
def allocation(a,c,delta):
    a=np.asarray(a,float);c=np.asarray(c,float)
    if delta>=sum(a):return np.ones_like(a)
    lo=0.;hi=1.
    def N(lam):return np.maximum(1,np.sqrt(lam*a/c))
    while np.sum(a/N(hi))>delta:hi*=4
    for _ in range(180):
        mid=(lo+hi)/2
        if np.sum(a/N(mid))>delta:lo=mid
        else:hi=mid
    return N(hi)
alloc=[]
for A in [.2,.1,.05]:
    for J in [2,3,5]:
        a=np.array([A*A]+[A**(j+1) for j in range(1,J+1)])
        c=np.array([1.]+[A**(-(j-1)) for j in range(1,J+1)])
        delta=A**J
        n=allocation(a,c,delta); nr=np.ceil(n)
        B=float(np.sum(np.sqrt(a*c)))
        easy=np.maximum(1,np.ceil((B/delta)*np.sqrt(a/c)))
        bill=float(c@easy);lower=max(float(sum(c)),B*B/delta)
        assert float(a@(1/nr))<=delta*(1+1e-10)
        assert float(a@(1/easy))<=delta*(1+1e-10)
        assert bill<=B*B/delta+sum(c)+1e-9
        assert bill<=2*lower+1e-9
        alloc.append(dict(A=A,J=J,delta=delta,continuous=n.tolist(),rounded=nr.tolist(),easy=easy.tolist(),law_bound=float(a@(1/nr)),easy_bill=bill,lower_bound=lower))

# Shared-keep Gaussian fixture: marginal buffered price vs actual source gap.
square_root_gap=[]
for A in [.2,.1,.05]:
    for N in [1,10,100,1000]:
        e=A/math.sqrt(N)
        marginal=math.sqrt(1+e*e)-1
        assert marginal<=e*e/2+1e-15
        square_root_gap.append(dict(A=A,N=N,actual_shared_keep_L2=e,marginal_W2=marginal))

exponents={
 'protected_v':str(F(19,15)), 'short_local_base':str(F(17,5)),
 'minimal_grid_offset':str(F(1,2)), 'minimal_grid_j_cost_offset':str(F(-1)),
 'minimal_pair_force_energy_j_offset':str(F(0)),
 'buffer_product_j_plus_cost_exponent':str(F(2)),
 'seven_offset_buffer_product_exponent':str(F(-9,2)),
 'known_force_physical_bias_J_offset':str(F(3,2)),
}
assert F(3,2)+F(3,2)*F(19,15)==F(17,5)
for j in range(1,11):
    jj=F(j)
    assert F(1)+jj+F(1,2)-F(3,2)==jj
    assert jj+F(1,2)-F(3,2)==jj-1
    assert (jj+1)-(jj-1)==2
    assert (jj+1)-(jj+7-F(3,2))==F(-9,2)
    assert F(1)+jj+F(1,2)==jj+F(3,2)

# The exact nested empirical coupling has a 1/sqrt(N) strong scale.
for N,M in [(1,2),(2,7),(13,101)]:
    v=np.full(M,-1/M,dtype=float);v[:N]+=1/N
    assert abs(float(v@v)-(1/N-1/M))<1e-13

# Bounded Hessians do not turn a small VALUE displacement into a small first.
for k in [10,100,1000,10000]:
    lam=.6;kap=.2;eps=math.pi/k
    value_gap=lam*eps+kap*math.sin(k*eps)/k
    derivative_gap=kap*(math.cos(k*eps)-1)
    assert abs(value_gap-lam*eps)<1e-13
    assert abs(derivative_gap+2*kap)<1e-13

out=dict(status='PASS',scope='Finite diagnostics, not complete hidden-family admission',num_quadratic_fixtures=len(fixtures),max_strong_ratio=max(x['strong_ratio'] for x in fixtures),max_law_ratio=max(x['law_ratio'] for x in fixtures),exponents=exponents,quadratic_fixtures=fixtures,coarse_grid_floor=coarse_floor,allocations=alloc,shared_keep_counterexample=square_root_gap)
(OUT/'strong_order_telescope_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ['quadratic_fixtures','coarse_grid_floor','allocations','shared_keep_counterexample']},indent=2))
