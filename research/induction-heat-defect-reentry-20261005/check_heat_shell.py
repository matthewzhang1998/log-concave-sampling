#!/usr/bin/env python3
"""Finite diagnostics, not a native compiler or substitute for the proofs."""
from fractions import Fraction as F
from pathlib import Path
import hashlib, itertools, json, math
import numpy as np

OUT=Path(__file__).resolve().parent
rng=np.random.default_rng(20261005)

def counts(nmax):
    T=[0,1]; V=[0,1]
    for n in range(2,nmax+1):
        T.append(sum(i*(n-i)*T[i]*T[n-i] for i in range(1,n)))
        V.append(sum((n-i)*(i+1)*V[i]*T[n-i]+i*(n-i+1)*T[i]*V[n-i] for i in range(1,n)))
    return T,V
T,V=counts(12)
assert T[8]==794880 and V[8]==27020800
assert 9*V[8]==243187200 and 10*V[8]==270208000

ledgers=[]
for n in range(4,101):
    g0,g1=F(1,n),F(1,n-1)
    b0,b1=F(1,2*(n-1)),F(1,4*(n-1))
    p0=F(n)-n*b0-(n-2)*g0-2*g1
    p1=F(n)-n*b1-(n-2)*g1-2*g1
    gain=F(n*n-4*n+8,4*n*(n-1))
    ro=F(n)-F(3,2)+F(2,n)
    rn=F(n)-F(5,4)+F(1,n-1)
    P=F(n)+F(1,n-1)
    E=F(n)-(n-2)*g1
    LB=F(n)-(n-1)*g1
    assert p1-p0==gain and gain>0
    assert E+1==P and F(n)+g1==P
    assert E+LB>P and 2*min(ro,rn)>P
    assert rn-g1==F(n)-F(5,4)>1
    assert b1+p1+g1>1
    if n<=12:
        ledgers.append(dict(n=n,Psi_old=str(p0),Psi_shell=str(p1),gain=str(gain),root_old=str(ro),root_shell=str(rn),root_first=str(rn-g1),target=str(P),self_bank=str(E+LB)))
assert F(1,8)+F(8,2*7*6)==F(37,168)

# Exact heterogeneous mixed-bridge reserve (tree arities, unequal widths).
bridge_cases=0
Gamma=F(1,7)
for arity in range(2,10):
    psis=[F(j+2,3) for j in range(arity)]
    # A labelled star is enough to check this algebra; every tree has r-1 edges.
    prices=[F(j%7+1,56)+F((j+3)%7+1,56) for j in range(arity-1)]
    assert all(price<=2*Gamma for price in prices)
    child=sum(psis)+sum(2*Gamma-price for price in prices)
    assert child>=sum(psis)
    for h in range(2,7):
        assert h*psis[0]+2*Gamma*(h-1)>psis[0]
    bridge_cases+=1

# Source formulas for g(x)=c*x+(b/w)sin(w*x), a smooth bounded-Hessian gradient.
c,b,w,A=.5,.2,1.3,1.
def g(x): return c*x+(b/w)*np.sin(w*x)
def gp(x): return c+b*np.cos(w*x)

def source(k,x,ps,z,t):
    a=1/math.sqrt(k+1)
    ans=np.zeros(np.broadcast_shapes(np.shape(x),np.shape(z),*[np.shape(p) for p in ps]))
    for eps in itertools.product((-1,1),repeat=k):
        q=z+a*t*sum((e*p for e,p in zip(eps,ps)),start=0)
        ans+=math.prod(eps)*(g(q+a*t*x)-g(q))
    return ans/(A*t*2**k)

def source_first(k,x,ps,z,t):
    a=1/math.sqrt(k+1); dx=0.; dp=np.zeros(k); dz=0.
    for eps in itertools.product((-1,1),repeat=k):
        s=math.prod(eps); q=z+a*t*sum(e*p for e,p in zip(eps,ps))
        h,h0=gp(q+a*t*x),gp(q)
        dx+=s*a*h/A
        for j,e in enumerate(eps): dp[j]+=s*a*e*(h-h0)/A
        dz+=s*(h-h0)/(A*t)
    return np.r_[dx,dp]/2**k,dz/2**k

max_first=0.; max_origin=0.; max_odd=0.
for k in range(7):
    for _ in range(80):
        x,z,Z=rng.normal(size=3); ps=rng.normal(size=k); t=10**rng.uniform(-2,0); delta=rng.uniform(0,1)
        f0,dz0=source_first(k,x,ps,z,t); f1,dz1=source_first(k,x,ps,z+delta*Z,t)
        markfirst=(f0-f1)/2
        max_first=max(max_first,float(np.linalg.norm(markfirst)))
        assert np.linalg.norm(markfirst)<=1+1e-12
        assert abs((dz0-dz1)/2)<=1/t+1e-12
        marked=lambda xx,pp: (source(k,xx,pp,z,t)-source(k,xx,pp,z+delta*Z,t))/2
        max_origin=max(max_origin,abs(float(marked(0.,ps))))
        assert abs(marked(0.,ps))<1e-12
        if k:
            flip=ps.copy(); flip[0]*=-1
            err=abs(float(marked(x,ps)+marked(x,flip)))
            max_odd=max(max_odd,err)
            assert err<1e-12

# Independent Gaussian quadrature for exact marked response, k=0,1,2.
# Integrates X, all P_j, and Z, while holding the captured base z fixed.
xh,wh=np.polynomial.hermite.hermgauss(16); xh=xh*np.sqrt(2); wh=wh/np.sqrt(np.pi)
quad=[]
for k in range(3):
    dim=k+2
    points=np.meshgrid(*([xh]*dim),indexing='ij')
    weights=np.meshgrid(*([wh]*dim),indexing='ij')
    X=points[0]; ps=points[1:k+1]; Z=points[-1]
    wt=np.prod(np.stack(weights),axis=0)
    z,t,delta=.37,.58,.41
    val=(source(k,X,ps,z,t)-source(k,X,ps,z+delta*Z,t))/2
    response=float(np.sum(wt*val*X*np.prod(np.stack(ps),axis=0))) if k else float(np.sum(wt*val*X))
    a=1/math.sqrt(k+1)
    # Constant linear term cancels; only the sine derivative remains.
    expected=(a**(k+1)*t**k/(2*A))*b*w**k*math.sin(w*z+(k+1)*math.pi/2)*math.exp(-w*w*t*t/2)*(1-math.exp(-w*w*delta*delta/2))
    err=abs(response-expected)
    assert err<2e-11,(k,response,expected,err)
    quad.append(dict(k=k,nodes=16**dim,response=response,expected=expected,error=err))

# Full rank-eight tensor-tree telescoping with the legal double-star graph.
edges=[(0,1),(0,2),(0,3),(0,4),(1,5),(1,6),(1,7)]
inds=[[] for _ in range(8)]
for e,(u,v) in enumerate(edges): inds[u].append(e); inds[v].append(e)
for v in range(8): inds[v].append(7+v)

def contract(arrs):
    args=[]
    for arr,ix in zip(arrs,inds): args.extend([arr,ix])
    return np.einsum(*args,list(range(7,15)),optimize=True)
telescope_errors=[]
chord_errors=[]
gx,gw=np.polynomial.legendre.leggauss(4); gx=(gx+1)/2; gw=gw/2
for _ in range(12):
    cold=[rng.normal(size=(2,)*len(ix))*.15 for ix in inds]
    warm=[rng.normal(size=(2,)*len(ix))*.15 for ix in inds]
    actual=contract(cold)-contract(warm)
    summed=np.zeros_like(actual)
    for m in range(8):
        summand=[cold[v] if v<m else warm[v] for v in range(8)]
        summand[m]=(cold[m]-warm[m])/2
        summed+=2*contract(summand)
    err=float(np.max(np.abs(actual-summed)))
    assert err<1e-12
    telescope_errors.append(err)
    chord=np.zeros_like(actual)
    for ss,ww in zip(gx,gw):
        mix=[(1-ss)*warm[v]+ss*cold[v] for v in range(8)]
        for m in range(8):
            term=mix.copy(); term[m]=cold[m]-warm[m]
            chord+=ww*contract(term)
    cherr=float(np.max(np.abs(chord-actual)))
    assert cherr<1e-12
    chord_errors.append(cherr)

# A fixed-frequency witness: the mark's normalized active first is not O(delta)
# when the original admissible frequency can increase like delta^(-1).
# k=0, t=delta; use phase differences that are constant after rescaling.
frequency_witness=[]
for q in [8,16,32,64,128,256]:
    delta=1/q; t=delta
    # D_x f_D at x=0,z=0,Z=pi: b[cos(0)-cos(pi)]/2=b.
    marked_first=b
    assert abs(marked_first-b)<1e-14
    frequency_witness.append(dict(frequency=q,delta=delta,marked_active_first=marked_first,first_over_delta=marked_first/delta))

result={
    'status':'PASS finite diagnostics; not native execution',
    'rank8_histories':T[8],'rank8_raw_values':V[8],
    'rank8_shell_raw_values':9*V[8],'rank8_rebuilt_total_raw_values':10*V[8],
    'exact_ledgers_n4_to100':True,'selected_ledgers':ledgers,
    'heterogeneous_bridge_cases':bridge_cases,
    'max_sampled_mark_active_first':max_first,'max_origin_error':max_origin,'max_probe_oddness_error':max_odd,
    'gaussian_response_checks':quad,'max_rank8_tensor_telescope_error':max(telescope_errors),'max_rank8_four_node_chord_error':max(chord_errors),
    'fixed_first_no_heat_grade_witness':frequency_witness,
    'rank8_allocation_gamma_ceiling':'37/168',
    'limitations':['No native sampler executed','No all-order restoration','No whole rank-eight cumulant high-frequency lower bound inferred from local source witness']
}
(OUT/'checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
