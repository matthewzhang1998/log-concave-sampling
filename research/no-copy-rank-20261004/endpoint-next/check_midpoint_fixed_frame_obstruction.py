#!/usr/bin/env python3
"""Linear-algebra diagnostics for the midpoint universal fixed-readout test."""
from pathlib import Path
import math,json
import numpy as np
OUT=Path(__file__).resolve().parent
rng=np.random.default_rng(77019)
a=1/16;h=.2;eps=.3;c=.6;s=.8
results=[]
for d in [2,3,4]:
    n=4*d
    rows=[]
    def H():
        A=rng.normal(size=(d,d));A=(A+A.T)/2
        return .5*np.eye(d)+.035*A/max(np.linalg.norm(A,2),1e-12)
    for _ in range(600):
        A0,A1,F,Fp,Tp,Tm=[H() for i in range(6)]
        DtS=np.eye(d)+(a*c/2)*(A0+A1)+(a*a/2)*A1@(F-Fp)
        DtU=(a*s/2)*(A0+A1)
        DtZ=-(a*a*eps/2)*A1@Fp
        J=np.concatenate([(Tp-Tm)@DtS,(Tp-Tm)@DtU,(Tp-Tm)@DtZ,h*(Tp+Tm)],axis=1)
        for i in range(d):
            for j in range(i+1,d):
                r=np.zeros((d,n));r[j,:]+=J[i,:];r[i,:]-=J[j,:]
                rows.append(r.ravel())
    L=np.array(rows)
    _,sv,vh=np.linalg.svd(L,full_matrices=False)
    expected=np.concatenate([np.zeros((d,3*d)),np.eye(d)],axis=1).ravel()
    expected/=np.linalg.norm(expected)
    null=vh[-1,:]
    alignment=abs(float(null@expected))
    assert alignment>1-1e-8,(d,alignment)
    assert sv[-1]<1e-11 and sv[-2]>1e-9,(d,sv[-2:])
    results.append(dict(d=d,equations=len(rows),unknowns=d*n,numerical_nullity=int(sum(sv<1e-11)),smallest_nonzero_singular=float(sv[-2]),null_singular=float(sv[-1]),alignment_with_V_scalar=alignment))

# One generic same-potential base record with all six jet sites distinct.
S=np.array([1.,2.]);U=np.array([3.,-1.]);Z=np.array([2.,1.]);V=np.array([-1.,2.]);m=.5
x0=c*S+s*U;x1=x0-a*m*eps*Z
tbar=S+a*m/2*(x0+x1)
qs=np.array([S,S+eps*Z,x0,x1,tbar+h*V,tbar-h*V])
separations=[np.linalg.norm(qs[i]-qs[j]) for i in range(6) for j in range(i)]
assert min(separations)>0
assert min(np.linalg.norm(qs,axis=1))>0

# Beta!=0 completion lower bound inherits 1/|beta|; beta=0 has nonzero N_VV.
kappa=.012;kk=4096;eta=1/80
betas=[-1.,-.7,-.1,.1,.7,1.]
beta_checks=[]
for beta in betas:
    lower=kappa*a*s*s*eta*kk/(1920*abs(beta)*h)
    beta_checks.append(dict(beta=beta,full_first_lower=lower))
# At the terminal bump center chi''(0)=-2 and zeta=h/480.
gpp_plus=(h/480)*(4/h)**2*(-2)
NVV=kappa*h/2*gpp_plus
assert abs(NVV+kappa/30)<1e-16
out=dict(symmetry_systems=results,six_distinct_jet_sites=qs.tolist(),minimum_jet_separation=min(separations),beta_checks=beta_checks,beta_zero_N_VV=NVV,
         note='Diagnostics only. Independent local-jet realizability and all-frame proof are analytic in MIDPOINT-FIXED-KNOWN-READOUT-NO-GO.md.')
(OUT/'midpoint_fixed_frame_obstruction_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
