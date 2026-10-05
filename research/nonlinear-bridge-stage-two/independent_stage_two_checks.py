#!/usr/bin/env python3
"""Independent finite checks for the new stage-two mathematical certificate.
These supplement the analytic audit; they are not a compiler or quadrature proof.
"""
from fractions import Fraction as F
from pathlib import Path
import json, math
import numpy as np
rng=np.random.default_rng(41052026)
checks=0

def ck(q):
    global checks
    assert bool(q)
    checks+=1

A=F(1,2); eta=A/2+A*A/4; q=A/2+A*A/2; r=A*A/4
cx=F(9,4)*(1+eta)+A*(F(13,4)+5*A/4)
cn=F(9,4)*q+A*(F(13,4)+3*A/2)
cm=F(9,4)*r+A*(F(5,2)+5*A/4)
cl=A*(F(5,2)+A)
sq=2*(F(3,4)*cx*cx+cn*cn+cm*cm+2*cl*cl)
ck(sq==F(547655,8192));ck(sq<81)
ck(F(169)*F(49,20)+357<25*32)
ck(F(49,20)**2>6);ck(357**2>126874)
ck(8*F(8)**4*F(2,5)**16<F(1,64))
ck(F(1038)*F(2,5)**16<F(1,64))

# A separate implementation of the normalized innovation covariance and its
# complex continuation. Stable expm1 expressions avoid artificial cancellation.
def dd(z): return np.sqrt(-np.expm1(-2*z))
def ff(z): return 2*z*np.exp(-z)/dd(z)
def RR(z,w):
    return np.array([[ff(z),ff(z)*(z-1)],
                     [np.exp(-z)*ff(w),np.exp(-z)*ff(w)*(2*z+w-1)]])
def qq(s): return np.exp(-s)*np.array([1.,2*s,2*s*(s-1)])

small=2.**-24
worst_g=worst_beta=worst_cov=0.
for a,h in np.exp(rng.uniform(-16,5,size=(1200,2))):
    R=RR(a,h); K=np.eye(2)-R@R.T
    ck(np.linalg.eigvalsh(K)[0]>1/16-1e-12)
    L=np.array([[dd(a),0],[np.exp(-h)*dd(a),dd(h)]])
    C=L@K@L.T
    target=np.array([[1-qq(a)@qq(a), np.exp(-h)-qq(a)@qq(a+h)],
                     [np.exp(-h)-qq(a)@qq(a+h),1-qq(a+h)@qq(a+h)]])
    er=np.linalg.norm(C-target);worst_cov=max(worst_cov,er);ck(er<2e-14)
    z=a*(1+small*np.exp(1j*rng.uniform(0,2*np.pi)))
    w=h*(1+small*np.exp(1j*rng.uniform(0,2*np.pi)))
    Dz=dd(z)/dd(a); Dw=dd(w)/dd(h)
    bottom=np.exp(-h)*np.expm1(-(w-h))*dd(z)/dd(h)
    D=np.array([[Dz,0],[bottom,Dw]])
    Rc=RR(z,w);Kc=np.eye(2)-Rc@Rc.T
    ev,O=np.linalg.eigh(K); inv=O@np.diag(ev**-.5)@O.T
    G=inv@D@Kc@D.T@inv
    er=np.linalg.norm(G-np.eye(2),2);worst_g=max(worst_g,er);ck(er<=1/128)
    # The normalized mean map from Y=(x,N,M) to the two innovations.
    xcol=np.array([np.exp(-a)*np.expm1(-(z-a))/dd(a),
                   np.exp(-z-h)*np.expm1(-(w-h))/dd(h)])
    M=np.column_stack([xcol,D@Rc-R]); bm=inv@M
    er=np.linalg.norm(bm,2);worst_beta=max(worst_beta,er);ck(er<=2**-13)

# Independent shifted cubic-jet identity, retaining arbitrary noncommuting
# first derivatives; the second derivative tensor is symmetric in its inputs.
worst_cubic=0.
for dim in [1,2,5]:
    for rep in range(100):
        hv,hw=rng.normal(size=(2,dim))
        Jx,J0=rng.normal(size=(2,dim,dim))
        H=rng.normal(size=(dim,dim,dim));H=(H+H.swapaxes(1,2))/2
        hh=lambda a,b: np.einsum('kij,i,j->k',H,a,b)
        weights=rng.dirichlet(np.ones(5)); hist=rng.normal(size=(5,dim))
        future=rng.normal(size=(5,5,dim));jt=rng.normal(size=(5,dim,dim))
        mean=weights@hist
        inner=sum(weights[j]*jt[j]@(weights@future[j]) for j in range(5))
        first=Jx@J0@hw-hh(hv,hv)/2+hh(hv,mean)
        ancestry=Jx@(inner-J0@hw)
        quadratic=sum(weights[j]*weights[k]*hh(hv-hist[j],hv-hist[k])/2
                      for j in range(5) for k in range(5))
        target=Jx@inner+hh(mean,mean)/2
        er=np.linalg.norm(first+ancestry+quadratic-target)
        worst_cubic=max(worst_cubic,er);ck(er<2e-12)

out={'status':'PASS','assertions':checks,
     'maximum_pair_covariance_identity_error':worst_cov,
     'maximum_complex_standardized_covariance_distance':worst_g,
     'maximum_complex_normalized_mean_map':worst_beta,
     'maximum_noncommuting_cubic_identity_error':worst_cubic,
     'exact_normalized_first_square_upper_bound':str(sq),
     'scope':'Independent identities and scalar constants. Analytic holomorphy, positive quadrature, full bias theorem and guarded compiler import are reviewed in the accompanying audit, not inferred from these samples.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
