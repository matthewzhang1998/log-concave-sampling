#!/usr/bin/env python3
"""Algebra/clock diagnostics, not a simulation of the imported native pair compiler."""
import hashlib,json,math
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
rng=np.random.default_rng(20261005)
I=np.eye(3)
def sym(M): return (M+M.T)/2
def psqrt(M):
    ev,Q=np.linalg.eigh(sym(M)); assert ev.min()>-1e-12
    return (Q*np.sqrt(np.maximum(ev,0)))@Q.T

# The finite 3-by-3 shield split, including clocks extremely near one.
max_shield_failure=0.; min_normalized_gap=100.
for k in range(2000):
    if k%2:
        q,s,r,t=np.sqrt(1-10**rng.uniform(-9,-.1,4))
    else: q,s,r,t=rng.uniform(.0001,.9999,4)
    cq,cs,cr,ct=np.sqrt(1-np.array([q,s,r,t])**2)
    L=np.array([[cq,0,0,0],[r*s*cq,r*cs,cr,0],[t*s*cq,t*cs,0,ct]])
    K=L@L.T; h=min(cq*cq,cr*cr,ct*ct); sig2=h/10
    gap=np.linalg.eigvalsh(K-sig2*I)[0]
    min_normalized_gap=min(min_normalized_gap,gap/sig2)
    max_shield_failure=max(max_shield_failure,(sig2-gap)/h)
    assert gap>=sig2-2e-15
    # Restoring independent site shields recovers covariance exactly.
    C=psqrt(K-sig2*I)
    assert np.max(np.abs(C@C.T+sig2*I-K))<3e-14

# Noncommuting PSD matrices in the actual leaf Y -> middle V -> root U order.
Hs=[]
for _ in range(3):
    Q,_=np.linalg.qr(rng.normal(size=(3,3)))
    Hs.append((Q*np.array([.1,.45,.9]))@Q.T)
HY,HV,HU=Hs
nu=.3; c=a=math.sqrt(nu/2); A=.07; w=.003
rhoY=rhoV=A; rhoU=-4*w*A/nu
BY,BV,BU=rhoY*HY,rhoV*HV,rhoU*HU
SY,SV,SU=[psqrt(I-B@B.T) for B in [BY,BV,BU]]
# Root as a linear map of independent p, xiY, xiV, xiU.
root=np.concatenate([BU@BV@BY,BU@BV@SY,BU@SV,SU],axis=1)
passive=np.concatenate([I,np.zeros((3,9))],axis=1)
readout=c*root+a*passive
actual=readout@readout.T
intended=nu*I-4*w*A**3*sym(HU@HV@HY)
wrong=nu*I-4*w*A**3*sym(HU@HY@HV)
assert np.max(np.abs(root@root.T-I))<3e-14
assert np.max(np.abs(actual-intended))<3e-14
orientation_gap=float(np.linalg.norm(intended-wrong))
assert orientation_gap>1e-9

# Finite positive endpoint dyadic Gauss rules, terminal midpoint included.
def rule(K,m=4):
    x,w=np.polynomial.legendre.leggauss(m)
    rr=[]; ww=[]
    for k in range(K):
        lo,hi=1-2.**(-k),1-2.**(-k-1)
        rr.extend((lo+hi)/2+(hi-lo)*x/2)
        ww.extend((hi-lo)*w/2)
    h=2.**(-K)
    rr.append(1-h/2); ww.append(h)
    return np.array(rr),np.array(ww)
clock_checks=[]
for K in [1,2,4,8,12]:
    r,beta=rule(K); mass=beta*r; cc=1-r*r
    h=np.minimum(np.minimum(cc[:,None,None],cc[None,:,None]),cc[None,None,:])
    sigma=np.sqrt(h/10)
    weights=mass[:,None,None]*mass[None,:,None]*mass[None,None,:]/2
    bill1=float(np.sum(weights/sigma)); bill2=float(np.sum(weights/sigma**2))
    assert abs(beta.sum()-1)<2e-15
    assert abs(mass.sum()-.5)<2e-15
    assert abs(weights.sum()-1/16)<2e-15
    assert bill1<2
    clock_checks.append(dict(panels=K,nodes=len(r),weight_mass=float(weights.sum()),sum_w_over_sigma=bill1,sum_w_over_sigma2=bill2))
# The exact quadratic target is -H^3/4, regardless of the finite clock rule.
H=Hs[0]*.08
quad=-4*clock_checks[-1]['weight_mass']*(H@H@H)
assert np.max(np.abs(quad+(H@H@H)/4))<1e-16

artifact=ROOT/'FINITE-C0-CUBIC-COVARIANCE-RESERVE.md'
result={
  'status':'PASS algebraic and finite-clock diagnostics',
  'scope':'Does not instantiate the imported fixed-order native pair compiler or prove its numerical implementation.',
  'shield_trials':2000,
  'min_shield_gap_div_sigma2':min_normalized_gap,
  'max_normalized_shield_bound_failure':max_shield_failure,
  'channel_covariance_error':float(np.max(np.abs(actual-intended))),
  'noncommuting_wrong_middle_order_gap':orientation_gap,
  'clock_checks':clock_checks,
  'quadratic_covariance_error':float(np.max(np.abs(quad+(H@H@H)/4))),
  'artifact_sha256':hashlib.sha256(artifact.read_bytes()).hexdigest(),
}
(ROOT/'finite_c0_cubic_covariance_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
