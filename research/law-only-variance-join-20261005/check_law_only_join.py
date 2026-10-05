#!/usr/bin/env python3
"""Independent-of-native-execution algebra and finite-clock diagnostics.
This does not simulate the imported LOW30 mean/pair compiler.
"""
import hashlib, json, math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss

checks=0
results={}
def ok(x,msg):
    global checks
    assert bool(x),msg
    checks+=1
rng=np.random.default_rng(5120506)
# Exact Gaussian disintegration using a direct block covariance calculation.
max_geom=0
for _ in range(1000):
    q=rng.uniform(.5,.99999); t=rng.uniform(0,.999999)
    s=1-q*q*t*t;c2=1-t*t;sig2=1-q*q
    v=c2*sig2/s
    joint=np.array([[c2,q*c2],[q*c2,s]])
    schur=joint[0,0]-joint[0,1]**2/joint[1,1]
    err=abs(v-schur);max_geom=max(max_geom,err)
    ok(err<1e-13,'conditional variance')
    ok(abs(1/v-(q*q/sig2+1/c2))<2e-7,'inverse variance')
    beta=q*c2/math.sqrt(s)
    ok(abs(beta*beta+v-c2)<1e-13,'carrier orthogonal rotation')
    O=np.array([[beta/math.sqrt(c2),math.sqrt(v/c2)],[-math.sqrt(v/c2),beta/math.sqrt(c2)]])
    ok(np.linalg.norm(O@O.T-np.eye(2))<2e-12,'known rotation orthogonal')
    if t<=q:
        ok(1-q-1e-14<=v<=2*(1-q)+1e-14,'bulk buffer bounds')
results['max_geometry_error']=max_geom
# Rebalanced complete ideal chronological path in noncommuting dimensions.
orientation=[]
for _ in range(150):
    D=3
    mats=[]
    for __ in range(3):
        B=rng.normal(size=(D,D));H=B@B.T;H/=np.linalg.norm(H,2)
        mats.append(H)
    HU,HV,HY=mats
    A=10**rng.uniform(-5,-1.3);u=A**.8;vp=u/2
    W=1/16;wj=rng.uniform(.0001,W);nu=vp*wj/W
    rho=(4*W*A**3/vp)**(1/3)
    Rs=[rho*HY,rho*HV,-rho*HU]
    # Each channel preserves a standard input covariance exactly.
    cov=np.eye(D);cross=np.eye(D)
    for R in Rs:
        innov=np.eye(D)-R@R.T
        ok(np.linalg.eigvalsh(innov).min()>0,'positive channel gap')
        cov=R@cov@R.T+innov;cross=R@cross
    actual=nu/2*(cov+np.eye(D)+cross+cross.T)
    M=HU@HV@HY
    target=nu*np.eye(D)-2*wj*A**3*(M+M.T)
    ok(np.linalg.norm(actual-target)<1e-13,'exact covariance calibration')
    ok(abs(nu*(-rho)*rho*rho+4*wj*A**3)<1e-16,'cubic amplitude balance')
    orientation.append(float(np.linalg.norm(M-HY@HV@HU)))
ok(max(orientation)>.01,'noncommuting orientation negative control')
results['max_orientation_reversal_difference']=max(orientation)
# All six raw-split variable-buffer exponents, before and after the terminal A.
x=F(4,5)
mean_exponents=[4-F(3,2)*x,4-x,4-x/2,4-x/2,F(9,2)-F(3,2)*x,5-F(3,2)*x]
ok(min(mean_exponents)==F(14,5),'variable-buffer LAW leading grade')
rows={'prefix':2+F(3,2)*x,'prefix_correction':3+x,'near':3+x/2,'bulk_gaussianization':4-x,'compiler_and_gram_mixture':5-F(3,2)*x,'balanced_pair_b5':6-F(7,6)*x,'mixed_covariance_mixture':7-F(3,2)*x,'covariance_target_restore':5-x/2,'final_mean':F(4)}
ok(min(rows.values())==F(16,5),'overall grade')
ok(1-x/2==F(3,5),'mean source radius power')
ok(F(3,2)-x==F(7,10),'self-reserve radius power')
ok(1-x/3==F(11,15),'pair radius power')
results['error_exponents']={k:str(v) for k,v in rows.items()}
# Positive endpoint-midpoint and dyadic-panel rule in correlation gap.
def clock(J,n,split=None):
    cuts=[0.0]+[2.0**(-j) for j in range(J, -1,-1)]
    if split and split not in cuts:cuts=sorted(cuts+[split])
    xs=[]; ws=[];g,gw=leggauss(n)
    for a,b in zip(cuts[:-1],cuts[1:]):
        if a==0:
            xx=np.array([(a+b)/2]);ww=np.array([b-a])
        else:
            xx=(a+b)/2+(b-a)/2*g;ww=(b-a)/2*gw
        xs.extend(1-xx);ws.extend(ww)
    return np.array(xs),np.array(ws)
near=[]
for w in [.2,.05,.01,.003,.0001]:
    ts,ww=clock(36,8,w);q=1-w;sig2=1-q*q
    vv=(1-ts*ts)*sig2/(1-q*q*ts*ts)
    ok(np.all(ts>0)&np.all(ts<1)&np.all(ww>0),'strict finite endpoint nodes')
    ok(abs(ww.sum()-1)<1e-13,'mass one')
    mask=ts>q
    debt=float(np.sum(ww[mask]/np.sqrt(vv[mask])))
    near.append({'w':w,'ratio_to_sqrt_w':debt/math.sqrt(w)})
    ok(debt<=4*math.sqrt(w),'finite near endpoint sum')
    ok(np.min(vv[~mask])>=w-1e-13,'all finite bulk buffers')
    for k in [0,1,2,7,31,128,1000,1000000]:
        moment=float(np.sum(ww*ts**k))
        ok(abs(moment-1/(k+1))<1e-7,'tested Hermite multiplier')
results['finite_near_branch']=near
# Factorized sqrt-weight shield bound, which avoids a N^4 materialization.
shield=[]
for J in [8,16,24,32]:
    ts,ww=clock(J,5);c=np.sqrt(1-ts*ts);N=len(ts)
    S=float(np.sum(np.sqrt(ww*ts)))
    R=float(np.sum(np.sqrt(ww*ts)/c))
    upper=3*math.sqrt(10)*R*S**3
    ok(R<=2*N,'one-clock sqrt-weight shield bound')
    ok(S<=math.sqrt(N),'one-clock sqrt weight count')
    ok(upper<=6*math.sqrt(10)*N**2.5,'factorized four-clock public-log bound')
    shield.append({'J':J,'nodes':N,'one_clock_shield_sum':R,'four_clock_bound':upper})
results['sqrt_weight_shields']=shield
results['checks']=checks
results['status']='PASS: algebra, noncommuting ideal channels, powers, finite singular sums; native compilers not executed'
root=Path(__file__).parent
(root/'check-results.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
