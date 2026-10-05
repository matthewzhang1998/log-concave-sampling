#!/usr/bin/env python3
"""Independent finite diagnostics for the guarded LAW-only variance join.

These test algebra, scalar quadrature, native substitutions and exact Gaussian
channel/row identities. They do not execute LOW30 or certify its guards.
"""
from pathlib import Path
import hashlib, json, math
import numpy as np
from numpy.polynomial.legendre import leggauss

HERE=Path(__file__).resolve().parent
LOW=Path('/workspace/shared/v9-curation-work/frozen/prerequisites/research-source/High Acc Ideas/ai-bucket/30_low_acc.tex')
PIN='7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8'
checks=0
maxima={}
def check(ok, label):
    global checks
    checks+=1
    if not bool(ok): raise AssertionError(label)
def close(a,b,label,tol=2e-11):
    check(np.max(np.abs(np.asarray(a)-np.asarray(b))) <= tol*(1+np.max(np.abs(np.asarray(b)))),label)
def note(name,x): maxima[name]=max(maxima.get(name,0.),float(x))
check(hashlib.sha256(LOW.read_bytes()).hexdigest()==PIN,'pinned LOW30 original')
rng=np.random.default_rng(621904)

# Gaussian disintegration and explicit source-zero carrier alignment.
for A in np.geomspace(1e-9,.08,73):
    w=A**.8; q=1-w; sig2=w*(2-w); sig=math.sqrt(sig2)
    for t in np.r_[np.linspace(0,q,17),q+(1-q)*np.array([.1,.5,.9,.999])]:
        c2=(1-t)*(1+t); c=math.sqrt(c2); d=c2*q*q+sig2
        v=c2*sig2/d; beta=q*c2/math.sqrt(d)
        close(1/v,q*q/sig2+1/c2,'reciprocal variance',tol=3e-10)
        close(beta*beta+v,c2,'conditional variance allocation')
        rot=np.array([[q*c/math.sqrt(d),sig/math.sqrt(d)],[-sig/math.sqrt(d),q*c/math.sqrt(d)]])
        close(rot@rot.T,np.eye(2),'orthogonal carrier pin')
        close(np.array([beta,-math.sqrt(v)])@rot,np.array([c,0.]),'same outer G')
        close(math.sqrt(d)*rot[0],np.array([q*c,sig]),'original Y row restored')
        if t<=q:
            check(v>=w*(1-1e-11) and v<=2*w*(1+1e-11),'bulk variance range')
        # Actual branch-carrier split retains independence before adding branches.
        shares=np.array([.5,.25,.25]); vec=np.sqrt(shares)
        Q,_=np.linalg.qr(np.column_stack([vec,rng.normal(size=(3,2))]))
        if Q[:,0]@vec<0: Q[:,0]*=-1
        close(Q.T@Q,np.eye(3),'independent complete branch carrier split')
        close(vec@Q,np.array([1.,0.,0.]),'sum of branch carriers')
    # Six raw-split terms at a conservative fixed share.
    u=w/8
    terms=np.array([A**4*u**(-1.5), A**4/u, A**4/math.sqrt(u), A**4/math.sqrt(u), A**4.5*u**(-1.5), A**5*u**(-1.5)])
    check(np.max(terms/(A**4*u**(-1.5)))<=1+1e-12,'raw-split law scaling')
    first=A+A**1.5/math.sqrt(u)
    check(first/A<=1+math.sqrt(8)*A**.1+1e-12,'raw-split actual first')
    close((A/math.sqrt(w))**2/math.sqrt(A),A**.7,'self-reserve guard power')
    close(A*w**(-1/3),A**(11/15),'pair guard power')
    # Leading powers after w=A^(4/5).
    close(A*A*w**1.5,A**(16/5),'prefix power')
    close(A**4/w,A**(16/5),'bulk power')
    close(A**3*math.sqrt(w),A**(17/5),'near power')
    close(A**5*w**(-1.5),A**(19/5),'mixture power')

# Finite positive endpoint rule and its integrable near singularity.
# Float64 diagnostic range keeps epsilon=A^3 above endpoint resolution;
# the theorem separately requires adequate scalar precision for smaller A.
def rule(eps):
    J=math.ceil(math.log2(4/eps)); n=max(1,math.ceil(math.log(16/eps,4)))
    x,g=leggauss(n); ts=[]; ws=[]; panels=[]
    for j in range(J):
        lo=1-2.**(-j); hi=1-2.**(-j-1)
        t=(lo+hi)/2+(hi-lo)*x/2; b=(hi-lo)*g/2
        panels.append((t,b,hi-lo)); ts.extend(t); ws.extend(b)
    h=2.**(-J); ts.append(1-h/2);ws.append(h)
    return np.array(ts),np.array(ws),panels,J,n,h
for A in np.geomspace(1e-4,.08,17):
    w=A**.8;q=1-w;s2=w*(2-w)
    t,omega,panels,J,n,h=rule(A**3)
    check(np.all((t>0)&(t<1)) and np.all(omega>0),'strict interior positive nodes')
    close(omega.sum(),1.,'outer mass')
    close(omega@t,.5,'outer first moment')
    check(h<w,'final midpoint in near branch')
    v=(1-t)*(1+t)*s2/(q*q*(1-t)*(1+t)+s2)
    near=float(np.sum(omega[t>q]/np.sqrt(v[t>q])))
    note('near_discrete_over_sqrt_w',near/math.sqrt(w))
    check(near/math.sqrt(w)<5,'discrete near integral numerical envelope')
    for deg in [0,1,2,3,10,100,1000,10000,1000000,1000000000]:
        check(abs(omega@(t**deg)-1/(deg+1))<=A**3+1e-13,'sampled OU Hermite multipliers')
    S0=float(np.sum(np.sqrt(omega*t)))
    S1=float(np.sum(np.sqrt(omega*t)/np.sqrt((1-t)*(1+t))))
    # Direct panel Cauchy-Schwarz proof ingredients.
    for tt,bb,width in panels:
        check(np.sqrt(bb).sum()<=math.sqrt(n*width)*(1+1e-12),'sqrt panel weight sum')
        check(np.sum(np.sqrt(bb)/np.sqrt((1-tt)*(1+tt)))<=math.sqrt(2*n)*(1+1e-12),'shielded sqrt panel sum')
    check(S1<=(J*math.sqrt(2*n)+2)*(1+1e-12),'one-clock shielded sqrt sum')
    # Four-clock sigma-min sum bounded by factorized majorant.
    # Smaller tensor subset evaluates actual fourfold expression directly.
    ix=np.linspace(0,len(t)-1,min(11,len(t)),dtype=int)
    tt=t[ix]; bb=omega[ix]; cc=np.sqrt((1-tt)*(1+tt))
    actual=0.
    for i in range(len(tt)):
      for j in range(len(tt)):
       for k in range(len(tt)):
        for l in range(len(tt)):
         weight=bb[i]*tt[i]*bb[j]*tt[j]*bb[k]*tt[k]*bb[l]*tt[l]
         sigma=min(cc[i],cc[k],cc[l])/math.sqrt(10)
         actual+=math.sqrt(weight)/sigma
    x0=np.sqrt(bb*tt).sum();x1=(np.sqrt(bb*tt)/cc).sum()
    check(actual<=3*math.sqrt(10)*x1*x0**3*(1+1e-12),'four-clock sqrt shield factorization')
    note('one_clock_shield_sum_over_Jsqrt_n',S1/(J*math.sqrt(n)))

# Exact rebalanced ordered three-channel covariance, genuinely noncommuting PSD factors.
for trial in range(180):
    D=3; A=float(10**rng.uniform(-6,-1.2));u=A**.8/8;vp=u/2;W=1/16
    weights=rng.uniform(.01,1,size=7);weights*=W/weights.sum();rho=(4*W*A**3/vp)**(1/3)
    mats=[]
    for _ in range(3):
        M=rng.normal(size=(D,D));H=M@M.T;H/=max(1.,np.linalg.norm(H,2));mats.append(H)
    HY,HV,HU=mats
    check(rho<.5,'sampled pair ideal gap')
    P=(-rho*HU)@(rho*HV)@(rho*HY)
    for weight in weights:
        nu=vp*weight/W
        close(nu*(-rho**3),-4*weight*A**3,'exact rho product calibration')
        cov=nu*np.eye(D)+(nu/2)*(P+P.T)
        target=nu*np.eye(D)-2*weight*A**3*(HU@HV@HY+(HU@HV@HY).T)
        close(cov,target,'ordered noncommuting covariance readout')
        check(np.linalg.eigvalsh(cov)[0]>=nu*(1-rho**3)-1e-14,'positive readout gap')
    prior=math.sqrt(u)*rho**5
    bound=(4*W*2)**(5/3)*A**5*u**(-7/6)
    close(prior,bound,'K prior scaling')
    check(A**5*u**(-7/6)<=A**4*u**(-1.5),'K prior below common law allowance')
    close(math.sqrt(u)*rho,(8*W)**(1/3)*A*u**(1/6),'K residual-first scaling')

# Root-count and final full-source accounting identities.
for _ in range(80):
    D=int(rng.integers(1,10));bulk=int(rng.integers(1,10));near=int(rng.integers(1,10))
    dims=D*rng.integers(2,100,size=bulk)
    aligned=D+int(np.sum(dims-D))+3*D*near
    unaligned=int(dims.sum())+4*D*near
    check(unaligned-aligned==D*(bulk+near-1),'common-carrier Gaussian dimension')

out={
 'status':'PASS', 'assertions':checks,
 'scope':'Finite algebra/geometry/quadrature/native-scale checks. The proof is in the independent audit; no native compiler execution or universal guard admission is claimed.',
 'low30_sha256':PIN,
 'draft_sha256_at_run':hashlib.sha256((HERE.parent/'LAW-ONLY-SHRINKING-BUFFER-JOIN.md').read_bytes()).hexdigest(),
 'maxima':maxima,
}
(HERE/'law_only_join_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
