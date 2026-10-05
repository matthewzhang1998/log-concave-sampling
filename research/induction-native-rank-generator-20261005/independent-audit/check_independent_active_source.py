"""Independent diagnostics: Fourier factorization, tensor cuts, full Jacobians,
radial geometry, exact grade/cost arithmetic, and audited input stability.
No author checker or its symbolic response routine is imported.
"""
from __future__ import annotations
import hashlib, importlib.util, itertools, json, math
from collections import Counter
from fractions import Fraction
from pathlib import Path
import numpy as np
from numpy.polynomial.hermite import hermgauss
HERE=Path(__file__).resolve().parent
BASE=HERE.parent
spec=importlib.util.spec_from_file_location('audited_active_source',BASE/'active_probe_source.py')
src=importlib.util.module_from_spec(spec); spec.loader.exec_module(src)
rng=np.random.default_rng(20261005)
counts=Counter(); extrema={}
def ck(ok,group,msg):
    counts[group]+=1
    if not bool(ok): raise AssertionError(group+': '+msg)
def update(name,x): extrema[name]=max(extrema.get(name,0.),float(x))
def tensorpower(u,n):
    out=np.array(1.)
    for _ in range(n): out=np.multiply.outer(out,u)
    return out

def gradient(q,dirs,freqs,bs):
    return .5*np.asarray(q)+sum(b/w*np.sin(w*np.dot(u,q))*u for u,w,b in zip(dirs,freqs,bs))
def exact_fourier_source(x,P,z,t,dirs,freqs,bs):
    k=len(P);a=1/math.sqrt(k+1)
    out=.5*a*x if k==0 else np.zeros_like(x)
    for u,w,b in zip(dirs,freqs,bs):
        v=np.exp(1j*w*np.dot(u,z))*np.expm1(1j*w*a*t*np.dot(u,x))
        for row in P: v*=1j*np.sin(w*a*t*np.dot(u,row))
        out+=b/w/t*v.imag*u
    return out

# A positive-Hessian admissible, non-coordinate-separable source.
for D in (1,3,11):
    dirs=rng.normal(size=(3,D));dirs/=np.linalg.norm(dirs,axis=1)[:,None]
    freqs=np.array([.7,1.6,3.1]);bs=np.array([.12,.09,.07])
    g=lambda q:gradient(q,dirs,freqs,bs)
    for k in range(10):
      for t in (.01,.2,1.7):
        x=rng.normal(size=D);P=rng.normal(size=(k,D));z=rng.normal(size=D)
        for R in (None,.7):
            y,J=(x,np.eye(D)) if R is None else src.radial_pullback(x,R)
            expected=J.T@exact_fourier_source(y,P,z,t,dirs,freqs,bs)
            got,ncalls=src.active_source(g,x,P,z,t,1.,R)
            err=np.linalg.norm(got-expected)
            update('fourier_source_absolute_error',err)
            ck(err<3e-12/t,'fourier','literal VALUE sum versus independent sine-product factorization')
            ck(ncalls==2**(k+1),'call_count','uncached 2^(k+1) VALUE count')
            zero,_=src.active_source(g,np.zeros(D),P,z,t,1.,R)
            ck(np.array_equal(zero,np.zeros(D)),'origin','same-center anchor')
            ck(np.linalg.norm(got)<=np.linalg.norm(x)/math.sqrt(k+1)+1e-11,'energy','pointwise one-energy envelope')
            if k:
              for j in (0,k-1):
                Pneg=P.copy();Pneg[j]*=-1
                other,_=src.active_source(g,x,Pneg,z,t,1.,R)
                ck(np.linalg.norm(got+other)<3e-12/t,'parity','odd in each tested actual probe')

# Width differentiation is certified from original VALUE geometry, roots frozen.
for k in range(7):
  D=3;dirs=rng.normal(size=(3,D));dirs/=np.linalg.norm(dirs,axis=1)[:,None]
  freqs=np.array([.7,1.6,3.1]);bs=np.array([.12,.09,.07]);g=lambda q:gradient(q,dirs,freqs,bs)
  for t in (.05,.4,1.2):
    x=rng.normal(size=D)*3;P=rng.normal(size=(k,D));z=rng.normal(size=D)
    for R in (None,1.):
      h=t*1e-4;fp,_=src.active_source(g,x,P,z,t+h,1.,R);fm,_=src.active_source(g,x,P,z,t-h,1.,R)
      y=x if R is None else src.radial_pullback(x,R)[0]
      bound=(2*np.linalg.norm(y)+np.linalg.norm(P))/math.sqrt(k+1)/t
      ck(np.linalg.norm((fp-fm)/(2*h))<=bound+1e-6,'width_first','actual frozen-root width derivative')

# Separately integrate Fourier response factors with Gaussian quadrature.
# This is deterministic integration of factorized independent banks, not MC.
raw,ww=hermgauss(96);ghx=raw*math.sqrt(2);ghw=ww/math.sqrt(math.pi)
for k in range(11):
  for t in (.15,.7,1.3):
    for w in (.4,1.7,3.2):
      a=1/math.sqrt(k+1);z=.37;c=w*a*t
      private=np.sum(ghw*ghx*np.expm1(1j*c*ghx))
      probe=np.sum(ghw*ghx*1j*np.sin(c*ghx))
      integrated=(np.exp(1j*w*z)*private*probe**k).imag/(w*t)
      derivative=a**(k+1)*t**k*w**k*math.sin(w*z+(k+1)*math.pi/2)*math.exp(-w*w*t*t/2)
      err=abs(integrated-derivative)
      update('gaussian_response_absolute_error',err)
      ck(err<3e-13,'gaussian_response','k+1 derivatives and total heat t^2')

# Origin zero and exposed-probe conditional means are different obligations.
z=.7;t=.9
k0_mean=.2*math.sin(z)*(math.exp(-t*t/2)-1)/t
k0_quadrature=sum(float(w)*src.active_source(lambda q:.5*q+.2*np.sin(q),np.array([x]),np.empty((0,1)),np.array([z]),t,1.)[0][0] for x,w in zip(ghx,ghw))
ck(abs(k0_mean-k0_quadrature)<1e-13 and abs(k0_mean)>.01,'conditional_means','k=0 private Gaussian mean is generally nonzero')
k=1;a=1/math.sqrt(2);P=np.array([[.8]])
exposed_mean=.2/t*(math.exp(-(a*t)**2/2)-1)*math.cos(z)*math.sin(a*t*P[0,0])
exposed_quad=sum(float(w)*src.active_source(lambda q:.5*q+.2*np.sin(q),np.array([x]),P,np.array([z]),t,1.)[0][0] for x,w in zip(ghx,ghw))
ck(abs(exposed_quad-exposed_mean)<1e-13 and abs(exposed_mean)>1e-3,'conditional_means','exposed active probe does not retain zero private mean')

# Actual tensor norms, all proper flattenings in non-coordinate directions.
for D in (2,3):
  dirs=rng.normal(size=(3,D));dirs/=np.linalg.norm(dirs,axis=1)[:,None]
  freqs=np.array([.7,1.6,3.1]);bs=np.array([.12,.09,.07])
  for k in range(6):
    a=1/math.sqrt(k+1);t=.63;z=rng.normal(size=D)
    T=.5*a*np.eye(D) if k==0 else np.zeros((D,)*(k+2))
    for u,w,b in zip(dirs,freqs,bs):
      amplitude=a**(k+1)*t**k*b*w**k*math.sin(w*np.dot(u,z)+(k+1)*math.pi/2)*math.exp(-w*w*t*t/2)
      T+=amplitude*tensorpower(u,k+2)
    ck(np.linalg.norm(T.ravel())<=a*math.sqrt(D)+1e-13,'tensor_energy','one vector energy, no extra dimension')
    for p in range(1,k+2):
      q=k+2-p
      independent=a**(k+1)*2**(k/2)*math.sqrt(math.factorial(p-1)*math.factorial(q-1))
      ck(abs(src.projection_cut_constant(k,p)-independent)<1e-14,'cut_constants','formula agreement')
      op=np.linalg.norm(T.reshape(D**p,D**q),2)
      ck(op<=independent+1e-13,'tensor_cuts','all proper tensor flattenings')
    ck(max(src.projection_cut_constant(k,p) for p in range(1,k+2)) <= a**(k+1)*2**(k/2)*math.sqrt(math.factorial(k))*(1+1e-14),'cut_constants','endpoint factorial maximizes every cut')

# Independent complete Jacobians, including caller after actual row substitution.
def jacobian(fn,v,h=2e-5):
    E=np.eye(len(v))*h
    return np.stack([(fn(v+e)-fn(v-e))/(2*h) for e in E],axis=1)
D=3
dirs=rng.normal(size=(3,D));dirs/=np.linalg.norm(dirs,axis=1)[:,None]
freqs=np.array([.7,1.6,3.1]);bs=np.array([.12,.09,.07]);g=lambda q:gradient(q,dirs,freqs,bs)
for k in (0,1,2,4):
  for radius in (.4,1.,1.35,1.999,2.,2.5):
    R=1.;x=np.array([radius,0.,0.]);P=rng.normal(size=(k,D));z=rng.normal(size=D);t=.37
    allv=np.r_[x,P.ravel(),z]
    for bounded in (False,True):
      RR=R if bounded else None
      fun=lambda v:src.active_source(g,v[:D],v[D:(k+1)*D].reshape(k,D),v[(k+1)*D:],t,1.,RR)[0]
      J=jacobian(fun,allv);Jx=J[:,:D];Jprivate=J[:,:D*(k+1)];Jz=J[:,D*(k+1):]
      update('private_curl_fd_error',np.linalg.norm(Jx-Jx.T,2))
      ck(np.linalg.norm(Jx-Jx.T,2)<3e-7,'firsts','private curl')
      ck(np.linalg.norm(Jx,2)<=(13 if bounded else 1)/math.sqrt(k+1)+1e-6,'firsts','private bound')
      ck(np.linalg.norm(Jprivate,2)<=(13 if bounded else 1)+1e-6,'firsts','complete active/private operator, no coordinatewise inference')
      ck(np.linalg.norm(Jz,2)<=1/t+1e-6,'firsts','caller bound 1/t')
      M=rng.normal(size=(D,4))
      joint=np.column_stack([Jprivate,Jz@M])
      guard=math.sqrt((13 if bounded else 1)**2+np.linalg.norm(M,2)**2/t**2)
      ck(np.linalg.norm(joint,2)<=guard+1e-6,'caller_rows','disjoint caller coordinates')
      N=rng.normal(size=(D,D*(k+1)))
      shared=Jprivate+Jz@N
      ck(np.linalg.norm(shared,2)<=(13 if bounded else 1)+np.linalg.norm(N,2)/t+1e-6,'caller_rows','shared coordinates require actual row sum')

# Analytic radial Hessian independent of implementation.
def radial_second(x,R,w,v):
    r=np.linalg.norm(x)
    if r<R:return np.zeros_like(x)
    u=x/r
    if r<2*R:
      s=r/R-1;h=R*(1+s-s**3+s**4/2);hp=1-3*s*s+2*s**3;hpp=(-6*s+6*s*s)/R
    else:h=1.5*R;hp=0.;hpp=0.
    b=(hp-h/r)/r; uw=u@w;uv=u@v
    return hpp*uw*uv*u+b*((w@v-uw*uv)*u+uv*(w-uw*u)+uw*(v-uv*u))
for R in (.17,1.,7.):
  for rratio in (.9,1.,1.00001,1.2,1.5,1.99999,2.,2.2,4.):
    for _ in range(8):
      x=rng.normal(size=3);x*=R*rratio/np.linalg.norm(x)
      w=rng.normal(size=3);w/=np.linalg.norm(w);v=rng.normal(size=3);v/=np.linalg.norm(v)
      analytic=radial_second(x,R,w,v);h=1e-5*R
      _,Jp=src.radial_pullback(x+h*w,R);_,Jm=src.radial_pullback(x-h*w,R)
      err=np.linalg.norm((Jp-Jm)@v/(2*h)-analytic)
      update('radial_second_scaled_fd_error',R*err)
      ck(R*err<7e-5,'radial','analytic bilinear Hessian versus implementation')
      ck(np.linalg.norm(analytic)<=8/R,'radial','conservative complete Hessian bound')
      y,J=src.radial_pullback(x,R)
      ck(np.linalg.norm(J,2)<=1+1e-12,'radial','global contraction')
      ck(np.linalg.norm(y)<=min(np.linalg.norm(x),1.5*R)+1e-12,'radial','radial range')

# Verify that imported B_n already can cover the half-heat split.
for e in range(40):
    # HS bound supplied by source Bessel: (e+1)^(e/2) * t^-e.
    # At t>=tau/sqrt(2), require [2(e+1)]^(e/2) <= 2^e e!.
    ck((2*(e+1))**e <= 4**e*math.factorial(e)**2,'smoothed_constants','one-HS local bound includes sqrt(2)^e')
    ck(2**e*math.factorial(e) <= 4**e*math.factorial(e)**2,'smoothed_constants','proper-cut squared local bound')

# Fraction arithmetic is separate from the author's recurrence code.
Ts=[0,1];Ss=[0,1];grades=[]
for n in range(2,41):
    Ts.append(sum(i*(n-i)*Ts[i]*Ts[n-i] for i in range(1,n)))
    Ss.append(sum(math.comb(n,i)*i*(n-i)*Ss[i]*Ss[n-i] for i in range(1,n)))
    if n<4:continue
    gamma=Fraction(1,n);beta=Fraction(1,2*(n-1));root=Fraction(n)-Fraction(1,2)-(n-2)*gamma
    target=Fraction(n)+gamma
    ck(root+(n-1)*beta==n-(n-2)*gamma,'grades','exact amplitude product with width loss')
    ck(2*root>target,'grades','conditional root-quadratic surplus')
    ck(n+1-(n-1)*gamma==target,'grades','one-hit versus heat balance')
    ck(root-gamma>1,'grades','complete root-caller first grade')
    ck(root+beta-gamma>1,'grades','nonroot ancestor caller grade')
    floor=target-n+(n-2)*gamma
    ck(floor==Fraction(n-1,n),'cost_floors','response and quadrature floor before width amplification')
    ck(floor+gamma==1,'cost_floors','VALUE absolute floor at minimum width, up to known constants')
    grades.append({'n':n,'tau':str(gamma),'beta':str(beta),'root':str(root),'root_caller':str(root-gamma),'nonroot_caller':str(root+beta-gamma),'conditional_quadratic_extra_width_ell_zero':str(2*root),'conditional_target':str(target),'main_response_floor_exponent':str(floor),'main_value_floor_exponent':str(floor+gamma)})
ck(Ts[8]==794880 and Ss[8]==32049561600,'census','independent rank-eight recurrences')
ck(2**6*math.factorial(6)==46080,'census','rank-eight B8')
ck(Ss[8]*46080==1476843798528000,'census','rank-eight C8')

# Independent degree-multiset histogram, rather than V_n recurrence.
hist={1:{(0,):1}};weighted={1:{(0,):1}}
for n in range(2,9):
  hh=Counter();ss=Counter()
  for i in range(1,n):
    for left,nleft in hist[i].items():
      for right,nright in hist[n-i].items():
        for dl,ml in Counter(left).items():
          for dr,mr in Counter(right).items():
            l=list(left);r=list(right);l[l.index(dl)]+=1;r[r.index(dr)]+=1
            degrees=tuple(sorted(l+r));hh[degrees]+=nleft*nright*ml*mr
            ss[degrees]+=math.comb(n,i)*weighted[i][left]*weighted[n-i][right]*ml*mr
  hist[n]=dict(hh);weighted[n]=dict(ss)
  ck(sum(hh.values())==Ts[n],'census','degree histogram counts unmerged histories')
  ck(sum(ss.values())==Ss[n],'census','degree histogram integer coefficient sum')
raw_v8=sum(multiplicity*sum(2**d for d in degrees) for degrees,multiplicity in hist[8].items())
ck(raw_v8==27020800,'census','raw VALUE census from complete degree distribution')
ck(Fraction(27,2)-Fraction(43,8)==Fraction(65,8),'grades','extra inverse-width exponent 43 exactly spends rank-eight own-return surplus')

# Finite VALUE perturbations: a deliberately worst-case scalar sign at k=0.
for t in (.1,1.,3.):
    delta=1e-7;x=np.array([2.]);z=np.array([0.]);P=np.empty((0,1))
    exact,_=src.active_source(lambda q:.5*q,x,P,z,t,1.)
    approximate,_=src.active_source(lambda q:.5*q+np.where(q>0,delta,-delta),x,P,z,t,1.)
    ck(abs(np.linalg.norm(approximate-exact)-2*delta/t)<1e-15,'value_precision','2 delta_g/(At) can be attained')

# Production source must not materialize the dense diagnostic radial Jacobian.
def forbidden_dense(*args,**kwargs):raise AssertionError('production allocated a dense radial matrix')
saved_eye,saved_outer=np.eye,np.outer
try:
    np.eye,np.outer=forbidden_dense,forbidden_dense
    for k in (0,3):
        answer,calls=src.active_source(lambda q:.5*q,np.ones(100),np.ones((k,100)),np.zeros(100),.7,1.,1.)
        ck(answer.shape==(100,) and calls==2**(k+1),'linear_storage','no np.eye or np.outer in actual source')
finally:np.eye,np.outer=saved_eye,saved_outer

hashes=json.loads((HERE/'AUDITED-INPUT-SHA256.json').read_text())
for path,expected in hashes.items():
    ck(hashlib.sha256(Path(path).read_bytes()).hexdigest()==expected,'input_hashes','audit input unchanged: '+path)
result={'status':'PASS','assertions':sum(counts.values()),'groups':dict(counts),'extrema':extrema,'rank8':next(r for r in grades if r['n']==8),'rank8_raw_VALUE_census':raw_v8,'rank8_degree_multisets':len(hist[8]),'scope':'Source formulas and bounds diagnosed independently. No marked-tree native return, full sampler, or cost exponent is proved by these tests.','grades':grades}
(HERE/'independent-checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='grades'},indent=2))
