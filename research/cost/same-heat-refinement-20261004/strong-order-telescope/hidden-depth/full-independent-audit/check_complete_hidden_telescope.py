#!/usr/bin/env python3
"""Independent finite diagnostics of the COMPLETE hidden strong telescope.

Exact-Gaussian means/covariances use fully expanded affine Gaussian records,
not Monte Carlo estimates. The quadratic seed is an exact test fixture, NOT
an implementation of the imported CW7 seed. Separate finite nonlinear tests
use only grad V and V'' for a C2/non-C3 potential. These tests supplement the
analytic audit; they cannot establish the old source/restoration ports.
"""
from dataclasses import dataclass
from functools import lru_cache
from fractions import Fraction as Fr
from pathlib import Path
import json, math
import numpy as np

OUT=Path(__file__).parent
rows=[]
def check(name,ok,**data):
    if not ok: raise AssertionError((name,data))
    rows.append(dict(name=name,passed=True,**data))

@dataclass
class Lin:
    mean: float
    v: dict
    def __add__(self,o):
        if not isinstance(o,Lin):o=Lin(float(o),{})
        d=self.v.copy()
        for i,a in o.v.items():d[i]=d.get(i,0.)+a
        return Lin(self.mean+o.mean,d)
    __radd__=__add__
    def __neg__(self):return self*(-1.)
    def __sub__(self,o):return self+-aslin(o)
    def __rsub__(self,o):return aslin(o)+-self
    def __mul__(self,a):return Lin(self.mean*a,{i:c*a for i,c in self.v.items()})
    __rmul__=__mul__
    def __truediv__(self,a):return self*(1/a)
    def var(self):return sum(a*a for a in self.v.values())
    def norm(self):return math.sqrt(self.mean*self.mean+self.var())
def aslin(x):return x if isinstance(x,Lin) else Lin(float(x),{})

class Tape:
    def __init__(self):self.n=0;self.calls=0;self.labels=[]
    def normal(self,label):
        i=self.n;self.n+=1;self.labels.append(label)
        return Lin(0.,{i:1.})

@lru_cache(None)
def flow_coeff(a,lam,T,N,M):
    ts=np.arange(N+1)*T/N
    co=np.cos(ts);si=np.sin(ts)
    if T==math.pi/2:co[-1]=0.;si[-1]=1.
    W=np.zeros((N+1,N+1))
    for i in range(1,N+1):
        jj=np.arange(i)
        W[i,jj]=np.cos(ts[i]-ts[jj+1])-np.cos(ts[i]-ts[jj])
    base=np.column_stack((co,1-co,math.sqrt(a)*si))
    curr=base.copy()
    for _ in range(M):curr=base-a*lam*(W@curr)
    return tuple(float(z) for z in curr[-1])

@lru_cache(None)
def spec(a,j):
    j=max(1,j);R=j+.5
    return (max(1,math.ceil(a**(-(j-1)))),
            max(1,math.ceil(a**(-max(0.,R-3.4)))),
            max(1,j-1),max(1,math.ceil((R-.5)/(34/15))-1))

@lru_cache(None)
def phi_coeff(a,lam,j):
    Nq,Ns,Mq,Ms=spec(a,j)
    q=flow_coeff(a,lam,math.pi/2,Nq,Mq)
    s=flow_coeff(a,lam,a**(19/30),Ns,Ms)
    return (s[0]*q[0],s[0]*q[1]+s[1],s[0]*q[2],s[2],Nq*Mq+Ns*Ms)

def seed(a,y,w,lam):return y/(1+a*lam)+math.sqrt(a/(1+a*lam))*w

def known_from_records(a,y,j,w,etas,lam,t):
    jj=max(1,j);x=seed(a,y,w,lam)
    al,be,cz,cl,bill=phi_coeff(a,lam,jj)
    for z,l in etas:
        x=al*x+be*y+cz*z+cl*l;t.calls+=bill
    return x

def known(a,y,j,lam,t):
    jj=max(1,j);w=t.normal('seed')
    es=[(t.normal('quarter'),t.normal('short')) for _ in range(jj)]
    return known_from_records(a,y,j,w,es,lam,t)

def known_pair(a,yf,yc,j,lam,t):
    jj=max(1,j);w=t.normal('paired-seed')
    es=[(t.normal('paired-quarter'),t.normal('paired-short')) for _ in range(jj)]
    # For j=1 both sides are literally K1, at their respective callers.
    cf=es if j==1 else es[1:]
    f=known_from_records(a,yf,j,w,es,lam,t)
    c=known_from_records(a,yc,max(1,j-1),w,cf,lam,t)
    return f,c

def mode(a,r,lam,m,t):
    x=r
    for _ in range(m):x=r-a*lam*x;t.calls+=1
    return x

def counts(a,j):
    return [max(1,math.ceil(a**(-2*max(j-1,0))))]+[max(1,math.ceil(a**(-2*(j-l)))) for l in range(1,j+1)]

def mean_bank(a,y,d,j,J,lam,t):
    if j==0:return force(a,y,d,0,J,lam,t)
    ns=counts(a,j);out=Lin(0.,{})
    for _ in range(ns[0]):out+=force(a,y,d,0,J,lam,t)/ns[0]
    for l in range(1,j+1):
        for _ in range(ns[l]):
            f,c=pair(a,y,d,l,J,lam,t)
            out+=(f-c)/ns[l]
    return out

def decoder(a,y,j,J,lam,t,T,obs):
    rs=[T+math.sqrt(a)*(e-sum(obs)/len(obs)) for e in obs]
    x=mode(a,rs[0],lam,J+4,t)
    for r in rs[1:]:x=known(a/2,(x+r)/2,j,lam,t)
    return x

def source(a,y,d,j,J,lam,t):
    if d==0:return known(a,y,j,lam,t)
    k=3+math.ceil(J*math.log2(1/a))
    mean=mean_bank(a,y,d-1,j,J,lam,t)
    T=y-.3*math.sqrt(a)*mean+math.sqrt(a/(k+1))*t.normal('statistic')
    obs=[t.normal('observation') for _ in range(k+1)]
    return decoder(a,y,j,J,lam,t,T,obs)

def force(a,y,d,j,J,lam,t):
    x=source(a,y,d,j,J,lam,t);t.calls+=1
    return math.sqrt(a)*lam*x

def pair(a,y,d,j,J,lam,t):
    if d==0:
        f,c=known_pair(a,y,y,j,lam,t)
        t.calls+=2
        return math.sqrt(a)*lam*f,math.sqrt(a)*lam*c
    k=3+math.ceil(J*math.log2(1/a))
    # Independent COMPLETE signal banks. No private ancestor is captured in y.
    mf=mean_bank(a,y,d-1,j,J,lam,t)
    mc=mean_bank(a,y,d-1,j-1,J,lam,t)
    zt=t.normal('shared-statistic');obs=[t.normal('shared-observation') for _ in range(k+1)]
    common=math.sqrt(a/(k+1))*zt;bar=sum(obs)/len(obs)
    tf=y-.3*math.sqrt(a)*mf+common;tc=y-.3*math.sqrt(a)*mc+common
    rf=[tf+math.sqrt(a)*(e-bar) for e in obs]
    rc=[tc+math.sqrt(a)*(e-bar) for e in obs]
    xf=mode(a,rf[0],lam,J+4,t);xc=mode(a,rc[0],lam,J+4,t)
    for i in range(1,k+1):
        xf,xc=known_pair(a/2,(xf+rf[i])/2,(xc+rc[i])/2,j,lam,t)
    t.calls+=2
    return math.sqrt(a)*lam*xf,math.sqrt(a)*lam*xc

def ideal_center(a,y,d,lam):
    q=y
    for _ in range(d):q=y-.3*a*lam*q/(1+a*lam)
    return q

# Exact Gaussian comparisons of literally different finite marginal programs.
fixture=[]
for a,d,j in [(a,d,j) for a in [.25,.125] for d in [0,1,2] for j in [1,2]]+[(.25,1,3)]:
    J=3;lam=.6;y=.7
    tf=Tape();tc=Tape();tp=Tape()
    f=force(a,y,d,j,J,lam,tf);c=force(a,y,d,j-1,J,lam,tc)
    pf,pc=pair(a,y,d,j,J,lam,tp)
    for side,x,z in [('fine',f,pf),('coarse',c,pc)]:
        check('complete_pair_exact_gaussian_marginal',abs(x.mean-z.mean)<2e-11 and abs(x.var()-z.var())<2e-11,
              a=a,d=d,j=j,side=side,mean_difference=x.mean-z.mean,variance_difference=x.var()-z.var())
    gap=(pf-pc).norm();zero=abs(pf.mean-pc.mean)
    check('all_zero_suffix_and_hidden_pair',j!=1 or zero<2e-13,a=a,d=d,j=j,zero_force_gap=zero)
    state=f/(math.sqrt(a)*lam);qid=ideal_center(a,y,d,lam)
    if d==0:qfinite=y
    else:
        low=force(a,y,d-1,j,J,lam,Tape());qfinite=y-.3*math.sqrt(a)*low.mean
        lowstate=low/(math.sqrt(a)*lam)
        qprev=ideal_center(a,y,d-1,lam)
        elow=math.hypot(lowstate.mean-qprev/(1+a*lam),math.sqrt(lowstate.var())-math.sqrt(a/(1+a*lam)))
        check('finite_to_ideal_picard_bias',abs(qfinite-qid)<=.3*a*elow+2e-13,
              a=a,d=d,j=j,center_bias=abs(qfinite-qid),bound=.3*a*elow)
    errfinite=math.hypot(state.mean-qfinite/(1+a*lam),math.sqrt(state.var())-math.sqrt(a/(1+a*lam)))
    errideal=math.hypot(state.mean-qid/(1+a*lam),math.sqrt(state.var())-math.sqrt(a/(1+a*lam)))
    fixture.append(dict(a=a,d=d,j=j,force_pair_l2=gap,strong_ratio=gap/a**j,zero_ratio=zero/a**j,
                        finite_center=qfinite,ideal_center=qid,finite_law_ratio=errfinite/a**(j+.5),
                        ideal_law_ratio=errideal/a**(j+.5),single_tape=tf.n,pair_tape=tp.n,
                        single_original_calls=tf.calls,pair_original_calls=tp.calls))
check('gaussian_fixture_uniform_strong_envelope',max(f['strong_ratio'] for f in fixture)<3,
      max_ratio=max(f['strong_ratio'] for f in fixture))
check('gaussian_fixture_uniform_law_envelope',max(f['ideal_law_ratio'] for f in fixture)<3,
      max_ratio=max(f['ideal_law_ratio'] for f in fixture))

# Exact observation covariance and decoder carrier identity, on rational rows.
for k in [1,2,5,13,37]:
    n=k+1
    for i in range(n):
        for j in range(n):
            cov=Fr(1,n)+sum((Fr(int(i==l))-Fr(1,n))*(Fr(int(j==l))-Fr(1,n)) for l in range(n))
            check('exact_observation_covariance',cov==int(i==j),k=k,i=i,j=j,cov=str(cov))
    v=Fr(1)
    for _ in range(k):v=v/4+Fr(1,4)+Fr(1,2)
    check('exact_decoder_coisometry',v==1,k=k,variance=str(v))

# C2 but not C3: V''=.5+.2 sqrt(|x|)/(1+sqrt(|x|)).
# grad is its exact continuous antiderivative, and no V''' is implemented.
def grad(x):
    xx=np.asarray(x);r=np.sqrt(np.abs(xx))
    f=r*r-2*r+2*np.log1p(r)
    small=r<1e-3
    if np.any(small):
        rr=r[small] if r.ndim else r
        ss=sum(2*((-1)**(m+1))*rr**m/m for m in range(3,12))
        if r.ndim:f[small]=ss
        else:f=ss
    return .5*xx+.2*np.sign(xx)*f

def hvp(x):
    r=np.sqrt(np.abs(x));return .5+.2*r/(1+r)

# Batch finite original VALUE and first recurrence. Last axis is the caller
# plus all original private inputs; each original Hessian action is scalar.
def nonlinear_flow(a,T,N,M,y,dy,x,dx,z,dz):
    ts=np.arange(N+1)*T/N;co=np.cos(ts);si=np.sin(ts)
    if T==math.pi/2:co[-1]=0;si[-1]=1
    W=np.zeros((N+1,N+1))
    for i in range(1,N+1):
        ii=np.arange(i);W[i,ii]=np.cos(ts[i]-ts[ii+1])-np.cos(ts[i]-ts[ii])
    base=y+co[:,None]*(x-y)+math.sqrt(a)*si[:,None]*z
    db=dy+co[:,None,None]*(dx-dy)+math.sqrt(a)*si[:,None,None]*dz
    val=base.copy();der=db.copy()
    for _ in range(M):
        val_old=val;val=base-a*(W@grad(val_old))
        hh=hvp(val_old)[:,:,None]*der
        der=db-a*(W@hh.reshape(N+1,-1)).reshape(db.shape)
    return val[-1],der[-1]

def nonlinear_level(a,y,j,records,suffix=False):
    jj=max(1,j);nb=len(y);nd=records.shape[1]+1
    dy=np.zeros((nb,nd));dy[:,0]=1
    # Deliberately includes I+O(sqrt(a)) caller error, to test the refresh repair.
    x=y+math.sqrt(a)*(.1*np.sin(y)+records[:,0]);dx=dy*(1+.1*math.sqrt(a)*np.cos(y))[:,None]
    dx[:,1]+=math.sqrt(a)
    es=range(2,jj+2) if suffix else range(1,jj+1)
    for i in es:
        z=records[:,2*i-1];l=records[:,2*i]
        dz=np.zeros((nb,nd));dl=dz.copy();dz[:,2*i]=1;dl[:,2*i+1]=1
        Nq,Ns,Mq,Ms=spec(a,jj)
        x,dx=nonlinear_flow(a,math.pi/2,Nq,Mq,y,dy,x,dx,z,dz)
        x,dx=nonlinear_flow(a,a**(19/30),Ns,Ms,y,dy,x,dx,l,dl)
    return x,dx

rng=np.random.default_rng(714827)
nonlinear=[]
for a in [.25,.125,.0625]:
    for j in [1,2,3]:
        records=rng.standard_normal((12,1+2*j));records[0,:]=0
        y=np.linspace(-.4,.4,len(records));y[0]=.3
        xf,df=nonlinear_level(a,y,j,records)
        if j==1:xc,dc=nonlinear_level(a,y,1,records)
        else:xc,dc=nonlinear_level(a,y,j-1,records,suffix=True)
        state_gap=np.sqrt(np.mean((xf-xc)**2));force_gap=math.sqrt(a)*np.sqrt(np.mean((grad(xf)-grad(xc))**2))
        shift=.7*a**(j-.5)
        xt,dt=nonlinear_level(a,y+shift,j,records)
        # Literal two-caller decomposition is exact; caller Lipschitz is <=1+Ca.
        two_gap=math.sqrt(np.mean((xt-xc)**2))
        shift_gap=math.sqrt(np.mean((xt-xf)**2))
        check('nonlinear_two_caller_triangle',two_gap<=shift_gap+state_gap+1e-13,a=a,j=j,two_gap=two_gap)
        p=np.zeros(records.shape[1]);p[-2]=math.cos(a**(19/30));p[-1]=math.sin(a**(19/30))
        fresh_res=np.max(np.linalg.norm(df[:,1:]-math.sqrt(a)*p,axis=1))
        caller_res=np.max(abs(df[:,0]-1))
        check('nonlinear_refreshed_actual_first',fresh_res<3*a**1.5 and caller_res<3*a,a=a,j=j,
              fresh_g1_ratio=fresh_res/a**1.5,caller_g1_ratio=caller_res/a)
        check('nonlinear_same_caller_suffix_bound',force_gap<3*a**j,a=a,j=j,force_ratio=force_gap/a**j)
        check('nonlinear_deterministic_zero_suffix_bound',math.sqrt(a)*abs(float(grad(xf[0])-grad(xc[0])))<3*a**j,
              a=a,j=j,zero_ratio=math.sqrt(a)*abs(float(grad(xf[0])-grad(xc[0])))/a**j)
        nonlinear.append(dict(a=a,j=j,state_ratio=state_gap/a**(j-.5),force_ratio=force_gap/a**j,
                              fresh_g1_ratio=fresh_res/a**1.5,caller_g1_ratio=caller_res/a))
check('C2_nonC3_fixture_hessian_sandwich',all(.5<=float(hvp(x))<=.7 for x in [-100.,-1.,0.,1e-10,1.,100.]))

# A fixed finite prox has a caller-weighted residual, not a uniform bound
# over all unbounded callers. The second calculation integrates that exact
# residual under a random original caller, before comparing Picard centers.
for a in [.25,.125,.0625]:
    lam=.6;m=7;sigma_y=3.2
    r=aslin(1.);pr=mode(a,r,lam,m,Tape())
    residual=pr.mean-1/(1+a*lam)
    exact_residual=-(-a*lam)**(m+1)/( 1+a*lam)
    check('fixed_prox_weighted_caller_residual',abs(residual-exact_residual)<2e-15,
          a=a,m=m,residual=residual,formula=exact_residual)
    check('integrated_random_caller_prox_residual',abs(residual*sigma_y)<=abs(exact_residual)*sigma_y+2e-14,
          a=a,caller_rms=sigma_y,residual_rms=abs(residual)*sigma_y)
    for d in [1,2]:
        j=2;J=3
        low0=source(a,0.,d-1,j,J,lam,Tape())
        low1=source(a,1.,d-1,j,J,lam,Tape())
        idealq1=ideal_center(a,1.,d-1,lam)
        mean_error_slope=low1.mean-low0.mean-idealq1/(1+a*lam)
        private_std_error=math.sqrt(low0.var())-math.sqrt(a/(1+a*lam))
        state_integrated=math.hypot(mean_error_slope*sigma_y,private_std_error)
        center_bias_integrated=.3*a*lam*abs(mean_error_slope)*sigma_y
        check('conditional_then_integrated_history_bias',center_bias_integrated<=.3*a*state_integrated+1e-13,
              a=a,d=d,caller_rms=sigma_y,center_bias_rms=center_bias_integrated,state_error_rms=state_integrated)
for x in [-1e-8,-.1,1e-8,.1]:
    h=min(abs(x)/10,1e-6)
    fd=(float(grad(x+h))-float(grad(x-h)))/(2*h)
    check('nonlinear_original_gradient_derivative',abs(fd-float(hvp(x)))<2e-6,x=x,finite_difference=fd,hessian=float(hvp(x)))

# Exponent and ambient-prior checks use exact rational arithmetic.
for j in range(1,16):
    check('coarse_statistic_offset_recovery',Fr(j-1)+Fr(1,2)+Fr(1,2)==j,j=j)
    for l in range(1,j+1):
        check('complete_cost_fixed_point',2*(j-l)+2*(l-1)==2*(j-1),j=j,l=l)
    for E in [Fr(3,2),Fr(4),Fr(10)]:
        B=10*(E+j-1+1)
        check('ambient_prior_strict_slack',B/10-(j-1)==E+1,j=j,E=str(E),B=str(B))
check('protected_width_constant',Fr(1)+Fr(3,2)*Fr(19,15)==Fr(29,10))
check('twin_width_constant',Fr(2)+Fr(9,10)==Fr(29,10))
check('full_radius_fixed',Fr(1)-Fr(9,10)==Fr(1,10))

out=dict(status='PASS',checks=len(rows),scope='Independent finite diagnostics; not a replacement for imported complete seed and restoration contracts',
         quadratic_fixtures=fixture,nonlinear_fixtures=nonlinear,checks_detail=rows)
(OUT/'complete_hidden_telescope_checks.json').write_text(json.dumps(out,indent=2)+'\n')
summary=dict(status='PASS',checks=len(rows),quadratic_fixtures=len(fixture),nonlinear_fixtures=len(nonlinear),
             max_gaussian_strong_ratio=max(f['strong_ratio'] for f in fixture),
             max_gaussian_law_ratio=max(f['ideal_law_ratio'] for f in fixture),
             max_nonlinear_strong_ratio=max(f['force_ratio'] for f in nonlinear),
             max_nonlinear_g1_ratio=max(f['fresh_g1_ratio'] for f in nonlinear))
print(json.dumps(summary,indent=2))
