#!/usr/bin/env python3
"""Independent diagnostics for the nested finite force-MEAN port.
No test below purports to establish a continuum or completed-consumer theorem.
"""
import json, math, pathlib
import numpy as np
from numpy.polynomial.legendre import leggauss
ROOT=pathlib.Path(__file__).resolve().parent
checks=0

def ck(test, description):
    global checks
    assert bool(test), description
    checks += 1

def rule(K,m):
    x,v=leggauss(m); rs=[]; ws=[]
    for j in range(K):
        a=2.**(-j-1)
        rs.extend(1-1.5*a-.5*a*x); ws.extend(.5*a*v)
    rs.append(1-2.**(-K-1)); ws.append(2.**(-K))
    return np.array(rs),np.array(ws)

rng=np.random.default_rng(941128)
r,w=rule(3,3); tau,v=rule(2,3)
q=np.sqrt(1-r*r); s=np.sqrt(1-tau*tau)
ck(np.all(w>0) and np.all(v>0),'positive mass')
ck(abs(w.sum()-1)<1e-14 and abs(v.sum()-1)<1e-14,'unit mass')
ck(abs(w@r-.5)<1e-14 and abs(v@tau-.5)<1e-14,'first moments')

# Conditional residuals are independent of the entire future; explicitly realize
# the finite triple (Z,X_r,{X_(r tau)}) with separate G and inner-bank columns.
genealogy=[]
for J in [3,7,17,41]:
    ts=(np.arange(J)+.5)/J
    R=np.minimum.outer(ts,ts)/np.maximum.outer(ts,ts)
    C=R-np.outer(ts,ts)
    eig,U=np.linalg.eigh(C)
    ck(eig.min()>-1e-13,'PSD conditional past covariance')
    L=U@np.diag(np.sqrt(np.maximum(eig,0)))
    ck(np.max(abs(L@L.T-C))<3e-13,'conditional residual factor')
    for rr in [.001,.2,.7,.99,1-1e-9]:
        qq=math.sqrt(1-rr*rr)
        rows=np.zeros((J+2,J+2)); rows[0,0]=1
        rows[1,:2]=[rr,qq]
        rows[2:,0]=ts*rr; rows[2:,1]=ts*qq; rows[2:,2:]=L
        times=np.r_[1,rr,rr*ts]
        target=np.minimum.outer(times,times)/np.maximum.outer(times,times)
        ck(np.max(abs(rows@rows.T-target))<4e-13,'full Markov OU covariance')
        ck(np.max(abs(rows[2:,1]-ts*qq))<1e-15,'outer G derivative coefficient q tau')
        ck(np.max(abs(np.sum(L*L,axis=1)-(1-ts*ts)))<4e-13,'row norm private past')
    energy=float(np.ones(J)@R@np.ones(J)/J**2)
    genealogy.append({'nodes':J,'inner_variance':energy,'residual_variance':energy-.25})
ck(abs(genealogy[-1]['inner_variance']-.5)<.002,'Riemann covariance convergence')

# A smooth gradient with noncommuting pointwise Hessians, within the C2 sandwich.
D=3
vectors=np.array([[1.,.3,-.2],[.2,1.,.4],[-.2,.4,1.]])
vectors/=np.linalg.norm(vectors,axis=1)[:,None]
A=.13

def g(x):
    return A*(.15*x+.25*np.tanh(vectors@x)@vectors)

def dg(x):
    a=vectors@x
    return A*(.15*np.eye(D)+.25*vectors.T@np.diag(1-np.tanh(a)**2)@vectors)

def actual(z,G,H):
    Bs=[]; Es=[]; Fz=np.zeros((D,D)); FG=Fz.copy(); FH=Fz.copy()
    Ez=Fz.copy(); EG=Fz.copy(); EH=Fz.copy()
    for wi,ri,qi in zip(w,r,q):
        x=ri*z+qi*G
        inn=[ti*x+si*H for ti,si in zip(tau,s)]
        h=sum(vj*g(y) for vj,y in zip(v,inn))
        T=sum(vj*ti*dg(y) for vj,ti,y in zip(v,tau,inn))
        S=sum(vj*si*dg(y) for vj,si,y in zip(v,s,inn))
        C=dg(x-h); B=dg(x)
        Bs.append(wi*g(x)); Es.append(wi*(g(x-h)-g(x)))
        Fz+=wi*ri*C@(np.eye(D)-T); FG+=wi*qi*C@(np.eye(D)-T)
        FH-=wi*C@S
        Ez+=wi*ri*((C-B)-C@T); EG+=wi*qi*((C-B)-C@T)
        EH-=wi*C@S
    b=sum(Bs); e=sum(Es)
    return b+e,b,e,np.hstack([Fz,FG,FH]),np.hstack([Ez,EG,EH])

max_fd=0.; max_curl=0.; max_comm=0.
for _ in range(24):
    z,G,H=rng.normal(size=(3,D))*rng.uniform(.1,4)
    f,b,e,df,de=actual(z,G,H)
    inp=np.r_[z,G,H]
    for index in range(3*D):
        unit=np.eye(3*D)[index]; eps=1e-5
        plus=actual(*(inp+eps*unit).reshape(3,D))
        minus=actual(*(inp-eps*unit).reshape(3,D))
        fdF=(plus[0]-minus[0])/(2*eps)
        fdE=(plus[2]-minus[2])/(2*eps)
        err=max(np.max(abs(fdF-df[:,index])),np.max(abs(fdE-de[:,index])))
        max_fd=max(max_fd,float(err)); ck(err<3e-10,'literal chain rule finite difference')
    depriv=de[:,D:]; dfpriv=df[:,D:]
    ck(np.linalg.norm(depriv,2)<=A*(float(w@q)+A)+1e-13,'complete private first E')
    ck(np.linalg.norm(dfpriv,2)<=A*(float(w@q)+A)+1e-13,'complete private first F')
    ck(np.linalg.norm(de[:,:D],2)<=A/2+A*A/4+1e-13,'captured carrier first')
    lift=np.vstack([depriv,np.zeros((D,2*D))]); skew=lift-lift.T
    c=float(np.linalg.norm(skew,2)); max_curl=max(max_curl,c)
    ck(c<=A*A*(float(w@q)+float(v@s))+1e-13,'complete square-lift curl')
    ori=actual(z,np.zeros(D),np.zeros(D))[2]
    ck(np.linalg.norm(ori)<=A*A*np.linalg.norm(z)/4+1e-13,'captured E origin amplitude')
    M=dg(z)@dg(G)-dg(G)@dg(z); max_comm=max(max_comm,float(np.linalg.norm(M)))
ck(max_comm>1e-6,'fixture truly noncommuting')
ck(np.max(abs(actual(np.zeros(D),np.zeros(D),np.zeros(D))[0]))==0,'all-zero graph exact zero')

# Quadratics give a fully explicit Euler-decoupling law, including translated z.
# The inner continuous residual variance is 1/4 and common-H variance beta^2.
U,_=np.linalg.qr(rng.normal(size=(D,D)))
lam=np.array([.1,.5,.9]); quadratic=[]
for amp in [.5,.1,.01]:
    B=amp*U@np.diag(lam)@U.T
    expect=.5*B-.25*B@B
    mean=sum(wi*ri*(B-B@B*float(v@tau)) for wi,ri in zip(w,r))
    ck(np.linalg.norm(mean-expect)<1e-14,'exact quadratic continuous force mean')
    for rr in [.0,.5,.99,1-1e-8]:
        qq=math.sqrt(1-rr*rr)
        for residual in [.25,float(v@s)**2]:
            b=amp*lam
            sd0=qq*(1-b/2)
            sd1=np.sqrt(sd0*sd0+residual*b*b)
            dist=float(np.linalg.norm(sd1-sd0))
            bound=amp*amp*math.sqrt(D)/qq
            ck(dist<=bound+1e-14,'physical quadratic decoupling q singularity')
            quadratic.append({'A':amp,'r':rr,'residual_variance':residual,'W2':dist,'scaled_ratio':dist/bound})

# Uniform Hermite quadrature applies to arbitrary L2 vector fields, so retain
# the literal multiplier check rather than assume the nonlinear Psi is smooth.
rr,ww=rule(7,7); delta=8*4.**(-7)+2.**(1-7)
max_moment=0.
for degree in list(range(100))+[1000,10000,100000,1000000]:
    err=abs(float(ww@(rr**degree))-1/(degree+1))
    max_moment=max(max_moment,err)
    ck(err<=delta,'uniform Hermite moment allowance')

report={'assertions':checks,'genealogy':genealogy,'max_derivative_finite_difference_error':max_fd,
        'max_square_curl':max_curl,'max_noncommuting_hessian_commutator':max_comm,
        'max_checked_hermite_error':max_moment,'quadratic_decoupling':quadratic,
        'scope':'Finite diagnostics only; continuum limit and mean-law compiler rely on the independent written proof and imported consumer guards.'}
(ROOT/'nested_force_mean_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'assertions':checks,'max_fd':max_fd,'max_curl':max_curl,'max_hermite':max_moment},indent=2))
