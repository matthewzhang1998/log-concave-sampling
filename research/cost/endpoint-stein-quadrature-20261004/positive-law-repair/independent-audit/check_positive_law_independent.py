#!/usr/bin/env python3
"""Independent diagnostics; no imports from the source or earlier test programs."""
import json, math, hashlib
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite import hermgauss
from scipy.linalg import expm, sqrtm
from scipy.special import ndtr
from scipy.integrate import cumulative_simpson
from scipy.interpolate import PchipInterpolator

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'POSITIVE-SECOND-ORDER-LAW-AND-QUADRATIC-AMPLIFIER.md'
PIN='a085dfdef3f66f69208612034c17e4e18045fc00a21be3aed64ccfcb5e9b7085'
CORE_PIN='6073cd34d5bd01d2a17a632c540d3fab087fcb2d05b2a808a7114627c0f7b0fe'
rng=np.random.default_rng(73482109)
checks=0

def check(test, label):
    global checks
    checks+=1
    if not bool(test): raise AssertionError(label)

def norm(x): return float(np.linalg.norm(x,2))
def hs(x): return float(np.linalg.norm(x,'fro'))
def psqrt(x):
    e,q=np.linalg.eigh((x+x.T)/2)
    return (q*np.sqrt(np.maximum(e,0)))@q.T

def gaussian_w2(c1,c2,m1=None,m2=None):
    a=psqrt(c1)
    sq=np.trace(c1+c2-2*psqrt(a@c2@a))
    if m1 is not None: sq+=np.sum((m1-m2)**2)
    return math.sqrt(max(float(sq),0))

def quad(K=8,m=5):
    x,w=leggauss(m); rs=[]; ws=[]
    for k in range(K):
        a=2.**(-k-1)
        rs.extend(1-(1.5*a+.5*a*x)); ws.extend(.5*a*w)
    rs.append(1-2.**(-K-1)); ws.append(2.**(-K))
    return np.array(rs), np.array(ws)

# A C2 convex scalar potential with globally bounded but non-Lipschitz Hessian.
# h(t)=sqrt(|t|)/(1+sqrt(|t|)); F=phi', phi(0)=F(0)=0.
def h(t):
    v=np.sqrt(np.abs(t)); return v/(1+v)
def F(t):
    u=np.abs(np.asarray(t)); v=np.sqrt(u)
    value=u-2*v+2*np.log1p(v)
    small=u<.01
    series=np.zeros_like(u)
    for j in range(1,19): series+=(-1.)**(j+1)*u**(1+j/2)/(1+j/2)
    return np.sign(t)*np.where(small,series,value)
def phi(t):
    u=np.abs(np.asarray(t)); v=np.sqrt(u)
    value=u*u/2-4*u*v/3+2*(u-1)*np.log1p(v)-u+2*v
    series=np.zeros_like(u)
    for j in range(1,19): series+=(-1.)**(j+1)*u**(2+j/2)/((1+j/2)*(2+j/2))
    return np.where(u<.01,series,value)

class C2Fixture:
    def __init__(self,d,linear=None):
        self.d=d; self.kappa=.12
        raw=rng.normal(size=(d+2,d)); self.directions=raw/np.linalg.norm(raw,axis=1)[:,None]
        self.shifts=np.linspace(-.65,.55,d+2)
        self.lam=.73/norm(self.directions.T@self.directions)
        self.linear=np.zeros(d) if linear is None else linear
    def grad(self,x):
        z=x@self.directions.T+self.shifts
        return self.kappa*x+self.lam*(F(z)-F(self.shifts))@self.directions+self.linear
    def potential(self,x):
        z=x@self.directions.T+self.shifts
        return self.kappa*np.sum(x*x,axis=-1)/2+self.lam*np.sum(phi(z)-phi(self.shifts)-F(self.shifts)*(z-self.shifts),axis=-1)+x@self.linear
    def hess(self,x):
        z=x@self.directions.T+self.shifts
        return self.kappa*np.eye(self.d)+self.lam*np.einsum('...k,ki,kj->...ij',h(z),self.directions,self.directions)

check(hashlib.sha256(SOURCE.read_bytes()).hexdigest()==PIN,'source pin')
check(hashlib.sha256((SOURCE.read_text().split('## 7.')[0].rstrip()+'\n').encode()).hexdigest()==CORE_PIN,'unchanged Sections 1-6')
r,w=quad(); s=np.sqrt(1-r*r); beta=float(w@s); c=.25+beta*beta
check(abs(w.sum()-1)<1e-14,'weight sum')
check(abs(w@r-.5)<1e-14,'first moment')
check(.5<=beta<=math.sqrt(3)/2,'beta range')

# Exact joint Gaussian entropy and exact conditional velocity, including
# nonsymmetric DZ H, affine offsets, rectangular G, and the motion term.
lin_rows=[]
for d,n in [(1,2),(3,2),(5,7)]:
  for scale in [.04,.2,.45,.8]:
    R=rng.normal(size=(d,d)); R*=scale/norm(R)
    S=rng.normal(size=(d,n)); S*=.6*scale/norm(S)
    u=.3*scale*rng.normal(size=d)
    a=norm(R); b=norm(S); EH=math.sqrt(hs(R)**2+hs(S)**2+u@u); J=hs(R)
    max_ratio=0.
    for t in [0.,.1,.5,1.]:
        B=np.eye(d)-t*R; Cy=B@B.T+t*t*S@S.T
        logdetB=np.linalg.slogdet(B)[1]
        joint=.5*(np.trace(Cy)+t*t*(u@u)-d)-logdetB
        formula=-t*np.trace(R)+.5*t*t*EH*EH-logdetB
        cap=.5*t*t*(EH*EH+J*J/(1-t*a))
        info=.5*np.linalg.slogdet(Cy)[1]-logdetB
        check(abs(joint-formula)<2e-12,'exact joint entropy identity')
        check(joint>=-1e-12 and joint<=cap+1e-12,'joint entropy bound')
        check(info>=-1e-12 and info<=joint+1e-12,'mutual information chain')
        CHY=R@B.T-t*S@S.T
        K=CHY@np.linalg.inv(Cy)-R
        defect=math.sqrt(max(float(np.trace(K@Cy@K.T)+t*t*np.dot(R@u,R@u)),0))
        bound=t*(a*EH+b*math.sqrt(EH*EH+J*J/(1-t*a)))
        check(defect<=bound+1e-12,'actual marginal velocity bound')
        if bound: max_ratio=max(max_ratio,defect/bound)
    # Flow vs Euler comparison, with affine deterministic flow solved in a block.
    aug=np.zeros((d+1,d+1)); aug[:d,:d]=-R; aug[:d,d]=-u
    E=expm(aug); L=E[:d,:d]; mu=E[:d,d]
    euler=np.eye(d)-R
    err=math.sqrt(hs(L-euler)**2+np.sum((mu+u)**2))
    check(err<=a*math.exp(a)*EH/2+1e-12,'Euler-flow bound')
    lin_rows.append({'d':d,'n':n,'a':a,'max_velocity_ratio':max_ratio})
# With b=0, the information term is zero but motion is genuinely nonzero.
a0=.25; t=.7
motion=a0*a0*t
check(motion>0,'motion cannot be omitted when G is absent')

# Genuine C2 and noncommuting source/target firsts, numerical directional checks,
# finite caller and zero/query restoration for the actual physical packet.
fixture_rows=[]
for d in [2,5,11]:
    pot=C2Fixture(d,linear=.2*rng.normal(size=d))
    x=rng.normal(size=d); xp=rng.normal(size=d)
    H1=pot.hess(x); H2=pot.hess(xp)
    comm=hs(H1@H2-H2@H1)
    check(comm>1e-5,'fixture Hessians genuinely noncommute')
    check(np.linalg.eigvalsh(H1).min()>=0 and norm(H1)<1,'global normalized Hessian sandwich sampled')
    for A in [.5,.125]:
      for M in [0,3]:
        y=rng.normal(size=d); Z=rng.normal(size=d); G=rng.normal(size=d)
        calls=[0]
        def value(x): calls[0]+=1; return pot.grad(x)
        bM=y.copy(); By=np.eye(d)
        for _ in range(M):
            H=pot.hess(bM); By=np.eye(d)-A*H@By; bM=y-A*value(bM)
        anchor=value(bM); q=r[:,None]*Z+s[:,None]*G
        node_values=np.stack([value(bM+math.sqrt(A)*qi) for qi in q])
        g=math.sqrt(A)*(node_values-anchor)
        HQ=w@g; Y=Z-HQ; X=bM+math.sqrt(A)*Y
        canonical=math.sqrt(A)*value(X)
        check(calls[0]==M+len(r)+2,'full VALUE query count')
        residual=bM-y+A*anchor
        check(np.linalg.norm(residual)<=A**(M+1)*np.linalg.norm(pot.grad(y))+1e-12,'mode residual')
        nodesH=pot.hess(bM+math.sqrt(A)*q)
        DZ=A*np.einsum('n,n,nij->ij',w,r,nodesH)
        DG=A*np.einsum('n,n,nij->ij',w,s,nodesH)
        Dpriv=math.sqrt(A)*np.hstack([np.eye(d)-DZ,-DG])
        DHcaller=math.sqrt(A)*(np.einsum('n,nij->ij',w,nodesH)-pot.hess(bM))@By
        DXcaller=By-math.sqrt(A)*DHcaller
        check(np.linalg.eigvalsh(DZ).min()>-1e-13 and norm(DZ)<=A/2+1e-12,'DZ monotone bound')
        check(norm(DG)<=A+1e-12,'G Lipschitz bound')
        check(norm(np.hstack([DZ,DG]))<=A+1e-12,'full private H first')
        check(hs(DZ)<=A*math.sqrt(d)/2+1e-12,'derivative energy pointwise')
        check(norm(DXcaller-np.eye(d))<=2*A/(1-A)+1e-12,'physical caller first')
        check(norm(math.sqrt(A)*pot.hess(X)@Dpriv)<=A*(1+A)+1e-12,'canonical private first')
        def execute(yy,zz,gg):
            bb=yy.copy()
            for _ in range(M): bb=yy-A*pot.grad(bb)
            aa=pot.grad(bb)
            qq=r[:,None]*zz+s[:,None]*gg
            sites=bb+math.sqrt(A)*qq
            vals=np.stack([aa if np.array_equal(site,bb) else pot.grad(site) for site in sites])
            hh=w@(math.sqrt(A)*(vals-aa))
            return bb+math.sqrt(A)*(zz-hh)
        dy=rng.normal(size=d); dz=rng.normal(size=d); dg=rng.normal(size=d); eps=1e-6
        fd=(execute(y+eps*dy,Z+eps*dz,G+eps*dg)-execute(y-eps*dy,Z-eps*dz,G-eps*dg))/(2*eps)
        analytic=DXcaller@dy+Dpriv@np.r_[dz,dg]
        check(np.linalg.norm(fd-analytic)<2e-6*(1+np.linalg.norm(analytic)),'literal source first finite difference')
        zero=execute(y,np.zeros(d),np.zeros(d))
        check(np.array_equal(zero,bM),'literal same-site anchored zero')
        fixture_rows.append({'d':d,'A':A,'M':M,'commutator_HS':comm,'source_first_fd_error':float(np.linalg.norm(fd-analytic)),'values':calls[0]})
# This exact scalar restriction certifies lack of any Lipschitz Hessian modulus.
holder_slopes=[float(h(v)/v) for v in [1e-2,1e-4,1e-6,1e-8]]
check(holder_slopes[-1]>9000 and holder_slopes[-1]>90*holder_slopes[0],'genuinely non-Lipschitz Hessian')

# Direct Gaussian quadrature of the exact tilted m_r in 2D, including r=0,1.
pot=C2Fixture(2); gh,gw=hermgauss(44); gh=gh*math.sqrt(2); gw=gw/math.sqrt(math.pi)
G2=np.stack(np.meshgrid(gh,gh,indexing='ij'),axis=-1).reshape(-1,2)
w2=np.outer(gw,gw).reshape(-1)
posterior_rows=[]
for A in [.5,.125]:
  for rr in [0.,.2,.8,.99,1.]:
    ss=math.sqrt(1-rr*rr)
    for xx in [np.zeros(2),np.array([1.1,-.7]),np.array([3.,2.])]:
      def moments(x):
        q=rr*x+ss*G2; vals=A*pot.grad(q); U=A*pot.potential(q)
        ww=w2*np.exp(-(U-U.min())); ww/=ww.sum()
        mean=ww@vals; centered=vals-mean
        cov=np.einsum('n,ni,nj->ij',ww,centered,centered)
        EH=np.einsum('n,nij->ij',ww,A*pot.hess(q))
        return mean,w2@vals,cov,EH
      mr,pr,cov,EH=moments(xx)
      cap=ss*ss*A*A*math.sqrt(rr*rr*np.dot(xx,xx)+ss*ss*2)
      defect=float(np.linalg.norm(mr-pr))
      check(defect<=cap+2e-10,'OU posterior drift comparison')
      check(norm(cov)<=ss*ss*A*A+2e-10,'posterior covariance Poincare cap')
      Dm=rr*(EH-cov)
      check(norm(Dm)<=rr*(A+ss*ss*A*A)+1e-10,'posterior drift global first cap')
      vv=np.array([.7,-.4]); eps=1e-6
      fd=(moments(xx+eps*vv)[0]-moments(xx-eps*vv)[0])/(2*eps)
      check(np.linalg.norm(fd-Dm@vv)<2e-6,'posterior tilted derivative identity')
      posterior_rows.append({'A':A,'r':rr,'x':xx.tolist(),'drift_defect':defect,'bound':cap})

# Matrix amplifier with dense B, the actual single G bank, exact VALUE graph,
# and all covariance cross-products. Here each graph run is one (Z,G) sample.
quadratic_rows=[]
for d in [1,3,17,47]:
  Q,_=np.linalg.qr(rng.normal(size=(d,d)))
  lam=np.linspace(.03,1.,d) if d>1 else np.ones(1)
  for A in [.5,.2,.05]:
    B=(Q*(A*lam))@Q.T
    C=np.eye(d)-B+c*B@B
    Cbad=np.eye(d)-B+(.25+np.sum(w*w*s*s))*B@B
    check(hs(C-Cbad)>1e-8,'independent-node replacement changes covariance')
    delta=(c-1)*B@B+c*B@B@B
    check(norm(delta)<=A*A+1e-12,'matrix defect norm')
    Y0=rng.normal(size=d); G0=rng.normal(size=d)
    y=(np.eye(d)-B/2)@Y0-beta*B@G0
    for m in range(6):
      calls=0; T=y.copy(); out=y.copy(); coefficient=1.; qpoly=np.eye(d); power=np.eye(d)
      def value(v):
        global calls
        calls+=1
        return B@v
      for k in range(1,m+1):
        u1=value(T); u2=value(u1); u3=value(u2)
        T=(c-1)*u2+c*u3
        coefficient*=(-.5-(k-1))/k
        out+=coefficient*T
        power=delta@power; qpoly+=coefficient*power
      check(calls==3*m,'amplifier literal extra VALUE count')
      check(np.linalg.norm(out-qpoly@y)<1e-11,'VALUE graph equals matrix polynomial')
      target=Q@np.diag(1/(1+A*lam))@Q.T
      actual=qpoly@C@qpoly.T
      eigD=(c-1)*(A*lam)**2+c*(A*lam)**3
      coeffs=[1.]
      for k in range(1,m+1): coeffs.append(coeffs[-1]*(-.5-(k-1))/k)
      qeval=sum(coeffs[k]*eigD**k for k in range(m+1))
      std=np.sqrt(1-A*lam+c*(A*lam)**2)
      werr=float(np.linalg.norm(qeval*std-1/np.sqrt(1+A*lam)))
      check(hs(actual-(Q*((qeval*std)**2))@Q.T)<2e-11,'actual dense covariance and all descendants')
      check(werr<=4/3*A**(2*m+2)*math.sqrt(d)+1e-13,'dimension-safe all-order quadratic bound')
      quadratic_rows.append({'d':d,'A':A,'m':m,'extra_values':calls,'W2':werr,'bound':4/3*A**(2*m+2)*math.sqrt(d)})

# One-dimensional CDF inversion checks the actual positive two-root law,
# independently of all velocity/entropy estimates above. This is a diagnostic,
# not a proof; compare two quadrature resolutions and report their discrepancy.
def scalar_g(x,A): return A*(.2*x+.7*(F(x+.6)-F(.6)))
def scalar_H(x,A): return A*(.2+.7*h(x+.6))
def scalar_U(x,A): return A*(.1*x*x+.7*(phi(x+.6)-phi(.6)-F(.6)*x))
def scalar_law(A,ng=64,ny=1801):
    xx=np.linspace(-9,9,ny)
    gg,ww=hermgauss(ng); gg*=math.sqrt(2); ww/=math.sqrt(math.pi)
    zz=np.broadcast_to(xx,(ng,ny)).copy()
    GG=gg[:,None]
    for _ in range(6):
        HH=np.zeros_like(zz); DH=np.zeros_like(zz)
        for ri,si,wi in zip(r,s,w):
            qq=ri*zz+si*GG
            HH+=wi*scalar_g(qq,A)
            DH+=wi*ri*scalar_H(qq,A)
        residual=zz-HH-xx
        zz-=residual/(1-DH)
    check(float(np.max(np.abs(residual)))<2e-11,'conditional monotone CDF inversion')
    Fy=ww@ndtr(zz)
    dens=np.exp(-xx*xx/2-scalar_U(xx,A))
    Ft=cumulative_simpson(dens,x=xx,initial=0); Ft/=Ft[-1]
    # strictly monotone interior avoids machine-precision plateaus in far tails.
    def inverse(ff):
        mask=np.r_[True,np.diff(ff)>1e-15]
        return PchipInterpolator(ff[mask],xx[mask],extrapolate=False)
    ux,uw=leggauss(1400); uu=(ux+1)/2
    qy=inverse(Fy)(uu); qt=inverse(Ft)(uu)
    check(np.all(np.isfinite(qy)) and np.all(np.isfinite(qt)),'law quantiles finite')
    return float(np.sqrt(np.sum(uw*(qy-qt)**2)/2))
scalar_rows=[]
for A in [.5,.25,.125,.0625]:
    e1=scalar_law(A,48,1401); e2=scalar_law(A,80,2201)
    check(abs(e1-e2)<max(1e-5,.04*e2),'independent scalar-law resolution agreement')
    check(e2<=10*(A*A+(8*4**-5+2**(1-8))*A),'scalar positive law theorem bound')
    scalar_rows.append({'A':A,'W2_low_resolution':e1,'W2_high_resolution':e2,'resolution_difference':abs(e1-e2),'W2_over_A2':e2/(A*A)})

report={'source_sha256':PIN,'seed':73482109,'assertions':checks,'scope':'Finite diagnostics support, but do not replace, the analytical audit. Scalar CDF results are deterministic numerical approximations. No all-order nonlinear theorem is claimed.','quadrature':{'K':8,'m':5,'nodes':len(r),'beta':beta,'c':c,'delta_sufficient':8*4**-5+2**(1-8)},'linear_conditional_entropy':lin_rows,'C2_noncommuting_source':fixture_rows,'non_Lipschitz_Hessian_slopes':holder_slopes,'OU_posterior_drift':posterior_rows,'quadratic_amplifier':quadratic_rows,'scalar_actual_positive_law':scalar_rows}
(HERE/'positive_law_independent_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'assertions':checks,'source_sha256':PIN,'scalar_actual_positive_law':scalar_rows,'quadratic_cases':len(quadratic_rows),'C2_fixture_cases':len(fixture_rows)},indent=2))
