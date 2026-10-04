#!/usr/bin/env python3
"""Author diagnostics. Finite checks supplement, and do not prove, the C2 theorem."""
import json, math
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite import hermgauss
from scipy.special import ndtr
from scipy.integrate import cumulative_simpson

OUT=Path(__file__).resolve().parent
checks=0

def check(ok,msg):
    global checks
    if not ok: raise AssertionError(msg)
    checks+=1

def rule(K,m):
    xx,ww=leggauss(m); rs=[]; ws=[]
    for k in range(K):
        a=2.**(-k-1)
        t=1.5*a+.5*a*xx
        rs.extend(1-t); ws.extend(.5*a*ww)
    rs.append(1-2.**(-K-1)); ws.append(2.**(-K))
    return np.array(rs),np.array(ws)

r,w=rule(12,8); s=np.sqrt(1-r*r); beta=w@s; c=.25+beta*beta
check(abs(w.sum()-1)<1e-14,'weight mass')
check(abs(w@r-.5)<1e-14,'first moment')
check(.5<=beta<=math.sqrt(3)/2,'beta limits')

# Quadratic, non-diagonal matrices: exact original-query graph and covariance.
rng=np.random.default_rng(914753)
quad=[]
for d in (1,2,5,13,31):
  for A in (.5,.25,.125,.0625):
    Q,_=np.linalg.qr(rng.normal(size=(d,d)))
    B=(Q*rng.uniform(.05*A,A,d))@Q.T
    Z=rng.normal(size=d); G=rng.normal(size=d)
    H=sum(wi*(B@(ri*Z+si*G)) for ri,si,wi in zip(r,s,w))
    Y=Z-H; C=np.eye(d)-B+c*(B@B)
    check(np.linalg.norm(Y-((np.eye(d)-B/2)@Z-beta*B@G))<2e-14,'literal shared root')
    Db=(c-1)*(B@B)+c*(B@B@B)
    check(np.linalg.norm((np.eye(d)+B)@C-np.eye(d)-Db)<2e-14,'covariance polynomial')
    check(np.linalg.norm(Db,2)<=A*A+1e-14,'small correction')
    T=Y.copy(); amp=Y.copy(); q=np.eye(d); coeff=1.; matpow=np.eye(d)
    vals=np.linalg.eigvalsh(B)
    for m in range(6):
      if m:
        u1=B@T; u2=B@u1; u3=B@u2
        T=(c-1)*u2+c*u3
        coeff*=(-.5-(m-1))/m
        amp+=coeff*T; matpow=Db@matpow; q+=coeff*matpow
      check(np.linalg.norm(amp-q@Y)<3e-13,'VALUE polynomial equality')
      target=1/(1+vals); cov=np.linalg.eigvalsh(q@C@q)
      w2=float(np.linalg.norm(np.sqrt(cov)-np.sqrt(target[::-1]))) # eigen order reverses target
      # spectral computation avoids small nonsimultaneous eigensolver sorting ambiguity
      db=(c-1)*vals**2+c*vals**3
      qq=np.ones(d); cp=np.ones(d); cc=1.
      for k in range(1,m+1):
        cp*=db; cc*=(-.5-(k-1))/k; qq+=cc*cp
      spectral=float(np.linalg.norm(qq*np.sqrt(1-vals+c*vals**2)-np.sqrt(target)))
      check(abs(w2-spectral)<3e-14,'full matrix Gaussian W2')
      check(spectral<=2*A**(2*m+2)*math.sqrt(d)+2e-14,'quadratic all-order bound')
      quad.append({'D':d,'A':A,'m':m,'W2':spectral})

# C2, non-C3 ridge potential. h'=sqrt(abs(t))/(1+sqrt(abs(t))).
def h(t):
    u=np.abs(t); v=np.sqrt(u)
    result=u-2*v+2*np.log1p(v)
    small=u<1e-4
    if np.any(small):
      a=u[small]; total=np.zeros_like(a)
      for k in range(1,19): total+=(-1)**(k+1)*a**(1+k/2)/(1+k/2)
      result=np.array(result); result[small]=total
    return np.sign(t)*result

def hp(t):
    v=np.sqrt(np.abs(t)); return v/(1+v)

def psi(t):
    u=np.abs(t); v=np.sqrt(u)
    result=.5*u*u-(4/3)*u*v+2*(u-1)*np.log1p(v)-u+2*v
    small=u<1e-3
    if np.any(small):
      a=u[small]; total=np.zeros_like(a)
      for k in range(1,21): total+=(-1)**(k+1)*a**(2+k/2)/((1+k/2)*(2+k/2))
      result=np.array(result); result[small]=total
    return result

axes=np.array([[1.,0.],[1/math.sqrt(2),1/math.sqrt(2)]])
def gv(x,A):
    dots=x@axes.T
    return A*(.2*x+.3*h(dots)@axes)

def hg(x,A):
    dots=x@axes.T
    ans=np.broadcast_to(.2*np.eye(2),x.shape[:-1]+(2,2)).copy()
    for j in range(2): ans+=.3*hp(dots[...,j])[...,None,None]*np.outer(axes[j],axes[j])
    return A*ans

# Tensor Gaussian integration on all four shared coordinates.
gh,gw=hermgauss(14); gh*=math.sqrt(2); gw/=math.sqrt(math.pi)
mesh=np.meshgrid(*([gh]*4),indexing='ij'); tape=np.stack(mesh,axis=-1).reshape(-1,4)
wm=np.meshgrid(*([gw]*4),indexing='ij'); weights=np.prod(np.stack(wm,axis=-1),axis=-1).reshape(-1)
z=tape[:,:2]; gg=tape[:,2:]
nonlinear=[]
for A in (.5,.25,.125):
    HH=np.zeros_like(z); J=np.zeros((len(z),2,2))
    for ri,si,wi in zip(r,s,w):
        q=ri*z+si*gg
        HH+=wi*gv(q,A); J+=wi*ri*hg(q,A)
    energy=float(weights@np.sum(HH*HH,axis=1))
    jenergy=float(weights@np.sum(J*J,axis=(1,2)))
    linear=float(weights@(np.sum(z*HH,axis=1)-np.trace(J,axis1=1,axis2=2)))
    # GH has a small error at the C2 cusp; use exact quadratic Gaussian IP separately above.
    comm=hg(np.array([[.1,2.],[2.,.1]]),A)
    commnorm=float(np.linalg.norm(comm[0]@comm[1]-comm[1]@comm[0]))
    check(commnorm>1e-6*A*A,'genuine noncommuting Hessians')
    check(np.max(np.linalg.eigvalsh(J))<=A/2+1e-14,'actual Z first')
    check(energy<=A*A*2,'one energy')
    check(jenergy<=A*A*2/4,'one Jacobian energy')
    for t in (.2,.6,1.):
      Y=z-t*HH
      ld=np.linalg.slogdet(np.eye(2)[None,:,:]-t*J)[1]
      kl=float(weights@(-t*np.sum(z*HH,axis=1)+t*t*np.sum(HH*HH,axis=1)/2-ld))
      klcancel=float(weights@(t*t*np.sum(HH*HH,axis=1)/2-ld-t*np.trace(J,axis1=1,axis2=2)))
      bound=.5*t*t*(energy+jenergy/(1-t*A/2))
      check(klcancel>=0,'positive remainder')
      check(klcancel<=bound+1e-14,'matrix logdet remainder bound')
      check(abs(kl-klcancel+t*linear)<1e-14,'KL linear cancellation ledger')
      nonlinear.append({'A':A,'t':t,'KL_integral':kl,'KL_cancelled':klcancel,'bound':bound,'IBP_quadrature_residual':linear,'Hessian_commutator_norm':commnorm})

# Scalar C2 target vs actual positive map, by monotone conditional CDF inversion.
def scalarW2(A):
    delta=A**3
    K=math.ceil(math.log2(4/delta)); m=math.ceil(math.log(16/delta,4))
    rr,ww=rule(K,m); ss=np.sqrt(1-rr*rr)
    xg,xw=hermgauss(64); xg*=math.sqrt(2); xw/=math.sqrt(math.pi)
    ys=np.linspace(-8.5,8.5,6001)
    FF=np.empty(len(ys))
    def force(x): return A*(.3*x+.5*h(x))
    def deriv(x): return A*(.3+.5*hp(x))
    for start in range(0,len(ys),100):
      y=ys[start:start+100,None]; zz=np.broadcast_to(y,(len(y),len(xg))).copy()
      for it in range(6):
        HH=np.zeros_like(zz); DD=np.zeros_like(zz)
        for ri,si,wi in zip(rr,ss,ww):
          q=ri*zz+si*xg[None,:]
          HH+=wi*force(q); DD+=wi*ri*deriv(q)
        step=(zz-HH-y)/(1-DD); zz-=step
      FF[start:start+len(y)]=ndtr(zz)@xw
    pot=.3*ys*ys/2+.5*psi(ys)
    dens=np.exp(-ys*ys/2-A*pot)
    target=cumulative_simpson(dens,x=ys,initial=0); target/=target[-1]
    pp,pw=leggauss(1000); pp=(pp+1)/2; pw/=2
    qy=np.interp(pp,FF,ys); qt=np.interp(pp,target,ys)
    error=float(np.sqrt(pw@((qy-qt)**2)))
    check(np.all(np.diff(FF)>=-1e-14),'actual CDF monotone')
    check(error<=10*(A*A+delta*A),'scalar theorem diagnostic')
    return {'A':A,'delta':delta,'nodes':len(rr),'W2_CDF_diagnostic':error,'W2_over_A2':error/A**2}

scalar=[scalarW2(A) for A in (.25,.125,.0625)]
report={'status':'PASS','assertions':checks,'quadrature_nodes':len(r),'beta':float(beta),'c':float(c),'quadratic':quad,'noncommuting_C2_KL':nonlinear,'scalar_C2_W2':scalar,'scope':'Finite floating-point checks support but do not prove the analytical theorem. CDF inversion uses original firsts only in diagnostics, not in the VALUE producer. GH IBP residual is reported, not silently set to zero.'}
(OUT/'positive_stein_law_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'status':report['status'],'assertions':checks,'scalar':scalar,'nonlinear_max_IBP_quad_error':max(abs(x['IBP_quadrature_residual']) for x in nonlinear)},indent=2))
