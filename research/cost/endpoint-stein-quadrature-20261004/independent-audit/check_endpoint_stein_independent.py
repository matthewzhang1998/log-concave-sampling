#!/usr/bin/env python3
"""Independent finite diagnostics; analytical proof is in the companion audit.
No imports from or calls to the author's diagnostic implementation.
"""
from pathlib import Path
import hashlib, json, math
import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite_e import hermegauss
from scipy.integrate import quad

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'ENDPOINT-STEIN-PACKET-AND-POLYLOG-QUADRATURE.md'
rng=np.random.default_rng(202610041935)
checks=0
records={}
def check(ok, msg):
    global checks
    checks+=1
    if not bool(ok): raise AssertionError(msg)
def close(a,b,atol=1e-11,rtol=1e-10,msg='close'):
    check(np.allclose(a,b,atol=atol,rtol=rtol),msg)
def op(a): return np.linalg.norm(a,2)
def quadrature(K,m):
    x,w=leggauss(m)
    ts=[]; ws=[]
    for k in range(K):
        a=2.0**(-k-1)
        ts.extend(a*(1.5+.5*x)); ws.extend(a*.5*w)
    h=2.0**(-K)
    ts.append(h/2); ws.append(h)
    return np.array(ts),np.array(ws)
def errors(t,w,degrees):
    # Compute from t, avoiding avoidable cancellation in forming 1-t.
    q=np.exp(np.asarray(degrees)[:,None]*np.log1p(-t)[None,:])@w
    return q-1/(np.asarray(degrees)+1)

# Ellipse geometry and moment integration, independently sampled through 10^12.
quad_rows=[]
for K in [1,2,4,8,12,20,28]:
  for m in [1,2,3,5,8,12]:
    t,w=quadrature(K,m)
    check(len(t)==K*m+1,'node count')
    check(np.all(w>0) and np.all((t>0)&(t<1)),'positive nodes and weights')
    close(w.sum(),1,atol=5e-15,rtol=0,msg='mass')
    close(w@(1-t),.5,atol=5e-15,rtol=0,msg='first moment')
    deg=np.unique(np.r_[np.arange(0,1200),np.round(np.geomspace(1200,1e12,650)),np.round(np.geomspace(.01,100,150)*2.0**K)]).astype(np.float64)
    err=errors(t,w,deg)
    bound=8*4.0**(-m)+2*2.0**(-K)
    for val in err:
      check(abs(val)<=bound+8e-15,'uniform-moment sampled bound')
      check(val<=8e-15,'positive Gauss plus midpoint underestimates integer moments')
    # Exact midpoint defect for the quadratic moment when panels integrate it.
    if m>=2:
      close(errors(t,w,[2])[0],-2.0**(-3*K)/12,atol=3e-15,rtol=1e-8,msg='terminal quadratic defect')
    quad_rows.append(dict(K=K,m=m,nodes=len(t),sampled_degrees=len(deg),max_degree=float(deg.max()),max_abs_error=float(abs(err).max()),worst_degree=float(deg[np.argmax(abs(err))]),certified_bound=bound))
for a in [.5,.25,2.0**-12,2.0**-50]:
  for theta in np.linspace(0,2*math.pi,1001):
    t=1.5*a+.625*a*math.cos(theta)+1j*.375*a*math.sin(theta)
    check(abs(1-t)<=1-7*a/8+4e-16,'rho=2 ellipse bound')
records['quadrature']=quad_rows

# Spectral multiplier acting on vector chaoses, with no dimension factor.
hermite_rows=[]
for D in [1,3,17,257]:
  for K,m in [(2,1),(5,3),(12,6)]:
    t,w=quadrature(K,m)
    degrees=np.arange(501,dtype=float)
    e=errors(t,w,degrees)
    coefficients=rng.normal(size=(len(degrees),D))
    lhs=np.sum((e[:,None]*coefficients)**2)
    rhs=float(np.max(abs(e))**2*np.sum(coefficients**2))
    check(lhs<=rhs*(1+1e-14),'vector Hermite contraction')
    target=np.zeros_like(coefficients); target[np.argmax(abs(e)),:]=rng.normal(size=D)
    close(np.sum((e[:,None]*target)**2),np.max(abs(e))**2*np.sum(target**2),atol=1e-18,rtol=3e-14,msg='spectral norm attainment')
    hermite_rows.append(dict(D=D,K=K,m=m,relative_squared_ratio=float(lhs/rhs)))
records['hermite']=hermite_rows

class Potential:
  """Noncommuting smooth Hessians, each lying in [0.15 I,0.85 I]."""
  def __init__(self,D):
    self.D=D
    q,_=np.linalg.qr(rng.normal(size=(D,D)))
    self.B=q@np.diag(np.linspace(.35,.65,D))@q.T
    self.u=rng.normal(size=(4,D)); self.u/=np.linalg.norm(self.u,axis=1)[:,None]
    self.c=np.array([.05,.05,.05,.05]); self.f=np.array([.13,1.7,11,83.])
    self.linear=rng.normal(size=D)*3
  def grad(self,x):
    return self.B@x+self.linear+((self.c/self.f)*np.sin(self.f*(self.u@x)))@self.u
  def hess(self,x):
    return self.B+self.u.T@np.diag(self.c*np.cos(self.f*(self.u@x)))@self.u

def packet(V,A,M,t,w,y,Z,G,with_jac=True):
  D=len(y); I=np.eye(D); b=y.copy(); B=I.copy(); sites=[]
  for _ in range(M):
    sites.append(b.copy())
    if with_jac: B=I-A*V.hess(b)@B
    b=y-A*V.grad(b)
  anchor=V.grad(b); sites.append(b.copy())
  H=np.zeros(D); JW=np.zeros((D,2*D)); JH_y=np.zeros((D,D))
  for ti,wi in zip(t,w):
    r=1-ti; s=math.sqrt(ti*(2-ti)); R=np.hstack((r*I,s*I)); q=r*Z+s*G
    x=b+math.sqrt(A)*q; sites.append(x.copy())
    H+=wi*math.sqrt(A)*(V.grad(x)-anchor)
    if with_jac:
      Hx=V.hess(x)
      JW+=wi*A*Hx@R
      JH_y+=wi*math.sqrt(A)*(Hx-V.hess(b))@B
  X=b+math.sqrt(A)*(Z-H); sites.append(X.copy())
  F=math.sqrt(A)*V.grad(X)
  result=dict(H=H,X=X,F=F,b=b,sites=sites)
  if with_jac:
    DXW=math.sqrt(A)*(np.hstack((I,np.zeros_like(I)))-JW)
    DXy=B-math.sqrt(A)*JH_y
    result.update(JHW=JW,DXW=DXW,DXy=DXy,DFW=math.sqrt(A)*V.hess(X)@DXW,DFy=math.sqrt(A)*V.hess(X)@DXy,B=B)
  return result

packet_rows=[]
for D in [1,3,11]:
  V=Potential(D)
  for A in [.5,.2,.03,.001]:
    for M in [0,1,4,9]:
      t,w=quadrature(6,4); n=len(t)
      for _ in range(3):
        y=rng.normal(size=D)*5; Z=rng.normal(size=D); G=rng.normal(size=D)
        p=packet(V,A,M,t,w,y,Z,G)
        check(len(p['sites'])==M+n+2,'complete original VALUE count')
        check(op(p['JHW'])<=A*(1+2e-14),'private H first')
        check(op(p['DXW']/math.sqrt(A)-np.hstack((np.eye(D),np.zeros((D,D)))))<=A*(1+2e-14),'leading private row')
        check(op(p['B'])<=1/(1-A)+1e-14,'finite mode caller')
        check(op(p['DXy']-np.eye(D))<=2*A/(1-A)+1e-14,'output caller')
        check(op(p['DFW'])<=A*(1+A)+1e-14,'force private first')
        check(op(p['DFy'])<=math.sqrt(A)*(1+2*A/(1-A))+1e-14,'force caller first')
        residual=p['b']-y+A*V.grad(p['b'])
        check(np.linalg.norm(residual)<=A**(M+1)*np.linalg.norm(V.grad(y))+3e-14,'finite mode residual')
        origin=packet(V,A,M,t,w,y,np.zeros(D),np.zeros(D))
        check(np.array_equal(origin['H'],np.zeros(D)),'literal anchored zero')
        check(np.array_equal(origin['X'],origin['b']),'physical origin')
        # All physical query positions have the caller profile, pathwise.
        b_bound=A/(1-A)*np.linalg.norm(V.grad(y))
        root=math.sqrt(np.dot(Z,Z)+np.dot(G,G))
        for x in p['sites']:
          check(np.linalg.norm(x-y)<=b_bound+math.sqrt(A)*(1+A)*root+1e-12,'query displacement profile')
        if M in [1,4] and A in [.2,.03]:
          dy=rng.normal(size=D); dy/=np.linalg.norm(dy)
          dz=rng.normal(size=D); dg=rng.normal(size=D)
          h=2e-7
          plus=packet(V,A,M,t,w,y+h*dy,Z,G,False)['X']
          minus=packet(V,A,M,t,w,y-h*dy,Z,G,False)['X']
          close((plus-minus)/(2*h),p['DXy']@dy,atol=2e-8,rtol=2e-7,msg='finite-difference caller')
          plus=packet(V,A,M,t,w,y,Z+h*dz,G+h*dg,False)['X']
          minus=packet(V,A,M,t,w,y,Z-h*dz,G-h*dg,False)['X']
          close((plus-minus)/(2*h),p['DXW']@np.r_[dz,dg],atol=2e-8,rtol=2e-7,msg='finite-difference private')
      packet_rows.append(dict(D=D,A=A,M=M,nodes=n,value_count=M+n+2))
records['finite_packet']=packet_rows

# A genuinely C2 potential with a non-Lipschitz Hessian at zero.
class C2Potential:
  def grad(self,x):
    u=abs(x)**.2
    return .5*x+.2*x*u/(1+u)+.7
  def hess(self,x):
    u=abs(x)**.2
    return np.diag(.5+.2*(u/(1+u)+.2*u/(1+u)**2))
V=C2Potential()
t,w=quadrature(7,3)
for yval in [0,1e-12,-1e-9,.2,-2.]:
  p=packet(V,.1,4,t,w,np.array([yval]),np.array([.7]),np.array([-.3]))
  check(op(p['JHW'])<=.1,'C2 only private first')
  check(op(p['DXy']-np.eye(1))<=2*.1/.9,'C2 only caller first')
records['c2_nonsmooth_hessian_cases']=5

# Explicit finite-mode target residual, exactly soluble for quadratics.
mode_rows=[]
for D in [1,3,23]:
  q,_=np.linalg.qr(rng.normal(size=(D,D))); H=q@np.diag(np.linspace(0,1,D))@q.T
  h=rng.normal(size=D); y=rng.normal(size=D)
  for A in [.5,.1,.001]:
    exact=np.linalg.solve(np.eye(D)+A*H,y-A*h)
    b=y.copy()
    for M in range(10):
      res=b-y+A*(H@b+h); ell=res/math.sqrt(A)
      displacement=np.linalg.solve(np.eye(D)+A*H,res)
      close(b-exact,displacement,atol=2e-14,rtol=1e-7,msg='mode target shift identity')
      check(np.linalg.norm(b-exact)<=math.sqrt(A)*np.linalg.norm(ell)+2e-14,'physical target residual bound')
      check(np.linalg.norm(res)<=A**(M+1)*np.linalg.norm(H@y+h)+2e-14,'mode residual power')
      mode_rows.append(dict(D=D,A=A,M=M,mode_distance=float(np.linalg.norm(b-exact)),physical_residual_bound=float(np.linalg.norm(res))))
      b=y-A*(H@b+h)
records['mode_residual']=mode_rows

# Gaussian covariance ledger retains all shared-root cross terms.
quadratic_rows=[]
for K,m in [(1,1),(3,2),(12,8),(24,12)]:
  t,w=quadrature(K,m); r=1-t; s=np.sqrt(t*(2-t)); beta=float(w@s)
  for D in [1,5,41]:
    q,_=np.linalg.qr(rng.normal(size=(D,D)))
    B=q@np.diag(np.linspace(.001,.2,D))@q.T; I=np.eye(D)
    actual=(I-B/2)@(I-B/2).T+beta**2*B@B
    ledger=I-B+(.25+beta**2)*B@B
    close(actual,ledger,atol=5e-15,rtol=1e-13,msg='same-root covariance')
    check(op(actual-np.linalg.inv(I+B))>1e-8,'uncancelled target covariance current')
    wrong=I-B+(.25+float((w*w)@(s*s)))*B@B
    check(op(actual-wrong)>1e-8,'independent node roots change law')
  quadratic_rows.append(dict(K=K,m=m,beta=beta,coefficient=.25+beta**2,wrong_independent_coefficient=.25+float((w*w)@(s*s))))
close(quadratic_rows[-1]['beta'],math.pi/4,atol=1e-10,rtol=0,msg='beta tends to pi/4')
# Weighted retained-coordinate factoring and the nonvanishing strong bridge.
weighted_rows=[]
for K,m in [(2,1),(8,4),(20,8)]:
  t,w=quadrature(K,m); r=1-t; s=np.sqrt(t*(2-t)); d=np.sqrt(w)
  stacked=np.column_stack((d*r,d*s))
  check(op(stacked)<=1+5e-15,'weighted full Gaussian stack')
  close(np.linalg.norm(d),1,atol=3e-15,rtol=0,msg='terminal weighted row')
  # Conservative Markov bound on Legendre derivative gives w_i >= 2^-K/m^4.
  check(w.min()>=2.0**(-K)/m**4*(1-1e-13),'finite inverse-width precision bound')
  beta=float(w@s); B=.1
  strong_bridge_l2=beta*B
  check(strong_bridge_l2>.06,'conditional-mean precision does not suppress bridge variance')
  weighted_rows.append(dict(K=K,m=m,stack_operator_norm=float(op(stacked)),max_inverse_width=float(1/d.min()),scalar_B=B,strong_bridge_l2=strong_bridge_l2))
records['weighted_graph']=weighted_rows
records['quadratic']=quadratic_rows
records['limit_quadratic_coefficient']=.25+math.pi**2/16

# Scalar nonlinear Stein and literal weak Taylor identity using independent GH
# integration. The quadrature difference is kept, rather than declared zero.
z,hg=hermegauss(72); hg=hg/math.sqrt(2*math.pi)
ZZ=z[:,None]; GG=z[None,:]; GW=hg[:,None]*hg[None,:]
sx,sw=leggauss(40); st=(sx+1)/2; sw=sw/2
weak_rows=[]
for K,m in [(2,1),(6,3)]:
  t,w=quadrature(K,m); r=1-t; s=np.sqrt(t*(2-t))
  A=.17; lam=.5; eta=.3; omega=1.4; kappa=.9
  Hval=np.zeros_like(GW)
  for ri,si,wi in zip(r,s,w):
    q=ri*ZZ+si*GG
    Hval+=wi*A*(lam*q+eta*np.sin(omega*q)/omega)
  cov=A*(-lam*kappa*kappa*math.exp(-kappa*kappa/2)/2-eta/(omega*omega)*(.5*(math.exp(-(omega-kappa)**2/2)+math.exp(-(omega+kappa)**2/2))-math.exp(-(omega*omega+kappa*kappa)/2)))
  integrand=lambda rr: -.5*A*eta*kappa/omega*(math.exp(-.5*(omega*omega+kappa*kappa)+omega*kappa*rr)-math.exp(-.5*(omega*omega+kappa*kappa)-omega*kappa*rr))
  reference=-A*lam*kappa*kappa*math.exp(-kappa*kappa/2)/2+quad(integrand,0,1,epsabs=1e-14,epsrel=1e-13)[0]
  close(reference,cov,atol=5e-15,rtol=1e-13,msg='analytical nonlinear Stein identity')
  observed_first=float(np.sum(GW*Hval*(-kappa*np.sin(kappa*ZZ))))
  q_first=-A*lam*kappa*kappa*math.exp(-kappa*kappa/2)/2+sum(wi*integrand(ri) for ri,wi in zip(r,w))
  close(observed_first,q_first,atol=5e-15,rtol=2e-13,msg='conditional Gaussian mean first term')
  remainder=sum(wi*(1-si)*float(np.sum(GW*Hval*Hval*(-kappa*kappa*np.cos(kappa*(ZZ-si*Hval))))) for si,wi in zip(st,sw))
  lhs=float(np.sum(GW*np.cos(kappa*(ZZ-Hval))))
  rhs=math.exp(-kappa*kappa/2)-cov-(q_first-reference)+remainder
  close(lhs,rhs,atol=2e-14,rtol=2e-13,msg='literal weak Taylor identity with quadrature defect')
  weak_rows.append(dict(K=K,m=m,covariance=cov,quadrature_first_defect=q_first-reference,uncancelled_taylor_current=remainder,weak_identity_residual=lhs-rhs))
records['nonlinear_weak_identity']=weak_rows

# Large-D radial limiting obstruction for one uniform time; finite 1D integrals.
radial=[]
for a in [.1,.03,.01,.003,.001]:
  tau=1/(1+a)
  value=quad(lambda rr:(math.sqrt(1-2*a*rr+a*a)-math.sqrt(tau))**2,0,1,epsabs=1e-17)[0]
  ratio=math.sqrt(value)/a
  check(.27<ratio<.31,'uniform-time radial O(a), not O(a^2), coefficient')
  radial.append(dict(a=a,limiting_W2_per_sqrtD=math.sqrt(value),ratio_to_a=ratio,ratio_to_a_squared=math.sqrt(value)/(a*a)))
records['uniform_time_radial_limit']=radial

# Common biased force errors sum by absolute weights. Independent numerical
# evaluations at the anchor are intentionally demonstrated not to cancel.
t,w=quadrature(8,4); bias=.003
close(w@np.full(len(w),bias),bias,atol=1e-16,rtol=0,msg='numerical bias no sqrt(n) gain')
check(abs(math.sqrt(.1)*(bias-(-bias)))>0,'independently biased zero is not exact')
records['numerical_zero']={'same_version_identical_site':'exact subtraction', 'different_bias_residual_example':math.sqrt(.1)*2*bias}

out={'status':'PASS','checks':checks,'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'seed':202610041935,'scope':'Finite diagnostics support, but do not replace, the analytical audit. No source-order or complete-cost theorem tested.', 'results':records}
(HERE/'endpoint_stein_independent_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'checks':checks,'source_sha256':out['source_sha256'],'max_nodes':max(x['nodes'] for x in quad_rows),'max_sampled_degree':max(x['max_degree'] for x in quad_rows),'limit_quadratic_coefficient':records['limit_quadratic_coefficient'],'weak_rows':weak_rows},indent=2))
