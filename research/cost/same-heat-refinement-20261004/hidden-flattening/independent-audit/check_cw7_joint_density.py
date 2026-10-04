#!/usr/bin/env python3
"""Deterministic diagnostics for the bounded independent CW7 audit."""
import hashlib, json, math
from pathlib import Path
import numpy as np
from scipy.integrate import quad

HERE=Path(__file__).resolve().parent
checks={}
def record(name, passed, **details):
    checks[name]={"passed":bool(passed), **details}
    if not passed: raise AssertionError(name)
def phi(v,b): return math.exp(-v*v/(2*b))/math.sqrt(2*math.pi*b)
def z(c,b,lam): return math.exp(-lam*c*c/(2*(1+b*lam)))/math.sqrt(1+b*lam)

# Actual CW7 common Gaussian reserve, and reference replacement.
n,a=7,0.1
C=a*np.eye(n)-a/(4*n)*np.ones((n,n))
record('actual_observation_covariance', np.allclose(np.linalg.eigvalsh(C), [3*a/4]+[a]*(n-1)), eigenvalues=np.linalg.eigvalsh(C).tolist())
Cref=a/n*np.ones((n,n))+a*(np.eye(n)-np.ones((n,n))/n)
record('reference_observations_independent', np.allclose(Cref,a*np.eye(n)))
r=np.linspace(-.3,.6,n); q=.2; rb=float(r.mean())
direct=float(np.prod([phi(float(t)-q,a) for t in r]))
formula=n**(-.5)*phi(rb-q,a/n)*(2*math.pi*a)**(-(n-1)/2)*math.exp(-float(np.sum((r-rb)**2))/(2*a))
record('observation_density_jacobian',abs(direct-formula)<1e-12,error=abs(direct-formula))

# Bayes product and stationary convolution, including normalizers.
lam,a,q,y,x=0.7,0.1,0.3,-0.1,0.2
lhs=phi(x-q,a)*phi(y-x,a)
rhs=phi(y-q,2*a)*phi(x-(q+y)/2,a/2)
record('gaussian_product_identity',abs(lhs-rhs)<1e-13)
conv=quad(lambda t: math.exp(-lam*t*t/2)*phi(t-q,a)/z(q,a,lam)*phi(y-t,a),-np.inf,np.inf,epsabs=1e-12)[0]
bayes=phi(y-q,2*a)*z((q+y)/2,a/2,lam)/z(q,a,lam)
record('stationary_normalizer_cancellation',abs(conv-bayes)<1e-11,convolution=conv,formula=bayes)
nonstat=phi(y, a)
record('fixed_predecessor_not_stationary_tilt',abs(nonstat-bayes)>1e-4,nonstationary_Y_density=nonstat,stationary_Y_density=bayes)

# Scores of normalized Gaussian convolution.
c,b,lam=.4,.07,.8
step=1e-5
fd=(math.log(z(c+step,b,lam))-math.log(z(c-step,b,lam)))/(2*step)
mu=c/(1+b*lam)
score=(mu-c)/b
record('partition_score_identity',abs(fd-score)<1e-10,finite_difference=fd,posterior_score=score)

# Distinct hidden laws. Rational fixture is intentionally readable; positive-heat formulas are general.
a=lam=tau2=1.
beta=1/(1+a*lam)
vars3=[a*beta,a*beta+beta**2*tau2,(a+tau2)/(1+lam*(a+tau2))]
record('posterior_of_mean_vs_mixture_vs_tilt',np.allclose(vars3,[.5,.75,2/3]),variances=vars3)

# Exact decoder recurrence, both readable and admissible-small-heat fixtures.
rows=[]
for a in [1.,.1,.001]:
  lam=1.; beta=1/(1+a*lam); alpha=1/(2+a*lam)
  s=a*beta**2
  for k in range(1,7):
    s=alpha**2*s+a*alpha**2+a*alpha
    exact=a*beta-a*a*lam*beta**2*alpha**(2*k)
    rows.append({'a':a,'k':k,'variance':s,'formula':exact,'target':a*beta})
    assert abs(s-exact)<1e-14
    assert s<a*beta
record('finite_decoder_nonstationarity',True,fixtures=rows)
record('small_heat_rational_fixture',abs(rows[6]['variance']-4751/53361)<1e-15)

# Omitting z reweights parent; integrate the one-step parent density as a separate check.
a=lam=1.
Z=quad(lambda r: phi(r,a)*z(r/2,a/2,lam),-np.inf,np.inf,epsabs=1e-12)[0]
varr=quad(lambda r:r*r*phi(r,a)*z(r/2,a/2,lam)/Z,-np.inf,np.inf,epsabs=1e-12)[0]
truex=a/(2+a*lam)+a/(2+a*lam)**2
naivex=3*a/(4+3*a*lam)
record('omitted_normalizer_parent_tilt',abs(varr-6/7)<1e-12 and abs(truex-4/9)<1e-12 and abs(naivex-3/7)<1e-12,parent_variance=varr,true_child_variance=truex,omitted_z_child_variance=naivex)

# Path base precision and ideal quadratic joint Hessian, no initial random state.
pathrows=[]
for k in [1,2,5,16,40]:
  a,L=.1,1.
  S=np.zeros((k,k))
  for i in range(1,k):S[i,i-1]=1.
  I=np.eye(k)
  # Coordinates w=(r,x); u=x-(Sx+r)/2, c=(Sx+r)/2.
  R=np.hstack([I,np.zeros((k,k))]); X=np.hstack([np.zeros((k,k)),I])
  U=np.hstack([-.5*I,I-.5*S]); C=np.hstack([.5*I,.5*S])
  B=(R.T@R+2*U.T@U)/a
  base=float(np.linalg.eigvalsh(B)[0])
  ideal=B+L*(X.T@X)-L/(1+a*L/2)*(C.T@C)
  mineig=float(np.linalg.eigvalsh(ideal)[0]); bound=1/(4*a)-L/2
  assert base>=1/(4*a)-1e-10 and mineig>=bound-1e-10
  assert np.linalg.norm(C,2)**2<=.5+1e-12
  pathrows.append({'k':k,'base_min_eigenvalue':base,'ideal_min_eigenvalue':mineig,'proved_lower_bound':bound})
record('affine_decoder_path_convexity',True,fixtures=pathrows)

# General matrix quadratic decoder and polynomial actions.
rng=np.random.default_rng(271828)
O,_=np.linalg.qr(rng.normal(size=(5,5))); A=O@np.diag([.1,.3,.5,.8,1.])@O.T
a=.08; I=np.eye(5); B=np.linalg.inv(I+a*A); D=np.linalg.inv(2*I+a*A)
S=a*B@B
for k in range(1,9):
  S=D@S@D+a*D@D+a*D
  Dk=np.linalg.matrix_power(D,k)
  formula=a*B+Dk@(a*B@B-a*B)@Dk
  assert np.linalg.norm(S-formula)<1e-14
record('matrix_quadratic_decoder_covariance',True,last_covariance=S.tolist())

b=.08; c=rng.normal(size=5); h=rng.normal(size=5); g=rng.normal(size=5); F=A@c+h
m=5; invpoly=np.eye(5); sqrtpoly=np.eye(5); power=np.eye(5); cj=1.
for j in range(1,m+1):
  power=power@(b*A)
  invpoly+=(-1)**j*power
  cj*= (2*j-1)/(2*j)
  sqrtpoly+=(-1)**j*cj*power
vals,vec=np.linalg.eigh(A); exactinv=np.linalg.inv(I+b*A); exactsqrt=(vec*(1+b*vals)**(-.5))@vec.T
meanapprox=c-b*invpoly@F; meanexact=exactinv@(c-b*h)
meanerr=float(np.linalg.norm(meanapprox-meanexact)); meanbound=b*b**(m+1)*np.linalg.norm(F)
noiseerr=float(np.linalg.norm(math.sqrt(b)*(sqrtpoly-exactsqrt),ord='fro')); noisebound=math.sqrt(5*b)*b**(m+1)/(1-b)
record('matrix_free_quadratic_polynomial_actions',meanerr<=meanbound+1e-14 and noiseerr<=noisebound+1e-14,degree=m,mean_error=meanerr,mean_bound=float(meanbound),gaussian_L2_noise_error=noiseerr,noise_bound=noisebound)

# Actual nonlinear-center one-edge curvature: smooth bounded-Hessian convex V.
# V(t)=lambda*t^2/2 + epsilon*(1-cos(t)); c(h)=-tau F(sqrt(a)*h).
lam,eps,tau,a,b=.6,.2,.05,.1,.05
hh=math.pi/(2*math.sqrt(a))
cc=-tau*(lam*math.sqrt(a)*hh+eps*math.sin(math.sqrt(a)*hh))
cp=-tau*math.sqrt(a)*(lam+eps*math.cos(math.sqrt(a)*hh))
cpp=tau*a*eps*math.sin(math.sqrt(a)*hh)
def weight(x):return math.exp(-lam*x*x/2-eps*(1-math.cos(x)))*phi(x-cc,b)
zz=quad(weight,-np.inf,np.inf,epsabs=1e-12)[0]
mu=quad(lambda x:x*weight(x),-np.inf,np.inf,epsabs=1e-12)[0]/zz
var=quad(lambda x:(x-mu)**2*weight(x),-np.inf,np.inf,epsabs=1e-12)[0]/zz
xx=mu+b/cpp*(2+cp*cp*var/(b*b))
curv=1+cp*cp*var/(b*b)+cpp*(mu-xx)/b
record('nonlinear_center_can_destroy_joint_convexity',curv<-.99,ancestor_h=hh,child_x=xx,second_directional_derivative=curv,center_second_derivative=cpp)

payload={'scope':'Deterministic identity checks, not a hidden-sampler or source-graph certificate.','all_passed':all(v['passed'] for v in checks.values()),'checks':checks}
(HERE/'cw7_joint_density_checks.json').write_text(json.dumps(payload,indent=2)+'\n')
print(json.dumps({'all_passed':payload['all_passed'],'checks':len(checks),'output':'cw7_joint_density_checks.json'},indent=2))
