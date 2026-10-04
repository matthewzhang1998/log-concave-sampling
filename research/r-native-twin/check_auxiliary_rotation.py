from pathlib import Path
import numpy as np,json
rng=np.random.default_rng(514377);checks=0;worst=0.
def ck(x):
 global checks
 checks+=1
 assert x,checks
def eq(x,y):
 global worst
 err=float(np.linalg.norm(np.asarray(x)-np.asarray(y)));worst=max(worst,err);ck(err<3e-12)
def root(A):
 vals,Q=np.linalg.eigh(A);return (Q*np.sqrt(np.maximum(vals,0)))@Q.T
for n in [1,2,4,7]:
 for rep in range(90):
  I=np.eye(n);O=np.zeros_like(I);a=.05;r=.07;sig=.4;eps=.3;c=.6;s=.8
  M=rng.normal(size=(n,n));M/=max(1,np.linalg.norm(M,2));
  if rep%3==0:M[:,n//2:]=0
  D=root(I-a*a*sig*sig*M@M.T);D2=root(I-a*a*sig*sig*M.T@M)
  inv=np.linalg.inv(I+D)
  rot=np.block([[D,a*sig*M],[-a*sig*M.T,D2]])
  eq(rot@rot.T,np.eye(2*n));eq(D@M,M@D2)
  S,U,Z,V=rng.normal(size=(4,n));SP=D@S+a*sig*M@V;VP=-a*sig*M.T@S+D2@V
  h=sig*V-a*sig*sig*M.T@inv@S
  eq(a*M@h,SP-S);eq(h,sig*VP+a*sig*sig*M.T@inv@SP)
  q1=rng.normal(size=n);q1/=np.linalg.norm(q1);q2=rng.normal(size=n);q2/=np.linalg.norm(q2)
  def g(x):return .5*x+.1*q1*np.sin(q1@x)+.1*q2*np.sin(q2@x)
  def H(x):return .5*I+.1*np.outer(q1,q1)*np.cos(q1@x)+.1*np.outer(q2,q2)*np.cos(q2@x)
  dold=a*(g(S)-g(S+eps*Z));dnew=a*(g(SP)-g(SP+eps*Z))
  x0=c*SP+s*U;x1=x0+M.T@dnew;p0=g(x0)+h;p1=g(x1)+h
  T0=SP+a*M@g(x0);T1=SP+a*M@g(x1);E=r*(g(T0)-g(T1))
  GS=r*(g(S+a*M@p0)-g(S+a*M@p1)+a*(g(S)-g(S+eps*Z))-dold)
  eq(GS,E);eq(r*(g(x0)-p0),-r*h);eq(r*(-g(x1)+p1),r*h)
  eq(r*(x1-x0-M.T@dold),r*M.T@(dnew-dold))
  # Regressed private lambda source is exactly gradient; original M may be singular.
  SS=D@SP-a*sig*M@VP;eq(SS,S)
  Jlam=r*a*a*sig*M.T@(H(SS)-H(SS+eps*Z))@M
  eq(Jlam,Jlam.T);ck(np.linalg.norm(Jlam,2)<=2*r*a*a*sig)
  # Exact fixed-row covariance regression for an arbitrary centered linear E=T S'.
  T=rng.normal(size=(n,n));HR=np.concatenate([-a*sig*sig*M.T@inv,sig*I],axis=1)
  ER=T@np.concatenate([D,a*sig*M],axis=1)
  eq(ER@HR.T,a*sig*sig*T@inv@M)
  # Anisotropic reversible OU Hermite multiplier bound, every tested multiindex.
  ds=np.linalg.eigvalsh(D)
  for _ in range(3):
   deg=rng.integers(0,8,n)
   ck(1-np.prod(ds**deg)<=np.max(1-ds)*np.sum(deg)+1e-12)
out={'status':'PASS','checks':checks,'max_exact_residual':worst,'scope':'Known noncommuting/singular-M rotation, exact graph and defect identities, lambda gradient, regression covariance and OU spectral frame.'}
Path(__file__).with_name('auxiliary_rotation_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
