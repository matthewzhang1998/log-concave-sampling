from pathlib import Path
import numpy as np,json
rng=np.random.default_rng(19101); checks=0; worst=0.0

def ck(x):
 global checks
 checks+=1
 if not x: raise AssertionError(checks)
def eq(x,y):
 global worst
 er=float(np.linalg.norm(x-y));worst=max(worst,er);ck(er<2e-12)
for d in [1,2,4,7]:
 for rep in range(35):
  I=np.eye(d);c=s=2**-.5; C=np.concatenate([I,np.zeros_like(I)],axis=1);P=np.concatenate([c*I,s*I],axis=1)
  eps=.3;a=.03;r=.07
  M=rng.normal(size=(d,d));M/=max(1,np.linalg.norm(M,2))
  v=rng.normal(size=d);v/=np.linalg.norm(v);u=rng.normal(size=d);u/=np.linalg.norm(u)
  def g(x):return .5*x+.1*v*np.sin(v@x)+.05*u*np.sin(u@x)
  def H(x):return .5*I+.1*np.outer(v,v)*np.cos(v@x)+.05*np.outer(u,u)*np.cos(u@x)
  a0=rng.normal(size=d);at=rng.normal(size=d)
  def g0(x):return g(x+a0)-g(a0)
  def gt(x):return g(x+at)-g(at)
  def H0(x):return H(x+a0)
  def Ht(x):return H(x+at)
  W=rng.normal(size=2*d);Z=rng.normal(size=d);S=P@W;X=C@W
  F=r*gt(S+a*M@g0(X)); p0p=np.zeros(d);ptp=np.zeros(d);p0m=np.zeros(d);ptm=np.zeros(d)
  vals=[]
  for k in range(1,9):
   p0p,ptp,p0m,ptm=g0(X+a*M.T@(ptp+ptm)),gt(S+a*M@p0p),-g0(X),-gt(S+a*M@p0p+eps*Z)
   vals.append(r*ptp)
   if k>=3:ck(np.linalg.norm(F-r*ptp)<=r*a*a*eps*np.linalg.norm(Z)+2e-12)
  eq(vals[1],F)
  Delta=gt(S)-gt(S+eps*Z);q0=X+a*M.T@Delta;qt=S+a*M@g0(q0)
  E=F-r*gt(qt);eq(vals[2],r*gt(qt))
  dq0=C+a*M.T@(Ht(S)-Ht(S+eps*Z))@P
  JH=r*Ht(qt)@(P+a*M@H0(q0)@dq0)
  JF=r*Ht(S+a*M@g0(X))@(P+a*M@H0(X)@C)
  JZE=r*a*a*eps*Ht(qt)@M@H0(q0)@M.T@Ht(S+eps*Z)
  JE=np.concatenate([JF-JH,JZE],axis=1)
  Prec=np.concatenate([P,np.zeros_like(I)],axis=1)
  curl=Prec.T@JE-JE.T@Prec
  ck(np.linalg.norm(JE,2)<=3*r);ck(np.linalg.norm(curl,2)<=5*r*a)
  ck(np.linalg.norm(JZE,2)<=r*a*a*eps)
  beta=(s**-2+eps**-2)**-.5
  B=beta*np.concatenate([np.zeros_like(I),I/s,-I/eps],axis=1)
  eq(B@B.T,I);eq(Prec@B.T,beta*I)
  Cdir=np.concatenate([C,P],axis=0);et=np.concatenate([np.zeros_like(I),I],axis=0)
  Ctw=np.block([[Cdir,np.zeros((2*d,d))],[Cdir,eps*et]])
  sel=np.zeros((d,4*d));sel[:,d:2*d]=I
  eq(B@Ctw.T,beta*sel)
  # Same-potential affine regression; no force readout sees Z in its forward P row.
  HH=.5*I+.1*np.outer(v,v)
  FF=r*HH@(S+a*M@HH@X)
  HH3=r*HH@(S+a*M@HH@(X+a*M.T@(HH@S-HH@(S+eps*Z))))
  Q=r*eps*HH@(a*M)@HH@(a*M).T@HH
  eq(FF-HH3,Q@Z)
  Jaff=np.concatenate([np.zeros((d,2*d)),Q],axis=1)
  eq(Jaff@Prec.T,np.zeros_like(I))
out={'status':'PASS','checks':checks,'max_identity_error':worst,'scope':'New literal K1-K8 VALUE iteration, same-potential nonlinear first/curl/column checks and affine copied-twin regression. Does not realize the missing feedback heat repair.'}
Path(__file__).with_name('two_node_third_twin_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
