from pathlib import Path
import numpy as np,json
rng=np.random.default_rng(740891)
checks=0;worst=0.
def ck(x):
 global checks
 checks+=1
 assert x,checks
def eq(x,y):
 global worst
 err=float(np.linalg.norm(x-y));worst=max(worst,err);ck(err<1e-12)
for n in [1,2,3,6]:
 for rep in range(100):
  I=np.eye(n); O=np.zeros((n,n)); a=.03;r=.07;eps=.4;c=.6;s=.8
  M=rng.normal(size=(n,n));M/=max(1,np.linalg.norm(M,2))
  v=rng.normal(size=n);v/=np.linalg.norm(v)
  u=rng.normal(size=n);u/=np.linalg.norm(u)
  def g(x):return .4*x+.1*v*np.sin(v@x)+.08*u*np.sin(u@x)
  def H(x):return .4*I+.1*np.outer(v,v)*np.cos(v@x)+.08*np.outer(u,u)*np.cos(u@x)
  anchor0=rng.normal(size=n);anchort=rng.normal(size=n)
  def g0(x):return g(x+anchor0)-g(anchor0)
  def gt(x):return g(x+anchort)-g(anchort)
  def H0(x):return H(x+anchor0)
  def Ht(x):return H(x+anchort)
  S,U,Z=rng.normal(size=(3,n));x=c*S+s*U;Delta=gt(S)-gt(S+eps*Z)
  x1=x+a*M.T@Delta;t0=S+a*M@g0(x);t1=S+a*M@g0(x1)
  K0=Ht(t0)@M@H0(x);K1=Ht(t1)@M@H0(x1);DK=K0-K1
  D=Ht(S)-Ht(S+eps*Z);Hp=Ht(S+eps*Z)
  J_S=r*(Ht(t0)-Ht(t1))+r*a*c*DK-r*a*a*K1@M.T@D
  J_U=r*a*s*DK;J_Z=r*a*a*eps*K1@M.T@Hp
  J=np.block([[J_S,J_U,J_Z],[O,O,O],[O,O,O]])
  lead=r*a*np.block([[c*(DK-DK.T),s*DK,O],[-s*DK.T,O,O],[O,O,O]])
  fb=r*a*a*np.block([[-K1@M.T@D+D@M@K1.T,O,eps*K1@M.T@Hp],[O,O,O],[-eps*Hp@M@K1.T,O,O]])
  eq(J-J.T,lead+fb);ck(np.linalg.norm(fb,2)<=5*r*a*a)
  d=(a/s)*M.T@Delta
  def Phi(z):return gt(S+a*M@g0(c*S+s*z))
  eq(r*(gt(t0)-gt(t1)),r*(Phi(U)-Phi(U+d)))
  # Joint-gradient lift, with its companion; compare its exact symmetric block.
  L1=np.concatenate([I,O],axis=1);L2=np.concatenate([I,eps*I],axis=1)
  HG=a*(L1.T@Ht(S)@L1-L2.T@Hp@L2)
  eq(HG,HG.T);ck(np.linalg.norm(HG,2)<=3*a)
  GD=a*np.concatenate([Delta,-eps*gt(S+eps*Z)])
  ck(np.linalg.norm(GD)<=a*eps*(np.linalg.norm(Z)+np.linalg.norm(S+eps*Z))+1e-12)
  eq(L1@L2.T,I);eq(L2@L2.T,(1+eps*eps)*I)
  # The first-only formula and a finite-difference check for this smooth fixture.
  W=np.concatenate([S,U,Z]);direction=rng.normal(size=3*n);hh=1e-5
  def Ef(w):
   ss,uu,zz=np.split(w,3);xx=c*ss+s*uu;dd=gt(ss)-gt(ss+eps*zz)
   return r*(gt(ss+a*M@g0(xx))-gt(ss+a*M@g0(xx+a*M.T@dd)))
  fd=(Ef(W+hh*direction)-Ef(W-hh*direction))/(2*hh)
  ck(np.linalg.norm(fd-np.concatenate([J_S,J_U,J_Z],axis=1)@direction)<1e-9)
out={'status':'PASS','checks':checks,'max_exact_residual':worst,'scope':'Exact K3 leading/subleading curl split and joint-gradient feedback source. No missing joint observer calibration is tested.'}
Path(__file__).with_name('leading_curl_split_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
