import numpy as np,json
from pathlib import Path
rng=np.random.default_rng(413145);checks=0

def ck(x):
 global checks
 checks+=1;assert bool(x),checks

def root(B):
 e,Q=np.linalg.eigh((B+B.T)/2);return (Q*np.sqrt(np.maximum(e,0)))@Q.T
for d in [1,2,3,5,8]:
 for rep in range(50):
  Q,_=np.linalg.qr(rng.normal(size=(d,d)));R,_=np.linalg.qr(rng.normal(size=(d,d)));sv=rng.uniform(0,1,d);sv[rng.random(d)<.45]=0;M=Q@np.diag(sv)@R.T
  a=r=rng.uniform(.01,.0625);sigma=rng.uniform(.01,1.5);eps=rng.uniform(.05,1);A=a*sigma*M;D=root(np.eye(d)-A@A.T);D2=root(np.eye(d)-A.T@A);R2=np.linalg.inv(np.eye(d)+D2);b=.5
  zero=np.zeros((d,d));C1=np.hstack([np.eye(d),zero,zero,zero,zero]);C2=np.hstack([np.eye(d),zero,eps*np.eye(d),zero,zero]);C3=np.hstack([D,zero,zero,-A,zero]);C4=np.hstack([D,zero,eps*np.eye(d),-A,zero]);Cs=[C1,C2,C3,C4];sg=[1,-1,-1,1]
  B0=b*np.hstack([M.T,zero,zero,-M.T@A@R2]);F=root(np.eye(d)-B0@B0.T);B=np.hstack([B0,F])
  ck(np.linalg.norm(A@R2@A.T-(np.eye(d)-D))<1e-12)
  ck(np.linalg.norm(B@B.T-np.eye(d))<1e-12)
  ck(np.linalg.norm(sum(s*C for s,C in zip(sg,Cs)))<1e-14)
  for C in Cs:ck(np.linalg.norm(B@C.T-b*M.T)<1e-12)
  W=rng.normal(size=5*d);Sp,U,Z,V,pad=np.split(W,5);S=D@Sp-A@V
  n=rng.normal(size=d);n/=np.linalg.norm(n);H=.5*np.eye(d);freq=4.7
  def g(x):return H@x+.1/freq*np.sin(freq*n@x)*n
  def jac(x):return H+.1*np.cos(freq*n@x)*np.outer(n,n)
  G=r*a/b*sum(s*C.T@g(C@W) for s,C in zip(sg,Cs));L=r*a*M.T@(g(Sp)-g(Sp+eps*Z)-g(S)+g(S+eps*Z))
  ck(np.linalg.norm(B@G-L)<1e-12)
  Hfull=r*a/b*sum(s*C.T@jac(C@W)@C for s,C in zip(sg,Cs))
  ck(np.linalg.norm(Hfull-Hfull.T)<1e-12)
  ck(np.linalg.norm(Hfull,2)<=12*r*a)
  blocks=np.array_split(G*b/(r*a),5)
  ck(np.linalg.norm(blocks[0]-(g(Sp)-g(Sp+eps*Z)-D@(g(S)-g(S+eps*Z))))<1e-12)
  ck(np.linalg.norm(blocks[2]-eps*(g(S+eps*Z)-g(Sp+eps*Z)))<1e-12)
  ck(np.linalg.norm(blocks[3]-A.T@(g(S)-g(S+eps*Z)))<1e-12)
  ck(np.linalg.norm(blocks[1])+np.linalg.norm(blocks[4])<1e-14)
result={'status':'PASS','checks':checks,'seed':413145,'scope':['singular and nonsymmetric M','known common coisometry and fixed-gap fill','exact raw lambda readout','genuine full-gradient Hessian and radius','all actual companion blocks and zero columns']}
Path(__file__).with_name('rotation_lambda_lift_checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
