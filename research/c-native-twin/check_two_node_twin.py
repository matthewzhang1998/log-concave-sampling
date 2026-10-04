# Publication copy; original/public SHA-256 mapping is in INVENTORY.json.
from pathlib import Path
import numpy as np,json
rng=np.random.default_rng(83406);checks=0;maxerr=0.
for d in range(1,7):
 for rep in range(30):
  Q0=np.linalg.qr(rng.normal(size=(d,d)))[0]; Qt=np.linalg.qr(rng.normal(size=(d,d)))[0]
  def G(Q,x):return .4*x+.2*Q.T@np.sin(Q@x)
  def H(Q,x):return .4*np.eye(d)+.2*Q.T@np.diag(np.cos(Q@x))@Q
  c=np.sqrt(.3);s=np.sqrt(.7);a=rng.uniform(.001,.0625);eps=rng.uniform(.001,1);r=.03
  M=rng.normal(size=(d,d));M/=max(1,np.linalg.norm(M,2))
  C=np.concatenate([np.eye(d),np.zeros((d,2*d))],axis=1)
  P=np.concatenate([c*np.eye(d),s*np.eye(d),np.zeros((d,d))],axis=1)
  Zrow=np.concatenate([np.zeros((d,2*d)),np.eye(d)],axis=1)
  V=rng.normal(size=3*d);x=C@V;S=P@V;z=Zrow@V
  p=np.zeros(4*d); Ct=np.vstack([C,P,C,P+eps*Zrow])
  Sd=np.zeros((4*d,4*d));Sd[:d,d:2*d]=a*M.T;Sd[:d,3*d:]=a*M.T;Sd[d:2*d,:d]=a*M;Sd[3*d:,:d]=a*M
  Jp=np.zeros((4*d,3*d))
  for K in range(3):
   args=Ct@V+Sd@p
   fs=[G(Q0,args[:d]),G(Qt,args[d:2*d]),-G(Q0,args[2*d:3*d]),-G(Qt,args[3*d:])]
   Hs=[H(Q0,args[:d]),H(Qt,args[d:2*d]),-H(Q0,args[2*d:3*d]),-H(Qt,args[3*d:])]
   Hblk=np.zeros((4*d,4*d))
   for j,Hj in enumerate(Hs):Hblk[j*d:(j+1)*d,j*d:(j+1)*d]=Hj
   p=np.concatenate(fs);Jp=Hblk@(Ct+Sd@Jp)
  delta=G(Qt,S)-G(Qt,S+eps*z);qfb=x+a*M.T@delta
  qF=S+a*M@G(Q0,x);qH=S+a*M@G(Q0,qfb)
  E=r*(G(Qt,qF)-G(Qt,qH));Eiter=r*G(Qt,qF)-r*p[d:2*d]
  err=np.linalg.norm(E-Eiter);maxerr=max(maxerr,err);assert err<1e-13;checks+=1
  assert np.linalg.norm(E)<=r*a*a*eps*np.linalg.norm(z)+1e-13;checks+=1
  FJ=r*H(Qt,qF)@(P+a*M@H(Q0,x)@C)
  EJ=FJ-r*Jp[d:2*d]
  expected=r*a*a*eps*H(Qt,qH)@M@H(Q0,qfb)@M.T@H(Qt,S+eps*z)
  err=np.linalg.norm(EJ@Zrow.T-expected);maxerr=max(maxerr,err);assert err<1e-13;checks+=1
  Curl=P.T@EJ-EJ.T@P;assert np.linalg.norm(Curl,2)<=5*r*a+1e-13;checks+=1
  beta=(s**-2+eps**-2)**-.5;B=beta*np.concatenate([np.zeros((d,d)),np.eye(d)/s,-np.eye(d)/eps],axis=1)
  Gfin=r/beta*Ct.T@p
  err=np.linalg.norm(B@Gfin-r*p[d:2*d]);maxerr=max(maxerr,err);assert err<1e-12;checks+=1
  assert np.linalg.norm(B@B.T-np.eye(d))<1e-12
  assert np.linalg.norm(P@B.T-beta*np.eye(d))<1e-12;checks+=2
  U0=s*x-c*(V[d:2*d]); kS=a*delta
  def Phi(u):return G(Qt,S+a*M@G(Q0,c*S+s*u))
  err=np.linalg.norm(E-r*(Phi(U0)-Phi(U0+M.T@kS/s)));maxerr=max(maxerr,err);assert err<1e-13;checks+=1
out={'checks':checks,'max_algebra_residual':maxerr,'scope':'Literal K3 VALUE twin, complete first algebra, protected/coisometric identities and exact conditional regrouping; no joint heat repair is certified.'}
print(json.dumps(out,indent=2));open(str(Path(__file__).with_name('two_node_twin_checks.json')),'w').write(json.dumps(out,indent=2)+'\n')
