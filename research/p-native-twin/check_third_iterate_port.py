from pathlib import Path
import numpy as np, json, hashlib
rng=np.random.default_rng(43104); count=0
for d in [1,2,4]:
 for rep in range(30):
  Q,_=np.linalg.qr(rng.normal(size=(d,d)))
  b0=rng.normal(size=d);bt=rng.normal(size=d)
  def G(x):return .35*x+.12*Q@np.sin(Q.T@x)
  def HG(x):return .35*np.eye(d)+.12*Q@np.diag(np.cos(Q.T@x))@Q.T
  g0=lambda x:G(b0+x)-G(b0)
  gt=lambda x:G(bt+x)-G(bt)
  h0=lambda x:HG(b0+x)
  ht=lambda x:HG(bt+x)
  M=rng.normal(size=(d,d));M/=max(1,np.linalg.norm(M,2))
  a=rng.uniform(.001,1/16);eps=rng.uniform(.001,1);r=.03
  c=.6;s=.8
  C=np.hstack((np.eye(d),np.zeros((d,d))))
  P=np.hstack((c*np.eye(d),s*np.eye(d)))
  Cd=np.vstack((C,P));e_t=np.vstack((np.zeros((d,d)),np.eye(d)))
  Ctw=np.block([[Cd,np.zeros((2*d,d))],[Cd,eps*e_t]])
  Ad=np.block([[np.zeros((d,d)),np.zeros((d,d))],[a*M,np.zeros((d,d))]])
  St=np.block([[Ad+Ad.T,Ad.T],[Ad,np.zeros_like(Ad)]])
  beta=1/np.sqrt(s**-2+eps**-2)
  B=beta*np.hstack((np.zeros((d,d)),np.eye(d)/s,-np.eye(d)/eps))
  Prec=np.hstack((P,np.zeros((d,d))))
  select=np.hstack((np.zeros((d,d)),np.eye(d),np.zeros((d,2*d))))
  assert np.linalg.norm(B@B.T-np.eye(d))<1e-12;count+=1
  assert np.linalg.norm(B@Ctw.T-beta*select)<1e-12;count+=1
  assert np.linalg.norm(Prec@B.T-beta*np.eye(d))<1e-12;count+=1
  w=rng.normal(size=3*d);W=w[:2*d];Z=w[2*d:];S=P@W;x=C@W
  pp=np.zeros(4*d);Jp=np.zeros((4*d,3*d)); hist=[];Jh=[]
  for K in range(1,9):
   args=Ctw@w+St@pp
   xx=[args[j*d:(j+1)*d] for j in range(4)]
   vals=[g0(xx[0]),gt(xx[1]),-g0(xx[2]),-gt(xx[3])]
   hs=[h0(xx[0]),ht(xx[1]),-h0(xx[2]),-ht(xx[3])]
   H=np.zeros((4*d,4*d))
   for j in range(4):H[j*d:(j+1)*d,j*d:(j+1)*d]=hs[j]
   pp=np.concatenate(vals);Jp=H@(Ctw+St@Jp)
   hist.append(r*pp[d:2*d]);Jh.append(r*Jp[d:2*d])
   GK=(r/beta)*Ctw.T@pp
   assert np.linalg.norm(B@GK-hist[-1])<1e-11;count+=1
  F=r*gt(S+a*M@g0(x));qfb=x+a*M.T@(gt(S)-gt(S+eps*Z))
  H3=r*gt(S+a*M@g0(qfb));E3=F-H3
  assert np.linalg.norm(hist[1]-F)<1e-12;count+=1
  assert np.linalg.norm(hist[2]-H3)<1e-12;count+=1
  for K in range(3,9):
   Ek=F-hist[K-1]
   assert np.linalg.norm(Ek)<=r*a*a*eps*np.linalg.norm(Z)+1e-13;count+=1
  expected=r*a*a*eps*ht(S+a*M@g0(qfb))@M@h0(qfb)@M.T@ht(S+eps*Z)
  assert np.linalg.norm(-Jh[2][:,2*d:]-expected)<1e-12;count+=1
  JF=r*ht(S+a*M@g0(x))@(Prec+a*M@h0(x)@np.hstack((C,np.zeros((d,d)))))
  for K in range(3,9):
   JE=JF-Jh[K-1]
   assert np.linalg.norm(JE,2)<=3*r+1e-12;count+=1
   curl=Prec.T@JE-JE.T@Prec
   assert np.linalg.norm(curl,2)<=5*r*a+1e-12;count+=1
  U=s*W[:d]-c*W[d:]
  assert np.linalg.norm(c*S+s*U-x)<1e-12;count+=1
  delta=a*M.T@(gt(S)-gt(S+eps*Z))
  Phi=lambda u:gt(S+a*M@g0(c*S+s*u))
  assert np.linalg.norm(r*(Phi(U)-Phi(U+delta/s))-E3)<1e-12;count+=1
out=Path(__file__).parent;f=out/'THIRD-ITERATE-TWO-NODE-WEAK-TWIN-PORT.md'
data={'status':'PASS','checks':count,'scope':'Actual finite nonlinear same-potential iterations, coisometry/column signs, same-Z mark, complete first/curl, and conditional-coordinate identities','source_sha256':hashlib.sha256(f.read_bytes()).hexdigest()}
(out/'third_iterate_port_checks.json').write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps(data,indent=2))
