from pathlib import Path
import hashlib,json
import numpy as np
rng=np.random.default_rng(41004)
checks=0
for n in [1,2,3,5,8]:
 for rep in range(40):
  X=rng.normal(size=(n,n)); Q,_=np.linalg.qr(X)
  D=Q@np.diag(rng.uniform(-1,1,n))@Q.T
  a=rng.uniform(.01,.49); B=np.eye(n)+a*D; Bi=np.linalg.inv(B)
  H=rng.normal(size=(n,n)); H=(H+H.T)/2; H/=max(1,np.linalg.norm(H,2))
  Hp=B@H@B
  C=np.sqrt(2)*np.linalg.norm(Bi,2)
  U=Bi@np.hstack((np.eye(n),-D))/C
  V=np.vstack((D,np.eye(n)))@Bi/C
  HD=np.block([[Hp,np.zeros_like(Hp)],[np.zeros_like(Hp),Hp]])
  W=H@D-D@H
  assert np.linalg.norm(U@HD@V-W/C**2)<1e-11;checks+=1
  assert np.linalg.norm(U,2)<=1+1e-12 and np.linalg.norm(V,2)<=1+1e-12;checks+=1
  s0=.125; rho=.05; AA=s0*rho*HD; K=U@AA@V
  cond=U@(np.eye(2*n)-AA@V@V.T@AA.T)@U.T+np.eye(n)-U@U.T
  assert np.linalg.norm(cond-(np.eye(n)-K@K.T))<1e-11;checks+=1
  J=rng.normal(size=(n,n)); J/=max(1,np.linalg.norm(J,2)); M=s0*J
  tau=1e-7; q=.3; b=.1; c=-q*tau*C*C/(2*s0*s0*rho*b)
  O=-tau*(J@W+(J@W).T)/2
  cross=b*c*(M@K.T+K@M.T)
  assert np.linalg.norm(cross+q*O)<1e-12;checks+=1
  assert b*b+c*c<.5;checks+=1
p=Path(__file__).parent
f=p/'ONE-PAIR-KNOWN-MATRIX-COMMUTATOR-ADAPTER.md'
data={'checks':checks,'status':'PASS','scope':'noncommuting matrix identities, contraction rows, exact reference covariance and cross sign; not nonlinear law tests','source_sha256':hashlib.sha256(f.read_bytes()).hexdigest()}
(p/'one_pair_commutator_checks.json').write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps(data,indent=2))
