import numpy as np,json
from pathlib import Path
rng=np.random.default_rng(409312)
def root(X):
 w,V=np.linalg.eigh((X+X.T)/2);assert w.min()>-1e-10;return (V*np.sqrt(np.maximum(w,0)))@V.T
def sym(n):
 X=rng.normal(size=(n,n));X=(X+X.T)/2;return X/max(1,np.linalg.norm(X,2))
cases=0;worst=0
for n in range(1,9):
 for _ in range(40):
  D=sym(n);H=sym(n);a=.3;B=np.eye(n)+a*D;Bi=np.linalg.inv(B);L0=np.linalg.norm(B,2)**2;Hp=B@H@B/L0
  Cn=np.sqrt(2)*np.linalg.norm(Bi,2)
  U=Bi@np.concatenate((np.eye(n),-D),axis=1)/Cn
  V=np.concatenate((D,np.eye(n)),axis=0)@Bi/Cn
  DD=np.block([[Hp,np.zeros((n,n))],[np.zeros((n,n)),Hp]])
  s0=.2;rho=.13;M=s0*rho*DD
  Fi=root(np.eye(2*n)-V@V.T);Fo=root(np.eye(n)-U@U.T);S=root(np.eye(2*n)-M@M.T)
  K=U@M@V
  Kexp=s0*rho*(H@D-D@H)/(L0*Cn**2)
  private=np.concatenate((U@M@Fi,U@S,Fo),axis=1)
  err=max(np.max(abs(K-Kexp)),np.max(abs(K+K.T)),np.max(abs(K@K.T+private@private.T-np.eye(n))))
  J=.05*rng.normal(size=(n,n))/max(1,np.linalg.norm(rng.normal(size=(n,n)),2));JE=s0*J
  tau=1e-6;b=.2;q=.1;c=-q*tau*Cn**2/(2*b*s0*s0*rho)
  W=(H@D-D@H)/L0;O=-(J@(tau*W)+(J@(tau*W)).T)/2
  cov=np.eye(n)+b*c*(JE@K.T+K@JE.T)
  err=max(err,np.max(abs(cov-(np.eye(n)-q*O))))
  worst=max(worst,float(err));assert err<1e-10;cases+=1
out={'status':'PASS','general_symmetric_D_cases':cases,'max_identity_residual':worst,'scope':'Known-contraction and stationary conditional covariance algebra, exact skew factor and fork readout; source/caller contracts still use the original finite pair proof.'}
Path(__file__).with_name('one_pair_affine_adapter_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
