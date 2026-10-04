import numpy as np,json
from pathlib import Path
rng=np.random.default_rng(681990)
def sym(n):
 X=rng.normal(size=(n,n));X=(X+X.T)/2;return X/max(1,np.linalg.norm(X,2))
def root(X):
 w,V=np.linalg.eigh((X+X.T)/2);assert w.min()>-1e-10;return (V*np.sqrt(np.maximum(w,0)))@V.T
cases=0;worst=0
for n in range(1,9):
 for _ in range(40):
  K0=rng.normal(size=(n,n));K0/=max(1,np.linalg.norm(K0,2));a=.3;B=np.eye(n)+a*K0;Bi=np.linalg.inv(B);S=-Bi@K0
  H=sym(n);JI=sym(n);JB=B.T@H@B;r=.002
  Jf=r*(Bi.T@JB-JI);Curl=Jf-Jf.T;word=S.T@JB-JB@S
  Cnorm=np.sqrt(1+np.linalg.norm(S,2)**2)
  U=np.concatenate((S.T,-np.eye(n)),axis=1)/Cnorm
  V=np.concatenate((np.eye(n),S),axis=0)/Cnorm
  DD=np.block([[JB,np.zeros((n,n))],[np.zeros((n,n)),JB]])
  s0=.2;rho=.1;M=s0*rho*DD;K=U@M@V
  Fi=root(np.eye(2*n)-V@V.T);Fo=root(np.eye(n)-U@U.T);Sp=root(np.eye(2*n)-M@M.T)
  priv=np.concatenate((U@M@Fi,U@Sp,Fo),axis=1)
  err=max(np.max(abs(Curl-r*a*word)),np.max(abs(K-s0*rho*word/Cnorm**2)),np.max(abs(K+K.T)),np.max(abs(K@K.T+priv@priv.T-np.eye(n))))
  q=.1;b=.2;c=-q*r*a*Cnorm**2/(2*b*s0*s0*rho);ME=s0*Jf;O=-(Jf@Curl+(Jf@Curl).T)/2
  err=max(err,np.max(abs(np.eye(n)+b*c*(ME@K.T+K@ME.T)-(np.eye(n)-q*O))))
  worst=max(worst,float(err));assert err<1e-10;cases+=1
out={'status':'PASS','known_nonsymmetric_B_cases':cases,'max_algebra_error':worst,'scope':'Known near-identity affine row, skew extraction, stationary covariance and fork readout. No random/nonlinear B or native nonlinear ancestor is inferred.'}
Path(__file__).with_name('nonsymmetric_affine_adapter_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
