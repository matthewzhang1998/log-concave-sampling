import numpy as np,json
from pathlib import Path
rng=np.random.default_rng(75144);count=0;maxerr=0.;minalign=1.;maxread=0.;maxfirst=0.
for n in [1,2,4,7]:
 for _ in range(450):
  L,_=np.linalg.qr(rng.normal(size=(n,n)));R,_=np.linalg.qr(rng.normal(size=(n,n)))
  ss=10**rng.uniform(-5,0,n);ss[rng.random(n)<.25]=0
  if not ss.any():ss[0]=.1
  M=L@np.diag(ss)@R.T;Q=M@M.T;q=np.linalg.norm(Q,'fro')
  H=[]
  for j in range(3):
   Z,_=np.linalg.qr(rng.normal(size=(n,n)))
   H.append((Z*rng.uniform(.4,.6,n))@Z.T)
  W=H[0]@M@H[1]@M.T@H[2]
  err=np.linalg.norm(W-.125*Q,'fro')/q;align=np.sum(W*Q)/(q*q)
  assert err<=.091+1e-12 and align>=.034-1e-12
  maxerr=max(maxerr,err);minalign=min(minalign,align);count+=2
  a=rng.uniform(.001,.15);sig=10**rng.uniform(-5,0);A=a*sig*M
  ev,U=np.linalg.eigh(np.eye(n)-A@A.T);D=(U*np.sqrt(ev))@U.T
  ev,U=np.linalg.eigh(np.eye(n)-A.T@A);D2=(U*np.sqrt(ev))@U.T
  beta=.5;B0=beta*np.concatenate([M.T,np.zeros((n,n)),np.zeros((n,n)),-M.T@A@np.linalg.inv(np.eye(n)+D2)],axis=1)
  eps=a**.9;O=np.zeros((n,n));I=np.eye(n)
  Cs=[np.concatenate(z,axis=1) for z in [(I,O,O,O),(I,O,eps*I,O),(D,O,O,-A),(D,O,eps*I,-A)]]
  for C in Cs:
   residual=np.linalg.norm(B0@C.T-beta*M.T)
   assert residual<1e-12;maxread=max(maxread,residual);count+=1
  assert np.linalg.norm(B0,2)**2<=5/16+1e-12;count+=1
  first=sum(np.linalg.norm(C,2)**2 for C in Cs)/beta
  assert first<=12+1e-12;maxfirst=max(maxfirst,first);count+=1
out={'checks':count,'max_relative_product_error':maxerr,'min_fixed_Q_pairing':minalign,'max_common_readout_residual':maxread,'max_gradient_first_multiplier':maxfirst,'scope':'Exact deterministic matrix diagnostics for the analytical rank-aware energy bound and projected-gradient offspring. No sample-estimated Gaussian lower bound is used.'}
Path(__file__).with_name('rank_aware_restoration_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
