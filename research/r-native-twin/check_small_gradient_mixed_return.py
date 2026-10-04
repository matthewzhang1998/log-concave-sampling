from pathlib import Path
from fractions import Fraction as F
import numpy as np,json
rng=np.random.default_rng(132892);checks=0;worst=0.
def ck(x):
 global checks
 checks+=1
 assert x,checks
def eq(x,y):
 global worst
 er=float(np.linalg.norm(np.asarray(x)-np.asarray(y)));worst=max(worst,er);ck(er<2e-12)
for d,n in [(1,3),(2,5),(3,8)]:
 for _ in range(160):
  A=.025;g=.5;kap=A**(1+g);rho=A**(g/4);alpha=rho/kap
  B=np.linalg.qr(rng.normal(size=(n,d)))[0].T
  G=rng.normal(size=(n,n));G=(G+G.T)/2;G*=kap/max(1,np.linalg.norm(G,2))
  JE=rng.normal(size=(d,n));JE*=A/max(1,np.linalg.norm(JE,2));JH=B@G
  weight=float(rng.uniform(.001,.1));bo=.15;co=.18;ze=.3
  C=alpha*JH.T;M=-weight*JE/(2*bo*co*alpha)
  ck(np.linalg.norm(C,2)<1);ck(np.linalg.norm(M,2)<1)
  # Exact Gaussian forward-chain covariance, with both private self completions.
  covQ=C@C.T+(np.eye(n)-C@C.T)
  covV=M@covQ@M.T+(np.eye(d)-M@M.T)
  got=(bo*bo+ze*ze)*np.eye(d)+co*co*covV+bo*co*(M@C+(M@C).T)
  want=(bo*bo+co*co+ze*ze)*np.eye(d)-weight*(JE@JH.T+JH@JE.T)/2
  eq(got,want);ck(np.linalg.eigvalsh(got).min()>.05)
for gi in range(1,201):
 g=F(gi,200);nu=g/4;ell=2+3*g/4;energy=1+3*g/4;target=3+2*g
 rows=[3+5*g/2,4+3*g/2,F(13,2)+3*g,7+3*g]
 for z in rows:ck(z>target)
 eq(float(ell*3+energy-F(1,2)),float(rows[2]));eq(float(ell*3+energy),float(rows[3]))
 ck(energy>1);ck(ell+nu>ell);ck(2*ell-F(1,2)>ell)
 # Explicit nominal profile and fixed-list child interpolation, not an actual-e claim.
 for j in [F(0),F(39,10),F(5)]:
  for p in [2,3,6]:
   needed=max(target+j-ell,j+energy-ell,F(0))+1
   B=int(np.ceil(float(needed*(p-1)/nu)))+1
   ck(ell+nu*B/(p-1)>target+j)
out={'status':'PASS','checks':checks,'max_reference_residual':worst,'scope':'Literal positive child/parent covariance, rational marked grades, complete first returns and declared nominal child-floor selection. Does not remove absolute finite-prior floors.'}
Path(__file__).with_name('small_gradient_mixed_return_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
