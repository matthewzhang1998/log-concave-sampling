import numpy as np, json
from pathlib import Path
rng=np.random.default_rng(73145)
la=.5; de=be=.025
checks=0; maxmark=0.; maxmixed=0.; maxsym=0.
def g(x):
 R=np.sqrt(1+x@x);e=np.zeros_like(x);e[0]=1
 return la*x+de*(e*(R-1)+x[0]*x/R)+be*np.sin(x[0])*e
def H(x):
 R=np.sqrt(1+x@x);e=np.zeros_like(x);e[0]=1
 return la*np.eye(len(x))+de*((np.outer(e,x)+np.outer(x,e))/R+x[0]*(np.eye(len(x))/R-np.outer(x,x)/R**3))+be*np.cos(x[0])*np.outer(e,e)
for d in [1,2,4,7]:
 for it in range(150):
  a=r=rng.uniform(.005,.06); eps=a**.9; sig=10**rng.uniform(-4,0)
  S,U,Z,V=rng.normal(size=(4,d)); x=(S+U)/np.sqrt(2)
  Q,_=np.linalg.qr(rng.normal(size=(d,d)));T,_=np.linalg.qr(rng.normal(size=(d,d)))
  M=Q@np.diag(rng.uniform(0,1,d))@T.T
  x1=x+a*M.T@(g(S)-g(S+eps*Z));t0=S+a*M@g(x);t1=S+a*M@g(x1)
  E=r*(g(t0)-g(t1));q=a*sig*M@V
  Es=r*(g(t0+q)-g(t1+q));ratio=np.linalg.norm(Es)/max(np.linalg.norm(E),1e-300)
  assert ratio<=1.5+1e-8;maxmark=max(maxmark,ratio); checks+=1
  bound=(.25/.4)*np.linalg.norm(q)*np.linalg.norm(E)
  mixed=np.linalg.norm(Es-E)/max(bound,1e-300)
  assert np.linalg.norm(Es-E)<=bound+2e-15*r;
  if bound>1e-13*r: maxmixed=max(maxmixed,mixed)
  checks+=1
  J=r*a*sig*M.T@(H(t0+q)-H(t1+q))@M
  sym=np.linalg.norm(J-J.T);assert sym<1e-12;maxsym=max(maxsym,sym);checks+=1
  direction=rng.normal(size=d);step=1e-4
  def active(v):
   qq=a*sig*M@v
   return M.T@(r*(g(t0+qq)-g(t1+qq)))
  fd=(active(V+step*direction)-active(V-step*direction))/(2*step)
  assert np.linalg.norm(fd-J@direction)<2e-10;checks+=1
  p0=g(x)+sig*V;p1=g(x1)+sig*V;feedback=a*(g(S)-g(S+eps*Z))
  GS=r*(g(S+a*M@p0)-g(S+a*M@p1)+a*(g(S)-g(S+eps*Z))-feedback)
  assert np.linalg.norm(GS-Es)<1e-13
  assert np.linalg.norm(r*(g(x)-p0)+r*sig*V)<1e-13
  assert np.linalg.norm(r*(-g(x1)+p1)-r*sig*V)<1e-13
  checks+=3
out={'checks':checks,'max_pointwise_mark_ratio':maxmark,'max_mixed_difference_ratio':maxmixed,'max_auxiliary_gradient_skew':maxsym,'scope':'Exact source identities and finite-dimensional diagnostics. The dimension-uniform drift restoration is proved in the companion note, not inferred from these samples.'}
path=Path(__file__).with_name('common_auxiliary_heat_checks.json');path.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
