import numpy as np,json
from pathlib import Path
rng=np.random.default_rng(412381);checks=0

def ck(x):
 global checks
 checks+=1;assert bool(x),checks

def gmake(d):
 Q,_=np.linalg.qr(rng.normal(size=(d,d)));T=Q@np.diag(rng.uniform(.45,.55,size=d))@Q.T
 v=rng.normal(size=d);v/=np.linalg.norm(v);k=rng.uniform(1,90)
 return lambda x:T@x+.05/k*np.sin(k*v@x)*v
for d in [1,2,4,7]:
 for rep in range(45):
  g0,gt=gmake(d),gmake(d)
  L,_=np.linalg.qr(rng.normal(size=(d,d)));R,_=np.linalg.qr(rng.normal(size=(d,d)));sv=rng.uniform(.3,.95,size=d);M=L@np.diag(sv)@R.T;sigma=min(sv)
  a=r=rng.uniform(.012,.0625);eb=a**.9;es=a**.95
  S,U,Z,Zp=rng.normal(size=(4,d));x=.6*S+.8*U;F=r*gt(S+a*M@g0(x));eta=.4-2*a*a/(1-2*a*a);cl=sigma*sigma*.4*.4*eta
  def calc(ep,z):
   u=np.zeros(d);vp=np.zeros(d);vm=np.zeros(d);ans={}
   for K in range(1,16):
    un=g0(x+a*M.T@(vp+vm));vpn=gt(S+a*M@u);vmn=-gt(S+a*M@u+ep*z)
    u,vp,vm=un,vpn,vmn
    ans[K]=F-r*vp
   return ans
  for ep in [eb,es]:
   A,B=calc(ep,Z),calc(ep,Zp)
   ck(np.linalg.norm(A[2])<1e-13)
   for K in range(3,16):
    ck(np.linalg.norm(A[K]-B[K])>=cl*r*a*a*ep*np.linalg.norm(Z-Zp)*(1-1e-9))
    ck(np.linalg.norm(A[K])<=r*a*a*ep*np.linalg.norm(Z)*(1+1e-9))
result={'status':'PASS','checks':checks,'scope':['literal simultaneous finite iterations K2 through K15','noncommuting strongly monotone original sources and nonsymmetric nondegenerate known M','uniform auxiliary inverse-Lipschitz and pointwise upper bounds'],'seed':412381}
Path(__file__).with_name('nondegenerate_relative_mark_checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
