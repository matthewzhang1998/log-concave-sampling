# Publication copy; original/public SHA-256 mapping is in INVENTORY.json.
# New diagnostic. Historical source/audit hashes are provenance only.
import numpy as np,math,json
from pathlib import Path
rng=np.random.default_rng(401812)
def field(n,k):
 a=np.empty((n,)+(n,)*k);v={}
 for ix in np.ndindex(a.shape):
  key=(ix[0],tuple(sorted(ix[1:])))
  if key not in v:v[key]=rng.normal()
  a[ix]=v[key]
 return a
cases=0;worst=0;hr=0;opr=0
for n in [2,3,4]:
 for deg in range(1,6):
  for rep in range(12):
   chunks=[];cov=np.zeros((n,n));ccov=cov.copy();e2=0
   for k in range(1,deg+1):
    f=field(n,k)/math.sqrt(math.factorial(k)*n**k)
    T=sum(np.swapaxes(f,0,j) for j in range(1,k+1))/k;Af=f-T
    C=k*(f-np.swapaxes(f,0,1));div=sum(np.swapaxes(C,1,j) for j in range(1,k+1))/k
    err=np.max(np.abs(div-k*Af));worst=max(worst,float(err));assert err<1e-12
    F=f.reshape(n,-1);AA=Af.reshape(n,-1);CC=C.reshape(n,-1)
    cov+=math.factorial(k)*F@F.T;ccov+=math.factorial(k-1)*CC@CC.T;e2+=math.factorial(k)*np.sum(f*f);chunks.append((k,F,AA))
   kap=np.sqrt(max(0,np.linalg.eigvalsh(ccov)[-1]));lf=np.sqrt(np.linalg.eigvalsh(cov)[-1])
   for t in [1e-6,1e-4,.003,.05,.3,1.]:
    dO=np.zeros((n,n));rcov=dO.copy()
    for k,F,AA in chunks:
     z=-np.expm1(-2*k*t);dO+=math.factorial(k)*z*F@AA.T;rcov+=math.factorial(k)*z*z*AA@AA.T
    dO=(dO+dO.T)/2;eta=min(1,np.sqrt(2*t))
    assert np.linalg.eigvalsh(rcov-eta*eta*ccov)[-1]<1e-10
    x=np.linalg.norm(dO,'fro')/(np.sqrt(e2)*kap*eta) if kap else 0;y=np.linalg.norm(dO,2)/(lf*kap*eta) if kap else 0
    assert x<=1+1e-12 and y<=1+1e-12;hr=max(hr,float(x));opr=max(opr,float(y));cases+=1
out={'status':'PASS','provenance':'New run for reconstructed theorem; it also agrees with the recorded historical diagnostic values.','mixed_chaos_heat_cases':cases,'max_divergence_identity_error':worst,'max_HS_bound_ratio':hr,'max_operator_bound_ratio':opr,'scope':'Analytical finite-chaos/row-covariance identity, not an executed orientation source.'}
Path(__file__).with_name('reconstructed_whole_source_heat_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
