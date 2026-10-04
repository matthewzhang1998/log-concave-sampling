# Publication copy; original/public SHA-256 mapping is in INVENTORY.json.
import json, math
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
rng=np.random.default_rng(903904)
checks=0; rows=[]
def ck(x):
 global checks
 checks+=1
 if not bool(x): raise AssertionError(checks)
for delta in [.25,.1,.03,.01,.003,.001,1e-5]:
 t0=delta*delta/128; T=.5*math.log(16/delta)
 J=math.ceil(math.log2(T/t0)); M=math.ceil(math.log(384/delta,4))
 z,w=leggauss(M); ts=[]; oms=[]; tsr=[]; omsr=[]
 eta=delta/2048
 for l in range(J):
  a=2**l*t0; t=a*(1.5+.5*z); om=.5*a*w
  tr=np.clip(t+rng.uniform(-eta*a,eta*a,size=M),a,2*a)
  wr=om*(1+rng.uniform(-eta/4,eta/4,size=M)); wr*=a/wr.sum()
  ck(np.all(om>0)); ck(np.all(wr>0)); ck(np.sum(np.abs(wr-om))<=eta*a)
  ck(np.max(np.abs(tr-t))<=eta*a*1.0001)
  ts.extend(t);oms.extend(om);tsr.extend(tr);omsr.extend(wr)
 ts=np.array(ts);oms=np.array(oms);tsr=np.array(tsr);omsr=np.array(omsr)
 W=2*oms*np.exp(-2*ts); V=-np.expm1(-2*ts)
 ck(W.sum()<=2); ck(np.all(W<=3*V*(1+1e-12)))
 ck((W/V).sum()<=J*(1+1e-12))
 ck((np.sqrt(W/V)).sum()<=J*math.sqrt(M)*(1+1e-12))
 ck((W/np.sqrt(V)).sum()<=4)
 ks=np.unique(np.concatenate([np.arange(1,1001),np.geomspace(1,1e14,2000)]))
 q=2*(np.exp(-2*ks[:,None]*ts[None,:])@oms)
 qr=2*(np.exp(-2*ks[:,None]*tsr[None,:])@omsr)
 err=np.abs(1-ks*q)/np.sqrt(ks)
 errr=np.abs(1-ks*qr)/np.sqrt(ks)
 for er,er2 in zip(err,errr): ck(er<=delta/2); ck(er2<=delta)
 rows.append({'delta':delta,'panels':J,'gauss_degree':M,'nodes':len(ts),'max_spectral_error':float(err.max()),'rounded_max_error':float(errr.max()),'weight_mass':float(W.sum()),'root_square_mass':float((W/V).sum()),'caller_mass':float(np.sqrt(W/V).sum())})
out={'status':'PASS','checks':checks,'rows':rows,'scope':'New positive scalar spectral clock checks and positive within-panel rounding. Exact Gaussian chaos factorization is proved in the source; native coefficient and root-observer admission are separate.'}
p=Path(__file__).with_name('positive_orientation_clock_checks.json');p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
