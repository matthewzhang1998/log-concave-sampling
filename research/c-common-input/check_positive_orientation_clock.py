import numpy as np,math,json
from numpy.polynomial.legendre import leggauss
from pathlib import Path
rng=np.random.default_rng(765210)
rows=[];total=0
for delta in [.25,.1,.03,.01,.003,.001,.0001]:
 t0=delta*delta/128;T=.5*np.log(16/delta);J=math.ceil(math.log2(T/t0));M=math.ceil(math.log(384/delta,4))
 z,gw=leggauss(M);tt=[];oo=[];tr=[];orr=[]
 for j in range(J):
  a=(2**j)*t0;t=1.5*a+.5*a*z;o=.5*a*gw
  eta=1e-4*delta;rt=t+rng.uniform(-1,1,M)*eta*a/4;ro=o*(1+rng.uniform(-1,1,M)*eta/4);ro*=a/ro.sum()
  tt.extend(t);oo.extend(o);tr.extend(rt);orr.extend(ro)
 tt=np.array(tt);oo=np.array(oo);tr=np.array(tr);orr=np.array(orr)
 k=np.unique(np.r_[np.arange(1,2001),np.round(np.geomspace(2001,max(2002,t0**-2),2000))])
 q=2*(np.exp(-2*k[:,None]*tt[None,:])@oo);qr=2*(np.exp(-2*k[:,None]*tr[None,:])@orr)
 err=np.abs(1-k*q)/np.sqrt(k);rerr=np.abs(1-k*qr)/np.sqrt(k)
 assert err.max()<delta and rerr.max()<delta
 w=2*oo*np.exp(-2*tt);v=np.sqrt(-np.expm1(-2*tt))
 assert np.max(w/(v*v))<=3+1e-12
 assert np.sum(w/(v*v))<=J+1e-10
 assert np.sum(np.sqrt(w)/v)<=J*np.sqrt(M)+1e-10
 rows.append({'delta':delta,'panels':J,'per_panel_degree':M,'nodes':len(tt),'spectral_tests':len(k),'max_relative_one_energy_multiplier':float(err.max()),'max_rounded_multiplier':float(rerr.max()),'sum_w':float(w.sum()),'sum_w_over_v':float(np.sum(w/v)),'sum_w_over_v_squared':float(np.sum(w/v**2)),'sqrt_weight_absolute_path_sum':float(np.sum(np.sqrt(w)/v))});total+=2*len(k)
out={'status':'PASS','spectral_node_rounding_cases':total,'rows':rows,'scope':'New independent scalar clock diagnostics. Uniform infinite-degree bounds and source-value encoding obligations are proved in the source, not inferred from this finite sweep.'}
Path(__file__).with_name('positive_orientation_clock_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
