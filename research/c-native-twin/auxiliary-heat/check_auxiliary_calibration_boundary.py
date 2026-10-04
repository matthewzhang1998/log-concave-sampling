import numpy as np,json
from pathlib import Path
rng=np.random.default_rng(521612);n=240000;m=.5;d=.1;p=.3;q=1.55
rows=[]
for a in [.02,.006,.002,.0006,.0002,.00006]:
 r=a;eps=a**.9;sig=a**p;k=a**(-q)
 S,U,Z=rng.normal(size=(3,n));x=(S+U)/np.sqrt(2)
 h=-a*m*eps*Z-(2*a*d/k)*np.cos(k*(S+eps*Z/2))*np.sin(k*eps*Z/2)
 gap=-a*(m*h+(2*d/k)*np.cos(k*(x+h/2))*np.sin(k*h/2))
 t0=S+a*(m*x+(d/k)*np.sin(k*x))
 osc=(2*r*d/k)*np.cos(k*(t0-gap/2))*np.sin(k*gap/2)
 E=r*m*gap+osc
 Ms=r*m*gap+np.exp(-.5*(k*a*sig)**2)*osc
 scale=r*a*a*eps;en=E/scale;mn=Ms/scale
 e=en.std();rel=np.sqrt(np.mean((mn-en)**2))/e
 leading=m*Z*(m+d*np.cos(k*x))*(m+d*np.cos(k*t0))
 row={'a':a,'sigma':sig,'k':k,'kaepsilon':k*a*eps,'kasigma':k*a*sig,'normalized_energy':e,'limit_energy':m*(m*m+d*d/2),'relative_value_change':float(rel),'limit_relative_value_change':d/np.sqrt(2*m*m+d*d),'L2_expansion_error':float(np.sqrt(np.mean((en-leading)**2))),'samples':n}
 assert abs(e-m*(m*m+d*d/2))<.003
 assert .12<rel<.16
 rows.append(row)
# Independent algebra tests of the orthogonal rotation and its two h formulas.
checks=0;maxorth=0.;maxh=0.
for dim in [1,2,4,8]:
 for j in range(120):
  M=rng.normal(size=(dim,dim));M/=max(1,np.linalg.norm(M,2));a=rng.uniform(.005,.2);s=rng.uniform(.001,1)
  ev,Q=np.linalg.eigh(np.eye(dim)-(a*s)**2*M@M.T);D=(Q*np.sqrt(ev))@Q.T
  ev,Q=np.linalg.eigh(np.eye(dim)-(a*s)**2*M.T@M);D2=(Q*np.sqrt(ev))@Q.T
  B=np.block([[D,a*s*M],[-a*s*M.T,D2]])
  err=np.linalg.norm(B@B.T-np.eye(2*dim));assert err<1e-12;maxorth=max(maxorth,err)
  S,V=rng.normal(size=(2,dim));Sp,Vp=(B@np.r_[S,V]).reshape(2,dim)
  C=np.linalg.inv(np.eye(dim)+D)
  h=s*V-a*s*s*M.T@C@S; hp=s*Vp+a*s*s*M.T@C@Sp
  err=np.linalg.norm(h-hp);assert err<1e-12;maxh=max(maxh,err)
  assert np.linalg.norm(a*M@h-(Sp-S))<1e-12
  checks+=3
out={'native_value_boundary':rows,'rotation_checks':checks,'max_rotation_orthogonality_error':maxorth,'max_h_regression_identity_error':maxh,'scope':'Native scalar common-heat source and exact known Gaussian algebra. Monte Carlo supports, but does not replace, the analytical limits.'}
Path(__file__).with_name('auxiliary_calibration_boundary_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
