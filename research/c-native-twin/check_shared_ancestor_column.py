# Publication copy; original/public SHA-256 mapping is in INVENTORY.json.
from pathlib import Path
import json,numpy as np
from numpy.polynomial.hermite import hermgauss
x,w=hermgauss(100); x=x*np.sqrt(2); w=w/np.sqrt(np.pi)
checks=0; max_ratio=0.; max_energy_ratio=0.
s0=np.sqrt(.7); c0=np.sqrt(.3)
for k0 in [.4,2.,7.]:
 for kt in [.3,3.,9.]:
  for S,Z in [(.3,1.2),(-.8,-.7)]:
   a=.05;r=.03;eps=.7;v=.6;mu=.2
   g0=lambda q:.3*q+.2*np.sin(k0*q)/k0
   H0=lambda q:.3+.2*np.cos(k0*q)
   gt=lambda q:.4*q+.2*np.sin(kt*q)/kt
   Ht=lambda q:.4+.2*np.cos(kt*q)
   shift=a*(gt(S)-gt(S+eps*Z))
   def E(u):
    q=c0*S+s0*u
    return r*(gt(S+a*g0(q))-gt(S+a*g0(q+shift)))
   q=c0*S+s0*(mu+v*x)
   je=np.dot(w,r*a*(Ht(S+a*g0(q))*H0(q)-Ht(S+a*g0(q+shift))*H0(q+shift)))
   vals=E(mu+v*x);ener=np.sqrt(np.dot(w,(vals-np.dot(w,vals))**2))
   for K in range(1,9):
    if K==1:t=np.array([.25])
    else:
     lo=1/(4*K*K); hi=.25
     t=(hi+lo)/2+(hi-lo)/2*np.cos(np.arange(K)*np.pi/(K-1))
    ds=np.array([np.prod([-t[l]/(t[j]-t[l]) for l in range(K) if l!=j]) for j in range(K)])
    assert abs(sum(ds)-1)<1e-10;assert sum(abs(ds))<=4+1e-8
    p=x[:,None];z=x[None,:];I=np.zeros((len(x),len(x)))
    for tj,dj in zip(t,ds):
     b=np.sqrt(tj);base=np.sqrt(1-tj)*z
     I+=dj*(E(mu+v*(base+b*p))-E(mu+v*(base-b*p)))/(2*b*s0*v)
    mean=I@w;err=np.sqrt(np.dot(w,(mean-je*x)**2))
    ratio=err/(ener*4**(-K)/(s0*v)+1e-25)
    max_ratio=max(max_ratio,ratio)
    assert ratio<4, (k0,kt,S,K,ratio)
    ie=np.sqrt(w@(I*I)@w); eratio=ie/(K*ener/(s0*v)+1e-25)
    max_energy_ratio=max(max_energy_ratio,eratio);assert eratio<4
    checks+=4
out={'checks':checks,'max_calibration_ratio_to_e_4minusK_over_s0v':max_ratio,'max_energy_ratio_to_K_e_over_s0v':max_energy_ratio,'scope':'Literal paired finite filter with one shared ancestor Gaussian and fixed terminal/twin feedback; no adjoint or stationary completion.'}
print(json.dumps(out,indent=2));open(str(Path(__file__).with_name('shared_ancestor_column_checks.json')),'w').write(json.dumps(out,indent=2)+'\n')
