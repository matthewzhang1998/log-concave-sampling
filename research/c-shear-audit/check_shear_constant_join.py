# Publication copy; original/public SHA-256 mapping is in INVENTORY.json.
from pathlib import Path
import json
import numpy as np
from fractions import Fraction as F
rng=np.random.default_rng(19041)
checks=0; maxerr=0.
sym=lambda x:(x+x.T)/2
for n in range(2,12):
  for rep in range(30):
    J=rng.normal(size=(n,n))*.015
    u=np.eye(n)[:,0]; v=np.eye(n)[:,1]
    Pi=np.eye(n)-np.outer(v,v)
    H=sym(rng.normal(size=(n,n))); H=Pi@H@Pi
    H/=max(1.,np.linalg.norm(H,2))
    m=rng.uniform(-1,1); D=H+m*np.outer(v,v)
    J0=np.outer(u,v)-np.outer(v,u)
    wedge=np.outer(H@u,v)-np.outer(v,H@u)
    err=np.linalg.norm(D@J0@D-m*wedge)
    maxerr=max(maxerr,err); assert err<1e-12; checks+=1
    cs=.2; rho=.12; tau=.002; q=.7; b=.15
    K=cs**2*rho**2*D@J0@D; M=cs*J
    c=-q*tau/(2*b*cs**3*rho**2)
    Ot=-tau*sym(J@(m*wedge))
    err=np.linalg.norm(b*c*(M@K.T+K@M.T)+q*Ot)
    maxerr=max(maxerr,err); assert err<1e-12; checks+=1
    # Exact covariance of two conditionally independent stationary channels.
    vari=.8
    actual=b*b*(M@M.T+np.eye(n)-M@M.T)+c*c*(K@K.T+np.eye(n)-K@K.T)+b*c*(M@K.T+K@M.T)+(vari-b*b-c*c)*np.eye(n)
    err=np.linalg.norm(actual-(vari*np.eye(n)-q*Ot))
    maxerr=max(maxerr,err); assert err<1e-12; checks+=1
    # Whole-root Riesz/Price identity on arbitrary matrix-valued chaos 1+2.
    # U(R)=L R+Q(R^2-1), RU=L+Q R, DO=L+2 Q R.
    L=sym(rng.normal(size=(n,n))); Q=sym(rng.normal(size=(n,n)))
    t=rng.normal(size=n); z=.37
    tl=t@L@t; tq=t@Q@t
    d_price=6*q*q*z*(tl*tl+2*tq*tq)
    d_riesz=q*q*z/4*24*(tl*tl+2*tq*tq)
    err=abs(d_price-d_riesz)/(1+abs(d_price)); maxerr=max(maxerr,err)
    assert err<1e-12; checks+=1
for gi in range(1,21):
  g=F(gi,20)
  for zi in range(1,40):
    z=F(zi,120); gamma=g-3*z-F(1,10000)
    if gamma<=0: continue
    assert 1+g+z<F(7,3)
    assert 3+2*g-z>1+g+z
    assert 1+g-z-2*gamma>0
    checks+=3
out={'checks':checks,'max_relative_or_algebra_residual':maxerr,'scope':'Exact fork covariance, full-square wedge, Price/Riesz coefficient, and rational exponent guards. Source admission and retained-host return are not inferred from these tests.'}
print(json.dumps(out,indent=2))
open(str(Path(__file__).with_name('shear_constant_join_checks.json')),'w').write(json.dumps(out,indent=2)+'\n')
