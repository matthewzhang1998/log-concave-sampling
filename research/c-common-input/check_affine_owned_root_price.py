import numpy as np,math,json
from pathlib import Path
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss
rng=np.random.default_rng(890671)
x,w=hermgauss(96);x=x*np.sqrt(2);w=w/np.sqrt(np.pi)
c,cw=leggauss(128);c=(c+1)/2;cw=cw/2
fac=np.exp(-(1-c*c)/2)
Rsin=(np.cos(x[:,None]*c[None,:])@(cw*fac))
Rcos=-(np.sin(x[:,None]*c[None,:])@(cw*fac))
sin=np.sin(x);cos=np.cos(x);mean_cos=np.exp(-.5)
def sym(n):
 X=rng.normal(size=(n,n));X=(X+X.T)/2;return X/max(1,np.linalg.norm(X,2))
def dfac(k):
 return math.prod(range(1,2*k,2))
cases=0;worst=0
for n in range(1,6):
 for rep in range(20):
  U=.015*sym(n);V=.01*sym(n);A=.02*sym(n);B=.01*sym(n)
  a=rng.normal(size=n);a/=np.linalg.norm(a)
  ur=(a@U@a)*sin+(a@V@a)*(cos-mean_cos)
  ru=(a@U@a)*Rsin+(a@V@a)*Rcos
  for z in [0.,.2,.7,1.]:
   q=.4
   var=1+(a@A@a)*sin+(a@B@a)*cos-q*z*ur
   dvar=(a@A@a)*cos-(a@B@a)*sin-q*z*((a@U@a)*cos-(a@V@a)*sin)
   assert var.min()>.9
   for k in [2,3,4]:
    direct=-q*k*dfac(k)*np.dot(w,ur*var**(k-1))
    d4=math.prod(range(2*k-3,2*k+1))*dfac(k-2)*var**(k-2)
    current=-q/4*np.dot(w,ru*dvar*d4)
    er=abs(direct-current);worst=max(worst,er);assert er<1e-11;cases+=1
out={'status':'PASS','same_root_Price_current_cases':cases,'worst_identity_residual':worst,'scope':'Independent Gaussian conditional-covariance tests of the exact −q/4 root/Price current, including shared sine/cosine root fields and noncommuting matrix coefficients. Does not prove a new native source law.'}
Path(__file__).with_name('affine_owned_root_price_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
