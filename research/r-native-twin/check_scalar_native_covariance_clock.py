from pathlib import Path
import numpy as np,json
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss
n=9;gh,gw=hermgauss(n);gh*=np.sqrt(2);gw/=np.sqrt(np.pi)
x=np.stack(np.meshgrid(gh,gh,gh,indexing='ij'),axis=-1).reshape(-1,3)
w=np.einsum('i,j,k->ijk',gw,gw,gw).reshape(-1)
a=.04;r=.06;eps=.3;c=.6;s=.8;scale=r*a*a*eps
# The same uniformly convex original potential, with its same actual query genealogy.
def g(z):return .5*z+.1*np.sin(.7*z)/.7
def H(z):return .5+.1*np.cos(.7*z)
def values_jac(z):
 S,U,Z=np.moveaxis(z,-1,0);xx=c*S+s*U;delta=g(S)-g(S+eps*Z);x1=xx+a*delta
 t0=S+a*g(xx);t1=S+a*g(x1);k0=H(t0)*H(xx);k1=H(t1)*H(x1);d=H(S)-H(S+eps*Z)
 val=(g(t0)-g(t1))/(a*a*eps)
 js=(H(t0)-H(t1)+a*c*(k0-k1)-a*a*k1*d)/(a*a*eps)
 ju=s*(k0-k1)/(a*eps);jz=k1*H(S+eps*Z)
 return val,np.stack([js,ju,jz],axis=-1)
val,j=values_jac(x);mu=w@val;var=w@((val-mu)**2)
t,tw=leggauss(16);t=(t+1)/2;tw/=2
cols=np.zeros(3)
for tj,wj in zip(t,tw):
 cv=np.sqrt(tj);vv=np.sqrt(1-tj);J=[]
 for start in range(0,len(x),60):
  q=cv*x[start:start+60,None,:]+vv*x[None,:,:]
  _,jac=values_jac(q)
  J.append(np.einsum('j,ijk->ik',w,jac))
 J=np.concatenate(J);cols+=wj*np.einsum('i,ij->j',w,J*J)
err=abs(var-cols.sum());assert err<2e-7
assert cols[1]>=0 and cols[2]>0
out={'status':'PASS','gauss_nodes_per_input':9,'clock_nodes':16,'native_source':'K3, same .5x+.1sin(.7x)/.7 primitive at both nodes','normalization':'all quantities divided by (r*a^2*epsilon)^2','variance':float(var),'S_forward_square':float(cols[0]),'U_true_square':float(cols[1]),'Z_true_square':float(cols[2]),'covariance_clock_residual':float(err),'scope':'Numerical quadrature of the actual same-source scalar K3 covariance identity; analytical theorem does not depend on quadrature accuracy.'}
Path(__file__).with_name('scalar_native_covariance_clock_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
