"""Checks of finite inverse / ambient Fenchel blocks; proof lives in the notes."""
from pathlib import Path
import math,json
import numpy as np
rng=np.random.default_rng(202610041736)
count=0;max_inv=0.; max_first=0.;max_grad=0.;max_common=0.;max_native=0.
eps=.05
v=np.array([1.,1.]); e=np.array([1.,0.])
def g(x): return .5*x+eps*np.sin(x[0])*e+(eps/2)*np.sin(x.sum())*v
def H(x): return .5*np.eye(2)+eps*np.cos(x[0])*np.outer(e,e)+(eps/2)*np.cos(x.sum())*np.outer(v,v)
def V(x): return .25*x@x-eps*np.cos(x[0])-(eps/2)*np.cos(x.sum())
def hinv(p,N):
 u=np.zeros(2);J=np.zeros((2,2))
 for j in range(N):
  J=(np.eye(2)-2*H(u))@J+2*np.eye(2)
  u=u+2*(p-g(u))
 return u,J
# Exact-reference numerical inverse is used only to inspect identities, not in the proposed VALUE code.
for k in range(100):
 p=rng.normal(size=2); h,_=hinv(p,40)
 for N in [1,2,3,5,10,15]:
  u,J=hinv(p,N)
  er=np.linalg.norm(u-h);bd=2.5*(.2**N)*np.linalg.norm(p)
  assert er<=bd+2e-15; max_inv=max(max_inv,er/max(bd,1e-30));count+=1
  first=np.linalg.norm(J,2);max_first=max(max_first,first);assert first<=2.5+1e-12;count+=1
# Explicit nonconservative finite inverse.
p=np.array([math.pi/4,math.pi/4]); _,J=hinv(p,3)
target=-4*eps**2*math.sin(2*eps)*np.array([[0.,1.],[-1.,0.]])
assert np.linalg.norm(J-J.T-target)<1e-14;count+=1
# Offgraph Fenchel Hessian for scalar g=m x + delta/k sin(kx).
B=np.array([[.4,-1],[-1,1/.6]])
assert np.linalg.det(B)<0;count+=1
# Full 8-block reference, with high-accuracy inverse in proof diagnostics.
def G(X,a,r,width,M,om):
 S,Z,x0,x1,p0,p1,d,lam=X
 t0=S+a*M@p0;t1=S+a*M@p1;D=g(S)-g(S+width*Z)
 h0=hinv(p0,45)[0];h1=hinv(p1,45)[0]
 return r*np.array([g(t0)-g(t1)+a*D-d,-a*width*g(S+width*Z),om*(g(x0)-p0)-lam,-om*(g(x1)-p1)+lam,a*M.T@g(t0)+om*(h0-x0),-a*M.T@g(t1)-om*(h1-x1),-S-M@lam,x1-x0-M.T@d])
def gradstar(p):
 h=hinv(p,45)[0];return p@h-V(h)
def potential(X,a,r,width,M,om):
 S,Z,x0,x1,p0,p1,d,lam=X
 t0=S+a*M@p0;t1=S+a*M@p1
 Q0=V(x0)+gradstar(p0)-x0@p0;Q1=V(x1)+gradstar(p1)-x1@p1
 return r*(V(t0)-V(t1)+om*(Q0-Q1)+a*(V(S)-V(S+width*Z))-S@d+lam@(x1-x0-M.T@d))
for k in range(30):
 a=rng.uniform(.01,1/16);r=a;width=a**.9;om=rng.choice([1.,a]); c=.6;s=.8
 M=rng.normal(size=(2,2));M/=max(1,np.linalg.norm(M,2))
 X=rng.normal(size=(8,2))
 GX=G(X,a,r,width,M,om)
 for j in range(16):
  dX=np.zeros((8,2));dX.flat[j]=1e-5
  fd=(potential(X+dX,a,r,width,M,om)-potential(X-dX,a,r,width,M,om))/(2e-5)
  max_grad=max(max_grad,abs(fd-GX.flat[j]));assert abs(fd-GX.flat[j])<1e-9;count+=1
 S,U,Z=rng.normal(size=(3,2));x0=c*S+s*U;delta=g(S)-g(S+width*Z);d=a*delta;x1=x0+M.T@d;p0=g(x0);p1=g(x1)
 X=np.array([S,Z,x0,x1,p0,p1,d,np.zeros(2)])
 GG=G(X,a,r,width,M,om);E=r*(g(S+a*M@p0)-g(S+a*M@p1))
 er=np.linalg.norm(GG[0]-E);max_native=max(max_native,er);assert er<1e-14;count+=1
 er=np.linalg.norm(GG[4]+GG[5]-a*M.T@E);max_common=max(max_common,er);assert er<1e-14;count+=1
 for ix in [2,3,7]: assert np.linalg.norm(GG[ix])<1e-14;count+=1
 # Same common heat: exact p-x constraints and complete physical value restoration bound.
 vv=rng.normal(size=2);sigma=a**2
 ps0=p0+sigma*vv;ps1=p1+sigma*vv
 Es=r*(g(S+a*M@ps0)-g(S+a*M@ps1))
 assert np.linalg.norm(Es-E)<=2*r*a*sigma*np.linalg.norm(M,2)*np.linalg.norm(vv)+1e-14;count+=1
# Quadratic opposite companion and false independent fiber variance.
for d in [1,2,5]:
 for a in [.01,.03,.0625]:
  r=a;wid=a**.9;m=.5;c=.6;s=.8;b=r*a*a*m**3*wid
  op_l2=math.sqrt(2)*r*a*m*math.sqrt(d*((1+a*m*c)**2+(a*m*s)**2+(a*a*m*m*wid/2)**2))
  assert op_l2>=math.sqrt(2)*r*a*m*math.sqrt(d);count+=1
  tau=a**4
  assert abs(2*r*r*a*a*m**3*tau-r*r*a*a*tau/4)<1e-30;count+=1
out={'status':'PASS','checks':count,'max_inverse_error_fraction_of_bound':max_inv,'max_inverse_first':max_first,'max_ambient_gradient_finite_difference_error':max_grad,'max_common_mode_identity_error':max_common,'max_native_physical_identity_error':max_native,'finite_inverse_skew':(J-J.T).tolist(),'scope':'Finite identity diagnostics; reference inverse in diagnostics only. No finite-output gradient or new graph law is inferred.'}
Path(__file__).with_name('fenchel_k3_lift_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
