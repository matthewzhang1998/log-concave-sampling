import numpy as np, math,json
from pathlib import Path
rng=np.random.default_rng(412221);checks=0

def ck(v):
 global checks
 checks+=1;assert bool(v),checks

def funcs(d):
 Q=rng.normal(size=(d,d));Q=(Q+Q.T)/2;Q=.5*np.eye(d)+.08*Q/max(1,np.linalg.norm(Q,2));n=rng.normal(size=d);n/=np.linalg.norm(n);k=2.7;delta=.07
 def g(x):return Q@x+delta/k*np.sin(k*(n@x))*n
 def V(x):return .5*x@Q@x-delta/k**2*np.cos(k*(n@x))
 return g,V
for d in [1,2,3,4]:
 for rep in range(30):
  g,V=funcs(d);g0,V0=g,V;gt,Vt=g,V
  M=rng.normal(size=(d,d));M/=max(1,np.linalg.norm(M,2));a=.04;r=.07;eps=.3
  X=rng.normal(size=(8,d))
  def grad(X):
   S,Z,x0,x1,p0,p1,dd,ll=X;t0=S+a*M@p0;t1=S+a*M@p1
   t0v,t1v,gs,gp,gx0,gx1=gt(t0),gt(t1),gt(S),gt(S+eps*Z),g0(x0),g0(x1)
   return r*np.array([t0v-t1v+a*(gs-gp)-dd,-a*eps*gp,gx0-p0-ll,-gx1+p1+ll,a*M.T@t0v-x0,-a*M.T@t1v+x1,-S-M@ll,x1-x0-M.T@dd])
  def pot(X):
   S,Z,x0,x1,p0,p1,dd,ll=X
   return r*(Vt(S+a*M@p0)-Vt(S+a*M@p1)+V0(x0)-V0(x1)-x0@p0+x1@p1+a*(Vt(S)-Vt(S+eps*Z))-S@dd+ll@(x1-x0-M.T@dd))
  hh=2e-5;fd=np.zeros_like(X);H=np.zeros((8*d,8*d))
  for j in range(8*d):
   shift=np.zeros(8*d);shift[j]=hh;shift=shift.reshape(8,d)
   fd.flat[j]=(pot(X+shift)-pot(X-shift))/(2*hh)
   H[:,j]=((grad(X+shift)-grad(X-shift))/(2*hh)).ravel()
  ck(np.linalg.norm(fd-grad(X))<2e-8)
  ck(np.linalg.norm(H-H.T)<2e-8)
  ck(np.linalg.norm(H,2)<=16*r)
  S,U,Z=rng.normal(size=(3,d));x0=.6*S+.8*U;dd=a*(gt(S)-gt(S+eps*Z));x1=x0+M.T@dd;p0=g0(x0);p1=g0(x1);GG=np.array([S,Z,x0,x1,p0,p1,dd,np.zeros(d)])
  ret=grad(GG);E=r*(gt(S+a*M@p0)-gt(S+a*M@p1))
  ck(np.linalg.norm(ret[0]-E)<1e-14)
  ck(np.linalg.norm(ret[2])+np.linalg.norm(ret[3])+np.linalg.norm(ret[7])<1e-14)
  ck(np.linalg.norm(ret[4]+ret[5]-a*M.T@E-r*M.T@dd)<1e-14)
  ck(np.linalg.norm(E)<=r*a*a*eps*np.linalg.norm(Z)*(1+1e-12))
  ck(np.linalg.norm(grad(np.zeros((8,d))))<1e-14)
for L in [8,12,20,40]:
 a=r=1/(2*L);eps=a**.9
 for N in [8,16,32,64,128,256,1024,4096]:
  k=2*math.pi*N;C=math.pi/4+.1
  g=lambda y:.5*y+.1/k*math.sin(k*y)
  H=lambda y:.5+.1*math.cos(k*y)
  S=1.;x=0.;Z=-L/(N*eps);de=g(S)-g(S+eps*Z);x1=x+a*de;t1=S+a*g(x1)
  ck(abs(x1-math.pi/(2*k))<1e-12)
  ck(abs(t1-(1+a*C/k))<1e-12)
  exact=r*a*(.1*k*g(t1)+a*(.6**3-.5**2*(.5+.1*math.cos(a*C))))
  ck(exact>=.05*r*a*k)
  ck(abs(r*(g(S)-g(t1)))<=r*a*a*eps*abs(Z)*(1+1e-8))
  pull=exact-r*(.1+.05*math.pi)
  ck(pull-exact<0)
  if a*k>=12:ck(pull>=.025*r*a*k)
result={'status':'PASS','checks':checks,'scope':['full eight-block ambient gradient matches scalar potential','ambient Hessian symmetry and small first','same-record graph restriction and zero blocks','exact extra p-channel identity','natural gradient and pullback unbounded-first fixture'],'seed':412221}
Path(__file__).with_name('joint_gradient_lifts_checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
