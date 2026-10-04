import numpy as np, json, math, hashlib
from pathlib import Path
rng=np.random.default_rng(410411)
checks=0; worst=0.
def ck(b):
 global checks
 checks+=1
 assert bool(b),checks
def sym(A): return (A+A.T)/2
for d in [2,3,5]:
 for _ in range(180):
  X=rng.normal(size=(d,d)); M=X/max(1,np.linalg.norm(X,2))
  hs=[]
  for i in range(4):
   Q=rng.normal(size=(d,d)); Q=sym(Q);hs.append(Q/max(1,np.linalg.norm(Q,2)))
  A0,H0,A1,H1=hs;c=.6;s=.8
  K0=A0@M@H0;K1=A1@M@H1;B=K0-K1
  P=np.hstack([np.eye(d),np.zeros((d,2*d))]);C=np.hstack([c*np.eye(d),s*np.eye(d),np.zeros((d,d))]);L=np.hstack([P.T,C.T]);J=np.block([[np.zeros((d,d)),M],[-M.T,np.zeros((d,d))]])
  def lift(A,H):
   D=np.block([[A,np.zeros((d,d))],[np.zeros((d,d)),H]])
   return L@D@J@D@L.T
  Clead=np.block([[c*(B-B.T),s*B,np.zeros((d,d))],[-s*B.T,np.zeros((d,d)),np.zeros((d,d))],[np.zeros((d,3*d))]])
  ck(np.linalg.norm(lift(A0,H0)-lift(A1,H1)-Clead)<1e-11)
  JS=rng.normal(size=(d,d));JZ=rng.normal(size=(d,d));row=np.hstack([JS,s*B,JZ])
  ck(np.linalg.norm(-sym((row@Clead)[:,:d])-(s*s*B@B.T-c*sym(JS@(B-B.T))))<1e-11)
  ck(np.linalg.norm(K0.T-H0@M.T@A0)<1e-12)
  # two independent finite fine banks: product of means, never mean of self-products
  probs=np.array([.2,.3,.5]);banks=[B,rng.normal(size=(d,d)),rng.normal(size=(d,d))]
  avg=sum(p*b for p,b in zip(probs,banks));prod=sum(probs[i]*probs[j]*(banks[i]@banks[j].T) for i in range(3) for j in range(3))
  ck(np.linalg.norm(prod-avg@avg.T)<1e-11)
T=np.array([[.5,.1],[.1,.5]]);N=np.diag([1.,0.]);e1=np.array([1.,0.]);delta=.1;c=.6;s=.8
for m in range(3,51):
 a=r=1/(2*math.pi*m);eps=a**.9;k=1/(a*eps)
 def g(y): return T@y+delta/k*np.sin(k*y[0])*e1
 def H(y): return T+delta*np.cos(k*y[0])*N
 S=np.zeros(2);U=np.zeros(2);Z=e1;x=c*S+s*U;Dlt=g(S)-g(S+eps*Z);x1=x+a*Dlt;t0=S+a*g(x);t1=S+a*g(x1)
 K0=H(t0)@H(x);K1=H(t1)@H(x1);D=H(S)-H(S+eps*Z)
 ds=.26+delta*math.sin(.5);u=math.cos(a*ds);v=math.cos(.5)
 expect=-delta*(u-v)*(N@T-T@N)
 ck(np.linalg.norm(K0-K1-(K0-K1).T-expect)<1e-12)
 ck(u-v>.12)
 JS=r*(H(t0)-H(t1))+r*a*c*(K0-K1)-r*a*a*K1@D
 JU=r*a*s*(K0-K1);JZ=r*a*a*eps*K1@H(S+eps*Z)
 ck(np.linalg.norm(JS-JS.T-r*a*c*expect)<1e-12)
 ck(np.linalg.norm(JU-JU.T-r*a*s*expect)<1e-12)
 ck(abs((JS-JS.T)[0,1])>=.0012*r*a*c)
 ck(np.linalg.norm(r*(g(t0)-g(t1)))<=r*a*a*eps*(1+1e-12))
 ck(np.linalg.norm(np.hstack([JS,JU,JZ]),2)<=3*r)
 ck(np.linalg.norm(D)<1e-12)
# Exact rational core of the physical transpose deficit.
B0=(T+delta*N)@N
ck(np.linalg.norm(B0@B0.T-sym(B0@B0)-np.array([[0.,.03],[.03,.01]]))<1e-14)
out={'status':'PASS','checks':checks,'random_seed':410411,'scope':['same-query reverse and skew-block algebra','complete physical leading current','independent conditional fine-bank Gram','48 exact same-potential native K3 frequency fixtures','actual mark and complete first at each fixture'],'not_claimed':['joint stationary source for nonlinear query jets','uniform covariance lower bound after clock cutoff','all-dimensional repair closure']}
p=Path(__file__).with_name('matrix_k3_word_checks.json');p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
