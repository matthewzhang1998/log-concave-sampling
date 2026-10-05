"""Diagnostics only. Not a native compiler implementation or numerical certificate."""
import json, math
from pathlib import Path
import numpy as np
RNG=np.random.default_rng(5102026)
r=1/math.sqrt(2); k=math.sqrt(3/8)
checks=0

def assert_(condition):
 global checks
 assert condition
 checks+=1

def jac(f,x,eps=1e-6):
 return np.stack([(f(x+eps*np.eye(len(x))[j])-f(x-eps*np.eye(len(x))[j]))/(2*eps) for j in range(len(x))],axis=1)

def normalized(g,x,N,M):
 v=(x+N)/2; w=x/4+N/2+M/4
 return g(x-r*g((v-k*g(w/k))/r))

def collective(g,x,Z,q):
 M=len(q); D=len(x); G=[]; y=x; prev=1.
 w=(2*np.arange(M,0,-1)-1)/M**2
 for j in range(M):
  t=q[j]/prev
  y=t*y+math.sqrt(1-t*t)*Z[j]
  G.append(y);prev=q[j]
 Y=[g(y) for y in G]; suffix=np.zeros(D); B=np.zeros(D)
 for j in range(M-1,-1,-1):
  lo=((M-j-1)/M)**2
  c=(q[j]-lo)/q[j]
  assert_(c>=0 and c<=1)
  assert_(abs(c+sum(w[j+1:])/q[j]-1)<1e-11)
  H=c*Y[j]+suffix/q[j]
  B+=w[j]*g(G[j]-H)
  suffix+=w[j]*Y[j]
 return g(x-B)

Avals=[.001,.03,.1,.5]
for A in Avals:
 for D in [1,2,4]:
  Q,_=np.linalg.qr(RNG.normal(size=(D,D)))
  K=Q@np.diag(A*RNG.random(D))@Q.T
  for _ in range(20):
   x,N,M=RNG.normal(size=(3,D));v=(x+N)/2;w=x/4+N/2+M/4
   got=normalized(lambda z:K@z,x,N,M)
   want=K@x-K@K@v+K@K@K@w
   assert_(np.linalg.norm(got-want)<1e-12)

D=2;a=np.array([1.,0.]);b=np.array([.6,.8])
for A in Avals:
 def g(y):return A*(y/2+(math.sin(a@y)*a+math.sin(b@y)*b)/8)
 for _ in range(12):
  x,N,M=RNG.normal(size=(3,D));z=RNG.normal(size=D)
  # One outer node t=1/2 has the exact first moment and positive mass.
  t=.5;c=math.sqrt(1-t*t)
  def en(R):
   G,N,M=R.reshape(3,D);x=t*z+c*G
   return normalized(g,x,N,M)-g(x)
  J=jac(en,RNG.normal(size=3*D))
  lift=np.zeros((3*D,3*D));lift[:D]=J
  assert_(np.linalg.norm(J,2)<1.2*A+1e-7)
  assert_(np.linalg.norm(lift-lift.T,2)<1.85*A*A+1e-7)
 for m in [1,2,3,7]:
  for _ in range(3):
   n=np.arange(m,0,-1);lo=((n-1)/m)**2;hi=(n/m)**2
   q=lo+(hi-lo)*RNG.uniform(.01,.99,size=m)
   x=RNG.normal(size=D);Z=RNG.normal(size=(m,D));count=[0]
   def counted(y):count[0]+=1;return g(y)
   collective(counted,x,Z,q)
   assert_(count[0]==2*m+1)
   j=jac(lambda zz:collective(g,x,zz.reshape(m,D),q),Z.ravel())
   assert_(np.linalg.norm(j,2)<=A*A*(1+A)+1e-7)

Sigma=np.array([[1,.5,.25],[.5,.5,.375],[.25,.375,.375]])
C=np.array([[-5.,12.,-8.],[-1.,1.,0.],[0.,-1.,1.]])
B=np.array([[-1.,0.,0.],[1.,-1.,0.],[0.,1.,-1.]])
assert_(np.linalg.norm(C@Sigma+Sigma@C.T+np.diag([2.,0.,0.]))<1e-12)
assert_(np.linalg.norm(C-Sigma@B.T@np.linalg.inv(Sigma))<1e-12)
assert_(np.linalg.eigvalsh(Sigma)[0]>0)
d=math.exp(-1)/4*(r*(r*math.cosh(r)-math.sinh(r))-(math.sinh(1)-2*math.cosh(1)+2))
assert_(d< -math.exp(-1)/2880)
out={'diagnostic_assertions':checks,'seed':5102026,'sine_raw_first_chaos':d,'native_compiler_executed':False,'claims':'Quadratic exactness, finite noncommuting first/curl diagnostics, convex weights, literal VALUE count, scalar Gaussian history identities.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
