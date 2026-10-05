"""Independent port and error-floor checks; not a compiler certification."""
import math
from fractions import Fraction as F
import numpy as np

# Exact endpoint certificate for the normalized private first.
f2=2*(F(3,4)*F(87,32)**2+F(21,16)**2+F(19,32)**2+F(1,2)**2)
assert f2==F(32231,2048) and f2<16
curl_coeff=math.sqrt(2)/16*(39*math.sqrt(3)+math.sqrt(2381))
assert curl_coeff<12

A=.5
a=np.array([1.,0.]); b=np.array([.6,.8]); I=np.eye(2)
def g(z): return A*(z/2+(np.sin(a@z)*a+np.sin(b@z)*b)/8)
def h(z): return A*(I/2+(np.cos(a@z)*np.outer(a,a)+np.cos(b@z)*np.outer(b,b))/8)
ss=[math.log(4),math.log(4/3)]
rows=[]
for s in ss:
 r=np.exp(-s); an=2*s*r; bn=2*s*(s-1)*r; c=np.sqrt(1-r*r-an*an-bn*bn)
 rows.append((r,an,bn,c))
assert abs(sum(row[0] for row in rows)/2-.5)<1e-14

def node(y):
 x,N,M,L=y.reshape(4,2)
 v=x/2+N/2; w=x/4+N/2+M/4
 gv=g(v); gw=g(w); gb=g(gw); S=x-gv+gb
 out=g(S)-g(x)
 H0=h(x); V=h(v); J=h(gw); W=h(w); H=h(S)
 Srows=[I-V/2+J@W/4,-V/2+J@W/2,J@W/4,np.zeros((2,2))]
 Hbar=H.copy(); Kterms=[np.zeros((2,2)) for _ in range(4)]
 for r,an,bn,c in rows:
  U=r*x+an*N+bn*M+c*L; D=gv-g(U)
  out+=(g(S+D)-g(S-D))/4
  Hp,Hm,HU=h(S+D),h(S-D),h(U)
  Hbar+=(Hp-Hm)/4
  K=(Hp+Hm)/2
  Drows=[V/2-r*HU,V/2-an*HU,-bn*HU,-c*HU]
  for j in range(4): Kterms[j]+=.5*K@Drows[j]
 jac=np.hstack([Hbar@Srows[j]+Kterms[j]-(H0 if j==0 else 0) for j in range(4)])
 return out,jac

def outer(y):
 z,G,N,M,L=y.reshape(5,2)
 out=np.zeros(2); jac=np.zeros((2,10))
 for t in [.2,.8]:
  c=math.sqrt(1-t*t); x=t*z+c*G
  f,D=node(np.array([x,N,M,L]).ravel())
  out+=f/2
  jac[:,0:2]+=t*D[:,0:2]/2
  jac[:,2:4]+=c*D[:,0:2]/2
  jac[:,4:]+=D[:,2:]/2
 return out,jac

rng=np.random.default_rng(1634); worst=0
for _ in range(100):
 y=rng.normal(size=10)*3
 f,J=outer(y); eps=1e-5
 num=np.column_stack([(outer(y+eps*np.eye(10)[j])[0]-outer(y-eps*np.eye(10)[j])[0])/(2*eps) for j in range(10)])
 err=np.linalg.norm(J-num,2); worst=max(worst,err)
 assert err<1e-8
 lift=np.zeros((8,8));lift[:2]=J[:,2:]
 assert math.sqrt(2)*np.linalg.norm(J[:,2:],2)<=4*A
 assert math.sqrt(2)*np.linalg.norm(lift-lift.T,2)<=12*A*A
print('Exact normalized first coefficient squared:',f2)
print('Normalized curl coefficient bound:',curl_coeff)
print('Worst noncommuting two-block finite-difference Jacobian discrepancy:',worst)
print('PASS. Native imported compilers were not tested or certified.')
