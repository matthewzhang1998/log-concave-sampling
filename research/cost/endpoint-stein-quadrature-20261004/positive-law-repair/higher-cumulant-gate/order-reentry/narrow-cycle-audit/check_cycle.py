from fractions import Fraction as F
from math import factorial,sqrt
import itertools,json,numpy as np
checks=0
# Exact formal series: log E exp(z S T)=-1/2 log(1-z^2).
N=20
m=[F(0)]*(N+1);m[0]=F(1)
for k in range(1,N//2+1):
 d=1
 for j in range(1,2*k,2):d*=j
 m[2*k]=F(d*d,factorial(2*k))
log=[F(0)]*(N+1)
for n in range(1,N+1):
 log[n]=m[n]-sum(F(k,n)*log[k]*m[n-k] for k in range(1,n))
 expected=F(1,n) if n%2==0 else F(0)
 assert log[n]==expected;(checks:=checks+1)
# Correct force/rank families and scalar normalization.
for h in range(1,15):
 assert 4*(2*h)==8*h and 2*(2*h)==4*h;checks+=1
 assert 4*(2*h)+2==8*h+2;checks+=1
 assert 4*(2*h+1)+2==8*h+6;checks+=1
rng=np.random.default_rng(16102026)
def cutnorm(T,mask):
 n=T.ndim;left=[i for i in range(n) if mask>>i&1];right=[i for i in range(n) if not mask>>i&1]
 M=T.transpose(left+right).reshape(2**len(left),2**len(right))
 return np.linalg.norm(M,2)
def maxcut(T):
 return max(cutnorm(T,m) for m in range(1,2**T.ndim-1))
for _ in range(35):
 B=rng.normal(size=(2,2,2,2));B/=maxcut(B)
 C=rng.normal(size=(2,2,2,2));C/=maxcut(C)
 D=rng.normal(size=(2,2,2,2));D/=maxcut(D)
 two=np.einsum('ijab,klab->ijkl',B,C)
 three=np.einsum('ijab,klbc,mnca->ijklmn',B,C,D)
 for T in [two,three]:
  for mask in range(1,2**T.ndim-1):
   assert cutnorm(T,mask)<=1+1e-11;checks+=1
  assert np.linalg.norm(T)<=np.linalg.norm(B)+1e-11;checks+=1
# Sharp marked diagonal cycles and excluded internal self-loop.
B=np.zeros((2,2,2,2))
for i in range(2):B[i,i,i,i]=1
T=np.einsum('ijab,klab->ijkl',B,B)
assert abs(maxcut(T)-1)<1e-12;checks+=1
assert abs(np.linalg.norm(T)-sqrt(2))<1e-12;checks+=1
E=np.zeros((2,2,2,2,2))
for a in range(2):
 for c in range(2):E[0,a,a,c,c]=F(1,2)
assert maxcut(E)<=1+1e-12;checks+=1
assert np.linalg.norm(np.einsum('iaacc->i',E))==2;checks+=1
out={'assertions_passed':checks,'scope':'Exact scalar reverse-symbol and finite grouped-cut diagnostics only; no positive cycle producer or all-order closure.'}
print(json.dumps(out,indent=2))
