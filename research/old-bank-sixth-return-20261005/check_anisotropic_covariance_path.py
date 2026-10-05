import numpy as np,math,json,itertools
from numpy.polynomial.hermite_e import hermegauss
rng=np.random.default_rng(6421)
def sym3(T):
 return sum(np.transpose(T,p) for p in itertools.permutations(range(3)))/6
K=sym3(rng.normal(size=(2,2,2)))/4
T=[sym3(rng.normal(size=(2,2,2)))/3 for _ in range(2)]
v=[(1-math.exp(-2))/2,(1+math.exp(-2))/2-math.exp(-1)]
C=np.array([[1.2,.23],[.23,.8]]);Q=np.linalg.inv(C);R=np.linalg.cholesky(C)
q,w=hermegauss(34); w=w/math.sqrt(2*math.pi)
pp=np.array(list(itertools.product(q,q)))@R.T
ww=np.array([a*b for a in w for b in w])
N=lambda A:np.einsum('iab,ac,bd,jcd->ij',A,Q,Q,A)
M=lambda A:np.einsum('iab,ac,jcd->ijbd',A,Q,A)
Nv=sum(vv*N(A) for vv,A in zip(v,T)); Mv=sum(vv*M(A) for vv,A in zip(v,T))
Ns=N(K);Ms=M(K)
allidx={n:list(itertools.product(range(2),repeat=n)) for n in range(6)}
D=[];Cov=[];Mean=[];Smean=[]
for p in pp:
 x=Q@p;cache={():1.}
 def H(I):
  if I in cache:return cache[I]
  i=I[0];J=I[1:]
  z=x[i]*H(J)-sum(Q[i,j]*H(J[:k]+J[k+1:]) for k,j in enumerate(J))
  cache[I]=z;return z
 H1=x;H2=np.array([H(I) for I in allidx[2]]).reshape(2,2)
 H3=np.array([H(I) for I in allidx[3]]).reshape(2,2,2)
 H5=np.array([H(I) for I in allidx[5]]).reshape((2,)*5)
 delta=-Nv@H1-2*np.einsum('ijbd,jbd->i',Mv,H3)
 for vv,A in zip(v,T):delta-=.5*vv*np.einsum('iab,jcd,abjcd->i',A,A,H5)
 ff=[np.einsum('iab,ab->i',A,H2) for A in T]
 Cov.append(sum(vv*np.outer(f,f) for vv,f in zip(v,ff)))
 D.append(delta);Mean.append(np.einsum('iab,ab->i',K,H2));Smean.append(-Ns@H1-2*np.einsum('ijbd,jbd->i',Ms,H3))
D=np.array(D);Cov=np.array(Cov);Mean=np.array(Mean);Smean=np.array(Smean)
checks=[]
for theta in [np.array([.3,.7]),np.array([-.8,.2]),np.array([1.,-.6])]:
 Ctheta=np.einsum('i,nij,j->n',theta,Cov,theta)
 J=1j*(D@theta)-.5*Ctheta
 leading=np.sum(ww*J*np.exp(1j*(pp@theta)))
 assert abs(leading)<2e-12
 for a in [.0,.4,1.]:
  row=[]
  for eps in [.12,.18,.27,.4]:
   shift=pp+eps**3*Mean+eps**6*Smean+a*eps**6*D
   phase=np.exp(1j*(shift@theta)-.5*(.55*np.dot(theta,theta)+a*eps**6*Ctheta))
   cur=eps**6*np.sum(ww*J*phase)
   row.append({'alpha':eps,'current_norm':float(abs(cur)),'normalized_alpha9':float(abs(cur)/eps**9)})
  checks.append({'theta':theta.tolist(),'path_a':a,'leading_error':float(abs(leading)),'cases':row})
assert max(x['normalized_alpha9'] for z in checks for x in z['cases'])<2
out={'physical_covariance':C.tolist(),'checks':checks,'max_leading_cancellation_error':max(z['leading_error'] for z in checks),'max_alpha9_normalized_current':max(x['normalized_alpha9'] for z in checks for x in z['cases']),'scope':'Two-dimensional non-diagonal covariance diagnostic for same-endpoint covariance/skew-square/rank-six cancellation. This does not execute the imported native VALUE programs.'}
open('/workspace/shared/old-bank-sixth-return-20261005/anisotropic_covariance_path_checks.json','w').write(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['max_leading_cancellation_error','max_alpha9_normalized_current']},indent=2))
