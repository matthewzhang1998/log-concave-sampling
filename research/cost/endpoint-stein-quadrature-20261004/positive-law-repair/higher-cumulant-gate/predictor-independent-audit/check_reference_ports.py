#!/usr/bin/env python3
"""Check finite graph ports and quadratic screen, never a finite law-error claim."""
from pathlib import Path
import json,math
import numpy as np
from numpy.polynomial.legendre import leggauss
out=Path(__file__).resolve().parent
checks=0
def ck(test,msg):
 global checks
 if not test:raise AssertionError(msg)
 checks+=1
rr,w=leggauss(5);r=(rr+1)/2;w=w/2
tt,v=leggauss(4);t=(tt+1)/2;v=v/2
T=sorted(set([1.]+list(r)+[ri*tj for ri in r for tj in t]),reverse=True)
N=len(T);idx={x:j for j,x in enumerate(T)}
G=np.array([[min(x,y)/max(x,y) for y in T] for x in T])
L=np.linalg.cholesky(G)
ck(np.linalg.eigvalsh(G).min()>0,'positive finite Markov Gram')
ck(np.max(abs(np.sum(L*L,axis=1)-1))<1e-14,'each actual Gaussian row has norm one')
# Independently produce the same Gram by successive OU transitions.
L2=np.zeros_like(L);L2[0,0]=1
for i in range(1,N):
 rho=T[i]/T[i-1];L2[i,:]=rho*L2[i-1,:];L2[i,i]=math.sqrt(1-rho*rho)
ck(np.max(abs(L2@L2.T-G))<1e-14,'literal sorted-time OU root construction')
outer=np.array([math.sqrt(wi)*L[idx[ri],:] for ri,wi in zip(r,w)])
inner=np.array([math.sqrt(wi*vj)*L[idx[ri*tj],:] for ri,wi in zip(r,w) for tj,vj in zip(t,v)])
ck(np.linalg.norm(outer,2)<=1+1e-14 and np.linalg.norm(inner,2)<=1+1e-14,'weighted stack norm, no node-count factor')
A=np.diag([.45,.57]);axis=np.array([1.,2.])/math.sqrt(5)
def f(x,amp):return amp*(A@x+.06*(math.cos(axis@x)-1)*axis+.03*np.array([math.cos(x[0])-1,0.]))
def Df(x,amp):return amp*(A-.06*math.sin(axis@x)*np.outer(axis,axis)-.03*np.diag([math.sin(x[0]),0.]))
def graph(R,amp):
 X=L@R;Y=X[0].copy();J=np.zeros((2,2*N));base=np.kron(L[0:1,:],np.eye(2))
 for ri,wi in zip(r,w):
  I=np.zeros(2);DI=np.zeros((2,2*N))
  for tj,vj in zip(t,v):
   ii=idx[ri*tj];I+=vj*f(X[ii],amp);DI+=vj*Df(X[ii],amp)@np.kron(L[ii:ii+1,:],np.eye(2))
  ii=idx[ri];point=X[ii]-I
  Y-=wi*f(point,amp)
  J-=wi*Df(point,amp)@(np.kron(L[ii:ii+1,:],np.eye(2))-DI)
 return Y,J,base
rng=np.random.default_rng(191923)
for amp in (.1,.2,.4):
 for trial in range(4):
  R=rng.normal(size=(N,2));Y,J,base=graph(R,amp)
  ck(np.linalg.norm(J,2)<=amp*(1+amp),'raw residual first A(1+A)')
  d=rng.normal(size=(N,2));d/=np.linalg.norm(d);h=2e-5
  fd=(graph(R+h*d,amp)[0]-graph(R-h*d,amp)[0])/(2*h)-base@d.reshape(-1)
  ck(np.max(abs(fd-J@d.reshape(-1)))<1e-9,'recorded HVP chain differentiates actual finite graph')
 ck(np.array_equal(graph(np.zeros((N,2)),amp)[0],np.zeros(2)),'exact all-root zero')
# The finite B^2 screen follows directly from actual row products.
H=sum(wi*L[idx[ri],:] for ri,wi in zip(r,w))
JJ=sum(wi*vj*L[idx[ri*tj],:] for ri,wi in zip(r,w) for tj,vj in zip(t,v))
S=float(H@H)
ck(abs(L[0]@H-.5)<1e-14,'Cov Z,H=1/2')
ck(abs(L[0]@JJ-.25)<1e-14,'Cov Z,J=1/4')
ck(abs(S+.5-1)>.001,'finite diagonal covariance screen is nonzero')
for amp in (.03,.07,.13):
 row=L[0]-amp*H+amp*amp*JJ
 actual=row@row
 polynomial=1-amp+(S+.5)*amp*amp-2*(H@JJ)*amp**3+(JJ@JJ)*amp**4
 ck(abs(actual-polynomial)<1e-14,'finite exact quadratic covariance')
# Analytical continuous two-substitution covariance, not a finite realization.
for amp in (.025,.1,.25,.5):
 cov=1-amp+amp*amp-.75*amp**3+.375*amp**4
 target=1/(1+amp)
 ck(abs(math.sqrt(cov)-math.sqrt(target))<=amp**3,'continuous quadratic W2 respects resolvent bound')
report=dict(status='PASS',assertions=checks,scope='Finite VALUE graph ports and quadratic screen; no finite-to-continuous order-three law comparison.',distinct_times=N,original_value_upper_bound=len(r)*len(t)+len(r),S_Q=S,B2_defect=S-.5,minimum_gram_eigenvalue=float(np.linalg.eigvalsh(G).min()))
(out/'reference_ports_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
