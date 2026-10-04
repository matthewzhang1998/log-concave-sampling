#!/usr/bin/env python3
"""Finite positive Gaussian clock calibration and smooth C2 diagonal gate.
All output formulas tested here are analytical amplitude coefficients of actual
VALUE predictors; no derivative-valued producer is executed.
"""
import hashlib,json,math,pathlib
import numpy as np
from numpy.polynomial.legendre import leggauss
ROOT=pathlib.Path(__file__).resolve().parent
checks=0

def ck(ok,msg):
 global checks
 assert ok,msg
 checks+=1

def rule(K,m):
 x,v=leggauss(m);r=[];w=[]
 for j in range(K):
  a=2.**(-j-1)
  r.extend(1-1.5*a-.5*a*x);w.extend(.5*a*v)
 r.append(1-2.**(-K-1));w.append(2.**(-K))
 return np.array(r),np.array(w)

def covariance_calibration(r,w):
 R=np.minimum.outer(r,r)/np.maximum.outer(r,r)
 R0=np.outer(r,r)+np.diag(1-r*r)
 q=np.sqrt(1-r*r);Rmax=np.outer(r,r)+np.outer(q,q)
 S=float(w@R@w);S0=float(w@R0@w);Sm=float(w@Rmax@w)
 ck(S0<.5<Sm,'calibration bracket')
 aux=R0 if S>=.5 else Rmax
 Sa=float(w@aux@w);eta=(.5-S)/(Sa-S) if S!=.5 else 0.
 C=(1-eta)*R+eta*aux
 ck(0<=eta<=1,'convex PSD coefficient')
 ck(np.max(abs(np.diag(C)-1))<3e-14,'unit marginals')
 ck(abs(w@C@w-.5)<4e-14,'exact calibrated quadratic energy')
 ck(np.linalg.eigvalsh(C-np.outer(r,r)).min()>-3e-13,'conditional covariance PSD')
 return C,{'raw_energy':S,'raw_energy_debt':S-.5,'independent_energy':S0,'mix_weight':eta}

def cos_cov(u,R):
 return .5*(np.exp(-u*(1-R))+np.exp(-u*(1+R))-2*np.exp(-u))

def inner_F(u,tau):
 return np.exp(-u/2)-.5*((1-tau)*np.exp(-u*(1-tau))+(1+tau)*np.exp(-u*(1+tau)))

def exact_I(u): return -math.expm1(-2*u)/(2*u)-math.exp(-u)
def exact_F(u): return math.exp(-u/2)+math.exp(-2*u)/u+math.expm1(-2*u)/(2*u*u)

ti,vi=rule(20,20)
reports=[]
for K,m in [(2,2),(4,4),(6,6),(8,8),(12,12),(20,20)]:
 r,w=rule(K,m);N=len(r);s2=float(w@w)
 ck(np.all(w>0),'positive clock weights')
 ck(abs(sum(w)-1)<3e-14 and abs(w@r-.5)<3e-14,'constant and linear exactness')
 C,rep=covariance_calibration(r,w)
 # A literal root factor: conditional Cholesky plus retained endpoint.
 L=np.linalg.cholesky(C-np.outer(r,r)+1e-16*np.eye(N))
 rows=np.column_stack((r,L))
 ck(np.max(abs(rows@rows.T-C))<1e-12,'actual finite Gaussian covariance draw')
 ck(np.linalg.norm((w[:,None]*rows).sum(axis=0))**2<.5+2e-13,'weighted Gaussian row norm')
 # Every inner pair uses a fresh Markov innovation at fixed outer root.
 for i in [0,N//2,N-1]:
  for tau in [.07,.5,.93]:
   row=np.r_[rows[i]*tau,math.sqrt(1-tau*tau)]
   ck(abs(row@row-1)<5e-13,'inner standard marginal')
   ck(abs(row[0]-r[i]*tau)<5e-14,'inner endpoint correlation')
   ck(abs(row[:-1]@rows[i]-tau)<5e-13,'inner outer correlation')
 u=4/s2;k=math.sqrt(u);I=exact_I(u)
 Kq=float(w@cos_cov(u,C)@w)
 Fq=float(vi@inner_F(u,ti));F=exact_F(u)
 ck(Kq>=s2*(1-math.exp(-u))**2/2-2e-14,'unavoidable even-chaos diagonal')
 ck((Kq-I)/u>=s2*s2/12-2e-14,'rank-two diagonal lower bound')
 ck(abs(Fq-F)<=s2*s2/48,'polylog inner clock error fits budget')
 defect=(Kq-I)/u+Fq-F
 ck(defect>=s2*s2/16-3e-14,'actual calibrated predictor nonlinear coefficient gate')
 # Hessian class is independent of the growing oscillation frequency.
 xx=np.linspace(-30,30,20001);a=.5;b=.25
 hh=a-b*np.sin(k*xx)
 ck(hh.min()>=.25-1e-14 and hh.max()<=.75+1e-14,'uniform C2 Hessian sandwich')
 # Quadratic matrix covariance through B²: the inner linear endpoint
 # coefficient is sum_i w_i r_i * sum_j v_j tau_j = 1/4.
 # Use independent local inner innovations; J's full row is explicitly built.
 for eig in [[.003],[.01,.035],[.005,.02,.06]]:
  dim=len(eig);rng=np.random.default_rng(dim+N);O,_=np.linalg.qr(rng.normal(size=(dim,dim)))
  B=O@np.diag(eig)@O.T
  # Aggregate outer row h and inner row ell are enough to execute linear graph.
  h=w@rows
  ell=np.r_[.5*h,np.sqrt(float((w*w).sum()*(vi*vi*(1-ti*ti)).sum()))]
  hhrow=np.r_[h,0.];zz=np.zeros_like(ell);zz[0]=1.
  blocks=[zz[j]*np.eye(dim)-hhrow[j]*B+ell[j]*(B@B) for j in range(len(ell))]
  CY=sum(T@T.T for T in blocks)
  cubic=-2*float(hhrow@ell)*(B@B@B)
  quartic=float(ell@ell)*np.linalg.matrix_power(B,4)
  ck(np.linalg.norm(CY-(np.eye(dim)-B+B@B+cubic+quartic))<3e-13,'literal matrix quadratic coefficient')
 rep.update({'K':K,'m':m,'outer_nodes':N,'inner_nodes':len(ti),'literal_value_count':N*(len(ti)+1),
             'weight_square_mass':s2,'frequency':k,'cosine_rank_two_debt':(Kq-I)/u,
             'inner_F_error':Fq-F,'normalized_variance_coefficient_debt':defect,
             'a_half_b_quarter_variance_coefficient_debt':defect/16,
             'proven_lower_bound':s2*s2/256})
 reports.append(rep)

# Tensor-product raw clocks fail even at degree 1, and at very high degree
# their moment retains diagonal mass. The triangular ratio rule does not.
moments=[]
for K,m in [(4,4),(8,8),(12,12)]:
 r,w=rule(K,m);R=np.minimum.outer(r,r)/np.maximum.outer(r,r)
 ds=[]
 for degree in [1,2,10,100,1000,10**6]:
  direct=float(w@(R**degree)@w);tri=float(w@(r**degree))
  target=1/(degree+1)
  ck(abs(tri-target)<=8*4**(-m)+2**(1-K)+5e-14,'triangle ratio law clock moment bound')
  if degree==10**6: ck(direct>=w@w-1e-14,'raw product diagonal retained')
  ds.append({'degree':degree,'product_error':direct-target,'triangle_error':tri-target})
 moments.append({'nodes':len(r),'moments':ds})

result={'status':'PASS','assertions':checks,'scope':'Positive finite quadratic calibration and a uniform-C2 obstruction for direct nested finite path predictors, not for all positive law compilers.', 'rules':reports,'moments':moments}
(ROOT/'clock_compression_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
