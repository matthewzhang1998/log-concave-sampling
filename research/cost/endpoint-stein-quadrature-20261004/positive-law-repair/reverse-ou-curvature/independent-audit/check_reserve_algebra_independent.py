"""Independent deterministic diagnostics for the reverse-OU curvature reserve.
Does not import author diagnostics or pretend to execute the finite LOW30 compiler.
"""
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
rng=np.random.default_rng(620041004)
checks=0
maxima={}
def check(ok,label):
    global checks
    checks+=1
    if not bool(ok): raise AssertionError(label)
def close(a,b,label,tol=1e-10):
    err=float(np.linalg.norm(a-b))
    maxima[label]=max(maxima.get(label,0),err)
    check(err<=tol*(1+np.linalg.norm(a)+np.linalg.norm(b)),label)
def op(a): return float(np.linalg.norm(a,2))
def quadrature(K,m):
    p,w=leggauss(m); ts=[]; ws=[]
    for k in range(K):
        lo=2.**(-k-1);hi=2.**(-k)
        ts.extend((lo+hi)/2+(hi-lo)*p/2);ws.extend((hi-lo)*w/2)
    ts.append(2.**(-K-1));ws.append(2.**(-K))
    return 1-np.array(ts),np.array(ws)
for K in [1,3,7]:
 for m in [1,2,4,9]:
    r,w=quadrature(K,m)
    close(np.array(w.sum()),np.array(1.),'quadrature_total')
    close(np.array(w@r),np.array(.5),'quadrature_first_moment')
    check(r.min()>=1/(4*m*m),'inverse_r_guard')
    check(w@(1/r)<=4*m*m,'lift_coefficient_guard')
r,w=quadrature(5,5);q=np.sqrt(1-r*r);n=len(r)
# Arbitrary anisotropic quadratic matrices, in an unrelated rotated basis.
for d in [2,5,13,31]:
 for A in [.5,.12,.007]:
  O,_=np.linalg.qr(rng.normal(size=(d,d)))
  B=(O*rng.uniform(0,A,d))@O.T
  for s in [1.,.47,.04,1e-4]:
    I=np.eye(d);P=np.concatenate([I,np.zeros((d,d))],axis=1)
    mats=[np.concatenate([rr*I,qq*I],axis=1) for rr,qq in zip(r,q)]
    KM=s*s*sum(ww*B@R for ww,R in zip(w,mats))
    FM=B@(P-KM); BM=B@P
    HM=sum(ww/rr*R.T@(s*s*B)@R for ww,rr,R in zip(w,r,mats))
    GM=P.T@BM
    Jp=s*GM+HM/s;Jm=s*GM-HM/s
    close(P@(Jp@Jp.T-Jm@Jm.T)@P.T,
          2*(KM@BM.T+BM@KM.T),'quadratic_polarization')
    Hstar=(FM@P.T+P@FM.T)/2-(KM@BM.T+BM@KM.T)/2
    close(Hstar,B-s*s*B@B,'quadratic_Hstar')
    H=B@np.linalg.inv(I+s*s*B)
    check(np.linalg.norm(Hstar-H)<=A**3*s**4*np.sqrt(d)+1e-14,
          'quadratic_third_order_bound')
    # An arbitrary retained-input pair: its cross covariance must be kept.
    Delta=min(s*s/2,.2);vD=(1-Delta)/2
    M=Delta/vD*P.T@FM
    joint_diff=2*np.eye(2*d)-M-M.T
    close(vD/2*P@joint_diff@P.T,
          vD*I-Delta*(FM@P.T+P@FM.T)/2,'retained_pair_readout')
# Exact one-energy estimate on arbitrary Gaussian couplings A_i=L_i Z,B_i=R_i Z.
for d in [1,2,11,64]:
 for k in range(12):
    L2=rng.normal(size=(d,d+4));R2=rng.normal(size=(d,d+4))
    L1=L2+.2*rng.normal(size=L2.shape);R1=R2+.2*rng.normal(size=R2.shape)
    lhs=np.linalg.norm(L1@R1.T-L2@R2.T)
    rhs=op(R1)*np.linalg.norm(L1-L2)+op(L2)*np.linalg.norm(R1-R2)
    check(lhs<=rhs*(1+1e-12),'one_energy_coupling')
# All exponents, including a very small s and Delta much smaller than s^2.
for A in [.5,.04,1e-6]:
 for s in [1.,.3,.001,1e-8]:
  for frac in [1.,.15,1e-5]:
    Delta=min(frac*s*s,.25);alpha=A*s*s;target=Delta*A**3*s**4
    terms=[Delta**2*A**2*alpha,Delta**4*A**4,
           Delta**4*A**4/np.sqrt(alpha),Delta*A**3*s**4,
           Delta**2*A**4*s**4,Delta**2*A**3.5*s**3]
    for val in terms: check(val<=target*(1+1e-12),'s_Delta_exponents')
# Genuine C2/non-C3 nonorthogonal ridge fixture, with no Hessian regularity rate.
V=np.array([[1.,0.],[.6,.8],[-.3,.7]])
V/=op(V)
def gp(t):
    a=np.abs(t)
    return np.sign(t)*np.where(a<=1,2/3*a**1.5,a-1/3)
def gh(t): return np.minimum(np.sqrt(np.abs(t)),1.)
def source(A):
    return (lambda y:A*V.T@gp(V@y),lambda y:A*V.T@np.diag(gh(V@y))@V)
comm=0.
for A in [.5,.11,.004]:
 g,hess=source(A)
 for s in [1.,.23,.02]:
  a=np.array([.4,-.25]);x=a.copy()
  for _ in range(4): x=a-s*s*g(x)
  anchor=g(x);I=np.eye(2);P=np.concatenate([I,np.zeros((2,2))],axis=1)
  mats=[np.concatenate([rr*I,qq*I],axis=1) for rr,qq in zip(r,q)]
  def raw(W,der=False):
    shifts=[x+s*R@W for R in mats]
    bars=[s*(g(y)-anchor) for y in shifts]
    K=sum(ww*b for ww,b in zip(w,bars));Q=x+s*(P@W-K)
    f=(g(Q)-anchor)/s;b=(g(x+s*P@W)-anchor)/s
    GH=sum(ww/rr*R.T@bb for ww,rr,R,bb in zip(w,r,mats,bars))
    Jp=s*P.T@b+GH/s;Jm=s*P.T@b-GH/s
    if not der:return f,Jp,Jm
    DK=sum(ww*s*s*hess(y)@R for ww,y,R in zip(w,shifts,mats))
    Df=hess(Q)@(P-DK)
    DGH=sum(ww/rr*s*s*R.T@hess(y)@R for ww,rr,R,y in zip(w,r,mats,shifts))
    Db=hess(x+s*P@W)@P
    return Df,s*P.T@Db+DGH/s,s*P.T@Db-DGH/s,DK
  close(np.concatenate(raw(np.zeros(4))),np.zeros(10),'C2_exact_zero')
  for _ in range(5):
    W=rng.normal(size=4);v=rng.normal(size=4);v/=np.linalg.norm(v)
    df,djp,djm,dk=raw(W,True)
    outp=raw(W+1e-5*v);outm=raw(W-1e-5*v)
    for mat,pv,mv in zip([df,djp,djm],outp,outm):
      close(mat@v,(pv-mv)/2e-5,'C2_value_directional_first',tol=3e-7)
    close(djp,djp.T,'C2_Jplus_full_gradient')
    close(djm,djm.T,'C2_Jminus_full_gradient')
    curl=P.T@df-df.T@P
    check(op(curl)<=3*A*A*s*s,'C2_square_curl')
    check(op(dk)<=A*s*s*(1+1e-12),'C2_K_first')
    h1=hess(x+s*P@W);h2=hess(x+s*(P@W-.1))
    comm=max(comm,np.linalg.norm(h1@h2-h2@h1))
check(comm>1e-7,'genuine_noncommuting_hessians')
# A scalar random-covariance mixture is not its mean Gaussian: fourth moment.
# Variance 1 +/- c equiprobably has E X^4=3(1+c^2), unlike 3.
c=.2;check(abs(3*(1+c*c)-3)>.1,'mixture_negative_control')
out={'assertions':checks,'seed':620041004,'maxima':maxima,
     'max_C2_commutator':comm,'quadrature_nodes':n,
     'scope':'Algebra/first-bound diagnostics only; no LOW30 compiler execution or numerical proof of general conditional law bound.'}
path=Path(__file__).with_name('independent_reserve_algebra_results.json')
path.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
