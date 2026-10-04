#!/usr/bin/env python3
"""Analytic Fourier and clock regressions for the held shear transfer.
No sampled Hessian is an emitted source. Nonlinear Fourier truncation is
explicitly bounded and used only for coefficient diagnostics.
"""
from pathlib import Path
import hashlib,json,math
import numpy as np
from scipy.special import jv
from numpy.polynomial.legendre import leggauss
import mpmath as mp
mp.mp.dps=70
ROOT=Path(__file__).resolve().parents[1];OUT=Path(__file__).resolve().parent
checks=[]
def check(name,ok,**kw):
 checks.append({'name':name,'pass':bool(ok),**kw})
 if not ok:raise AssertionError((name,kw))
def norm2(x):return float(np.dot(x,x))
def sin_cov(w,v):return .5*(math.exp(-.5*norm2(w-v))-math.exp(-.5*norm2(w+v)))
def canonical(pair):
 a,b=pair
 return (-a,-b) if a<0 or (a==0 and b<0) else (a,b)
def nonlinear_source(a,r,h,d,lam,k,s,N=14):
 kap=a*r;den=lam*lam+k*k;d0=d*k*k/den;t=lam*a*d0/k
 records=[((0,1),np.array([kap*h*d0/k,0.]))]
 wave=np.array([lam,k])
 for m in range(-N,N+1):
  coef=float(mp.besselj(0,mp.mpf(t))-1) if m==1 else float(jv(m-1,t))
  records.append(((1,m),r*d*wave/den*coef))
 e2=0.
 for fi,vi in records:
  wi=np.array([fi[0]*lam,fi[1]*k])
  for fj,vj in records:
   wj=np.array([fj[0]*lam,fj[1]*k])
   e2+=float(vi@vj)*sin_cov(wi,wj)
 e=math.sqrt(max(e2,0.))
 J=[]
 for f,vec in records:
  w=np.array([f[0]*lam,f[1]*k]);mat=np.outer(vec,w)*math.exp(-s*s*norm2(w)/2)
  if np.any(mat):J.append((f,mat))
 O={}
 for fi,Ji in J:
  for fj,Jj in J:
   X=Ji@(Jj-Jj.T);mat=-(X+X.T)/4
   for pair in ((fi[0]+fj[0],fi[1]+fj[1]),(fi[0]-fj[0],fi[1]-fj[1])):
    pair=canonical(pair);O[pair]=O.get(pair,np.zeros((2,2)))+mat
 c=math.sqrt(1-s*s);o2=0.
 for fi,Oi in O.items():
  wi=c*np.array([fi[0]*lam,fi[1]*k])
  for fj,Oj in O.items():
   wj=c*np.array([fj[0]*lam,fj[1]*k])
   ec=.5*(math.exp(-norm2(wi-wj)/2)+math.exp(-norm2(wi+wj)/2))
   o2+=float(np.sum(Oi*Oj))*ec
 # |J_n(t)| <= exp(t²/4)(|t|/2)^n/n!, omitted n>=N conservatively.
 x=abs(t)/2;tail=2*math.exp(x*x)*x**N/math.factorial(N)/(1-x/(N+1))
 return e,math.sqrt(max(o2,0)),d0,t,tail

nonlinear=[]
a=.2;r=.15;h=.5;d=.1;kap=a*r
for s in (.03,.1,.3,.7):
 for b in (.03,.1,.3):
  eta=b*s*s
  for ls in (.25,1.,4.):
   lam=ls/s
   for ke in (.3,1.,3.):
    k=ke/eta
    e,onorm,d0,t,tail=nonlinear_source(a,r,h,d,lam,k,s)
    g0norm=d0/k*math.sqrt(-math.expm1(-2*k*k)/2)
    check('same_potential_source_energy_sandwich',kap*(h-d)*g0norm*(1-1e-9)<=e<=kap*(h+d)*g0norm*(1+1e-9))
    check('nonlinear_Bessel_tail_certificate',tail<1e-40)
    err=-math.expm1(-k*k*eta*eta/2)*onorm
    allowed=2*kap*eta*e/(s*s)
    check('nonlinear_matched_transfer',err<=allowed*(1+1e-9))
    # Exact elementary z-Lipschitz bound of the transverse heat Hessian.
    Dz=d*lam*lam/(lam*lam+k*k)*math.exp(-lam*lam*s*s/2)*(abs(k)+abs(lam*a*d0))
    check('transverse_heat_first_no_ancestor_frequency_loss',Dz<=(1+a)/s)
    nonlinear.append({'s':s,'b':b,'eta':eta,'lambda_s':ls,'k_eta':ke,'e':e,'error':err,'error_over_kappa_eta_e_s2':err/(kap*eta*e/(s*s)),'Bessel_tail_bound':tail})

quadratic=[]
for s in (.05,.2,.7):
 for b in (.02,.1,.4):
  eta=b*s*s
  for regime in ('k_s_one','k_eta_one','k_eta_ten'):
   k={'k_s_one':1/s,'k_eta_one':1/eta,'k_eta_ten':10/eta}[regime]
   amp=kap*h*d;e=amp/k*math.sqrt(-math.expm1(-2*k*k)/2)
   c=math.sqrt(1-s*s);cos4=(3+4*math.exp(-2*k*k*c*c)+math.exp(-8*k*k*c*c))/8
   err=amp*amp*math.exp(-k*k*s*s)*(-math.expm1(-k*k*eta*eta/2))*math.sqrt(cos4)
   check('exact_quadratic_terminal_transfer',err<=2*kap*eta*e/(s*s))
   quadratic.append({'s':s,'b':b,'regime':regime,'normalized_error':err/(kap*eta*e/(s*s))})
# Same source, but deliberately INVALID replacement of conditional J by raw Df.
s=.25;eta=.001;k=1/eta;amp=kap*h*d
e=amp/k*math.sqrt(-math.expm1(-2*k*k)/2)
wrong=amp*amp*(-math.expm1(-k*k*eta*eta/2))*(1+math.exp(-2*k*k))/2
wrong_ratio=wrong/(kap*eta*e/(s*s))
check('raw_same_Z_mark_is_a_real_counterexample_to_illegal_genealogy',wrong_ratio>100)
# A retained ancestor observer likewise destroys the unqualified weak transfer.
observer_loss=(-math.expm1(-.5))*(1+math.exp(-2*k*k*s*s))/2
check('unpriced_fast_ancestor_observer_not_covered',observer_loss>.19)

clock=[]
for delta in (.2,.1,.03,.01,.003):
 t0=delta*delta/128;T=.5*math.log(16/delta)
 J=math.ceil(math.log2(T/t0));M=math.ceil(math.log(384/delta,4))
 x,ww=leggauss(M);ts=[];ws=[]
 for panel in range(J):
  aa=2**panel*t0
  tt=1.5*aa+.5*aa*x;weights=.5*aa*ww
  ts.extend(tt);ws.extend(2*weights*np.exp(-2*tt))
 ts=np.array(ts);ws=np.array(ws);ss=np.sqrt(-np.expm1(-2*ts))
 check('clock_positive',np.all(ws>0))
 check('clock_total_mass',sum(ws)<=1+1e-12)
 check('clock_w_over_s2',sum(ws/(ss*ss))<=J*(1+1e-10))
 check('clock_individual_weight',np.all(ws<=3*ss*ss))
 # Actual sqrt-weight source path bound, with its declared logarithmic factor.
 path=sum(np.sqrt(ws)/(ss*ss));scaled=path*min(ss)/math.sqrt(M)
 check('clock_inverse_eta_path_no_extra_power',scaled<10)
 b=delta;rho=.2
 weighted=kap/(b*rho)*sum(np.sqrt(ws)/(b*ss*ss))
 formula=kap/(b*b*rho)*path
 check('literal_weighted_side_path_identity',abs(weighted-formula)<=1e-12*(1+formula))
 calibration=sum(ws*2*kap*(b*ss*ss)/(ss*ss))
 check('literal_weighted_transfer_identity',abs(calibration-2*kap*b*sum(ws))<=1e-12)
 clock.append({'delta':delta,'nodes':len(ws),'panels':J,'M':M,'s_min':float(min(ss)), 'sum_w':float(sum(ws)),'sum_w_over_s2':float(sum(ws/(ss*ss))), 'sum_sqrtw_over_s2':float(path),'path_times_smin_over_sqrtM':float(scaled)})
# Parameter dependence: smaller clock tolerance really changes inverse-width cost.
p01=next(x for x in clock if x['delta']==.1)['sum_sqrtw_over_s2']
p001=next(x for x in clock if x['delta']==.01)['sum_sqrtw_over_s2']
check('smaller_delta_has_real_extra_path_cost',p001/p01>8)
# Exact genealogy/sign identities, numerically at unrelated nonsymmetric marks.
rng=np.random.default_rng(419)
for dim in (2,3,9,32):
 for trial in range(8):
  uv=np.linalg.qr(rng.normal(size=(dim,2)))[0];u=uv[:,0];v=uv[:,1]
  Pi=np.eye(dim)-np.outer(v,v);Z=rng.normal(size=(dim,dim))
  H=Pi@((Z+Z.T)/2)@Pi;H/=max(1,np.linalg.norm(H,2))
  scalar=float(rng.uniform(-1,1));D=H+scalar*np.outer(v,v)
  knownJ=np.outer(u,v)-np.outer(v,u);fill=np.eye(dim)-np.outer(u,u)-np.outer(v,v)
  W=scalar*(np.outer(H@u,v)-np.outer(v,u@H))
  check('full_square_DJD_skew_identity',np.linalg.norm(D@knownJ@D-W)<1e-12)
  check('known_middle_fill_covariance',np.linalg.norm(knownJ@knownJ.T+fill@fill.T-np.eye(dim))<1e-12)
  marker=rng.normal(size=(dim,dim));tau=.03;rho=.12;csel=.2
  O=-tau*(marker@W+(marker@W).T)/2
  M=csel*marker;K=csel*csel*rho*rho*W
  cross=M@K.T+K@M.T
  check('same_passive_cross_orientation_and_transposes',np.linalg.norm(cross-2*csel**3*rho*rho*O/tau)<1e-12)
  q=.7;bb=.1;cc=-q*tau/(2*csel**3*rho*rho*bb)
  check('negative_fork_readout_sign',np.linalg.norm(bb*cc*cross+q*O)<1e-12)
  # This unscaled random case checks algebra, not the positivity budget.

paths=['p-shear-heat/NONLINEAR-SHEAR-CONDITIONAL-HEAT-CANDIDATE.md','r-true-orientation/POSITIVE-COMMON-INPUT-ORIENTATION-CLOCK.md','p-shear-heat/FULL-SQUARE-SHEAR-SOURCE-GENEALOGY-AND-COST.md']
out={'all_pass':all(x['pass'] for x in checks),'count':len(checks),'checks':checks,
 'nonlinear_fourier_cases':nonlinear,'quadratic_exact_formula_cases':quadratic,'clock_cases':clock,
 'invalid_raw_mark_ratio':wrong_ratio,'invalid_retained_ancestor_loss':observer_loss,
 'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
 'scope':'Coefficient-transfer audit only. Fourier truncation diagnostics have explicit tail bounds; executable retained-host closure is not inferred.'}
(OUT/'shear_transfer_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'all_pass':out['all_pass'],'checks':len(checks),'nonlinear_cases':len(nonlinear),'quadratic_cases':len(quadratic),'clock_cases':len(clock),'max_nonlinear_normalized_error':max(x['error_over_kappa_eta_e_s2'] for x in nonlinear),'invalid_raw_mark_ratio':wrong_ratio},indent=2))
