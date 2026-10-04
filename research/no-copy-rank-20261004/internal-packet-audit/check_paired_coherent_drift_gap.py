"""Analytic limiting formulas plus exact native radial sufficient-statistic diagnostics."""
import math,json
from decimal import Decimal,localcontext,ROUND_HALF_EVEN
from pathlib import Path
import numpy as np
rng=np.random.default_rng(2026100431)
lam=.5;delta=beta=.025;c=s=2**-.5;t=.75;count=80000
checks=0
def ck(b):
 global checks
 checks+=1
 assert bool(b),checks

def gap(t):return 2*beta**2*math.exp(-1)*(math.cosh(c)-math.cosh(t*c))
true_lim=2*beta*lam*math.exp(-.5)+2*beta**2*math.exp(-1)*math.cosh(c)
split_lim=lambda t:2*beta*lam*math.exp(-.5)+2*beta**2*math.exp(-1)*math.cosh(t*c)
for tt in np.linspace(0,1,101):
 ck(abs(true_lim-split_lim(float(tt))-gap(float(tt)))<1e-17)
 ck(gap(float(tt))>=0)
 if tt<1:ck(gap(float(tt))>0)
clock_integral=2*beta**2*math.exp(-1)*(math.cosh(c)-math.sinh(c)/c)
z,w=np.polynomial.legendre.leggauss(40)
ck(abs(sum(ww*gap((zz+1)/2)/2 for zz,ww in zip(z,w))-clock_integral)<1e-18)
ck(clock_integral>0)

# Entire two-bank standard Gaussian covariance after common-root averaging.
Sigma=np.block([[np.eye(3),t*np.eye(3)],[t*np.eye(3),np.eye(3)]])
chol=np.linalg.cholesky(Sigma)
ck(np.linalg.norm(chol@chol.T-Sigma)<1e-14)
rows=[]
for NN in [200,1000,5000,20000]:
 with localcontext() as ctx:
  ctx.prec=70
  pi=Decimal('3.1415926535897932384626433832795028841971693993751058209749445923078164')
  aa=Decimal(NN)**(-Decimal(5)/Decimal(9))
  dd=(2*pi*Decimal(NN)/(Decimal('0.025')*aa))**2
  D=int(dd.to_integral_value(rounding=ROUND_HALF_EVEN));a=float(aa)
 eps=NN**(-.5);df=float(D-1)
 first=chol@rng.normal(size=(6,count))
 # Exact Wishart_6(D-1,Sigma) Bartlett sufficient statistics, in floating arithmetic.
 bart=np.zeros((count,6,6))
 for j in range(6):
  bart[:,j,j]=np.sqrt(rng.chisquare(df-j,size=count))
  for ell in range(j):bart[:,j,ell]=rng.normal(size=count)
 root=np.einsum('ij,njk->nik',chol,bart)
 gram=np.einsum('nik,njk->nij',root,root)
 def bank(off):
  S,U,Z=first[off:off+3];x=c*S+s*U
  SS=gram[:,off,off];UU=gram[:,off+1,off+1];ZZ=gram[:,off+2,off+2]
  SU=gram[:,off,off+1];SZ=gram[:,off,off+2]
  xx=c*c*SS+s*s*UU+2*c*s*SU
  RS=np.sqrt(1+S*S+SS);RSp=np.sqrt(1+(S+eps*Z)**2+SS+2*eps*SZ+eps*eps*ZZ)
  rad_diff=(2*eps*(S*Z+SZ)+eps*eps*(Z*Z+ZZ))/(RSp+RS)
  hh=a*(-lam*eps*Z-delta*rad_diff+delta*(S*S/RS-(S+eps*Z)**2/RSp)+beta*(np.sin(S)-np.sin(S+eps*Z)))
  RX=np.sqrt(1+x*x+xx);RXh=np.sqrt(1+(x+hh)**2+xx)
  rdiff=-hh*(2*x+hh)/(RX+RXh)
  d0=-lam*hh+delta*(rdiff-hh*(2*x+hh)/RX-(x+hh)**2*rdiff/(RX*RXh))+beta*(np.sin(x)-np.sin(x+hh))
  K=a*delta*math.sqrt(float(D))
  phase=(K-2*math.pi*NN)+S+a*(lam*x+delta*(RX-math.sqrt(float(D))-1+x*x/RX)+beta*np.sin(x))
  t0=2*math.pi*NN+phase;t1=t0-a*d0
  R0=np.sqrt(1+SS+t0*t0);R1=np.sqrt(1+SS+t1*t1)
  A0=lam+delta*(3*x/RX-x**3/RX**3)+beta*np.cos(x)
  A1=lam+delta*(3*(x+hh)/RXh-(x+hh)**3/RXh**3)+beta*np.cos(x+hh)
  T0=lam+delta*(3*t0/R0-t0**3/R0**3)+beta*np.cos(phase)
  T1=lam+delta*(3*t1/R1-t1**3/R1**3)+beta*np.cos(phase-a*d0)
  for H in [A0,A1,T0,T1]:ck(float(H.min())>=.4 and float(H.max())<=.6)
  return dict(S=S,x=x,h=hh,d0=d0,A0=A0,A1=A1,T0=T0,T1=T1)
 ba=bank(0);bt=bank(3)
 true=ba['T0']*ba['A0']-ba['T1']*ba['A1']
 split=bt['T0']*ba['A0']-bt['T1']*ba['A1']
 diff=true-split
 # Two algebraically identical paired expressions; do not independently regenerate i=0,1.
 ck(np.max(np.abs(true-(ba['T0']*(ba['A0']-ba['A1'])+(ba['T0']-ba['T1'])*ba['A1'])))<1e-15)
 ck(np.max(np.abs(split-(bt['T0']*(ba['A0']-ba['A1'])+(bt['T0']-bt['T1'])*ba['A1'])))<1e-15)
 true_sample_lim=2*beta*(lam+beta*np.cos(ba['S']))*np.cos(ba['x'])
 split_sample_lim=2*beta*(lam+beta*np.cos(bt['S']))*np.cos(ba['x'])
 row={
 'N':NN,'a':a,'dimension':str(D),'samples':count,
 'shift_L2_error':float(np.sqrt(np.mean((ba['h']+math.pi)**2))),
 'ancestor_difference_L2_error':float(np.sqrt(np.mean((ba['A0']-ba['A1']-2*beta*np.cos(ba['x']))**2))),
 'terminal0_L2_error':float(np.sqrt(np.mean((ba['T0']-lam-beta*np.cos(ba['S']))**2))),
 'terminal1_L2_error':float(np.sqrt(np.mean((ba['T1']-lam-beta*np.cos(ba['S']))**2))),
 'true_word_mean':float(true.mean()),'split_word_mean':float(split.mean()),
 'gap_sample_mean':float(diff.mean()),'gap_sample_standard_error':float(diff.std(ddof=1)/math.sqrt(count)),
 'gap_limit':gap(t),
 'true_word_to_limit_L2_error':float(np.sqrt(np.mean((true-true_sample_lim)**2))),
 'split_word_to_limit_L2_error':float(np.sqrt(np.mean((split-split_sample_lim)**2))),
 'paired_gap_to_coupled_limit_L2_error':float(np.sqrt(np.mean((diff-(true_sample_lim-split_sample_lim))**2)))}
 ck(row['shift_L2_error']<.01)
 rows.append(row)
ck(rows[-1]['gap_sample_mean']>0)
ck(rows[-1]['paired_gap_to_coupled_limit_L2_error']<rows[0]['paired_gap_to_coupled_limit_L2_error'])
out={'status':'PASS','checks':checks,'seed':2026100431,'clock_t':t,'fine_width':math.sqrt(1-t),
 'analytic_true_limit':true_lim,'analytic_split_limit':split_lim(t),'analytic_paired_gap':gap(t),
 'analytic_unit_clock_integral_gap':clock_integral,'native_diagnostics':rows,
 'scope':'Exact explicit Hessian-value limits prove the theorem; finite Wishart samples only diagnose actual native sufficient statistics. No certified finite N0, orientation lower bound from the word gap, or universal producer impossibility is claimed.'}
Path(__file__).with_name('paired_coherent_drift_gap_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
