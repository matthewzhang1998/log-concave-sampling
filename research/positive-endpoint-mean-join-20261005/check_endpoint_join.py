#!/usr/bin/env python3
"""Independent algebraic diagnostics of the stated join; no native compiler is simulated."""
from pathlib import Path
from fractions import Fraction as F
import json,math,hashlib,random
import numpy as np
P=Path(__file__).resolve().parent
counts={}
def ck(c,key):
    if not c: raise AssertionError(key)
    counts[key]=counts.get(key,0)+1
p,q=F(39,10),F(59,15)
h=F(1,2); rho=F(3,4)
ck(q-p==F(1,30),'residual_margin')
ck(2*p-1==F(34,5),'local_grade')
ck(2*q-1==F(103,15),'local_grade')
ck(2*(p-2)/5==F(19,25),'stop_exponent')
ck(p+F(1,2)==F(22,5),'physical_readout')
ck(sum([h*h,h*h*F(1,4),h*h*F(1,4),F(1,8),F(1,4)])==1-h*h,'positive_budget')
ck(h*h+h*h*F(1,2)+F(1,8)==F(1,2),'visible_covariance')
ck(F(1,4)/(1-h*h)==F(1,3),'keep_fraction')
ck(h*h/F(1,8)==2,'cubic_radius')
ck(-(h**3)/6==-F(1,48),'cubic_sign')
for i in range(0,201):
    # parameter r^2 lets all bridge covariance identities remain rational.
    r2=F(i,202); s2=1-r2; D=s2/4; t2=r2+D
    v0=D/t2; v=(1-t2)*D/(t2*s2)
    ck(v==v0*(1-h*h),'bridge_noise')
    ck(D**2/(t2*s2)==v0*h*h,'bridge_signal_squared')
    # t*(a+b*r)=r checked after factoring r, valid including r=0.
    ck((1-t2)/s2+D/s2==1,'exact_contraction')
    ck(s2>0 and D>0 and v>0,'positive_bridge_gaps')
    ck((D/s2)*s2**4==D*s2**3,'mode_readout_power')
# Positive finite schedule and actual telescoping weights.
for exponent in np.linspace(.01,50,400):
    A=10.0**(-float(exponent)); cutoff=A**float(F(19,25)); s2=1.;j=0;alpha_min=A
    ksum={k:0. for k in [float(p),float(q),4.]};states=[]
    while s2>cutoff:
        prev=s2; s2*=.75;j+=1;alpha_min=min(alpha_min,A*prev)
        d=prev/4;r=math.sqrt(1-prev);t=math.sqrt(1-s2)
        states.append((r,t))
        for k in ksum: ksum[k]+=d*prev**((2*k-1)/2)
    ck(s2<=cutoff and s2>.75*cutoff,'minimal_terminal_cutoff')
    ck(s2**2.5<=A**float(p-2)*(1+1e-13),'terminal_grade')
    ck(j<=math.ceil(float(F(19,25))*math.log(1/A)/math.log(4/3))+1,'stage_count')
    ck(math.log(1/alpha_min)<=float(F(44,25))*math.log(1/A)*(1+1e-12),'uniform_local_log_scale')
    for k,sm in ksum.items():
        ck(sm <=1/(4*(1-.75**(k+.5)))*(1+1e-13),'weighted_geometric_sum')
    for k0 in [0,len(states)//2,len(states)-1]:
        weight=states[-1][1]
        for k in range(k0+1,len(states)):weight*=states[k][0]/states[k][1]
        ck(abs(weight-states[k0][1])<=2e-13,'propagated_exact_weight')
# Noncommuting anisotropic conditional Gaussian regression identity.
rng=np.random.default_rng(163741)
for dim in range(1,9):
  for _ in range(40):
    B=rng.normal(size=(dim,dim));Sig=.01*(B@B.T)/max(1,dim)
    C=.5*np.eye(dim)+.25*Sig;Ci=np.linalg.inv(C)
    beta2=.03+rng.uniform(0,.08)
    S=rng.normal(size=dim)
    mean=beta2*Ci@S
    var=beta2*np.eye(dim)-beta2**2*Ci
    lhs=(np.outer(mean,mean)+var)/beta2**2-np.eye(dim)/beta2
    rhs=np.outer(Ci@S,Ci@S)-Ci
    ck(np.allclose(lhs,rhs,rtol=1e-11,atol=1e-11),'anisotropic_Hermite_regression')
    ck(np.linalg.eigvalsh(var).min()>0,'private_regression_covariance')
    ck(np.linalg.eigvalsh(C).min()>=.5-1e-13,'visible_positive_gap')
    # One-dimensional output-slot asymmetry has no defect; higher dimensions test tensor symmetrization.
    N=rng.normal(size=(dim,dim,dim));N=(N+N.transpose(0,2,1))/2
    import itertools
    Sym=sum(N.transpose(perm) for perm in itertools.permutations(range(3)))/6
    D=N-Sym
    x=rng.normal(size=dim)
    ck(abs(np.einsum('ijk,i,j,k',D,x,x,x))<1e-9,'zero_symmetric_leading_current')
# Exact OU resolvent moments: split the double integral at the time diagonal.
def conditional_ou_moment(m,n):
    return (F(math.factorial(m+n+1),2**(m+n+2)*math.factorial(m)*math.factorial(n))
            *(F(1,n+1)+F(1,m+1))-F(1,2**(m+n+2)))
expected={(0,0):F(1,4),(0,1):F(1,4),(0,2):F(3,16),
          (1,1):F(5,16),(1,2):F(9,32),(2,2):F(19,64)}
for indices,value in expected.items():
    ck(conditional_ou_moment(*indices)==value,'exact_conditional_history_moment')
# Quadratic same-history target negative control: canonical mean is not posterior-force mean.
for a in [F(i,100) for i in range(1,26)]:
    m3=a/2-a*a/4+a**3/8
    post=a*F(1,2)/(1+a*F(3,4))
    ck(m3!=post,'canonical_vs_posterior_negative_control')
    cov2=a*a/4-a**3/2+F(5,16)*a**4
    cov3=cov2+2*conditional_ou_moment(0,2)*a**4-2*conditional_ou_moment(1,2)*a**5+conditional_ou_moment(2,2)*a**6
    # Exact conditional Gaussian moments of H,J,K yield difference beginning at degree four.
    ck(abs(float(cov3-cov2))<=float(a**4),'covariance_target_grade')
# Source pins prove lineage unchanged, not that imported compilers have numerically run.
pins=json.loads((P/'INPUT-PINS.json').read_text())['inputs']
for pin in pins:
    ck(hashlib.sha256(Path(pin['path']).read_bytes()).hexdigest()==pin['sha256'],'input_pin')
result={'status':'PASS','assertions':sum(counts.values()),'checks':counts,'native_compilers_executed':False,'scope':'exact bridge/budget/exponent algebra, finite schedule and anisotropic Gaussian reference diagnostics; mathematical proof and literal native admission required'}
(P/'endpoint_join_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
