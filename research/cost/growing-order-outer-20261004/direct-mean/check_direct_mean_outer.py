#!/usr/bin/env python3
"""Exact accounting plus diagnostic finite maps; not a replacement for proof."""
from fractions import Fraction as F
import json, math, hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
checks=[]
def chk(name, ok, details=None):
    assert bool(ok), (name, details)
    checks.append({'name':name,'pass':True,**({'details':details} if details is not None else {})})
# Entire direct producer grammar. No coefficient/root rows are generated.
for j in range(2,61):
    P=2*j-1; source_c=2*j-2
    m,v,g,z=F(0),F(2*j-1,2),F(1,2),F(0)
    b=3*m+2*v+g+z
    I=b-F(1,2); prefix=I; final=b+2*m+2*g+F(1,2)
    chk(f'j{j}:literal_direct_grades',(I,prefix,final)==(P,P,P+2))
    chk(f'j{j}:reference_capacity',2*j>P)
    chk(f'j{j}:reset_capacity',2*j+1==P+2)
    chk(f'j{j}:phase_cost',source_c==P-1)
    chk(f'j{j}:total_cost',1+2*source_c==2*P-1)
    for l in range(1,j+1):
        chk(f'j{j}:bank{l}:complete_fine',2*(j-l)+2*(l-1)==source_c)
        chk(f'j{j}:bank{l}:complete_coarse',2*(j-l)+2*max(l-2,0)<=source_c)
        chk(f'j{j}:bank{l}:strong',l+(j-l)==j)
    # Base variance scale A and allocation A^-2(j-1).
    chk(f'j{j}:base_variance',1+(j-1)==j)
    dim_exp=F(2*P-1,2*P); kap_exp=1+F(2*P-1,P)
    chk(f'j{j}:global_rate',dim_exp==1-F(1,2*P) and kap_exp==3-F(1,P))
# Fixed depth recurrence, including complete independent lower-depth sources.
for maxj in (2,3,5,9,17):
    c=[0]+[2*(j-1) for j in range(1,maxj+1)]
    for depth in range(1,3*maxj+1):
        nxt=[0]
        for j in range(1,maxj+1):
            bill=[j-1,2*(j-1)+c[0]]
            for l in range(1,j+1): bill += [2*(j-l)+c[l],2*(j-l)+c[l-1]]
            nxt.append(max(bill))
        chk(f'cost_depth{depth}_maxj{maxj}',nxt==c)
        c=nxt
# All actual epochs share comparable heat; reset contraction guard.
for kap in (1,2,10,100):
    for L in (10,20,100):
        T=L**-2; H=1e-3/(kap*L**3); n=math.ceil(T/H); h=T/n
        A0=L**4*H**2; Afinal=A0+1e-3*n*h**3
        chk(f'epoch:k{kap}:L{L}:heat',Afinal/A0<=1.25)
        chk(f'epoch:k{kap}:L{L}:row',h*h/A0<=L**-4*(1+1e-12))
        chk(f'epoch:k{kap}:L{L}:reset',A0<T/(100*kap))
# Interpolation includes factorial smoothing loss, not only Picard defect.
for L in (10,20,40,80):
    J=math.ceil(2*L*L); Ca=12; ratio=L**(-Ca/2)
    log_err=(J+1)*math.log(4)+.5*math.lgamma(J+1)+J*math.log(ratio)
    chk(f'interpolation:L{L}',log_err < -10*L, {'log_bound':log_err})
# Scalar finite quarter, differentiating actual recorded queries.
def grad(x): return .9*x-.4*np.sign(x)*np.log1p(np.abs(x))
def hess(x): return .9-.4/(1+np.abs(x))
def mode(a,y,steps=40):
    x=float(y); dy=1.
    for _ in range(steps): dy=1-a*hess(x)*dy; x=y-a*grad(x)
    return x,dy

def quarter(a,y,x0,z,N,M,dx0=1.,dy0=0.):
    t=np.linspace(0,math.pi/2,N+1)
    co=np.cos(t); si=np.sin(t); co[-1]=0.; si[-1]=1.
    W=np.zeros((N+1,N+1))
    for i in range(1,N+1):
        q=np.arange(i); W[i,:i]=np.cos(t[i]-t[q+1])-np.cos(t[i]-t[q])
    base=y+co*(x0-y)+math.sqrt(a)*si*z
    x=base.copy(); dx=co*dx0; dy=1-co+co*dy0
    for _ in range(M):
        H=hess(x)
        dx=co*dx0-a*W@(H*dx)
        dy=1-co+co*dy0-a*W@(H*dy)
        x=base-a*W@grad(x)
    return float(x[-1]),float(dx[-1]),float(dy[-1])

rng=np.random.default_rng(2601004)
max_input_ratio=0.; max_center_ratio=0.; seed_residual_ratio=0.
for a in (.005,.02,.08):
    for y in (-20.,-.3,0.,.6,30.):
        b,db=mode(a,y)
        chk(f'mode:a{a}:y{y}:residual',abs(b-y+a*grad(b))<1e-10*(1+abs(y)))
        chk(f'mode:a{a}:y{y}:caller',abs(db-1)<=a/(1-a)+1e-12)
        for N in (3,7,17):
            for M in (2,3,5):
                z0,z=rng.normal(size=2); x0=b+math.sqrt(a)*z0
                x,dx,dy=quarter(a,y,x0,z,N,M)
                max_input_ratio=max(max_input_ratio,abs(dx)/(a/(1-a)))
                chk(f'quarter:a{a}:y{y}:N{N}:M{M}:input',abs(dx)<=a/(1-a)+1e-12)
                # Caller derivative includes seed derivative.
                _,_,dtotal=quarter(a,y,x0,z,N,M,0.,db)
                max_center_ratio=max(max_center_ratio,abs(dtotal-1)/a)
                chk(f'quarter:a{a}:y{y}:N{N}:M{M}:caller',abs(dtotal-1)<=3*a)
                eps=1e-5
                xp=quarter(a,y,x0+eps,z,N,M)[0]; xm=quarter(a,y,x0-eps,z,N,M)[0]
                chk(f'quarter:a{a}:y{y}:N{N}:M{M}:finite_diff',abs((xp-xm)/(2*eps)-dx)<2e-7)
                # Unanchored finite Picard layer has a real deterministic caller residual.
                zz=quarter(a,y,b,0.,N,M)[0]
                chk(f'quarter:a{a}:y{y}:N{N}:M{M}:zero_residual',abs(zz-b)<=a**(M+1)*abs(grad(b))/(1-a)+1e-12*(1+abs(y)))
# Complete bank Gaussian fixture at huge ambient sample dimension: output variance depends on physical dimension.
for j in range(2,8):
    a=F(1,4); n0=a**(-2*(j-1))
    vari=a*a/n0
    for l in range(1,j+1):
        n=a**(-2*(j-l)); vari+=(a**(2*l))/n
    chk(f'bank_variance_j{j}',vari==(j+1)*a**(2*j))
# Hidden decoder reference covariance exact, any finite length.
for k in (1,2,5,17,100):
    Q=np.eye(k+1)-np.ones((k+1,k+1))/(k+1)
    C=np.ones((k+1,k+1))/(k+1)+Q
    chk(f'decoder_covariance_k{k}',np.max(np.abs(C-np.eye(k+1)))<1e-14)
report={'claim':'Exact accounting and finite C2/non-C3 quarter/mean diagnostics only; not a full sampler execution',
        'assertions':len(checks),'all_pass':True,'max_quarter_input_contraction_ratio':max_input_ratio,
        'max_quarter_center_residual_over_A':max_center_ratio,'checks':checks}
(ROOT/'direct_mean_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='checks'},indent=2))
