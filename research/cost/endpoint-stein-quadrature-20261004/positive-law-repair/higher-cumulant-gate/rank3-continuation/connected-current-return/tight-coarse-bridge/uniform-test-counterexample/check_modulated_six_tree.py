#!/usr/bin/env python3
"""Finite diagnostics for a printed analytic counterexample, not a proof by sampling."""
import hashlib, itertools, json, math
from pathlib import Path
import numpy as np
import sympy as sp


# Publication-only source locator and original/public pin mapping.
import sys
_ARCHIVE=next(p for p in Path(__file__).resolve().parents if (p/'publication_sources.py').is_file())
sys.path.insert(0,str(_ARCHIVE))
from publication_sources import source_path, source_label, verify_source_pin, pin_report

HERE=Path(__file__).resolve().parent
checks=[]
def check(label, ok, **data):
    if not bool(ok):
        raise AssertionError((label,data))
    checks.append({'name':label,**data})

kappa=.5; eta=.125
# Independently differentiate the original scalar potential in a symbolic model.
x0,x1,x2,sig=sp.symbols('x0 x1 x2 sig', real=True, nonzero=True)
xs=[x0,x1,x2]
u=sp.Rational(1,4)*sum(x*x for x in xs)+sp.Rational(1,8)*sig**2*sp.cos(x0)*(sp.cos(x1/sig)+sp.cos(x2/sig))
check('symbolic Hessian ii', sp.simplify(sp.diff(u,x1,2)-(sp.Rational(1,2)-sp.Rational(1,8)*sp.cos(x0)*sp.cos(x1/sig)))==0)
check('symbolic Hessian 0i', sp.simplify(sp.diff(u,x0,x1)-sp.Rational(1,8)*sig*sp.sin(x0)*sp.sin(x1/sig))==0)
check('symbolic C2 iiii', sp.simplify(sig**2*sp.diff(u,x1,4)-sp.Rational(1,8)*sp.cos(x0)*sp.cos(x1/sig))==0)
check('symbolic C2 iii0', sp.simplify(sig**2*sp.diff(u,x1,3,x0)+sp.Rational(1,8)*sig*sp.sin(x0)*sp.sin(x1/sig))==0)
check('symbolic forbidden mixed center index',sp.simplify(sp.diff(u,x1,3,x2))==0)

rng=np.random.default_rng(721365)
def hess(n,x):
    s=n**-.5
    E=np.zeros((n+1,n+1))
    E[0,0]=-eta*s*s*np.cos(x[0])*np.cos(x[1:]/s).sum()
    E[np.arange(1,n+1),np.arange(1,n+1)]=-eta*np.cos(x[0])*np.cos(x[1:]/s)
    E[0,1:]=E[1:,0]=eta*s*np.sin(x[0])*np.sin(x[1:]/s)
    return kappa*np.eye(n+1)+E

def center(n,x):
    s=n**-.5; D=n+1; d=math.exp(-(1+s*s)/2)
    T=np.zeros((D,)*4)
    for i in range(1,D):
        for sign in [-1,1]:
            v=np.zeros(D); v[0]=1; v[i]=sign/s
            T += eta*s**4*d/2*np.cos(x[0]+sign*x[i]/s)*np.einsum('a,b,c,d->abcd',v,v,v,v)
    return T

def cuts(T):
    r=T.ndim
    for mask in range(1,2**r-1):
        # A and complement have equal singular values; inspect one representative.
        if not mask&1: continue
        a=[i for i in range(r) if mask>>i&1]; b=[i for i in range(r) if not(mask>>i&1)]
        yield mask,np.linalg.norm(T.transpose(a+b).reshape(T.shape[0]**len(a),-1),2)

for n in [3,8,16,64]:
    for trial in range(8):
        x=rng.normal(size=n+1)*rng.uniform(.1,20)
        eig=np.linalg.eigvalsh(hess(n,x))
        check(f'global Hessian fixture n={n} sample={trial}',eig.min()>=.25-1e-12 and eig.max()<=.75+1e-12)

for n in [3,5,8]:
    x=rng.normal(size=n+1)
    T=center(n,x)
    for mask,norm in cuts(T): check(f'C2 proper cut n={n} mask={mask}',norm<=2+1e-11)
    d=math.exp(-(1+1/n)/2)
    for i in range(1,n+1):
        check(f'C2 diagonal n={n} i={i}',abs(T[i,i,i,i]-eta*d*np.cos(x[0])*np.cos(x[i]*np.sqrt(n)))<1e-12)
        check(f'C2 cross n={n} i={i}',abs(T[i,i,i,0]+eta*d/np.sqrt(n)*np.sin(x[0])*np.sin(x[i]*np.sqrt(n)))<1e-12)

# Full finite tensor, including the actual four heat-smoothed leaves.
n=4; D=n+1
roots=rng.normal(size=(6,D))
T,Tp=center(n,roots[0]),center(n,roots[1])
leaf_d=math.exp(-(n+1)/8)
Bs=[kappa*np.eye(D)+leaf_d*(hess(n,r)-kappa*np.eye(D)) for r in roots[2:]]
gamma=.75*(1-2/n)
def tree(bs):
    b1,b3,b1p,b3p=bs
    return gamma*np.einsum('au,cv,buvh,dw,fz,ewzh->abcdef',b1,b3,T,b1p,b3p,Tp,optimize=True)
J=tree(Bs); J0=tree([kappa*np.eye(D)]*4)
H=np.zeros((D,)*6)
for i in range(1,D): H[(i,)*6]=1/np.sqrt(n)
check('HS unit diagonal six-test',abs(np.linalg.norm(H.ravel())-1)<1e-14)
manual=0
for i in range(1,D):
    manual += kappa**4*gamma*sum(T[i,i,i,h]*Tp[i,i,i,h] for h in range(D))/np.sqrt(n)
check('complete six-tree diagonal contraction',abs(np.vdot(H,J0)-manual)<1e-13)
for mask,norm in cuts(J): check(f'six-tree proper cut mask={mask}',norm<=4*gamma+1e-10)
for mask,norm in cuts(J-J0): check(f'leaf difference proper cut mask={mask}',norm<=4*leaf_d+1e-10)
check('leaf difference HS error',np.linalg.norm((J-J0).ravel())<=4*leaf_d*np.sqrt(D)+1e-10)

# Exact Gaussian moments, with eps=1-t carried separately to avoid rounding near one.
def moments(v,eps):
    ep=math.exp(-v*eps); em=math.exp(-v*(2-eps))
    e2=math.exp(-2*v)
    e4p=math.exp(-4*v*eps); e4m=math.exp(-4*v*(2-eps))
    meanM=(ep+em)/2; meanN=(ep-em)/2
    M2=(1+2*e2+(e4p+e4m)/2)/4
    N2=(1-2*e2+(e4p+e4m)/2)/4
    MN=(e4p-e4m)/8
    return meanM,meanN,M2,N2,MN

records=[]
for n in [16,32,64,128,256,512,1024,4096,16384,65536,10**6]:
    sigma2=1/n; v=15/16-7/(8*n); gamma=.75*(1-2/n)
    c=gamma*kappa**4*eta**2*math.exp(-(1+sigma2))
    for panel_fraction in [0,.125,.5,1,2]:
        eps=panel_fraction/n**3
        M,N,M2,N2,MN=moments(v,eps)
        mu,nu,X2,Y2,XY=moments(v/sigma2,eps)
        a=mu; b=sigma2*nu
        var_condmean=a*a*M2+b*b*N2+2*a*b*MN-(a*M+b*N)**2
        noise=M2*(X2-mu**2)+sigma2**2*N2*(Y2-nu**2)+2*sigma2*MN*(XY-mu*nu)
        var_ideal=c*c*(n*var_condmean+noise)
        check(f'positive exact conditional variance n={n} f={panel_fraction}',var_condmean>=0 and var_ideal>0)
        check(f'printed sd coefficient bound n={n} f={panel_fraction}',var_condmean>=1/32**2 and c>=2**-14)
        error=4*math.exp(-(n+1)/8)*math.sqrt(n+1)
        robust=max(0,math.sqrt(var_ideal)-error)**2
        if n>=128:
            check(f'robust actual-leaf variance lower n={n} f={panel_fraction}',robust>=n/2**40,lower=robust,target=n/2**40)
        records.append({'n':n,'price_panel_fraction':panel_fraction,'var_ideal_over_n':var_ideal/n,'robust_actual_lower_over_n':robust/n})

v=15/16; gamma=.75
limit=(gamma*kappa**4*eta**2/math.e)**2*(1-math.exp(-4*v))**2/32
check('nonzero analytic limiting variance',limit>0,value=limit)
for rec in records[-5:]: check('large-n exact moment agrees with limit '+str(rec['price_panel_fraction']),abs(rec['var_ideal_over_n']/limit-1)<1e-4)
check('explicit leaf-error constant at n=128',6*math.exp(-129/8)<=2**-20)
# An arbitrary dimension multiplier dominates any fixed public log power eventually;
# the mathematical statement is elementary, not inferred from these fixtures.
sources=[
    HERE.parent/'SEALED-MANIFEST.json',
    HERE.parent/'COVARIANCE-TEST-PORT-REDUCTION.md',
    HERE.parent/'TIGHT-BRIDGE-LEDGER-AND-REROOT-GATE.md',
    HERE.parents[2]/'EXACT-FIVE-CLOCK-NATIVE-SOURCE-ADAPTER.md',
    HERE.parents[2]/'grouped-kernel-response'/'PAIRED-ROTATION-COVARIANCE-ALL-CUT-LEMMA.md',
]
pins={source_label(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
check('sealed source manifest unchanged',pins[source_label(sources[0])]=='c42ce2f9f086ae94abe59fe744af90dc01127221459614fe868c118314448179')
report={'status':'PASS','assertions':len(checks),'limiting_variance_over_n':limit,'scope':'Diagnostics supporting the printed original-gradient counterexample. Not a positive producer or general-law no-go.','source_pins':pins,'checks':checks,'exact_moment_records':records}
(HERE/'checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'status':report['status'],'assertions':report['assertions'],'limiting_variance_over_n':limit,'report':str(HERE/'checks.json')}))
