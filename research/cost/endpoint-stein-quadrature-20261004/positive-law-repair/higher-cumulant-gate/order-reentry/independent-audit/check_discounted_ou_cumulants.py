#!/usr/bin/env python3
"""Independent polynomial and true-Markov checks of the analytical hierarchy."""
from pathlib import Path
import hashlib
import itertools as it
import json
import math
import numpy as np
import sympy as s


# Publication-only source locator and original/public pin mapping.
import sys
_ARCHIVE=next(p for p in Path(__file__).resolve().parents if (p/'publication_sources.py').is_file())
sys.path.insert(0,str(_ARCHIVE))
from publication_sources import source_path, source_label, verify_source_pin, pin_report

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'DISCOUNTED-OU-CONNECTED-CUMULANT-RECURRENCE.md'
PINS={str(SOURCE):'cd37a782bbcef3cdf5ec39a33be4874df56ff06af75170485b182c04b361d05e',
      str(HERE.parent/'FIRST-RESOLVENT-CLOCK-QUADRATURE.md'):
          '56968a4ff1880fee1dff80cb6b8f0044c27ee8108dca4c0d8a3ec8c45ee08926'}
count=0
def req(v,msg):
    global count
    count+=1
    if not bool(v): raise AssertionError(msg)
def eq(a,b,msg): req(s.expand(a-b)==0,msg)
for name,pin in PINS.items():
    req(hashlib.sha256(Path(name).read_bytes()).hexdigest()==pin,'audited source pin')
def generator(p,vs):
    return s.expand(sum(s.diff(p,x,2)-x*s.diff(p,x) for x in vs))
def resolvent(p,k,vs):
    ans=0
    for powers,c in s.Poly(s.expand(p),*vs).terms():
        choices=[]
        for n,z in zip(powers,vs):
            choices.append([(n-2*j,s.factorial(n)/(2**j*s.factorial(j)*s.factorial(n-2*j))
                              *s.hermite_prob(n-2*j,z)) for j in range(n//2+1)])
        for tup in it.product(*choices):
            ans+=c*s.prod(t[1] for t in tup)/(k+sum(t[0] for t in tup))
    ans=s.expand(ans)
    eq(k*ans-generator(ans,vs),p,'Hermite resolvent inverse')
    return ans
def hierarchy(g,vs,maxn):
    # Raw moment hierarchy follows directly from the linear MGF equation.
    moments=[s.Integer(1)]
    kappas=[s.Integer(0)]
    for n in range(1,maxn+1):
        moments.append(resolvent(n*g*moments[n-1],n,vs))
        kappa=s.expand(moments[n]-sum(s.binomial(n-1,j-1)*kappas[j]*moments[n-j]
                                     for j in range(1,n)))
        kappas.append(kappa)
        if n==1:
            eq(kappa-generator(kappa,vs),g,'first connected resolvent')
        else:
            rhs=sum(s.binomial(n,j)*sum(s.diff(kappas[j],x)*s.diff(kappas[n-j],x)
                                      for x in vs) for j in range(1,n))
            eq(n*kappa-generator(kappa,vs),rhs,'connected binomial hierarchy versus raw moments')
    return kappas

x=s.symbols('x')
k_linear=hierarchy(x,[x],6)
eq(k_linear[1],x/2,'linear mean')
eq(k_linear[2],s.Rational(1,4),'linear variance')
for n in range(3,7): eq(k_linear[n],0,'linear higher connected cumulants vanish')
k_quad=hierarchy(x*x-1,[x],6)
eq(k_quad[1],(x*x-1)/3,'quadratic mean')
eq(k_quad[2],s.Rational(2,9)*(x*x+1),'quadratic covariance')
eq(k_quad[3],s.Rational(16,135)*(3*x*x+2),'quadratic third cumulant')
hierarchy(x+(x*x-1)/4,[x],6)

# A nontrivially coupled two-output polynomial gradient tests tensor factors.
x,y,l1,l2=s.symbols('x y l1 l2')
vs=[x,y]
gs=[x*x+x*y/3,y*y+x*x/6]
eq(s.diff(gs[0],y),s.diff(gs[1],x),'polynomial gradient fixture')
projected=l1*gs[0]+l2*gs[1]
ks=hierarchy(projected,vs,3)
vv=[resolvent(g,1,vs) for g in gs]
B=s.Matrix([[s.diff(v,z) for z in vs] for v in vv])
eq(B[0,1],B[1,0],'symmetric Dv for gradient source')
C=s.Matrix(2,2,lambda i,j:2*resolvent((B*B.T)[i,j],2,vs))
ls=[l1,l2]
eq(ks[2],sum(ls[i]*ls[j]*C[i,j] for i,j in it.product(range(2),repeat=2)),
   'tensor k2 equals 2 R2 BBstar')
third_source=sum(ls[i]*ls[j]*ls[k]*sum(B[i,a]*s.diff(C[j,k],vs[a]) for a in range(2))
                 for i,j,k in it.product(range(2),repeat=3))
eq(ks[3],6*resolvent(third_source,3,vs),'tensor k3 equals 6 R3 Sym B dot DC')

# Independent Wick-integral check from the actual conditioned Markov covariance.
r,ss,t=s.symbols('r ss t',positive=True)
cov_rs=r/ss-r*ss  # ordered r <= ss
variance=4*s.integrate(s.integrate(cov_rs**2,(r,0,ss)),(ss,0,1))
eq(variance,s.Rational(2,9),'true-Markov quadratic variance at z=0')
triple=(r/ss-r*ss)*(ss/t-ss*t)*(r/t-r*t)
third=48*s.integrate(s.integrate(s.integrate(triple,(r,0,ss)),(ss,0,t)),(t,0,1))
eq(third,s.Rational(32,135),'true-Markov quadratic third cumulant at z=0')
q=s.symbols('q',positive=True)
for k in range(1,7):
    for n in range(12):
        eq(s.integrate(q**(k-1+n),(q,0,1)),s.Rational(1,k+n),'positive-clock Hermite multiplier')

# Unsymmetrized covariance tensors cannot acquire a Sym on only one side.
# g(x)=x in dimension two gives B=I/2. Take A_field(x)=x_1 e2 tensor e2.
cov_value=s.integrate(s.exp(-ss)*s.exp(-(t-ss)),(ss,0,t))
req(s.simplify(cov_value-t*s.exp(-t))==0,'linear martingale matrix covariance')
req(s.simplify(cov_value-cov_value/3)!=0,'unilateral full symmetrization changes covariance tensor')

# One-derivative positive quadrature: actual dyadic interior nodes and terminal midpoint.
max_moment_ratio=0.
max_derivative_ratio=0.
for panels,gauss_order in ((4,3),(8,5),(12,7)):
    gl,gw=np.polynomial.legendre.leggauss(gauss_order)
    nodes=[];weights=[]
    for j in range(1,panels+1):
        a=2.**(-j)
        nodes.extend((1-1.5*a+.5*a*gl).tolist())
        weights.extend((.5*a*gw).tolist())
    h=2.**(-panels)
    nodes=np.array(nodes+[1-h/2]);weights=np.array(weights+[h])
    req(abs(weights.sum()-1)<1e-14,'positive rule total mass')
    req(np.all(weights>0) and np.all((nodes>0)&(nodes<1)),'positive interior rule')
    exponents=list(range(121))+[1000,10000,100000,1000000]
    for degree in exponents:
        value=float(weights@np.exp(degree*np.log(nodes)))
        ratio=value*(degree+1)
        max_moment_ratio=max(max_moment_ratio,ratio)
        req(ratio<=6+1e-12,'absolute positive moment bound')
    for k in range(1,7):
        ns=np.array(exponents,dtype=float)
        errs=np.array([float(weights@np.exp((n+k-1)*np.log(nodes)))-1/(n+k) for n in ns])
        delta_sample=float(np.max(np.abs(errs)))
        for n,e in zip(ns,errs):
            req(abs(e)<=7/(n+1)+1e-14,'decaying multiplier error bound')
            ratio=n*e*e/(7*delta_sample)
            max_derivative_ratio=max(max_derivative_ratio,float(ratio))
            req(ratio<=1+1e-10,'one-derivative interpolation bound on tested chaoses')
    n=10**9;k=3
    e=float(weights@np.exp((n+k-1)*np.log(nodes)))-1/(n+k)
    req(abs(math.sqrt(n*(n-1))*abs(e)-1)<1e-7,'twice differentiated high-chaos obstruction')

report={'status':'PASS','assertions':count,
        'pins':{source_label(p):v for p,v in PINS.items()},
        'quadratic_conditional_cumulants':{'kappa_1':'(z^2-1)/3','kappa_2':'2(z^2+1)/9',
                                         'kappa_3':'16(3z^2+2)/135'},
        'scope':'Independent finite polynomial checks of raw-moment versus connected recurrence, tensor symmetrization factors, actual conditional-Markov Wick integrals, and Hermite resolvents. Polynomial fixtures test identities, not global Lipschitz bounds.',
        'max_positive_moment_times_degree_plus_one':max_moment_ratio,
        'max_tested_derivative_interpolation_ratio':max_derivative_ratio,
        'display_correction_verified':'The initial unilateral Sym error in section 3 was removed at the audited pin; the exact covariance identity is now unsymmetrized.'}
(HERE/'discounted_ou_cumulant_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
