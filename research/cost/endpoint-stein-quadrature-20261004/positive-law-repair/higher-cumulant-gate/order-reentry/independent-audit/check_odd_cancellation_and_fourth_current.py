#!/usr/bin/env python3
"""Independent exact cumulant, Markov, Riesz and buffer-current diagnostics."""
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
BASE=HERE.parent.parent
PINS={
    str(HERE.parent/'ODD-CUMULANT-CANCELLATION-AND-TARGET-IDENTITY.md'):
        'ce0d7e8e262a8168aac8818a7830fdb720443fc75b8ad4b1c64343287ef0d9c9',
    str(BASE/'rank3-continuation/FOURTH-ORDER-GAUSSIANIZATION-WHEN-THIRD-CUMULANT-VANISHES.md'):
        '06ab222f8f739490c0142582a565ebaf3c9605a7ef0963385502bab32105c880',
}
count=0
def req(value,message):
    global count
    count+=1
    if not bool(value): raise AssertionError(message)
def eq(a,b,message): req(s.expand(a-b)==0,message)
for name,pin in PINS.items():
    req(hashlib.sha256(Path(name).read_bytes()).hexdigest()==pin,'source pin')

weights=[s.Rational(1,6)]*8+[s.Rational(-1,3)]
for k,val in enumerate([1,s.Rational(1,3),0,s.Rational(1,54),s.Rational(-1,324)],1):
    eq(sum(w**k for w in weights),val,'nine-copy cumulant coefficient')
vandermonde=[]
for m in range(1,8):
    bs=list(range(1,m+2))
    matrix=s.Matrix([[b**(2*j+3) for b in bs] for j in range(m)])
    vec=matrix.nullspace()[0]
    lcm=s.ilcm(*[v.q for v in vec])
    multiplicities=[int(lcm*v) for v in vec]
    mean=sum(d*b for d,b in zip(multiplicities,bs))
    req(mean!=0,'finite prescribed odd cancellations can preserve mean after scaling')
    for j in range(m):
        eq(sum(d*s.Rational(b,mean)**(2*j+3) for d,b in zip(multiplicities,bs)),0,
           'exact finite odd cancellation via signed multiplicities')
    eq(sum(d*s.Rational(b,mean) for d,b in zip(multiplicities,bs)),1,'normalized coefficient mean')
    # The d values count signed copies, not arbitrary probability weights.
    det=s.Matrix([[b**3*(b*b)**j for b in bs] for j in range(m+1)]).det()
    req(det!=0,'all-odd cancellation Vandermonde invertibility')
    vandermonde.append({'odd_ranks_through':2*m+1,'distinct_magnitudes':m+1,
                         'complete_copies':sum(abs(d) for d in multiplicities)})

r,x=s.symbols('r x',positive=True)
covL=s.integrate(x/r-x*r,(x,0,r))+s.integrate(r/x-r*x,(x,r,1))
req(s.simplify(covL+r*s.log(r))==0,'true conditional Markov covariance')
eq(s.integrate(-r*s.log(r),(r,0,1)),s.Rational(1,4),'Var L')
eq(s.integrate(r*r*s.log(r)**2,(r,0,1)),s.Rational(2,27),'exact squared covariance integral')
eq(s.Rational(1,8)*s.Rational(2,27),s.Rational(1,108),'mixed cubic lower bound')
eq(8*s.Rational(15,27),s.Rational(40,9),'centered J third absolute moment bound')
eps=s.Rational(1,20)
num=3*eps/s.Integer(108)-eps**3*s.Rational(40,9)
eq(num,s.Rational(1,1200),'cubic numerator floor')
eq(num/(1+eps)**3,s.Rational(20,27783),'explicit positive third cumulant floor')
req(5/(3*math.sqrt(2*math.pi))>.5,'elementary normal central probability lower bound')
req(math.cosh(1)<2,'sech squared floor')

# Exact polynomial identities test the Riesz algebra, not the global-Lipschitz theorem.
v=s.symbols('v0:2')
out=s.symbols('y0:2')
def moment(p):
    total=0
    for powers,c in s.Poly(s.expand(p),*v).terms():
        prod=c
        for k in powers:
            prod*=0 if k%2 else (s.factorial2(k-1) if k else 1)
        total+=prod
    return s.expand(total)
def riesz(p):
    eq(moment(p),0,'Riesz field centered')
    inv=0
    for powers,c in s.Poly(s.expand(p),*v).terms():
        choices=[]
        for n,z in zip(powers,v):
            choices.append([(n-2*k,s.factorial(n)/(2**k*s.factorial(k)*s.factorial(n-2*k))
                             *s.hermite_prob(n-2*k,z)) for k in range(n//2+1)])
        for tup in it.product(*choices):
            degree=sum(t[0] for t in tup)
            if degree: inv+=c*s.prod(t[1] for t in tup)/degree
    return [s.expand(s.diff(inv,z)) for z in v]
def sym(tensor,indices):
    perms=list(it.permutations(indices))
    return s.expand(sum(tensor[p] for p in perms)/len(perms))

nonzero_third=False
for f in (
    [v[0]+(v[0]**2-1)/3, v[1]+v[0]*v[1]/4+(v[0]**2-1)/5],
    [v[0]+(v[0]**3-3*v[0])/7, v[1]+v[0]**2*v[1]/5],
):
    f=[s.expand(p-moment(p)) for p in f]
    J=[[s.diff(p,z) for z in v] for p in f]
    rf=[riesz(p) for p in f]
    Sigma={(i,j):moment(f[i]*f[j]) for i,j in it.product(range(2),repeat=2)}
    B={(i,j):s.expand(sum(rf[i][a]*J[j][a] for a in range(2))-Sigma[i,j])
       for i,j in it.product(range(2),repeat=2)}
    rb={ij:riesz(p) for ij,p in B.items()}
    rawT={(i,j,k):s.expand(sum(rb[i,j][a]*J[k][a] for a in range(2)))
          for i,j,k in it.product(range(2),repeat=3)}
    T={ijk:sym(rawT,ijk) for ijk in rawT}
    M={ijk:moment(p) for ijk,p in T.items()}
    for i,j,k in T:
        eq(2*M[i,j,k],moment(f[i]*f[j]*f[k]),'second Stein tensor mean is kappa3/2')
        nonzero_third|=bool(M[i,j,k]!=0)
    rt={ijk:riesz(p-M[ijk]) for ijk,p in T.items()}
    rawQ={(i,j,k,l):s.expand(sum(rt[i,j,k][a]*J[l][a] for a in range(2)))
          for i,j,k,l in it.product(range(2),repeat=4)}
    Q={ijkl:sym(rawQ,ijkl) for ijkl in rawQ}
    for phi in (out[0]**2,out[0]**2*out[1],out[0]*out[1]**3,
                out[0]**3*out[1]**2,out[0]**2*out[1]**4):
        def deriv(indices):
            return s.diff(phi,*[out[i] for i in indices]).subs(dict(zip(out,f)),simultaneous=True)
        left=sum(moment(B[ij]*deriv(ij)) for ij in B)
        mid=sum(moment(T[ijk]*deriv(ijk)) for ijk in T)
        right=sum(M[ijk]*moment(deriv(ijk)) for ijk in T)+sum(moment(Q[ijkl]*deriv(ijkl)) for ijkl in Q)
        eq(left,mid,'same-endpoint second-to-third current')
        eq(mid,right,'same-endpoint centered third-to-fourth current')
req(nonzero_third,'general skew fixture tests nonzero constant third term')

# Exact third-Hermite vector isometry with a non-diagonal inverse buffer root.
H=s.Matrix([[1,s.Rational(1,4)],[s.Rational(1,4),s.Rational(1,2)]])
Q={idx:s.Rational(1+sum(idx),7) for idx in it.product(range(2),repeat=4)}
C={idx:sum(Q[(idx[0],j,k,l)]*H[j,idx[1]]*H[k,idx[2]]*H[l,idx[3]]
           for j,k,l in it.product(range(2),repeat=3)) for idx in Q}
def h3(a,b,c):
    return v[a]*v[b]*v[c]-(int(a==b)*v[c]+int(a==c)*v[b]+int(b==c)*v[a])
velocity=[sum(C[i,a,b,c]*h3(a,b,c) for a,b,c in it.product(range(2),repeat=3)) for i in range(2)]
eq(sum(moment(p*p) for p in velocity),6*sum(p*p for p in C.values()),'factor sqrt(3!) in buffer current')
normH=np.linalg.norm(np.array(H).astype(float),2)
req(float(sum(p*p for p in C.values()))<=normH**6*float(sum(p*p for p in Q.values()))*(1+1e-12),
    'three inverse-root factors operator bound')
eq(s.integrate(r**3,(r,0,1)),s.Rational(1,4),'fourth-order dynamic path length factor')

report={'status':'PASS','assertions':count,'pins':{source_label(p):v for p,v in PINS.items()},'finite_odd_patterns':vandermonde,
        'markov_third_cumulant_lower_bound':'20/27783',
        'buffer_constant':'sqrt(6)/(4 sigma^3)',
        'scope':'Exact coefficient and Markov identities; independent finite-polynomial Riesz/Stein identity tests and Hermite isometry. Polynomial fixtures test algebra only, not the globally Lipschitz theorem.'}
(HERE/'odd_cancellation_and_fourth_current_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
