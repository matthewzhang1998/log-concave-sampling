#!/usr/bin/env python3
"""Independent native rank-three target, covariance-root and current checks."""
from pathlib import Path
import hashlib,itertools as it,json,math
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
str(HERE.parent/'THIRD-CUMULANT-POSITIVE-CLOCK-COMPRESSION.md'):'c3ed7cdf2145707391c0e16bfe675f9f154beedf795ef1ab7b49dc1116737e53',
str(BASE/'rank3-continuation/POSITIVE-RANK3-PACKET-FEEDBACK-LEMMA.md'):'c689f0a947bcf6358f38bc6470ede05f564dcec11ff5e4171df8250e4af78193',
str(BASE/'rank3-continuation/POSITIVE-QUADRATIC-REFERENCE-FOR-SKEW-GAUSSIANIZATION.md'):'a797723ef51ae41a28bdd47c99bb7d5fccee49de6599d07ec8b085578c907c41',
str(BASE/'rank3-continuation/EXACT-FIVE-CLOCK-NATIVE-SOURCE-ADAPTER.md'):'dd6707cae6dabe52f2172202ce0cbcd087c7a069ba1dd47800a9c89300216c26',
}
count=0
def req(v,msg):
    global count
    count+=1
    if not bool(v):raise AssertionError(msg)
def eq(a,b,msg):req(s.expand(a-b)==0,msg)
for name,pin in PINS.items():req(hashlib.sha256(Path(name).read_bytes()).hexdigest()==pin,'source pin')
def moment(p,vs):
    ans=0
    for powers,c in s.Poly(s.expand(p),*vs).terms():
        ans+=c*s.prod(0 if n%2 else s.factorial2(n-1) if n else 1 for n in powers)
    return s.expand(ans)
def resolvent(p,k,vs):
    ans=0
    for powers,c in s.Poly(s.expand(p),*vs).terms():
        choices=[[(n-2*j,s.factorial(n)/(2**j*s.factorial(j)*s.factorial(n-2*j))*s.hermite_prob(n-2*j,z))
                   for j in range(n//2+1)] for n,z in zip(powers,vs)]
        for choices_j in it.product(*choices):
            degree=sum(a for a,b in choices_j)
            if degree+k:ans+=c*s.prod(b for a,b in choices_j)/(degree+k)
    return s.expand(ans)
def riesz(p,vs):
    eq(moment(p,vs),0,'private Riesz input centered')
    return [s.diff(resolvent(p,0,vs),v) for v in vs]

# Gaussian regression of the native reverse polynomial, a=b=c=1, beta^2=3.
S,U,V,A,q,t,z=s.symbols('S U V A q t z')
P=S/3+U/s.sqrt(2)+V/s.sqrt(6)
Y=S/3-U/s.sqrt(2)+V/s.sqrt(6)
xi=S/3-2*V/s.sqrt(6)
eq(P+Y+xi,S,'visible Gaussian readout')
C=A*P*Y
Cbar=A*(S*S/9-s.Rational(1,3))
eq(moment(C,[U,V]),Cbar,'exact conditional Hermite regression')
E=s.expand(C-Cbar)
RE=riesz(E,[U,V])
current=s.expand(sum(a*s.diff(E,v) for a,v in zip(RE,[U,V])))
W=S+q*Cbar+t*q*E+z
for degree in range(1,5):
    lhs=s.diff(moment(W**degree,[U,V]),t)
    rhs=t*q*q*degree*(degree-1)*moment(current*W**(degree-2),[U,V]) if degree>=2 else 0
    eq(lhs,rhs,'first same-endpoint private-Gaussian Riesz current')

# Coarse-root decoupling algebra; linear A(Q) is an identity fixture only.
Q,a0,a1=s.symbols('Q a0 a1')
H=S*S/9-s.Rational(1,3)
Wq=S+q*a0*H+t*q*a1*Q*H+z
for degree in range(1,6):
    lhs=s.diff(moment(Wq**degree,[Q]),t)
    rhs=t*q*q*a1*a1*H*H*degree*(degree-1)*moment(Wq**(degree-2),[Q]) if degree>=2 else 0
    eq(lhs,rhs,'owned-coarse-root same-endpoint current')

G=s.symbols('G')
quad=s.sqrt(3)*G+q*A*(G*G-1)/3
eq(moment(quad,[G]),0,'quadratic reference centered')
eq(s.expand(moment(quad**3,[G])).coeff(q,1),6*A,'native leading cumulant factor six')
w,aa,bb,cc,amp,sigma=s.symbols('w aa bb cc amp sigma',nonzero=True)
rho_root=4*w*amp/(aa*bb*cc*sigma)
eq(aa*bb*cc*rho_root*amp**2*sigma/amp**3,4*w,'five-clock root logarithmic coefficient')
eq(6*4*w,24*w,'five-clock connected third normalization')

# Exact coupled-gradient compressed identity and an explicit wrong-slot detector.
x,y,l1,l2=s.symbols('x y l1 l2');vs=[x,y];ls=[l1,l2]
g=[x*x+x*y/3,y*y+x*x/6]
v=[resolvent(p,1,vs) for p in g]
B=s.Matrix([[s.diff(p,z) for z in vs] for p in v])
T={(j,b,a):s.diff(B[j,b],vs[a]) for j,b,a in it.product(range(2),repeat=3)}
J={(j,k,a):s.expand(sum(T[j,b,a]*B[k,b] for b in range(2)))
   for j,k,a in it.product(range(2),repeat=3)}
JU={ijk:resolvent(p,3,vs) for ijk,p in J.items()}
correct=24*resolvent(sum(ls[i]*ls[j]*ls[k]*B[i,a]*JU[j,k,a]
                        for i,j,k,a in it.product(range(2),repeat=4)),3,vs)
Cmat=s.Matrix(2,2,lambda j,k:2*resolvent((B*B.T)[j,k],2,vs))
recurrence=6*resolvent(sum(ls[i]*ls[j]*ls[k]*B[i,a]*s.diff(Cmat[j,k],vs[a])
                          for i,j,k,a in it.product(range(2),repeat=4)),3,vs)
eq(correct,recurrence,'raw inner derivative slot yields exact compressed target')
Jwrong={idx:s.expand(sum(J[p] for p in it.permutations(idx))/6) for idx in J}
wrong=24*resolvent(sum(ls[i]*ls[j]*ls[k]*B[i,a]*resolvent(Jwrong[j,k,a],3,vs)
                      for i,j,k,a in it.product(range(2),repeat=4)),3,vs)
req(s.expand(correct-wrong)!=0,'full inner three-slot symmetrization changes target')

# Noncommuting Gaussian channel collapse and actual covariance-root ledger.
rng=np.random.default_rng(231003)
max_root_ratio=0.
for n in (2,3,5):
    for _ in range(30):
        fields=[]
        for j in range(3):
            F=rng.normal(size=(n,n));F*=.12/np.linalg.norm(F,2);fields.append(F)
        cross=fields[0]@fields[1]@fields[2]
        wrong_cross=fields[2]@fields[1]@fields[0]
        req(np.linalg.norm(cross-wrong_cross)>1e-7,'chronological matrices do not commute')
        cov=np.eye(n)
        fullcross=np.eye(n)
        for F in reversed(fields):
            cov=F@cov@F.T+np.eye(n)-F@F.T
            fullcross=F@fullcross
        req(np.linalg.norm(cov-np.eye(n))<1e-12,'all channel marginals remain standard')
        req(np.linalg.norm(fullcross-cross)<1e-12,'ordered path cross matrix')
        # Reverse law R=M^T Y+sqrt(I-M^T M)xi retains the same terminal pair.
        evals,evecs=np.linalg.eigh(np.eye(n)-cross.T@cross)
        root=(evecs*np.sqrt(evals))@evecs.T
        req(np.linalg.norm(cross.T@cross+root@root.T-np.eye(n))<1e-12,'reverse terminal marginal covariance')
        change=np.linalg.norm(root-np.eye(n),'fro')
        bound=np.linalg.norm(cross,2)*np.linalg.norm(cross,'fro')
        ratio=change/bound
        max_root_ratio=max(max_root_ratio,ratio)
        req(ratio<=1+1e-7,'strong covariance-root feedback product bound')

max_dilation_ratio=0.
for panels,order in ((4,3),(8,5),(12,7)):
    gl,gw=np.polynomial.legendre.leggauss(order)
    rr=[];ww=[]
    for j in range(1,panels+1):
        a=2.**(-j);rr.extend((1-1.5*a+.5*a*gl).tolist());ww.extend((.5*a*gw).tolist())
    h=2.**(-panels);rr=np.array(rr+[1-h/2]);ww=np.array(ww+[h])
    rp=np.sqrt((1+rr*rr)/2)
    req(np.all((1-rp)>=(1-rr)/4),'half-shield endpoint gaps comparable')
    for degree in list(range(101))+[1000,10000,100000]:
        value=float((ww*rr)@np.exp(degree*np.log(rp)))
        max_dilation_ratio=max(max_dilation_ratio,(degree+2)*value)
        req((degree+2)*value<=24,'dilated positive-clock multiplier envelope')

report={'status':'PASS','assertions':count,'pins':{source_label(p):v for p,v in PINS.items()},
 'max_root_feedback_to_operator_times_HS_ratio':max_root_ratio,
 'max_dilated_moment_times_degree_plus_two':max_dilation_ratio,
 'wrong_inner_symmetrization_difference':str(s.factor(correct-wrong)),
 'scope':'Independent exact native regression/current identities, coupled-gradient five-clock ordering, amplitude normalization, and noncommuting Gaussian-channel/root diagnostics. Actual imported pair compilers are not run.'}
(HERE/'native_rank3_feedback_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
