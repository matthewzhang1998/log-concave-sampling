#!/usr/bin/env python3
"""Exact Price/Riesz/Hermite factors, proper cuts and four-path diagnostics."""
from pathlib import Path
from fractions import Fraction
from itertools import combinations, permutations
import hashlib, json, math
import numpy as np
import sympy as S
from scipy.special import roots_hermitenorm
HERE=Path(__file__).resolve().parent
checks=0

def check(ok,msg):
    global checks
    checks+=1
    if not ok: raise AssertionError(msg)

# Exact Gaussian polynomial calculation. C=a*G^2 is used ONLY to test identities;
# it does not have the theorem's bounded first/operator interface.
g,p,t,a,b,s=S.symbols('g p t a b s',real=True)
K=2*a*a
Q=-K/(8*s**3)
h=Q*(p**3-3*p)
h1=S.diff(h,p);h2=S.diff(h,p,2);h3=S.diff(h,p,3)
J=s+t*t*h1
mean=s*p+t*t*h
variance=b+a*(1-t)+t*a*g*g
T2=4*a**3*g*g

def normal_moment(k):
    if k%2:return S.Integer(0)
    return S.Integer(1) if k==0 else S.factorial2(k-1)

def egp(poly):
    terms=S.Poly(S.expand(poly),g,p)
    out=S.Integer(0)
    for (i,j),coef in terms.terms():
        out+=coef*normal_moment(i)*normal_moment(j)
    return S.expand(out)

conditional={}
for n in range(9):
    conditional[n]=S.expand(sum(S.binomial(n,2*j)*normal_moment(2*j)*variance**j*mean**(n-2*j) for j in range(n//2+1)))

def dtest(n,r,coef=1):
    if n<r:return S.Integer(0)
    return S.factorial(n)/S.factorial(n-r)*egp(coef*conditional[n-r])

scalar=[]
for n in (2,4,6,8):
    direct=S.diff(egp(conditional[n]),t)
    fourth=t*K*dtest(n,4)/4
    sixth=t*t*dtest(n,6,T2)/8
    packet=2*t*dtest(n,1,h)
    check(S.expand(direct-fourth-sixth-packet)==0,f'exact Price plus Riesz order {n}')
    ibp=2*t*Q*(dtest(n,4,J**3)+3*t*t*dtest(n,3,J*h2)+t*t*dtest(n,2,h3))
    check(S.expand(packet-ibp)==0,f'exact triple Hermite IBP order {n}')
    feedback=2*t*Q*(dtest(n,4,J**3-s**3)+3*t*t*dtest(n,3,J*h2)+t*t*dtest(n,2,h3))
    check(S.expand(direct-sixth-feedback)==0,f'exact same-law cancellation order {n}')
    scalar.append({'test_degree':n,'identities':'passed'})
check(S.expand(t*K/4+2*t*Q*s**3)==0,'exact negative one-eighth factor')
# Cubic correction cannot erase its positive second-order covariance bill.
second=egp(conditional[2])
check(S.expand(second-(b+a+s**2+6*t**4*Q**2))==0,'positive cubic covariance feedback')
# Centered C has R(C-EC)=aG, T1=2a^2G^2, R(T1-ET1)=2a^2G.
check(S.expand(g*(a*g)-S.diff(a*g,g)-a*(g*g-1))==0,'first Gaussian divergence')
check(S.expand(g*(2*a*a*g)-S.diff(2*a*a*g,g)-2*a*a*(g*g-1))==0,'second Gaussian divergence')

# Every allocation is exact, including the internal packet keep.
for theta in (Fraction(1,2),Fraction(3,5),Fraction(4,5)):
    ug,up,uk=Fraction(1,2),Fraction(1,4),Fraction(1,4)
    visible=theta*up
    keep=uk+(1-theta)*up
    check(ug+visible+keep==1,'no internal keep double count')
    check(visible>0 and keep>0,'positive visible and untouched shares')

# Proper cuts of a genuine bounded smooth positive Gaussian covariance field.
# C(G)=alpha^2[I+sum_j tanh(G_j) H_j], sum ||H_j||op<=1/4.
rng=np.random.default_rng(204205)
gh,gw=roots_hermitenorm(128)
v=float((gw/math.sqrt(2*math.pi))@(np.tanh(gh)**2))
cut_results=[]
for d in (2,3,4,6):
    for rep in range(4):
        m=3
        mats=[]
        for _ in range(m):
            z=rng.normal(size=(d,d));hmat=(z+z.T)/2
            hmat/=4*m*np.linalg.norm(hmat,2)
            mats.append(hmat)
        alpha=.2
        radius=alpha**2*sum(np.linalg.norm(x,2) for x in mats)
        # Conservative Hilbert operator first bound, sufficient for all checks.
        lip=alpha**2*math.sqrt(sum(np.linalg.norm(x,'fro')**2 for x in mats))
        energy=alpha**2*math.sqrt(v*sum(np.linalg.norm(x,'fro')**2 for x in mats))
        k0=v*alpha**4*sum(np.einsum('ij,kl->ijkl',hmat,hmat) for hmat in mats)
        k=sum(np.transpose(k0,perm) for perm in permutations(range(4)))/24
        hs=np.linalg.norm(k)
        check(hs<=lip*energy*(1+1e-12),'one-energy K HS')
        maxcuts={1:0.,2:0.,3:0.}
        for size in (1,2,3):
            for left in combinations(range(4),size):
                right=tuple(j for j in range(4) if j not in left)
                flat=np.transpose(k,left+right).reshape(d**size,d**(4-size))
                op=np.linalg.norm(flat,2)
                bound=lip*radius if size!=2 else max(lip*lip,radius*radius)
                check(op<=bound*(1+2e-12),'source-qualified physical proper cut')
                maxcuts[size]=max(maxcuts[size],float(op))
        q1=k.reshape(d,d**3)
        q2=k.reshape(d*d,d*d)
        for label,flat in [('three-index',q1),('two-index',q2)]:
            gram=flat@flat.T
            check(np.linalg.norm(gram,'fro')<=np.linalg.norm(flat,2)*hs*(1+2e-12),label+' feedback Gram ideal')
        for _ in range(3):
            c=alpha**2*(np.eye(d)+sum(math.tanh(float(rng.normal()))*z for z in mats))
            check(np.linalg.eigvalsh(c)[0]>=.75*alpha**2*(1-1e-12),'positive smooth covariance sample')
        cut_results.append({'dimension':d,'replicate':rep,'HS':float(hs),'energy_bound':lip*energy,'max_propercuts':maxcuts})

# Four original derivative placements have the exact path contraction.
for d in (2,3,4):
    Bs=[];Ts=[]
    for _ in range(4):
        z=rng.normal(size=(d,d)); Bs.append((z+z.T)/2)
        z=rng.normal(size=(d,d,d));Ts.append(sum(np.transpose(z,perm) for perm in permutations(range(3)))/6)
    left=[np.einsum('iem,je->ijm',Ts[0],Bs[1]),np.einsum('ie,jem->ijm',Bs[0],Ts[1])]
    right=[np.einsum('kfm,lf->klm',Ts[2],Bs[3]),np.einsum('kf,lfm->klm',Bs[2],Ts[3])]
    direct=np.einsum('ijm,klm->ijkl',left[0]+left[1],right[0]+right[1])
    expanded=np.zeros_like(direct)
    for hl in (0,1):
        for hr in (0,1):
            edges=[(0,1),(2,3),(hl,2+hr)]
            degrees=[sum(i in edge for edge in edges) for i in range(4)]
            check(sorted(degrees)==[1,1,2,2],'literal four-path topology')
            check(degrees[hl]==degrees[2+hr]==2,'hit vertices are path centers')
            expanded+=np.einsum('ijm,klm->ijkl',left[hl],right[hr])
    check(np.allclose(expanded,direct,rtol=1e-12,atol=1e-12),'four derivative placements equal covariance gradient product')

# Exact scaling exponents for the substantive correction terms.
# sqrt(u)*alpha^r = A^r u^((1-r)/2).
for r,expected in [(5,Fraction(-2)),(6,Fraction(-5,2)),(8,Fraction(-7,2))]:
    check(Fraction(1-r,2)==expected,'physical correction scaling')
for A in np.logspace(-8,-2,7):
    for beta in (.4,.8,1.,1.5):
        u=A**beta
        alpha=A/math.sqrt(u)
        main=A**5/u**2
        for term in (A**5/u**1.5,A**6/u**2.5,A**8/u**3.5):
            check(term<=main*(1+1e-12),'lower-order correction absorption')

# Finite positive dyadic shield sums used by the unheated path specialization.
# Positive within-panel weights need no lower individual-weight bound.
endpoint_diagnostics=[]
for panels in (4,8,16,32,64,128,256):
    sums=np.zeros(4)
    for j in range(1,panels+1):
        delta=math.ldexp(1.,-j)
        fractions=rng.uniform(.1,1.,size=7);fractions/=fractions.sum()
        weights=delta*fractions
        shields=np.sqrt(delta*rng.uniform(1.,2.,size=7))
        sums+=np.array([sum(weights/shields),sum(weights/shields**2),
                        sum(weights**2/shields**2),sum(weights**2/shields**3)])
    check(sums[0]<=1/(math.sqrt(2)-1),'unheated first-inverse sum')
    check(sums[1]<=panels,'unheated second-inverse logarithmic sum')
    check(sums[2]<=1,'unheated squared second-inverse sum')
    check(sums[3]<=1/(math.sqrt(2)-1),'unheated squared third-inverse sum')
    endpoint_diagnostics.append({'panels':panels,'sums':sums.tolist()})
for order in range(5,21):
    beta=Fraction(2*(order-4),order-3)
    check(Fraction(order-4)-beta*Fraction(order-3,2)==0,'general cutoff order equality')
    check(beta<2,'fixed cutoff remains below two')

sources=[HERE/'INDEPENDENT-POSITIVE-GRAM-CORRECTION-AUDIT.md',Path(__file__),
Path('/workspace/shared/retuned-gram-native-20261005/POSITIVE-CORRECTED-GRAM-MIXTURE.md'),
Path('/workspace/shared/retuned-gram-native-20261005/GRAM-NOISE-FOUR-PATH-TARGET-AND-SHARPNESS.md'),
Path('/workspace/shared/fourth-cumulant-return-20261005/POSITIVE-SMOOTHED-FOURTH-CUMULANT-PACKET.md'),
Path('/workspace/shared/fourth-cumulant-return-20261005/independent-audit/AUDIT.md'),
Path('/workspace/shared/fourth-cumulant-return-20261005/independent-audit/MANIFEST.json')]
result={'status':'passed','assertions':checks,
'scope':'Exact polynomial Price/Riesz/Hermite identities, source-qualified tensor cut diagnostics, exact variance allocations, original derivative placements and scaling; not native compiler execution.',
'scalar_identity_checks':scalar,'proper_cut_diagnostics':cut_results,
'endpoint_diagnostics':endpoint_diagnostics,
'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
(HERE/'positive_gram_correction_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['status','assertions','scope']},indent=2))
