"""Independent finite diagnostics for the bounded tight-bridge audit.
No native compiler, higher-current producer, or all-order theorem is executed.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib
import json
import math
import platform
import numpy as np
from numpy.polynomial.hermite import hermgauss


# Publication-only source locator and original/public pin mapping.
import sys
_ARCHIVE=next(p for p in Path(__file__).resolve().parents if (p/'publication_sources.py').is_file())
sys.path.insert(0,str(_ARCHIVE))
from publication_sources import source_path, source_label, verify_source_pin, pin_report

HERE = Path(__file__).resolve().parent
BRIDGE = HERE.parent
RANK3 = BRIDGE.parent.parent
GATE = RANK3.parent
NATIVE = Path('external:LOW30')
PINNED = {
    BRIDGE/'TIGHT-BRIDGE-LEDGER-AND-REROOT-GATE.md': 'd36ca9949386d3203ecde7a60d83ced13ac094b89cf1262bfb7de9f6acc8cc6a',
    HERE/'REVIEWED-TIGHT-BRIDGE-SOURCE.md': 'd36ca9949386d3203ecde7a60d83ced13ac094b89cf1262bfb7de9f6acc8cc6a',
    BRIDGE/'check_tight_bridge.py': '743db9d0b096f7eb9618e354488202d73339af75f17925deda19e15de7340910',
    RANK3/'grouped-kernel-response/PAIRED-ROTATION-COVARIANCE-ALL-CUT-LEMMA.md': 'b8383590bdea35e5783c4fe86fdd68551d3d817a1ce81ee1c9004de8964ca668',
    RANK3/'grouped-kernel-response/INTERMEDIATE-SCALE-B-MARK-COROLLARY.md': '070cac5569e9c6a1386bc26c5f98146bc8e25b5e00850cd146950c769e2235c7',
    RANK3/'EXACT-FIVE-CLOCK-NATIVE-SOURCE-ADAPTER.md': 'dd6707cae6dabe52f2172202ce0cbcd087c7a069ba1dd47800a9c89300216c26',
    BRIDGE.parent/'CONDITIONAL-QUARTIC-PACKET-AND-CYCLE-RETURN.md': 'ccaae83e8308c578921cc3a0f85239758bdaafd44c209a9a75cd39f9fd02e02b',
    NATIVE: '7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8',
}
counts = Counter()

def check(group, condition):
    if not bool(condition):
        raise AssertionError(f'Failed {group} assertion #{counts[group]+1}')
    counts[group] += 1

for path, expected in PINNED.items():
    status = verify_source_pin(path, expected)
    if status is not None: check('source_pins', status)

# Six-force/six-mark source topology; not an implementation of its C2 filters.
edges = [(0,1),(1,2),(2,3),(1,4),(2,5)]
degree = [sum(v in e for e in edges) for v in range(6)]
check('native_topology', degree == [1,3,3,1,1,1])
check('native_topology', [d+1-2 for d in degree] == [0,2,2,0,0,0])
check('native_topology', sum(degree)-6 == 4)  # sum(j_v-1)=M-2 for tree.
check('native_topology', 2*3+4-2*4+1 == 3)  # two 4-vertex spines, four edges.

# Enumerate all nonroot inverse-delta allocations at three exact fractional p's.
for p in (F(1,4),F(1,2),F(3,4)):
    for flags in product((0,1), repeat=5):
        root = 1+2*p+p*sum(flags)
        others = [1-p*f for f in flags]
        spine = root+sum(others[:3])
        check('amplitude_ledger', root+sum(others) == 6+2*p)
        check('amplitude_ledger', spine == 4+2*p+p*sum(flags[3:]))
        check('amplitude_ledger', min([root]+others)>0)
        if p == F(1,2):
            check('amplitude_ledger', 2*spine == 10+sum(flags[3:]))

# Exact Gaussian-polynomial cumulants of four separately linear probe slots.
# The monomial exponent records the shifts in colors P,Q,S,T, in that order.
def gm(k):
    return 0 if k%2 else math.prod(range(1,k,2))

def mul(a,b):
    out=Counter()
    for ka,va in a.items():
        for kb,vb in b.items():
            out[tuple(x+y for x,y in zip(ka,kb))] += va*vb
    return dict(out)

moments={0:{(0,0,0,0):1}}
cumulants={}
for n in range(1,7):
    univ={k:math.comb(n,k)*gm(n-k) for k in range(n+1) if (n-k)%2==0}
    mu={ks:math.prod(univ[k] for k in ks) for ks in product(univ,repeat=4)}
    moments[n]=mu
    kap=dict(mu)
    for k in range(1,n):
        for term,value in mul(cumulants[k],moments[n-k]).items():
            kap[term]=kap.get(term,0)-math.comb(n-1,k-1)*value
    kap={k:v for k,v in kap.items() if v}
    cumulants[n]=kap
    for shifts,value in kap.items():
        check('wick_exact_cumulants', all(k%2==n%2 for k in shifts))
        check('wick_exact_cumulants', value>0)
        for sides in product((0,1),repeat=2):
            spine=F(5)+F(sum(sides),2)
            grade=n*spine+sum(k*(1-F(f,2)) for k,f in zip(shifts[2:],sides))
            if n>=2:
                check('wick_exact_cumulants',grade>=2*spine)
            check('wick_exact_cumulants', 4*n+sum(shifts[2:])>=4*n)
check('wick_exact_cumulants',cumulants[1]=={(1,1,1,1):1})
check('wick_exact_cumulants',cumulants[2][(0,0,0,0)]==1)
check('wick_exact_cumulants',(0,0,0,0) not in cumulants[3])

# Dense noncommuting Gaussian loading calculation retains the original R.
rng=np.random.default_rng(20261005)
def sqrt_psd(a):
    eig,vec=np.linalg.eigh(a)
    return (vec*np.sqrt(eig))@vec.T
for d in (1,2,3,7):
    eye=np.eye(d)
    mp=rng.normal(size=(d,d)); mp*=.42/np.linalg.norm(mp,2)
    mm=rng.normal(size=(d,d)); mm*=.37/np.linalg.norm(mm,2)
    sp=sqrt_psd(eye-mp@mp.T); sm=sqrt_psd(eye-mm@mm.T)
    # Four independent input blocks: original R, Z0, Z+, Z-.
    lp=np.hstack((mp/math.sqrt(2),mp/math.sqrt(2),sp,np.zeros((d,d))))
    lm=np.hstack((mm/math.sqrt(2),-mm/math.sqrt(2),np.zeros((d,d)),sm))
    lr=np.hstack((eye,np.zeros((d,3*d))))
    ly=(lp-lm)/math.sqrt(2)
    for got,want in ((lp@lp.T,eye),(lm@lm.T,eye),(lp@lm.T,np.zeros((d,d))),
                     (ly@ly.T,eye),(ly@lr.T,(mp-mm)/2),(lr@lr.T,eye)):
        check('retained_R_channel',np.max(np.abs(got-want))<2e-13)
    # At equal crosses, law is independent but its noise loading still changes.
    equal_y=np.hstack((np.zeros((d,d)),mp,sp/math.sqrt(2),-sp/math.sqrt(2)))
    check('retained_R_channel',np.max(np.abs(equal_y@equal_y.T-eye))<2e-13)
    check('retained_R_channel',np.linalg.norm(equal_y[:,d:2*d])>0)

# Exact bounded-source sine witness and Hermite conditional expansion.
x,w=hermgauss(96); x*=math.sqrt(2); w/=math.sqrt(math.pi)
def hermite(n,r):
    a,b=1.,r
    if n==0:return a
    for k in range(1,n): a,b=b,r*b-k*a
    return b
for delta,eps in product((.04,.2,.6),(.08,.23)):
    m=eps*np.sin(delta*x)
    v=eps**2*(-math.expm1(-2*delta*delta))/2
    check('bounded_reroot',abs(w@m)<1e-14)
    check('bounded_reroot',abs(w@(m*m)-v)<1e-14)
    for r in (0.,.5,1.,1.5,-1.):
        m2=float(w@((r*m)**2+1-m*m))
        check('bounded_reroot',abs(m2-(1+(r*r-1)*v))<1e-13)
        # Coefficient of t^n in conditional MGF from direct Gaussian moments.
        for n in range(0,11):
            direct=0.
            for j in range(0,n+1,2):
                direct += float(w@(math.comb(n,j)*(r*m)**(n-j)*(1-m*m)**(j/2)*gm(j)))
            series=sum(math.factorial(n)*hermite(k,r)*float(w@(m**k)) /
                       (math.factorial(k)*2**((n-k)//2)*math.factorial((n-k)//2))
                       for k in range(n+1) if (n-k)%2==0)
            check('bounded_reroot',abs(direct-series)<2e-9)
        if r*r==1:
            fourth=float(w@((r*m)**4+6*(r*m)**2*(1-m*m)+3*(1-m*m)**2))
            check('bounded_reroot',abs(fourth-(3-2*float(w@(m**4))))<1e-12)
            check('bounded_reroot',fourth<3)
        else: check('bounded_reroot',abs(m2-1)>1e-7)
    a,c=.3,.4
    k4=3*float(w@((a*a+c*c+2*a*c*m)**2))-3*(a*a+c*c)**2
    check('bounded_reroot',abs(k4/24-(a*c)**2*v/2)<1e-14)

# Literal original gradient is anchored and globally strictly increasing.
A=.02; kappa=.45; eta=.17; sigma=.41; t=1/math.sqrt(3)
g=lambda q:A*(kappa*q+eta*(1-np.cos(q)))
dg=lambda q:A*(kappa+eta*np.sin(q))
check('original_gradient_witness',g(0)==0)
check('original_gradient_witness',0<kappa-eta<kappa+eta<1)
# Direct C1 selected-Jacobian mean, before any finite first-chaos filter.
# Averaging x,z gives variance 2 in the argument; mean is cos(q) times sine(P).
for q,p in product((-.2,math.pi/2,1.4),(.3,.8,-.6)):
    direct=float(w@((dg(q+sigma*t*(math.sqrt(2)*x+p))-dg(q+sigma*t*(math.sqrt(2)*x-p)))/(2*A)))
    target=eta*math.exp(-sigma*sigma*t*t)*math.cos(q)*math.sin(sigma*t*p)
    check('original_gradient_witness',abs(direct-target)<2e-15)
    # External OU preparation is still a bounded odd sine times cos(q).
    u=.7; attenuation=math.exp(-.5*(sigma*t)**2*(1-u*u))
    prepared=float(w@np.sin(sigma*t*(u*p+math.sqrt(1-u*u)*x)))
    check('original_gradient_witness',abs(prepared-attenuation*math.sin(sigma*t*u*p))<2e-15)
for z in np.linspace(-2,2,15):
    check('original_gradient_witness',abs(math.cos(math.pi/2+.21*z)-math.cos(math.pi/2-.21*z)+2*math.sin(.21*z))<1e-14)

# All-proper-cut and one-HS grouped-network diagnostics, no self-contractions.
def cut_norms(a):
    n=a.ndim
    out=[]
    for mask in range(1,2**n-1):
        left=[i for i in range(n) if mask>>i&1]
        right=[i for i in range(n) if not(mask>>i&1)]
        matrix=a.transpose(left+right).reshape(math.prod(a.shape[i] for i in left),-1)
        out.append(np.linalg.norm(matrix,2))
    return out
for trial in range(4):
    blocks=[]
    for j in range(3):
        b=rng.normal(size=(2,)*6)
        b/=max(cut_norms(b)); blocks.append(b)
    two=np.einsum('abijkl,cdijkl->abcd',blocks[0],blocks[1])
    three=np.einsum('abijkl,cdijmn,efklmn->abcdef',*blocks)
    for out in (two,three):
        for cut in cut_norms(out):check('grouped_tensor_cuts',cut<=1+1e-12)
        check('grouped_tensor_cuts',np.linalg.norm(out)<=math.sqrt(2)+1e-12)
    check('grouped_tensor_cuts',np.linalg.norm(three)<=min(np.linalg.norm(b) for b in blocks)+1e-12)

# Tight endpoint mass diagnoses the chosen upper-bound ledger only.
for h in range(2,9):
    for k in (2,6,12):
        sig=F(1,2**k); delta_mass=sig**2
        beta=delta_mass**2/sig**4
        check('width_ledger',beta==1)
        check('width_ledger',beta**h*sig**(-2*(h-1))==sig**(-2*(h-1)))
        check('width_ledger',delta_mass**4/sig**10==sig**-2)
    for bp in (F(1),F(7,2),F(7),F(10)):
        grade=7*h-2*bp*(h-1)
        check('width_ledger',grade==(7-2*bp)*h+2*bp)
        if bp>=F(7,2):check('width_ledger',grade<=(7-2*bp)*(h-1)+2*bp)
check('width_ledger',14-2*F(7)==0)
# Demonstrate expressly that a divergent derivative bound is not a variance lower bound.
for k in (4,16,64):
    exact_variance=(1-math.exp(-2*k*k))/2
    check('width_ledger',exact_variance<=.5 and k*k>exact_variance)

# Joint score: invert the correlated two-bank covariance without using a marginal.
for delta in (.05,.2,.5):
    tau=1-2*delta*delta
    cov=np.array([[1,tau],[tau,1]])
    B=.013
    displacement=np.array([B,-B])
    energy=float(displacement@np.linalg.solve(cov,displacement))
    check('joint_score',abs(energy-(B/delta)**2)<1e-12)
    check('joint_score',abs(displacement[0]**2-B*B)<1e-15)
for N in (7,10,12):
    for j in range(1,9):check('joint_score',F(N)+j*(1-F(1,2))==F(N)+F(j,2))

report={
    'status':'PASS',
    'scope':'Independent finite diagnostics and pin verification for a bounded candidate/failure audit; not a native compiler or closure proof.',
    'checks':sum(counts.values()), 'groups':dict(counts),
    'source_sha256':PINNED[BRIDGE/'TIGHT-BRIDGE-LEDGER-AND-REROOT-GATE.md'],
    'source_pins':{source_label(p):sha for p,sha in PINNED.items()},
    'publication_pin_verification':pin_report(),
    'versions':{'python':platform.python_version(),'numpy':np.__version__},
    'open':['tight whole-coarse positive current compiler','all-proper-cut higher joint physical-R score return','source-qualified all-order closure'],
}
(HERE/'tight_bridge_independent_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='source_pins'},indent=2))
