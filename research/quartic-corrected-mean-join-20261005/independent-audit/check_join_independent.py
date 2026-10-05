#!/usr/bin/env python3
"""Independent scalar/algebraic diagnostics; never executes native compilers."""
from fractions import Fraction as F
from pathlib import Path
import hashlib, json, math
import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
checks = []

def check(label, condition, detail=None):
    if not bool(condition):
        raise AssertionError((label, detail))
    checks.append({"label": label, "passed": True, "detail": detail})

alpha, beta, gamma, grade = F(14,15), F(19,10), F(19,10), F(39,10)
check("bridge optimal h exponent", gamma == (1+3*alpha)/2)
check("bridge balanced grade", 2+gamma == 3+3*alpha-gamma == grade)
check("positive cutoff ordering", 0 < alpha < beta < 2)
check("native mean order 25 domination", F(25-4)-beta*F(25-3,2)==F(1,10))
check("native radius margin", 1-beta/2==F(1,20))
check("retuned reserve margin", F(3,2)-3*beta/4==F(3,40))
check("mixed K radius margin", 1-beta/3==F(11,30))
check("mixed K covariance gap", 3-beta>0)

nonbuffer = {
    "smoothing commutator":2+gamma,
    "smoothing drift":1+2*gamma,
    "prefix omission and coherent shift":3+alpha,
    "bridge mean-only decoupling":3+3*alpha-gamma,
    "bridge quadrature epsilon=A^2":4+alpha,
    "RAW square-root cutoff":3+beta/2,
    "RAW cutoff over sqrt(w)":3+beta-alpha/2,
    "RAW mean mismatch":4+beta,
    "outer fixed target":F(4),
    "force mean target":F(4),
    "optional own-mean compiler":F(4),
}
for name, power in nonbuffer.items():
    check("concrete nonbuffer: "+name, power>=grade, str(power))
check("smallest strict nonbuffer exponent", min(x for x in nonbuffer.values() if x>grade)==F(59,15))

# Each row is the power AFTER the terminal A-Lipschitz force.
rows = [
    ("retuned mean and corrected Gram native",F(5),F(1)),
    ("mixed K native prior",F(6),F(7,6)),
    ("covariance target restoration",F(5),F(1,2)),
    ("mixed K mixture",F(7),F(3,2)),
    ("cubic target mismatch",F(5),F(1)),
    ("cubic order-five prior",F(6),F(2)),
    ("cubic feedback and joint pure-cubic feedback",F(7),F(5,2)),
    ("quartic native and smoothing target",F(6),F(2)),
    ("quartic F2 to H target",F(6),F(3,2)),
    ("corrected Gram sixth current and order-six prior",F(7),F(5,2)),
    ("Gram j-to-g target before absorption",F(6),F(3,2)),
    ("joint cubic-quartic feedback",F(8),F(3)),
    ("joint pure-quartic feedback",F(9),F(7,2)),
    ("true quartic consumer leading",F(6),F(2)),
]
ledger = []
for name, a, p in rows:
    bulk=a-p*alpha
    endpoint=a-(p-1)*beta if p>1 else a
    check("concrete bulk: "+name, bulk>grade, str(bulk))
    check("concrete endpoint: "+name, endpoint>grade, str(endpoint))
    ledger.append({"name":name,"a":str(a),"p":str(p),"bulk":str(bulk),"endpoint":str(endpoint)})
check("local-buffer rows strictly above four", min(min(F(x["bulk"]),F(x["endpoint"])) for x in ledger)>4)

# Rational grid independently tests the bounded-order family; the proof uses
# the affine differences recorded in the audit rather than this grid alone.
for n in range(1,100):
    eps=F(n,100)
    al=1-2*eps/3; be=ga=2-eps; gr=4-eps
    b=math.floor(3+2/eps)+1
    if F(b)<=3+2/eps: b+=1
    check(f"family {n}: balance",2+ga==3+3*al-ga==gr)
    check(f"family {n}: native order",b-4-be*F(b-3,2)>0)
    check(f"family {n}: residual first",F(3,2)-be/4>1)
    for name,a,p in rows:
        check(f"family {n}: bulk {name}", a-p*al>gr)
        check(f"family {n}: endpoint {name}", a-(p-1)*be>gr if p>1 else a>gr)
    for x in [1+2*ga,3+al,4+al,3+be/2,3+be-al/2,4+be,F(4)]:
        check(f"family {n}: other forcing {x}",x>gr)

# Exact geometry and positive dyadic-rule fixtures, including the truncated
# final panel. This verifies the executing-sum estimates, not an integral proxy.
for A in [2.0**(-k) for k in [5,10,15,20]]:
    w=A**float(alpha); eta=A**float(beta); q=1-w
    boundaries=sorted(set([0.,q,1-eta]+[1-2.**(-j) for j in range(1,100) if 0<1-2.**(-j)<1-eta]))
    nodes=[((a+b)/2,b-a) for a,b in zip(boundaries,boundaries[1:])]
    for t,omega in nodes:
        sig2=w*(2-w); c2=(1-t)*(1+t); s=sig2+q*q*c2; v=c2*sig2/s
        check("reciprocal conditional variance",abs(1/v-(q*q/sig2+1/c2)) <= 1e-10/v)
        check("standard Y marginal",abs((q*t)**2+s-1)<1e-14)
        check("conditional buffer lower bound",v>=eta/2*(1-1e-7))
        # Conditional mean coefficient a(z,Y), Y=qt z+sqrt(s) GY.
        az=t*sig2/s; ay=q*c2/s
        check("source-zero z row",abs(az+ay*q*t-t)<1e-12)
        check("source-zero total private variance",abs(ay*ay*s+v-c2)<1e-12)
    for p in [0.5,1,7/6,1.5,2,2.5,3,3.5]:
        total=sum(omega*(((1-t)*(1+t)*w*(2-w)/(w*(2-w)+q*q*(1-t)*(1+t)))**(-p)) for t,omega in nodes)
        bound=w**(-p)+(eta**(1-p) if p>1 else math.log(1/eta) if p==1 else 1)
        check(f"finite positive sum p={p}",total/bound<6, total/bound)

# Actual variance budget: five outer services and half keep; corrected Gram
# has its own three shares entirely inside its outer allocation.
u=F(1,10)
check("outer variance allocation",5*u+F(1,2)==1)
check("corrected Gram nested allocation",u/2+u/4+u/4==u)
for theta in [F(1,5),F(1,2),F(4,5)]:
    check("Gram visible/internal/external keep budget",u/2+theta*u/4+(1-theta)*u/4+u/4==u)

# Independent exact simultaneous degree-two/degree-three Gaussian regression
# and conditional-centering identity. This detects omission of mixed feedback.
S,U,V,t=sp.symbols("S U V t", real=True)
x1=S/sp.sqrt(3)+U/sp.sqrt(2)+V/sp.sqrt(6)
x2=S/sp.sqrt(3)-U/sp.sqrt(2)+V/sp.sqrt(6)
P=(x1*x1-1)/10+(x2**3-3*x2)/17
def normal_moment(n):
    return 0 if n%2 else sp.factorial2(n-1) if n else 1
def mean_uv(poly):
    ans=0
    for (i,j),coef in sp.Poly(sp.expand(poly),U,V).terms():
        ans+=coef*normal_moment(i)*normal_moment(j)
    return sp.simplify(ans)
projection=mean_uv(P)
check("simultaneous Hermite regression",sp.simplify(projection-((S*S-1)/30+(S**3-3*S)/(51*sp.sqrt(3))))==0)
E=sp.expand(P-projection)
check("regressed complement exactly centered",mean_uv(E)==0)
H=[sp.Integer(1),U,U*U-1,U**3-3*U]
Hv=[x.subs(U,V) for x in H]
coef={}
# Monomial-to-Hermite basis for degrees at most three in each variable.
mono={0:[(0,1)],1:[(1,1)],2:[(2,1),(0,1)],3:[(3,1),(1,3)]}
for (i,j),c in sp.Poly(E,U,V).terms():
    for ii,ci in mono[i]:
        for jj,cj in mono[j]:
            coef[ii,jj]=sp.simplify(coef.get((ii,jj),0)+c*ci*cj)
check("zero complement chaos absent",sp.simplify(coef.get((0,0),0))==0)
potential=sp.expand(sum(c*H[i]*Hv[j]/(i+j) for (i,j),c in coef.items() if i+j))
REgrad=sp.expand(sp.diff(potential,U)*sp.diff(E,U)+sp.diff(potential,V)*sp.diff(E,V))
base=S+projection
W=base+t*E
for degree in range(1,6):
    lhs=mean_uv(E*degree*W**(degree-1))
    rhs=0 if degree==1 else t*mean_uv(REgrad*degree*(degree-1)*W**(degree-2))
    check(f"positive simultaneous centering identity degree {degree}",sp.simplify(lhs-rhs)==0)

# Orthogonal known-row alignment preserves each node's entire original
# Gaussian law; shared roots need not preserve cross-node independence.
rng=np.random.default_rng(9052026)
for n in [3,7,19]:
    row=rng.normal(size=n); row/=np.linalg.norm(row)
    target=np.zeros(n); target[0]=1
    d=row-target
    Hm=np.eye(n)-2*np.outer(d,d)/(d@d)
    check(f"known row rotation {n}: orthogonal",np.linalg.norm(Hm@Hm.T-np.eye(n))<1e-12)
    check(f"known row rotation {n}: exact active row",np.linalg.norm(row@Hm-target)<1e-12)

pins={}
for name in ["fourth-cumulant-return-20261005","retuned-gram-native-20261005","combined-bridge-skew-20261005","combined-mean-reentry-20261005"]:
    path=Path("/workspace/shared")/name/"MANIFEST.json"
    pins[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
out={"status":"PASS","assertions":len(checks),"scope":"Independent scalar/algebraic diagnostics; native compiler not numerically executed.","reviewed_input_manifests":pins,"concrete_buffer_ledger":ledger,"checks":checks}
(HERE/"independent_checks.json").write_text(json.dumps(out,indent=2)+"\n")
print(f"PASS: {len(checks)} independent assertions")
