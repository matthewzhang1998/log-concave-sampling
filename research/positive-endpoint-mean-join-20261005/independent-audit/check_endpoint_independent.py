#!/usr/bin/env python3
"""Independent algebra/geometry diagnostics; not a native compiler execution."""
from pathlib import Path
from fractions import Fraction as F
import hashlib, json, math
import sympy as s

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
checks=[]
def check(name, condition, evidence=None):
    if not bool(condition):
        raise AssertionError(name)
    checks.append({'name':name, **({'evidence':str(evidence)} if evidence is not None else {})})

def zero(x): return s.expand(x)==0
p=F(39,10); q=F(59,15); rho=F(3,4); h=F(1,2)
# 1. Exact allocations, exponents, and the source-zero row.
check('conditional source-zero buffer', h*h+h*h*F(1,2)+F(1,8)+F(1,4)==1-h*h)
check('full source-zero variance', h*h+h*h+h*h*F(1,2)+F(1,8)+F(1,4)==1)
check('regressed visible baseline', h*h+h*h*F(1,2)+F(1,8)==F(1,2))
check('untouched keep fraction', F(1,4)/(1-h*h)==F(1,3))
check('strict public-log margin', q-p==F(1,30))
check('physical grade', p+F(1,2)==F(22,5))
check('terminal squared-clock exponent', F(2,5)*(p-2)==F(19,25))
check('smallest alpha exponent', 1+F(2,5)*(p-2)==F(44,25))
check('first bridge spatial exponent', 2*p-1==F(34,5))
check('second bridge spatial exponent', 2*q-1==F(103,15))
check('native skew normalized radius', s.Rational(1,2)/s.sqrt(s.Rational(1,8))==s.sqrt(2))
check('negative skew tensor scale', -h**3/F(6)==-F(1,48))

# 2. Exact coherent linear F3 example, using the OU history itself.
a=s.symbols('a', real=True)
Cfull=s.zeros(3); Ccond=s.zeros(3)
means=s.Matrix([s.Rational(1,2**(i+1)) for i in range(3)])
for i in range(3):
    for j in range(3):
        Cfull[i,j]=(s.factorial(i+j+1)/(2**(i+j+2)*s.factorial(i)*s.factorial(j))
                    *(s.Rational(1,i+1)+s.Rational(1,j+1)))
        Ccond[i,j]=Cfull[i,j]-means[i]*means[j]
        check(f'OU conditional covariance identity {i},{j}',zero(Cfull[i,j]-Ccond[i,j]-means[i]*means[j]))
c=s.Matrix([a,-a*a,a**3]); m=(c.T*means)[0]
check('canonical linear m3',zero(m-(a/2-a*a/4+a**3/8)))
varF=(c.T*Cfull*c)[0]; condF=(c.T*Ccond*c)[0]
varP=1-2*m+varF
varJoined=s.Rational(1,4)*(1-m)**2+s.Rational(3,4)+s.Rational(1,4)*condF
check('full mean/covariance Gaussian target matches buffered P3',zero(varJoined-(s.Rational(1,4)*varP+s.Rational(3,4))))
check('P3 Gaussian covariance agrees with target through cubic order',zero(s.series(varP-1/(1+a),a,0,4).removeO()))
check('P3 quartic discrepancy is generally nonzero',s.expand(s.series(varP-1/(1+a),a,0,5).removeO()).coeff(a,4)!=0)
check('posterior force is a different target',not zero(m-(a*s.Rational(1,2)/(1+s.Rational(3,4)*a))))
check('F3-F2 covariance difference begins at fourth order',min(term[0][0] for term in s.Poly(s.expand(condF-(s.Matrix([a,-a*a,0]).T*Ccond*s.Matrix([a,-a*a,0]))[0]),a).terms())==4)

# 3. Anisotropic Hermite regression, including trace term and physical packet fills.
x=s.Matrix(s.symbols('x0:3'))
B=s.Matrix([[1,2,0],[0,1,-1],[2,0,1]])
Sigma=B*B.T/s.Integer(200)
C=s.eye(3)/2+Sigma/4
Ci=C.inv()
packet_variances=[s.Rational(1,24),s.Rational(1,12)]
beta2=[3*v/4 for v in packet_variances]
fill=sum(v/4 for v in packet_variances)
J=s.eye(3)*(s.Rational(1,4)+s.Rational(1,8)+fill)+Sigma/4
check('packet fill stays inside one-eighth share',sum(beta2)+fill==s.Rational(1,8))
check('anisotropic total reference covariance',J+sum(beta2)*s.eye(3)==C)
for n,b in enumerate(beta2):
    conditional_mean=b*Ci*x
    conditional_cov=b*s.eye(3)-b*b*Ci
    projected=(conditional_cov+conditional_mean*conditional_mean.T)/(b*b)-s.eye(3)/b
    target=Ci*x*x.T*Ci-Ci
    for i in range(3):
        for j in range(3):
            check(f'anisotropic full Hermite projection packet {n}, entry {i},{j}',zero(projected[i,j]-target[i,j]))
    wrong=conditional_mean*conditional_mean.T/(b*b)-s.eye(3)/b
    check(f'trace-free regression shortcut is detected packet {n}',any(not zero(wrong[i,j]-target[i,j]) for i in range(3) for j in range(3)))

# 4. Output-slot asymmetry has zero leading contraction but a real quadratic map.
from itertools import permutations, product
N={idx:s.Rational(1+idx[0]*7+(idx[1]+idx[2])*3+idx[1]*idx[2],97) for idx in product(range(3),repeat=3)}
Ns={idx:sum(N[tuple(j)] for j in permutations(idx))/6 for idx in N}
D={idx:N[idx]-Ns[idx] for idx in N}
check('nontrivial output-slot defect',any(v!=0 for v in D.values()))
for idx in N:
    check(f'full defect symmetrization {idx}',sum(D[tuple(j)] for j in permutations(idx))==0)
phi=(x[0]+2*x[1]-3*x[2])**5+x[0]*x[1]*x[2]
contract=sum(v*s.diff(phi,x[i],x[j],x[k]) for (i,j,k),v in D.items())
check('leading output-slot current vanishes exactly',zero(contract))
poly_defect=s.Matrix([sum(D[(i,j,k)]*(x[j]*x[k]-int(j==k)) for j in range(3) for k in range(3)) for i in range(3)])
check('same leading cumulant does not mean identical positive map',any(not zero(z) for z in poly_defect))

# 5. Exact same-endpoint positive complement interpolation, no removed roots observed.
z,u,k,t=s.symbols('z u k t')
S=3*z/2
E=S*u/17+(u*u-1)/29
RE=S/17+u/29
DE=s.diff(E,u)
W=S+(S*S-s.Rational(9,4))/31+t*E+k/2

def gm(n):
    return s.Integer(0) if n%2 else s.factorial2(n-1) if n else s.Integer(1)
def expectation(poly, variables):
    total=0
    for powers,coef in s.Poly(s.expand(poly),*variables).terms():
        total+=coef*math.prod(gm(power) for power in powers)
    return s.expand(total)
check('complement is conditionally centered',expectation(E,[u])==0)
for degree in range(1,6):
    lhs=s.diff(expectation(W**degree,[z,u,k]),t)
    rhs=0 if degree==1 else t*degree*(degree-1)*expectation(RE*DE*W**(degree-2),[z,u,k])
    check(f'exact positive Riesz current degree {degree}',zero(lhs-rhs))
check('quadratic skew map adds order-six covariance',zero(expectation((z+a**3*(z*z-1))**2,[z])-1-2*a**6))

# 6. Rational bridge identities and numerical schedule bounds for many scales.
for s2 in [F(1),F(3,4),F(9,16),F(1,16),F(1,1024)]:
    r2=1-s2; delta=s2/4; t2=r2+delta
    # Squared identities avoid introducing irrational fixtures.
    b2=delta*delta/(t2*s2*s2); v0=delta/t2
    v=(1-t2)*delta/(t2*s2)
    check(f'exact bridge mean square s2={s2}',b2*s2==v0*h*h)
    check(f'exact bridge variance s2={s2}',v==v0*(1-h*h))
    # a+b*r=r/t after removing common r/t.
    check(f'exact reference contraction s2={s2}',(1-t2)/s2+delta/s2==1)
for log2A in [1,2,4,8,16,32,64,128,256,512]:
    A=2.0**(-log2A)
    stop=A**float(F(19,25)); J=math.ceil(math.log(stop)/math.log(float(rho)))
    # Float first-index correction, so exact-power near-integers cannot shift J.
    while float(rho)**J>stop: J+=1
    while J>0 and float(rho)**(J-1)<=stop:J-=1
    sJ2=float(rho)**J
    check(f'first terminal threshold A=2^-{log2A}',float(rho)*stop*(1-1e-13)<sJ2<=stop*(1+1e-13))
    check(f'terminal debt reaches grade p A=2^-{log2A}',sJ2**2.5<=A**float(p-2)*(1+1e-12))
    for exponent in [p,q,F(4)]:
        powers=[float(rho)**((j-1)*float(exponent+F(1,2)))/4 for j in range(1,J+1)]
        upper=1/(4*(1-float(rho)**float(exponent+F(1,2))))
        check(f'exact-kernel weighted geometric error k={exponent}, A=2^-{log2A}',sum(powers)<=upper*(1+1e-13))
    # Local alpha logarithm never exceeds 44/25 times base log; fixed factor for terminal.
    check(f'buffered alpha log bound A=2^-{log2A}',1+(J-1)*(-math.log(float(rho)))/(log2A*math.log(2))<=float(F(44,25))+1e-13)

# 7. Source pins. These guard audit identity, not native numerics.
sources=[
 ROOT/'MEAN-COVARIANCE-SKEW-ENDPOINT.md',
 Path('/workspace/shared/quartic-corrected-mean-join-20261005/ACTUAL-QUARTIC-CORRECTED-MEAN-JOIN.md'),
 Path('/workspace/shared/quartic-corrected-mean-join-20261005/independent-audit/INDEPENDENT-JOIN-AUDIT.md'),
 Path('/workspace/shared/shrinking-buffer-skew-join-20261005/POSITIVE-SKEW-SHRINKING-BUFFER-JOIN.md'),
 Path('/workspace/shared/fourth-cumulant-return-20261005/comparison/FOURTH-CONDITIONAL-CUMULANT-COMPARISON.md'),
]
base=Path('/workspace/shared/recovery-20261004/cost/endpoint-stein-quadrature-20261004/positive-law-repair/higher-cumulant-gate')
sources += [base/x for x in [
 'THIRD-ORDER-FULL-REVERSE-OU-LAW-AND-ENDPOINT.md',
 'same-carrier-feedback/SAME-CARRIER-P3-MEAN-AND-FULL-COVARIANCE-GATE.md',
 'same-carrier-feedback/FULL-COVARIANCE-CLOSURE-AND-RESUMMED-MEAN-FRONTIER.md',
 'rank3-continuation/POSITIVE-RANK3-PACKET-FEEDBACK-LEMMA.md',
 'rank3-continuation/POSITIVE-QUADRATIC-REFERENCE-FOR-SKEW-GAUSSIANIZATION.md',
 'order-reentry/independent-audit/INDEPENDENT-NATIVE-RANK3-PACKET-AUDIT.md',
]]
pins=[]
for path in sources:
    data=path.read_bytes()
    pins.append({'path':str(path),'sha256':hashlib.sha256(data).hexdigest()})
    check(f'readable nonempty source {path.name}',len(data)>0)
# The exact graph and imported pin ledgers are independently checked as data.
graph=json.loads((ROOT/'ENDPOINT-GRAPH.json').read_text())
check('graph baseline total',sum(F(b['baseline_after_readout']) for b in graph['banks'])==F(3,4))
check('graph mean is complete LAW',next(b for b in graph['banks'] if b['id']=='M')['input_type']=='completed own-mean LAW')
check('graph component identities',set(b['id'] for b in graph['banks'])=={'M','H','K','V3','keep'})
for ledger in [ROOT/'INPUT-PINS.json',Path('/workspace/shared/quartic-corrected-mean-join-20261005/INPUT-PINS.json')]:
    for n,pin in enumerate(json.loads(ledger.read_text())['inputs']):
        check(f'imported hash {ledger.parent.name} {n}: {Path(pin['path']).name}',hashlib.sha256(Path(pin['path']).read_bytes()).hexdigest()==pin['sha256'])
for filename in ['ENDPOINT-GRAPH.json','INPUT-PINS.json']:
    path=ROOT/filename
    pins.append({'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
result={'status':'PASS','assertions':len(checks),'scope':'Exact algebra, target examples, reference regression/current, bridge and schedule diagnostics. Imported native compilers are not numerically executed.','checks':checks,'sources':pins}
(HERE/'independent_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':result['status'],'assertions':result['assertions'],'native_compilers_executed':False},indent=2))
