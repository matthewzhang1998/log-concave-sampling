#!/usr/bin/env python3
"""Independent exact checks. PASS is scoped; tests do not establish analytic theorems."""
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
import hashlib, itertools, json

OUT = Path(__file__).resolve().parent
SOURCE = Path('/workspace/shared/connected-value-cumulant-source-20261005')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

integrity = []
for line in (SOURCE/'SHA256SUMS').read_text().splitlines():
    expected, name = line.split(maxsplit=1)
    actual = sha(SOURCE/name)
    assert actual == expected, name
    integrity.append({'path':str(SOURCE/name),'sha256':actual})
for name, expected in json.loads((SOURCE/'INPUT-PINS.json').read_text()).items():
    actual = sha(Path(name))
    assert actual == expected, name
    integrity.append({'path':name,'sha256':actual})

# Independent of the source's partition/Stirling enumeration: count surjections
# by inclusion-exclusion, then choose occupied replica labels.
variance = []
for n in range(1,9):
    total = F(0)
    total_maps = 0
    terms = []
    for k in range(1,n+1):
        onto = sum((-1)**j * comb(k,j) * (k-j)**n for j in range(k+1))
        maps = comb(n,k) * onto
        coeff = F((-1)**(k-1)*factorial(k-1), factorial(n)//factorial(n-k))
        total += maps*coeff**2
        total_maps += maps
        terms.append({'occupied_labels':k,'maps':maps,'coefficient':str(coeff)})
    assert total_maps == n**n
    variance.append({'rank':n,'variance_factor':str(total),'terms':terms})
assert variance[3]['variance_factor'] == '10/3'
assert variance[7]['variance_factor'] == '46212/35'

# Coefficientwise moment-cumulant equality. Integer partitions of the physical
# slot count describe the multiplicities of nonempty replica labels. The
# kernel's occupancy count divided by (n)_k must equal the formal log-MGF
# coefficient. This verifies every formal moment monomial through rank eight.
def integer_partitions(n, minimum=1):
    if n == 0:
        yield ()
    for r in range(minimum,n+1):
        for rest in integer_partitions(n-r,r):
            yield (r,)+rest

def log_mgf_coefficient(parts):
    k=len(parts)
    counts={r:parts.count(r) for r in set(parts)}
    return F((-1)**(k-1)*factorial(k-1)*factorial(sum(parts)),
             prod(factorial(r) for r in parts)*prod(factorial(m) for m in counts.values()))

def prod(values):
    answer=1
    for x in values: answer*=x
    return answer

formal_checks=[]
for n in range(1,9):
    checked=0
    for parts in integer_partitions(n):
        k=len(parts)
        multiplicities={r:parts.count(r) for r in set(parts)}
        label_assignments=F(factorial(n),factorial(n-k)*prod(factorial(m) for m in multiplicities.values()))
        ordered_slots=F(factorial(n),prod(factorial(r) for r in parts))
        kernel=label_assignments*ordered_slots*F((-1)**(k-1)*factorial(k-1),factorial(n)//factorial(n-k))
        expected=log_mgf_coefficient(parts)
        assert kernel==expected
        checked+=1
    formal_checks.append({'rank':n,'formal_moment_monomials':checked})

# Independent formal log series on a noncentered three-point scalar law.
# Exact occupancy expansion of the replica kernel is compared with the
# coefficient of log(sum E[B^j] t^j/j!).
points=[(F(-2),F(1,6)),(F(1),F(1,2)),(F(4),F(1,3))]
mom=[sum(p*x**j for x,p in points) for j in range(9)]
series=[mom[j]/factorial(j) for j in range(9)]
def poly_mul(a,b,N):
    out=[F(0)]*(N+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<=N: out[i+j]+=x*y
    return out

u=series[:]; u[0]=0
log=[F(0)]*9
power=[F(1)]+[F(0)]*8
for k in range(1,9):
    power=poly_mul(power,u,8)
    for j in range(9): log[j]+=F((-1)**(k-1),k)*power[j]
scalar=[]
for n in range(1,9):
    ek=sum(log_mgf_coefficient(parts)*prod(mom[r] for r in parts) for parts in integer_partitions(n))
    assert ek==factorial(n)*log[n]
    scalar.append({'rank':n,'cumulant':str(ek)})

# Eight orthogonal all-distinct centered-rank-four monomials.
centered={tuple([0]*4):F(1,2),tuple([1]*4):F(1,2)}
for chosen in itertools.combinations(range(4),2):
    centered[tuple(int(i in chosen) for i in range(4))]=F(-1,2)
assert len(centered)==8 and sum(x*x for x in centered.values())==2

# Exact sharpness examples for both rank inequalities, including biased ones.
rank_checks=[]
for D in [2,3,5,8]:
    for M in range(1,D+1):
        for c,eta in [(F(3,2),F(0)),(F(3,2),F(1,4))]:
            amplitude=F(D,M)*(c-eta)
            subsets=list(itertools.combinations(range(D),M))
            mean=[sum(amplitude for S in subsets if i in S)/len(subsets) for i in range(D)]
            assert mean==[c-eta]*D
            hs2=M*amplitude**2
            op2=amplitude**2
            assert hs2==F(D*D,M)*(c-eta)**2
            assert op2==F(D*D,M*M)*(c-eta)**2
            assert sum((x-c)**2 for x in mean)==D*eta**2
            rank_checks.append({'D':D,'M':M,'c':str(c),'eta':str(eta),'HS_squared':str(hs2),'operator_squared':str(op2)})

# Exact conservative explicit epsilon certificate. All mixed cumulant entries
# in aY, S, C have magnitude <= n^n (n-1)! 2^n. There are <=3^n terms.
epsilon=F(1,2**80)
witness=[]
for n in [4,8]:
    derivative_lower=F(factorial(n),2**n*n**n) # uses exp(-1/2)>1/2
    remainder_majorant=(6*n)**n*factorial(n-1)
    assert epsilon<=F(1,8)
    assert remainder_majorant*epsilon < derivative_lower/2
    witness.append({'rank':n,'epsilon':str(epsilon),'derivative_magnitude_lower':str(derivative_lower),'remainder_coefficient_upper':remainder_majorant,'epsilon_times_remainder_divided_by_derivative_bound':str(epsilon*remainder_majorant/derivative_lower),'kappa_upper_strictly_below':str(-epsilon*derivative_lower/2)})
assert F(1,2)-F(1,8)-F(1,8)**2==F(23,64)
assert F(1,2)+F(1,8)+F(1,8)**2==F(41,64)

# Finite OU covariance witness at t_j = j log 2. Covariance e^{-|t-s|}
# minus e^{-(t+s)} is then rational; positivity is checked exactly.
times=[1,2,3]
weights=[F(1,3),F(1,4),F(5,12)]
cov=[[F(1,2**abs(j-k))-F(1,2**(j+k)) for k in times] for j in times]
yq_cov=[sum(weights[k]*cov[k][j] for k in range(3)) for j in range(3)]
assert all(x>0 for x in yq_cov)
variance_yq=sum(weights[j]*yq_cov[j] for j in range(3))
assert variance_yq>0

# Quantifier tests: these are counterexamples to overextensions, not to the
# correctly scoped theorem. For a standard caller, c(y)=y^2-1 has c(0)=-1
# and mean zero. T=0 therefore matches the unconditional averaged target,
# though its conditional error is nonzero. The first moments are exact.
normal2=F(1); normal4=F(3)
assert normal2-1==0
assert normal4-2*normal2+1==2

result={'status':'PASS_SCOPED','scope':'Exact finite algebra, source integrity, sharpness examples and an analytic-witness arithmetic certificate. Infinite-dimensional OU identities and nuclear-norm proofs are audited in the companion text; diagnostics alone are not proofs. No active native-source construction or impossibility is certified.',
 'source_integrity':integrity,'variance':variance,'formal_moment_coefficients':formal_checks,
 'noncentered_exact_cumulants':scalar,'centered_rank4_variance_factor':2,
 'sharp_rank_examples':rank_checks,'explicit_witness':witness,
 'finite_ou_covariances':[str(x) for x in yq_cov],'finite_ou_variance':str(variance_yq),
 'unconditional_bias_counterexample':{'c':'y^2 - 1','E_c':'0','E_c_squared':'2','T':'0','meaning':'Unconditional averaging does not preserve the m=E|c| bound.'}}
(OUT/'independent-checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['variance','sharp_rank_examples']},indent=2))
