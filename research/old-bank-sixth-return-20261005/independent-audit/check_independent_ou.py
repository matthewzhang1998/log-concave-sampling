"""Independent finite checks. These do not certify the pending old-law remainder."""
from pathlib import Path
import hashlib, json, math
import numpy as np

ROOT=Path(__file__).resolve().parent
rng=np.random.default_rng(613440)
ratios=[]
cov_errors=[]
for _ in range(12000):
    # Include points near all five original endpoints independently.
    gaps=10.0**rng.uniform(-11,0,5)
    q,s,r1,r2,r3=np.sqrt(1-gaps)
    sig=np.sqrt((1-np.array([r1,r2,r3])**2)/2)
    A=np.array([
        [r1*np.sqrt(1-q*q),0,sig[0],0,0],
        [r2*s*np.sqrt(1-q*q),r2*np.sqrt(1-s*s),0,sig[1],0],
        [r3*s*np.sqrt(1-q*q),r3*np.sqrt(1-s*s),0,0,sig[2]],
    ])
    Gamma=A@A.T
    m=min(sig[0]**2,sig[2]**2,1-s*s)
    ratios.append(float(np.linalg.eigvalsh(Gamma)[0]/m))
    c=float(rng.random())
    h2=1-c
    t2=sig**2+h2*m/32
    residual=h2*(Gamma-m*np.eye(3)/32)
    assert np.linalg.eigvalsh(residual)[0]>=-1e-13
    cov_errors.append(float(np.max(np.abs(residual+np.diag(t2)-(h2*Gamma+np.diag(sig**2))))))
assert min(ratios)>1/16-1e-5
assert max(cov_errors)<1e-14

# The standard correlation c has r=sqrt(c), h²=1-c in its symmetric form.
from numpy.polynomial.legendre import leggauss
rules=[]
for J,n in [(4,4),(8,6),(16,10),(24,14)]:
    z,v=leggauss(n)
    nodes=[];weights=[]
    for j in range(J):
        lo=1-2.**(-j);hi=1-2.**(-j-1)
        nodes.extend((lo+hi)/2+(hi-lo)*z/2)
        weights.extend((hi-lo)*v/2)
    nodes=np.array(nodes);weights=np.array(weights)
    max_error=max(abs(float(weights@(nodes**k))-1/(k+1)) for k in list(range(100))+[200,1000,10000,1000000])
    envelope=float(np.sum(weights/np.sqrt(1-nodes)))
    assert np.all(weights>0)
    assert envelope<=1/(math.sqrt(2)-1)+1e-12
    assert max_error<=2**(-J)+10*4**(-n)
    rules.append(dict(J=J,n=n,count=J*n,error=max_error,bridge_mass_envelope=envelope))

# Normalization of all nine derivative-hit trees.
census=[]
for a in range(3):
    for b in range(3):
        orders=[0,1,0,0,1,0]
        orders[a]+=1;orders[3+b]+=1
        # A longest path whose endpoints are physical leaves of the original copies.
        spine=6-(a==1)-(b==1)
        census.append(dict(hit_pair=[a,b],C0=orders.count(0),C1=orders.count(1),C2=orders.count(2),spine=spine))
assert [sum(x[k] for x in census) for k in ('C0','C1','C2')]==[24,24,6]
assert [sum(x['spine']==k for x in census) for k in (4,5,6)]==[1,4,4]

# Generic trace shortcut remains false even if every proper tensor cut is one.
# F(q)=sin(q_0)*sum_i e_i^⊗6 has every physical proper cut <=1,
# while an HS-unit diagonal test has variance n*(1-e^-2)/2.
# This is a generic cut-only diagnostic, not a newly source-qualified witness.
trace=[dict(n=n,variance=n*(1-math.exp(-2))/2,derivative_root_to_tensor_HS=math.sqrt(n)) for n in [1,4,16,64,256]]

# A fixed-Q matrix-chaos bound cannot be promoted to sup_Q before integration.
# Over sign directions v in {+-1/sqrt(D)}^D, max_v(v.P)^2=||P||_1²/D.
sup_counter=[dict(D=d,fixed_direction_L2=math.sqrt(3),
                  expected_maximum=1+(d-1)*2/math.pi)
             for d in [1,4,16,64,256]]

snapshot=ROOT/'REVIEWED-OU-SNAPSHOT.md'
out=dict(
    reviewed_sha256=hashlib.sha256(snapshot.read_bytes()).hexdigest(),
    random_genealogies=12000,
    minimum_Gamma_over_m=min(ratios),
    maximum_covariance_error=max(cov_errors),
    positive_bridge_rules=rules,
    hit_census=census,
    derivative_adapter_totals=dict(C0=24,C1=24,C2=6),
    generic_cut_only_countertest=trace,
    conditional_supremum_countertest=sup_counter,
    verdict='Exact OU/heat/clock and conditional polynomial-current layers only. The full old-law same-endpoint remainder is not certified by this checker.'
)
(ROOT/'independent_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['reviewed_sha256','minimum_Gamma_over_m','maximum_covariance_error','positive_bridge_rules','verdict']},indent=2))
