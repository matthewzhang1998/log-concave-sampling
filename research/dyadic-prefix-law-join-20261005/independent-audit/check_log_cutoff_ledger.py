from fractions import Fraction as F
from pathlib import Path
import json, hashlib

out = Path(__file__).parent
# Bounds are actual positive outer sums: p>1 gives w^-p+eta^(1-p).
# w=A, eta=A^2 K^2. Logarithmic K-powers do not change A-grade.
rows = {
    'order_six_mean_prior': (F(7), F(5,2)),
    'native_cubic_and_calibration': (F(5), F(1)),
    'native_curl': (F(5), F(1,2)),
    'native_centered_reserve': (F(11,2), F(5,4)),
    'native_quadratic_reserve': (F(6), F(3,2)),
    'mixed_K_native': (F(6), F(7,6)),
    'quartic_consumer_packet': (F(6), F(2)),
    'corrected_Gram_sixth_and_cubic_feedback': (F(7), F(5,2)),
    'fourth_target': (F(6), F(3,2)),
    'mixed_K_mixture': (F(7), F(3,2)),
}
ledger = {}
for name, (a,p) in rows.items():
    grades = [a-p]
    kpowers = [F(0)]
    if p>1:
        grades.append(a+2*(1-p)); kpowers.append(2*(1-p))
    elif p==1:
        grades.append(a); kpowers.append('log(1/eta)')
    else:
        grades.append(a); kpowers.append(F(0))
    assert min(grades)>=4,(name,grades)
    ledger[name]={'A_grades':list(map(str,grades)), 'K_powers':list(map(str,kpowers))}

alpha=F(1); gamma=F(2); d=F(1,2); B=4
prefix={
    'native_mean': 1+B*(1+F(3,2)*alpha)-(B-1)*gamma,
    'mean_quadrature_eps_A2':2+alpha+2,
    'covariance_calibration':3+3*alpha+d-gamma,
    'centered_action':5+6*alpha-3*gamma-d/2,
    'covariance_feedback':5+6*alpha-3*gamma,
    'covariance_clock_eps_A':3+3*alpha-gamma+1,
    'true_prefix_third':4+F(9,2)*alpha-2*gamma,
    'covariance_caller_shift':4+F(5,2)*alpha-gamma-d/2,
    'mean_caller_shift_and_nested_omission':3+alpha,
    'smoothing_commutator':2+gamma,
    'smoothing_drift':1+2*gamma,
}
assert min(prefix.values())==4
# Native half-step source exponents are positive through fixed order six.
b=F(2); step=F(1,16); least=F(1)
while b<6:
    e=step/b
    assert 0<e<1
    least=min(least,1-e)
    b+=step
report={'status':'PASS', 'outer_rows':ledger, 'prefix_grades':{k:str(v) for k,v in prefix.items()}, 'minimum_native_halfstep_positive_exponent':str(least), 'scope':'Exact scalar ledger only; not execution of native compilers.'}
(out/'log_cutoff_ledger_checks.json').write_text(json.dumps(report,indent=2)+'\n')
inputs=[
 '/workspace/shared/dyadic-bridge-square-20261005/DYADIC-POSITIVE-BRIDGE-COVARIANCE.md',
 '/workspace/shared/multilevel-covariance-independent-audit-20261005/INDEPENDENT-AUDIT.md',
 '/workspace/shared/quartic-corrected-mean-join-20261005/ACTUAL-QUARTIC-CORRECTED-MEAN-JOIN.md',
 '/workspace/shared/v9-curation-work/frozen/prerequisites/research-source/High Acc Ideas/ai-bucket/30_low_acc.tex']
pins={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in inputs}
(out/'INPUT-PINS.json').write_text(json.dumps(pins,indent=2)+'\n')
print(json.dumps(report,indent=2))
