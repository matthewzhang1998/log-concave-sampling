#!/usr/bin/env python3
"""Independent exact exponent/count checks; does not run any native compiler."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parent.parent
n = 0

def check(b):
    global n
    assert b
    n += 1

# Physical-readout exponents of the six LOW30 raw:split rows.
for j in range(1, 101):
    alpha = F(j, 100)
    powers = [4-F(3,2)*alpha, 4-alpha, 4-alpha/2,
              4-alpha/2, F(9,2)-F(3,2)*alpha, 5-F(3,2)*alpha]
    check(min(powers) == 4-F(3,2)*alpha)
    check(min(1, F(3,2)-alpha/2) == 1)
    check(1-alpha/2 > 0)
    check(F(3,2)-alpha > 0)
    check(1-alpha/3 > 0)
    check(1+alpha/6 >= 1)

for alpha, wanted in [(F(4,5), (F(3,5), F(7,10), F(11,15))),
                       (F(1), (F(1,2), F(1,2), F(2,3)))]:
    check((1-alpha/2, F(3,2)-alpha, 1-alpha/3) == wanted)

# Before alignment a bulk node owns Y plus its complete service banks;
# a near node owns X, Y and the two F_Q roots. Share exactly one row.
for bulk in range(1, 9):
    for near in range(0, 9):
        service_dims = sum(19+3*j for j in range(bulk))
        before = bulk + service_dims + 4*near
        after = before - (bulk+near-1)
        check(after == 1+service_dims+3*near)

# Check the affine error recurrence without assuming equal forcing terms.
for depth in range(3, 21):
    A = F(1, 7)
    e2 = F(2, 13)
    forcing = {j: F(j, 37+j) for j in range(3, depth+1)}
    e = e2
    for j in range(3, depth+1):
        e = A*e+forcing[j]
    check(e == A**(depth-2)*e2 + sum(A**(depth-j)*forcing[j]
                                        for j in range(3, depth+1)))

result = {
    'status': 'PASS', 'checks': n,
    'reviewed_sha256': hashlib.sha256((root/'REUSABLE-PORT-AND-DEPTH-RECURRENCE.md').read_bytes()).hexdigest(),
    'scope': 'Independent exact rational checks of native scale exponents, physical source closure, root-count algebra and finite affine recurrence only. No native VALUE/compiler execution or nonlinear theorem is numerically tested.'
}
print(json.dumps(result, indent=2))
Path(__file__).with_name('depth-port-checks.json').write_text(json.dumps(result, indent=2)+'\n')
