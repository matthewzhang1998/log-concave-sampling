from fractions import Fraction as F
from pathlib import Path
import json,hashlib
checks=0
for gi in range(1,101):
 g=F(gi,100)
 for zi in range(1,20):
  z=g*F(zi,100); gamma=(g-z)/2
  assert 2*z<gamma<g-3*z;checks+=1
  assert gamma-2*z>0 and g-3*z-gamma>0;checks+=1
  assert 1+g-z-2*gamma==1;checks+=1
  assert 1+g+z<F(7,3);checks+=1
 for z,gamma in [(g/10,F(9,20)*g)]:
  assert gamma-2*z==g/4 and g-3*z-gamma==g/4;checks+=1
p=Path(__file__).parent
f=p/'INDEPENDENT-SCALAR-SHEAR-RETENTION-BOUNDARY-AUDIT.md'
x={'status':'PASS','checks':checks,'scope':'exact rational law/first/cross-width/variance exponent checks; source graph proof is in the text','audit_sha256':hashlib.sha256(f.read_bytes()).hexdigest()}
(p/'retention_guard_checks.json').write_text(json.dumps(x,indent=2)+'\n')
print(json.dumps(x,indent=2))
