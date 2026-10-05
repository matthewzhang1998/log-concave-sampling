#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
R=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
pins=json.loads((R/'INPUT-PINS.json').read_text())['inputs']
for item in pins:assert sha(Path(item['path']))==item['sha256']
main='NEXT-HISTORY-M4-AND-MARKED-CURRENT-PORT.md'
assert sha(R/main)=='b85af0391c9e6064e7acd1386da0e0d7c339d092772ec6064ac1502005c6bf5a'
assert sha(R/main) in (R/'independent-audit/AUDIT-REPORT.md').read_text()
files={str(p.relative_to(R)):sha(p) for p in sorted(R.rglob('*')) if p.is_file() and p.name not in ['MANIFEST.json','SHA256SUMS'] and '__pycache__' not in str(p)}
manifest={'date':'2026-10-05','status':'Independent PASS within exact analytical target-isolation scope', 'target':'m4 and F3 conditional mean/covariance/skew/quartic marked currents; subsequent F4 covariance sign obstruction', 'constructive_status':'Finite positive original-VALUE marked-history supplier and accuracy-raising step remain open','new_results':['Exact coherent history-difference targets and dimension-safe ranks 1 through 4','Ordered shifted four-force response chain','Positive joint block-Gram and buffered marked current','Actual same-g negative next-history covariance increment rules out repeatable independent additive covariance banks'], 'checks':{'author_exact_assertions':1691,'independent_exact_assertions':704,'input_pins':len(pins)},'files':files,'inputs':pins,'nonclaims':['New native compiler admission','Accuracy-raising VALUE producer','RAW promotion of m3 mean graph','Arbitrary-order/end-to-end cost theorem','Endpoint distribution theorem']}
(R/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
entries={**files,'MANIFEST.json':sha(R/'MANIFEST.json')}
(R/'SHA256SUMS').write_text(''.join(h+'  '+p+'\n' for p,h in sorted(entries.items())))
print('Manifest SHA256 '+sha(R/'MANIFEST.json'))
print('Files '+str(len(files))+'; inputs '+str(len(pins)))
