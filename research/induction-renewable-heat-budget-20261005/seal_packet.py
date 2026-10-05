#!/usr/bin/env python3
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
pins=json.loads((ROOT/'INPUT-PINS.json').read_text())
for item in pins['inputs']:
    assert sha(Path(item['path']))==item['sha256'],item['path']
assert (ROOT/'independent-audit'/'INDEPENDENT-CLOCK-MASS-AUDIT.md').exists(),'Independent report required before sealing'
assert (ROOT/'FINAL-STATUS.md').exists()
files=[p for p in sorted(ROOT.rglob('*')) if p.is_file() and '__pycache__' not in p.parts and p not in {ROOT/'MANIFEST.json',ROOT/'SHA256SUMS'}]
manifest={
 'date':'2026-10-05','status':'two-step source/positive-current theorem; finite clock credit; not arbitrary-order closure',
 'result':{'rank':8,'heat_exponents':['1/8','1/7','1/5'],'common_Gamma':'1/5','Psi':['879/140','226/35','228/35'],'gains':['5/28','2/35'],'final_current_grade':'41/5','outside_old_first':'O(alpha)','new_local_inverse_alpha_query_exponent':0,'raw_values_per_clock_tuple':513395200},
 'qualification':'Relative to sealed native interfaces and actual radius/first/gap/readout/precision admission. Common H is the total source-zero row from independently owned groups. No raw shared-carrier replacement or completed higher-order endpoint is asserted.',
 'independent_audit':'independent-audit/INDEPENDENT-CLOCK-MASS-AUDIT.md',
 'files':[{'path':str(p.relative_to(ROOT)),'sha256':sha(p),'bytes':p.stat().st_size} for p in files]
}
(ROOT/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
files.append(ROOT/'MANIFEST.json')
(ROOT/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.relative_to(ROOT)}\n' for p in sorted(files)))
print('Verified all',len(pins['inputs']),'input pins; sealed',len(files),'files.')
