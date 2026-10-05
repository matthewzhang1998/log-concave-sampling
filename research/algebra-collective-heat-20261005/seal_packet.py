#!/usr/bin/env python3
"""Pin completed Task B artifacts and recheck read-only source hashes."""
from pathlib import Path
import hashlib,json,datetime
ROOT=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
pins=json.loads((ROOT/'INPUT-PINS.json').read_text())
checks=[]
for row in pins['sources']:
 p=Path(row['path'])
 checks.append({'path':str(p),'expected_sha256':row['sha256'],'actual_sha256':sha(p),'matches':sha(p)==row['sha256']})
assert all(row['matches'] for row in checks),'A pinned input changed; do not silently repin it.'
author=json.loads((ROOT/'checks.json').read_text())
assert author['status']=='PASS' and all(c['passed'] for c in author['checks'])
audit=json.loads((ROOT/'independent-audit/AUDIT-RESULT.json').read_text())
assert audit['verdict']=='PASS_FOR_REVISED_SCOPED_CLAIMS' and not audit['blocking_findings']
audited=json.loads((ROOT/'independent-audit/AUDITED-FILES.json').read_text())
for row in audited['files']:
 assert sha(ROOT/row['path'])==row['sha256'], 'Audited author file changed: '+row['path']
verification={'sealed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'source_inputs_unchanged':True,'source_count':len(checks),'source_checks':checks,
 'author_checks':author['test_count'],'author_status':'PASS',
 'independent_audit':audit['verdict'],'independent_checks':audit['independent_check_count'],'audited_versions_unchanged':True,
 'native_sampler_executed':False,'external_sharing':False,
 'scope':'Finite positive-covariance heat regrouping and guarded leading-current realization; no universal saturated-clock reserve gain; no exact all-orders native right inverse.'}
(ROOT/'SEAL-VERIFICATION.json').write_text(json.dumps(verification,indent=2)+'\n')
files=[p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name not in ['MANIFEST.json','SHA256SUMS']]
files.sort()
manifest={'files':[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in files],
 'source_inputs_unchanged':True,'native_sampler_executed':False}
(ROOT/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
files.append(ROOT/'MANIFEST.json');files.sort()
(ROOT/'SHA256SUMS').write_text(''.join(sha(p)+'  '+str(p.relative_to(ROOT))+'\n' for p in files))
print(json.dumps({'sealed_files':len(files),'source_count':len(checks),'author_checks':author['test_count']},indent=2))
