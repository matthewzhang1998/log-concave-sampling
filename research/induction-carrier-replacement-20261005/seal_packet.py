#!/usr/bin/env python3
"""Verify the source and audit pins and create a deterministic package manifest."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
inputs=json.loads((ROOT/'INPUT-PINS.json').read_text())
for row in inputs['inputs']:
 p=Path(row['path']); assert sha(p)==row['sha256'],str(p)
pins=json.loads((ROOT/'independent-audit/REVIEW-PIN.json').read_text())
for name,h in pins['source_sha256'].items():
 assert sha(ROOT/name)==h, f'audit pin mismatch: {name}'
checks=json.loads((ROOT/'independent-audit/independent-check-results.json').read_text())
for name,h in checks['source_sha256'].items():
 assert sha(ROOT/name)==h, f'independent check source mismatch: {name}'
files=[]
for p in sorted(ROOT.rglob('*')):
 if not p.is_file() or '__pycache__' in p.parts or p.name in {'MANIFEST.json','SHA256SUMS'}: continue
 files.append({'path':str(p.relative_to(ROOT)),'sha256':sha(p),'bytes':p.stat().st_size})
manifest={'date':'2026-10-05','status':'sealed; independent audit passed under stated imported interfaces',
 'scope':['Actual original shared-carrier replacement with explicit finite-dimensional smallness guard.',
          'Independent variance-partition native alternative with explicit count/share and shrinking-eta eligibility.',
          'No public-log-only dimensional shared-carrier theorem; no all-order or sublinear-exponent claim.'],
 'audit':'independent-audit/CARRIER-REPLACEMENT-AUDIT.md','files':files}
(ROOT/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
allfiles=files+[{'path':'MANIFEST.json','sha256':sha(ROOT/'MANIFEST.json')}]
(ROOT/'SHA256SUMS').write_text(''.join(f"{r['sha256']}  {r['path']}\n" for r in sorted(allfiles,key=lambda r:r['path'])))
print(json.dumps({'files':len(allfiles),'manifest_sha256':sha(ROOT/'MANIFEST.json'),'source_and_audit_pins_verified':True},indent=2))
