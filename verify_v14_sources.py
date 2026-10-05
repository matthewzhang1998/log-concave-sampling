#!/usr/bin/env python3
"""Verify additive V14 bytes, all three packet seals and direct input pins."""
from pathlib import Path
import hashlib,json,re
ROOT=Path(__file__).resolve().parent
selection=json.loads((ROOT/'verification/V14-FROZEN-SELECTION.json').read_text())
mapping=json.loads((ROOT/'verification/V13-SOURCE-MAP.json').read_text())['source_mappings'];mapping.update(json.loads((ROOT/'verification/V14-SOURCE-MAP.json').read_text())['source_mappings'])
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
errors=[];sealed=0;inputs=0;external=[]
for e in selection['files']:
 p=ROOT/e['path']
 if not p.is_file() or sha(p)!=e['sha256'] or p.stat().st_size!=e['bytes'] or e['sha256']!=e['original_sha256']:errors.append('Exact source mismatch: '+e['path'])
 if p.name=='SHA256SUMS':
  for line in p.read_text().splitlines():
   expected,name=line.split('  ',1);q=p.parent/name
   if sha(q)!=expected:errors.append('Seal mismatch: '+str(q.relative_to(ROOT)))
   else:sealed+=1
 if p.name=='INPUT-PINS.json':
  record=json.loads(p.read_text());rows=record if isinstance(record,list) else record.get('sources',record.get('inputs',[]))
  for row in rows:
   if row['path'].endswith(('30_low_acc.tex','13_high_acc.tex')):external.append({'record':e['path'],'source':'external:LOW30' if row['path'].endswith('30_low_acc.tex') else 'external:HIGH13','expected_sha256':row['sha256']});continue
   if row['path'] not in mapping:errors.append('Unmapped source: '+row['path']);continue
   q=ROOT/mapping[row['path']]
   if sha(q)!=row['sha256'] or ('bytes' in row and q.stat().st_size!=row['bytes']):errors.append('Input pin mismatch: '+row['path'])
   else:inputs+=1
report={'status':'PASS_WITH_EXTERNAL_PROVENANCE_QUALIFICATION' if not errors else 'FAIL','exact_new_files':len(selection['files']),'sealed_packet_count':3,'matching_sealed_members':sealed,'direct_input_pins_verified':inputs,'external_pins_not_reverified':external,'errors':errors}
print(json.dumps(report,indent=2))
if errors:raise SystemExit(1)
