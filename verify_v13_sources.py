#!/usr/bin/env python3
"""Verify exact V13 copies, source mappings, seals, and direct input pins.

Historical mismatches are explicit exceptions, never successful hash checks.
The full external LOW30/HIGH13 manuscripts are not read by this validator.
"""
from pathlib import Path
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parent
selection=json.loads((ROOT/'verification/V13-FROZEN-SELECTION.json').read_text())
source_map=json.loads((ROOT/'verification/V13-SOURCE-MAP.json').read_text())['source_mappings']
historical=json.loads((ROOT/'verification/V13-HISTORICAL-PINS.json').read_text())
files={e['path']:e for e in selection['files']}
errors=[];seals=0;sealed_members=0;inputs=0;external=[];resolved=[];exceptions=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def source_path(name,base):
    name=str(name)
    if name.endswith(('30_low_acc.tex','13_high_acc.tex')) or name.startswith('external:'):return None
    if name in source_map:return ROOT/source_map[name]
    if name.startswith('/workspace/shared/'):
        return ROOT/'research'/name[len('/workspace/shared/'):]
    p=Path(name)
    return p if p.is_absolute() else base/p
for rel,e in files.items():
    p=ROOT/rel
    if not p.is_file() or sha(p)!=e['sha256'] or p.stat().st_size!=e['bytes'] or e['original_sha256']!=e['sha256']:
        errors.append('Frozen source mismatch: '+rel)
for source,rel in source_map.items():
    if rel not in files:errors.append('Unselected mapping: '+rel)
    if '/agent_notes/' in source or '/user_notes/' in source:errors.append('Invalid source locator')
known_path='research/positive-endpoint-mean-join-20261005/independent-audit/INDEPENDENT-ENDPOINT-JOIN-AUDIT.md'
known_old='cc4f5bc8d50672a8da492cef001f498f0c9275c837261d59734470a7354f89ef'
known_actual='962eb398c771a25ee4fa1ddeb4a30123fef46a1ea71ab0269031c8c4664c8b10'
for rel in files:
    if Path(rel).name!='SHA256SUMS':continue
    p=ROOT/rel;seals+=1
    for n,line in enumerate(p.read_text().splitlines(),1):
        if not line.strip():continue
        match=re.fullmatch(r'([a-f0-9]{64})\s+\*?(.+)',line)
        if not match:errors.append('Unrecognized seal line: '+rel+':'+str(n));continue
        expected,name=match.groups();q=source_path(name,p.parent)
        if q is None:external.append({'record':rel,'source':name,'expected':expected});continue
        if not q.is_file():errors.append('Missing sealed member: '+str(q));continue
        actual=sha(q);relative=str(q.relative_to(ROOT))
        if expected!=actual:
            if (relative,expected,actual)==(known_path,known_old,known_actual):
                exceptions.append({'record':rel,'path':relative,'recorded':expected,'actual':actual,'status':'documented_historical_mismatch'})
            else:errors.append('Seal mismatch: '+relative)
        else:sealed_members+=1

def pin_records(x):
    if isinstance(x,dict):
        if isinstance(x.get('path'),str) and isinstance(x.get('sha256'),str) and re.fullmatch('[0-9a-f]{64}',x['sha256']):
            yield x['path'],x['sha256']
        for key,value in x.items():
            if isinstance(value,str) and re.fullmatch('[0-9a-f]{64}',value) and ('/' in key or '.' in key):yield key,value
            elif isinstance(value,(list,dict)):yield from pin_records(value)
    elif isinstance(x,list):
        for item in x:yield from pin_records(item)
for rel in files:
    if not re.fullmatch(r'(?:.*-)?(?:INPUT|SOURCE)-PINS\.json',Path(rel).name):continue
    p=ROOT/rel
    for name,expected in pin_records(json.loads(p.read_text())):
        q=source_path(name,p.parent)
        if q is None:external.append({'record':rel,'source':name,'expected':expected});continue
        if q.is_file() and sha(q)==expected:inputs+=1;continue
        match=next((h for h in historical['historical_pin_resolutions'] if h['from']==rel and h['sha256']==expected and h.get('resolved_path')),None)
        if match:
            reviewed=ROOT/match['resolved_path']
            if reviewed.is_file() and sha(reviewed)==expected:
                resolved.append({'record':rel,'source':name,'exact_reviewed_snapshot':match['resolved_path'],'sha256':expected});continue
        errors.append('Direct input pin mismatch or missing: '+rel+' -> '+name)
if len(exceptions)!=1:errors.append('Expected exactly one explicitly documented checksum-seal discrepancy')
result={'status':'PASS_WITH_EXPLICIT_HISTORICAL_AND_EXTERNAL_QUALIFICATIONS' if not errors else 'FAIL','exact_new_files_verified':len(files),'source_mappings_verified':len(source_map),'checksum_lists_verified':seals,'matching_checksum_members':sealed_members,'direct_input_pins_verified':inputs,'exact_reviewed_snapshot_resolutions':resolved,'historical_seal_exceptions':exceptions,'external_pins_not_reverified':external,'other_historical_hash_references':'See V13-HISTORICAL-PINS.json; earlier audit output hashes are not silently recertified against current sources.','errors':errors}
print(json.dumps(result,indent=2))
if errors:raise SystemExit(1)
