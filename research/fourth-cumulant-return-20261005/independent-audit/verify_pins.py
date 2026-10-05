#!/usr/bin/env python3
"""Verify the reviewed input bytes and independent deliverables, without writes."""
from pathlib import Path
import hashlib,json
p=Path(__file__).resolve().parent
m=json.loads((p/'MANIFEST.json').read_text())
count=0
for item in m['inputs']+m['outputs']:
    path=Path(item['path'])
    got=hashlib.sha256(path.read_bytes()).hexdigest()
    if got!=item['sha256']:raise SystemExit(f"HASH MISMATCH: {path}: {got}")
    count+=1
print(json.dumps({'pin_status':'PASS','files_checked':count,'mathematical_verdict':m['mathematical_verdict'],'scope':m['scope']},indent=2))
