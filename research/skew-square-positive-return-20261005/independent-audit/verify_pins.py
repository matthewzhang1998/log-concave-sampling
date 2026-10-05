from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parent
manifest=json.loads((root/'MANIFEST.json').read_text())
for entry in manifest['inputs']+manifest['outputs']:
    actual=hashlib.sha256(Path(entry['path']).read_bytes()).hexdigest()
    if actual!=entry['sha256']:raise SystemExit('HASH MISMATCH: '+entry['path'])
print(json.dumps({'pin_status':'PASS','files_checked':len(manifest['inputs'])+len(manifest['outputs']),'verdict':manifest['verdict'],'scope':manifest['scope']},indent=2))
