#!/usr/bin/env python3
"""Verify the frozen matrix-free packet; optionally export public provenance.

Publication adaptation of the historical author regeneration utility. Does not
rewrite sealed manifests, historical results, source pins, or archive files.
"""
from pathlib import Path
import argparse,hashlib,json
HERE=Path(__file__).resolve().parent
ROOT=next(p for p in HERE.parents if (p/'INVENTORY.json').is_file())
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.output and (args.output.resolve()==ROOT or ROOT in args.output.resolve().parents):
        parser.error('--output must resolve outside the immutable archive directory')
    entries={e['path']:e for e in json.loads((ROOT/'INVENTORY.json').read_text())['files']}
    manifest=json.loads((HERE/'MANIFEST.json').read_text())
    files=[]
    for item in manifest['inventory']+[{'path':'MANIFEST.json'},{'path':'SHA256SUMS'}]:
        path=(HERE/item['path']).resolve()
        if HERE not in path.parents: raise AssertionError('File outside frozen package')
        entry=entries[str(path.relative_to(ROOT))]
        assert hashlib.sha256(path.read_bytes()).hexdigest()==entry['sha256']
        assert path.stat().st_size==entry['bytes']
        if 'sha256' in item: assert item['sha256']==entry['original_sha256']
        files.append(entry)
    actual={str(p.relative_to(HERE)) for p in HERE.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix not in {'.pyc','.pyo'}}
    assert actual=={i['path'] for i in manifest['inventory']}|{'MANIFEST.json','SHA256SUMS'}
    incoming=manifest['incoming_manifest']
    dep=(HERE/incoming['path']).resolve();de=entries[str(dep.relative_to(ROOT))]
    assert de['original_sha256']==incoming['sha256']
    assert hashlib.sha256(dep.read_bytes()).hexdigest()==de['sha256']
    for item in json.loads(dep.read_text())['inventory']:
        path=(dep.parent/item['path']).resolve();entry=entries[str(path.relative_to(ROOT))]
        assert entry['original_sha256']==item['sha256']
        assert hashlib.sha256(path.read_bytes()).hexdigest()==entry['sha256']
    report={'status':'PASS','files':files,'incoming_original_manifest_sha256':incoming['sha256'],'scope':manifest['scope'],'historical_source_records':'Preserved unchanged; historical README source record maps to the exact separately bundled pre-final README.'}
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'PASS','packet_files':len(files),'incoming_inventory_pins_verified':len(json.loads(dep.read_text())['inventory'])}))
if __name__=='__main__':main()
