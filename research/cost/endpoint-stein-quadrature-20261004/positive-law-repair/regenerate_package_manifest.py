#!/usr/bin/env python3
"""Validate and export the frozen publication manifest without editing the archive.

Publication adaptation of the original regeneration utility. Historical source
pins remain in MANIFEST.json; current public hashes and original/public mappings
come from the root INVENTORY.json. LOW30 remains external and is not reverified.
The exact four-component inclusion list is retained in MANIFEST.json.
"""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ARCHIVE = ROOT.parents[3]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='Optional destination outside the archive; by default print summary only.')
    args = parser.parse_args()
    if args.output:
        target = args.output.resolve()
        if target == ARCHIVE or ARCHIVE in target.parents:
            parser.error('--output must resolve outside the immutable archive directory')
    inventory = json.loads((ARCHIVE/'INVENTORY.json').read_text())
    entries = {e['path']:e for e in inventory['files']}
    frozen = json.loads((ROOT/'MANIFEST.json').read_text())
    files = []
    for item in [*frozen['files'], {'path':'MANIFEST.json'}]:
        path = ROOT/item['path']
        if path.parent != ROOT and ROOT not in path.resolve().parents:
            raise RuntimeError('File outside frozen package')
        relative = str(path.relative_to(ARCHIVE))
        entry = entries[relative]
        if sha(path) != entry['sha256'] or path.stat().st_size != entry['bytes']:
            raise RuntimeError('Publication integrity mismatch: '+relative)
        if 'sha256' in item and item['sha256'] != entry['original_sha256']:
            raise RuntimeError('Original source mapping mismatch: '+relative)
        files.append({k:entry[k] for k in ('path','sha256','bytes','original_sha256','original_bytes','provenance')})
    actual = {str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix not in ('.pyc','.pyo')}
    expected = {e['path'] for e in frozen['files']} | {'MANIFEST.json'}
    # V10 adds exactly two separately frozen packages beside the V9 selection.
    for package in ('higher-cumulant-gate', 'reverse-ou-curvature'):
        manifest_path = ROOT/package/'MANIFEST.json'
        manifest = json.loads(manifest_path.read_text())
        items = manifest['files']
        if isinstance(items, dict):
            items = [dict(path=k, **v) for k, v in items.items()]
        for item in [*items, {'path':'MANIFEST.json'}]:
            local = Path(package)/item['path']
            path = ROOT/local
            if ROOT not in path.resolve().parents:
                raise RuntimeError('File outside frozen package')
            entry = entries[str(path.relative_to(ARCHIVE))]
            if sha(path) != entry['sha256'] or path.stat().st_size != entry['bytes']:
                raise RuntimeError('V10 publication integrity mismatch: '+str(local))
            if 'sha256' in item and item['sha256'] != entry['original_sha256']:
                raise RuntimeError('V10 original source mapping mismatch: '+str(local))
            expected.add(str(local))
    if actual != expected:
        raise RuntimeError('Frozen package file-set mismatch')
    dependencies = []
    for dep in frozen['external_prerequisites']:
        path = 'research/'+dep['path']
        if path in entries:
            entry = entries[path]
            if entry['original_sha256'] != dep['sha256'] or sha(ARCHIVE/path) != entry['sha256']:
                raise RuntimeError('Prerequisite mapping mismatch: '+path)
            dependencies.append({'path':path,'original_sha256':dep['sha256'],'public_sha256':entry['sha256'],'status':'public_copy_verified'})
        else:
            dependencies.append({'path':'external:LOW30','original_sha256':dep['sha256'],'status':'external_not_reverified'})
    report = {'package':frozen['package'],'components':len(frozen['components']),'files':files,'dependencies':dependencies,
        'historical_independent_assertions':4901,'public_independent_assertions':4900,
        'scope':'Four bounded components only; no general nonlinear all-order or sublinear-cost theorem.'}
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'PASS','files':len(files),'components':len(frozen['components']),'v10_frozen_additional_files':63,'external_prerequisites_not_reverified':1},indent=2))

if __name__ == '__main__':
    main()
