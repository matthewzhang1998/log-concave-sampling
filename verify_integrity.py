#!/usr/bin/env python3
"""Verify research payload bytes, provenance mappings, and bundle checksums."""
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parent
IGNORED_PARTS = {'.git', '.venv', '__pycache__', '.pytest_cache'}
errors = []

def included(path):
    relative = path.relative_to(ROOT)
    return (path.is_file() and not any(part in IGNORED_PARTS for part in relative.parts)
            and path.name != '.DS_Store' and path.suffix not in {'.pyc', '.pyo'})

def check(relative, wanted, size=None):
    path = ROOT / relative
    if Path(relative).is_absolute() or '..' in Path(relative).parts:
        errors.append('Invalid path: '+relative)
        return
    if not path.is_file():
        errors.append('Missing: '+relative)
        return
    if hashlib.sha256(path.read_bytes()).hexdigest() != wanted:
        errors.append('Hash mismatch: '+relative)
    if size is not None and path.stat().st_size != size:
        errors.append('Size mismatch: '+relative)

inventory = json.loads((ROOT / 'INVENTORY.json').read_text())
entries = {entry['path']: entry for entry in inventory['files']}
if len(entries) != len(inventory['files']):
    errors.append('Duplicate inventory path')
for entry in inventory['files']:
    check(entry['path'], entry['sha256'], entry['bytes'])
    changed = entry['sha256'] != entry['original_sha256']
    if changed != (entry['provenance'] == 'sanitized_publication_copy'):
        errors.append('Inconsistent public/original mapping: '+entry['path'])
actual_payload = {str(p.relative_to(ROOT)) for p in (ROOT / 'research').rglob('*') if included(p)}
if actual_payload != set(entries):
    errors.append('Research file-set mismatch')
if inventory['payload_files'] != len(entries) or inventory['payload_bytes'] != sum(e['bytes'] for e in entries.values()):
    errors.append('Inventory totals mismatch')
provenance = json.loads((ROOT / 'PROVENANCE-CHECKS.json').read_text())
for binding in provenance['bindings']:
    entry = entries[binding['path']]
    if binding['public_sha256'] != entry['sha256'] or binding['original_sha256'] != entry['original_sha256']:
        errors.append('Provenance mapping mismatch: '+binding['path'])
listed = []
for line in (ROOT / 'SHA256SUMS').read_text().splitlines():
    wanted, relative = line.split('  ', 1)
    listed.append(relative)
    check(relative, wanted)
if len(listed) != len(set(listed)):
    errors.append('Duplicate checksum path')
actual = {str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if included(p) and p.name != 'SHA256SUMS'}
if actual != set(listed):
    errors.append('Bundle file-set mismatch: '+str(sorted(actual.symmetric_difference(listed))))
if errors:
    print('\n'.join(errors))
    raise SystemExit(1)
print(f"PASS: {len(entries)} research files, {len(provenance['bindings'])} provenance mappings, "
      f"{len(listed)} bundle checksums, {sum(p.stat().st_size for p in ROOT.rglob('*') if included(p))} total bytes.")
