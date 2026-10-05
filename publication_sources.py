#!/usr/bin/env python3
"""Verify original source pins against archived public bytes and explicit mappings.

Sanitized files are never represented as original bytes. Missing LOW30 is
reported separately and is never counted as a verified source assertion.
"""
from pathlib import Path
import argparse, hashlib, json
ROOT=Path(__file__).resolve().parent
ENTRIES={e['path']:e for e in json.loads((ROOT/'INVENTORY.json').read_text())['files']}
LOW30_PIN='7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8'
parser=argparse.ArgumentParser(add_help=False)
parser.add_argument('--low30-source',type=Path)
ARGS,_=parser.parse_known_args()
SKIPPED_EXTERNAL=[]
LOW30_VERIFIED=False
def source_path(name,base=None):
    value=str(name)
    if value=='external:LOW30':return Path(value)
    if value.startswith('research/'):return ROOT/value
    p=Path(name)
    return p if p.is_absolute() else (Path(base)/p if base is not None else ROOT/p)
def source_label(path):
    path=source_path(path)
    if str(path)=='external:LOW30':return str(path)
    return str(path.resolve().relative_to(ROOT))
def verify_source_pin(path,pin):
    global LOW30_VERIFIED
    path=source_path(path)
    if str(path)=='external:LOW30':
        if pin!=LOW30_PIN:raise AssertionError('Unexpected external source pin')
        if ARGS.low30_source is None:
            SKIPPED_EXTERNAL.append('external:LOW30')
            return None
        if hashlib.sha256(ARGS.low30_source.read_bytes()).hexdigest()!=pin:
            raise AssertionError('LOW30 source pin mismatch')
        LOW30_VERIFIED=True
        return True
    relative=source_label(path)
    entry=ENTRIES[relative]
    actual=hashlib.sha256(path.read_bytes()).hexdigest()
    if actual!=entry['sha256'] or path.stat().st_size!=entry['bytes']:
        raise AssertionError('Public byte integrity mismatch: '+relative)
    if pin!=entry['original_sha256']:
        raise AssertionError('Original/public pin mapping mismatch: '+relative)
    return True
def pin_report():
    return {'external_hash_assertions_not_run':len(SKIPPED_EXTERNAL),
            'LOW30_status':'verified_original_bytes' if LOW30_VERIFIED else 'external_not_reverified',
            'verification':'Bundled source assertions verify public bytes and their explicit original-pin mappings; mathematical assertions are unchanged.'}
