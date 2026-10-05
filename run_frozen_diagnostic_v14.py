#!/usr/bin/env python3
"""Run one exact V14 source in a disposable public archive.

Only historical research locators are redirected. Hash functions, assertions,
source files and mathematical code are unchanged. Explicitly identified
LOW30 provenance-only reads are optional, with missing checks reported.
This file is an execution helper for run_diagnostics.py, not a sandbox.
"""
import builtins, hashlib, importlib.util, io, json, os, pathlib, sys

ROOT=pathlib.Path(__file__).resolve().parent
SCRIPT=pathlib.Path(sys.argv[1]).resolve()
MAPPING=json.loads((ROOT/'verification/V13-SOURCE-MAP.json').read_text())['source_mappings']
MAPPING.update(json.loads((ROOT/'verification/V14-SOURCE-MAP.json').read_text())['source_mappings'])
LOW30_PIN='7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8'
LOW30=os.environ.get('SAMPLING_LOW30_SOURCE')
SKIPPED=[]
ADAPTATIONS=[]
_INIT=pathlib.Path.__init__;_OPEN=builtins.open;_IO_OPEN=io.open
_SPEC=importlib.util.spec_from_file_location

def is_external(value):
    return str(value).endswith('/30_low_acc.tex') or str(value)=='external:LOW30'

def location(value):
    if not isinstance(value,(str,bytes,os.PathLike)):return value
    if isinstance(value,bytes):return value
    text=os.fspath(value)
    if is_external(text):
        return LOW30 if LOW30 else str(ROOT/'external-NOT-BUNDLED/30_low_acc.tex')
    if text in MAPPING:return str(ROOT/MAPPING[text])
    if text=='/workspace/shared':return str(ROOT/'research')
    old='/workspace/shared/recovery-20261004/'
    if text.startswith(old):
        rel=text[len(old):]
        exact=ROOT/'research/v13-source-dependencies/recovery'/rel
        return str(exact if exact.exists() else ROOT/'research'/rel)
    if text.startswith('/workspace/shared/'):
        return str(ROOT/'research'/text[len('/workspace/shared/'):])
    if text.startswith('/workspace/'):
        raise FileNotFoundError('Unmapped historical research locator: '+text)
    return value

def path_init(self,*args,**kwargs):
    if args:args=(location(args[0]),)+args[1:]
    _INIT(self,*args,**kwargs)
def mapped_open(file,*args,**kwargs):return _OPEN(location(file),*args,**kwargs)
def mapped_io_open(file,*args,**kwargs):return _IO_OPEN(location(file),*args,**kwargs)
def spec(name,location_arg=None,*args,**kwargs):
    return _SPEC(name,location(location_arg),*args,**kwargs)

def publication_hash(path):
    if is_external(path):
        if not LOW30:
            SKIPPED.append('external:LOW30');return None
        with _OPEN(LOW30,'rb') as stream:
            value=hashlib.sha256(stream.read()).hexdigest()
        if value!=LOW30_PIN:raise AssertionError('Supplied LOW30 hash mismatch')
        return value
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()

def publication_pin(path,expected):
    value=publication_hash(path)
    if value is None:return None
    if value!=expected:raise AssertionError('Source hash mismatch: '+str(path))
    return True

source=SCRIPT.read_text()
relative=str(SCRIPT.relative_to(ROOT/'research'))
# Explicit provenance-only execution adapters. Exact scripts remain on disk.
replacements={
'dyadic-prefix-law-join-20261005/independent-audit/check_log_cutoff_ledger.py': [('hashlib.sha256(Path(p).read_bytes()).hexdigest()','publication_hash(p)')],
'induction-shared-variables-20261005/independent-audit/check_decorated_independently.py': [('hashlib.sha256(path.read_bytes()).hexdigest()','publication_hash(path)')],
'law-only-variance-join-20261005/independent-audit/check_law_only_join.py': [("check(hashlib.sha256(LOW.read_bytes()).hexdigest()==PIN,'pinned LOW30 original')","external_status=publication_pin(LOW,PIN)\nif external_status is not None:check(external_status,'pinned LOW30 original')")],
'marked-shifted-local-current-20261005/check_shifted_atom.py': [("for path,digest in pins.items():ok('sealed_input_pin',hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest)","for path,digest in pins.items():\n status=publication_pin(path,digest)\n if status is not None:ok('sealed_input_pin',status)")],
'retuned-gram-native-20261005/independent-audit/check_retuned_gram.py': [('hashlib.sha256(p.read_bytes()).hexdigest()','publication_hash(p)')],
 'two-cycle-positive-return-20261005/independent-audit/check_two_cycle_audit.py': [('hashlib.sha256(p.read_bytes()).hexdigest()','publication_hash(p)')],
}
replacements.update({'combined-bridge-skew-20261005/independent-audit/check_combined_join.py': [('        check(digest(member) == expected_hash, f"pinned imported source {member}")\n        verified_imports += 1', '        status=publication_pin(member,expected_hash)\n        if status is not None:\n            check(status, f"pinned imported source {member}")\n            verified_imports += 1')], 'positive-endpoint-mean-join-20261005/independent-audit/check_endpoint_independent.py': [("    check(f'imported hash {ledger.parent.name} {n}: {Path(pin['path']).name}',hashlib.sha256(Path(pin['path']).read_bytes()).hexdigest()==pin['sha256'])", '    status=publication_pin(pin[\'path\'],pin[\'sha256\'])\n    if status is not None:check(f\'imported hash {ledger.parent.name} {n}: {Path(pin["path"]).name}\',status)')], 'positive-fourth-order-endpoint-20261005/check_endpoint.py': [(" ck(hashlib.sha256(Path(pin['path']).read_bytes()).hexdigest()==pin['sha256'],'input_pin')", " status=publication_pin(pin['path'],pin['sha256'])\n if status is not None:ck(status,'input_pin')")], 'positive-fourth-order-endpoint-20261005/independent-audit/check_fourth_endpoint_independent.py': [("    check(f'new input pin {i}: '+Path(pin['path']).name,sha(Path(pin['path']))==pin['sha256'])", "    status=publication_pin(pin['path'],pin['sha256'])\n    if status is not None:check(f'new input pin {i}: '+Path(pin['path']).name,status)")], 'public-readout-cycle-contract-20261005/independent-audit/check_grouped_contract_independently.py': [("        check(sha(f) == entry['sha256'], f'direct prerequisite pin: {f}')\n        pins.append({'path': str(f), 'sha256': sha(f)})", "        status=publication_pin(f,entry['sha256'])\n        if status is not None:check(status,f'direct prerequisite pin: {f}')\n        pins.append({'path':str(f),'sha256':publication_hash(f)})")], 'quartic-corrected-mean-join-20261005/check_join.py': [(" check(p.is_file(),f'input exists {p}')\n check(hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256'],f'pin {p}')", " status=publication_pin(p,item['sha256'])\n if status is not None:\n  check(p.is_file(),f'input exists {p}')\n  check(status,f'pin {p}')")]})
replacements.update({'fourth-cumulant-return-20261005/independent-audit/verify_pins.py': [('    got=hashlib.sha256(path.read_bytes()).hexdigest()\n    if got!=item[\'sha256\']:raise SystemExit(f"HASH MISMATCH: {path}: {got}")\n    count+=1', "    status=publication_pin(path,item['sha256'])\n    if status is not None:count+=1")], 'skew-square-positive-return-20261005/independent-audit/verify_pins.py': [("for entry in manifest['inputs']+manifest['outputs']:\n    actual=hashlib.sha256(Path(entry['path']).read_bytes()).hexdigest()\n    if actual!=entry['sha256']:raise SystemExit('HASH MISMATCH: '+entry['path'])", "verified=0\nfor entry in manifest['inputs']+manifest['outputs']:\n    status=publication_pin(entry['path'],entry['sha256'])\n    if status is not None:verified+=1"), ("'files_checked':len(manifest['inputs'])+len(manifest['outputs'])", "'files_checked':verified")]})
replacements.update({'algebra-current-image-20261005/independent-audit/independent_checks.py': [("    b=Path(x['path']).read_bytes()\n    good=len(b)==x['bytes'] and hashlib.sha256(b).hexdigest()==x['sha256']\n    check(good,'pin:'+x['path'])\n    pins.append(dict(path=x['path'],sha256=hashlib.sha256(b).hexdigest(),bytes=len(b),match=good))", "    status=publication_pin(x['path'],x['sha256'])\n    if status is not None:\n        b=Path(x['path']).read_bytes()\n        good=len(b)==x['bytes'] and hashlib.sha256(b).hexdigest()==x['sha256']\n        check(good,'pin:'+x['path'])\n        pins.append(dict(path=x['path'],sha256=hashlib.sha256(b).hexdigest(),bytes=len(b),match=good))\n    else:\n        pins.append(dict(path=x['path'],sha256=None,bytes=None,match=None,status='external_not_reverified'))")]})
for old,new in replacements.get(relative,[]):
    if source.count(old)!=1:raise AssertionError('Execution adapter source drift: '+relative)
    source=source.replace(old,new);ADAPTATIONS.append('Optional external LOW30 provenance read; no mathematical assertion changed')
pathlib.Path.__init__=path_init
builtins.open=mapped_open
io.open=mapped_io_open
importlib.util.spec_from_file_location=spec
sys.path.insert(0,str(SCRIPT.parent))
sys.argv=[str(SCRIPT)]
namespace={'__name__':'__main__','__file__':str(SCRIPT),'__package__':None,
           'publication_hash':publication_hash,'publication_pin':publication_pin}
try:
    exec(compile(source,str(SCRIPT),'exec'),namespace)
finally:
    print('\nPUBLICATION_EXECUTION '+json.dumps({'script':relative,'source_sha256':hashlib.sha256(SCRIPT.read_bytes()).hexdigest(),'locator_mapping':'exact public research copies','external_hash_reads_not_run':len(SKIPPED),'external_sources_not_reverified':sorted(set(SKIPPED)),'execution_adapters':ADAPTATIONS}),file=sys.stderr)
