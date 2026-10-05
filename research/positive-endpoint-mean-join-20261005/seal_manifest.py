from pathlib import Path
import hashlib,json
P=Path(__file__).resolve().parent
report=P/'independent-audit/INDEPENDENT-ENDPOINT-JOIN-AUDIT.md'
assert report.is_file(), 'Independent report must exist before sealing'
checks=json.loads((P/'independent-audit/independent_checks.json').read_text())
assert checks['status']=='PASS'
main=P/'MEAN-COVARIANCE-SKEW-ENDPOINT.md'
current=hashlib.sha256(main.read_bytes()).hexdigest()
assert any(x['path']==str(main) and x['sha256']==current for x in checks['sources']), 'Independent reviewed main pin changed'
files={str(f.relative_to(P)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(P.rglob('*')) if f.is_file() and f.name not in {'MANIFEST.json','SHA256SUMS'} and '__pycache__' not in f.parts}
manifest={
 'date':'2026-10-05','status':'independent guarded PASS',
 'theorem':'positive finite canonical-mean/full-covariance/negative-skew reverse-OU endpoint',
 'standardized_bound':'C A^(39/10) sqrt(D)+Lambda A^(59/15) sqrt(D)+absolute floors',
 'log_absorption_window':'Lambda A^(1/30)<=1',
 'physical_leading_grade':'22/5 under inherited external sqrt(A) readout',
 'parameters':{'bridge_h':'1/2','covariance_unscaled_buffer':'1/2','cubic_buffer':'1/8','untouched_keep':'1/4','rho':'3/4','terminal_squared_clock_exponent':'19/25','mean_native_order':25},
 'outer_composition':'explicit local reverse-OU bridge, exact reference contraction, geometric sum, finite-mode moment closure, separately priced positive Stein terminal',
 'native_compilers_numerically_executed':False,
 'author_assertions':json.loads((P/'endpoint_join_checks.json').read_text())['assertions'],
 'independent_assertions':checks['assertions'],
 'limitations':['literal imported native numerical/caller/dimension/root guards required','no grade four theorem','no terminal LAW to RAW conversion','no arbitrary-order or order-uniform cost theorem'],
 'files':files,'inputs':json.loads((P/'INPUT-PINS.json').read_text())['inputs']}
(P/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
entries=files|{'MANIFEST.json':hashlib.sha256((P/'MANIFEST.json').read_bytes()).hexdigest()}
(P/'SHA256SUMS').write_text(''.join(f'{v}  {k}\n' for k,v in sorted(entries.items())))
print(entries['MANIFEST.json'])
