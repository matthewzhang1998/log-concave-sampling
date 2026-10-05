from pathlib import Path
import hashlib,json
P=Path(__file__).resolve().parent
checks=json.loads((P/'independent-audit/independent_checks.json').read_text())
assert checks['status']=='PASS'
assert not checks['native_compilers_executed']
main=P/'FOURTH-ORDER-POSITIVE-ENDPOINT.md'
# Independent source pins may be stored in the check or separate reviewed-source file.
review=json.loads((P/'independent-audit/reviewed-source.json').read_text())
source_pins=checks.get('sources',[])
if isinstance(review,list):source_pins+=review
elif isinstance(review,dict):source_pins+=review.get('sources',review.get('files',[]))
if isinstance(source_pins,dict):source_pins=[{'path':k,'sha256':v} for k,v in source_pins.items()]
for name in ['FOURTH-ORDER-POSITIVE-ENDPOINT.md','ENDPOINT-GRAPH.json','INPUT-PINS.json','PROVENANCE-NOTE.md']:
 f=P/name;sha=hashlib.sha256(f.read_bytes()).hexdigest()
 assert any(x.get('path') in [str(f),name] and x.get('sha256')==sha for x in source_pins), 'Exact independent review pin missing: '+name
for pin in json.loads((P/'INPUT-PINS.json').read_text())['inputs']:
 assert hashlib.sha256(Path(pin['path']).read_bytes()).hexdigest()==pin['sha256']
author=json.loads((P/'endpoint_checks.json').read_text());assert author['status']=='PASS'
files={str(f.relative_to(P)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(P.rglob('*')) if f.is_file() and f.name not in {'MANIFEST.json','SHA256SUMS'} and '__pycache__' not in f.parts}
manifest={
 'date':'2026-10-05','status':'independent guarded PASS',
 'theorem':'finite positive general-C2 canonical-m3/full-Cov(F3)/negative-skew reverse-OU log-grade-four endpoint',
 'normalized_bound':'Lambda_ep A^4 sqrt(D)+absolute floors',
 'physical_bound':'Lambda_ep A^(9/2) sqrt(D)+restored physical mode/numerical floors, under inherited sqrt(A) readout',
 'leading_public_logs_retained':True,
 'parameters':{'bridge_h':'1/2','covariance_unscaled_buffer':'1/2','cubic_buffer':'1/8','untouched_keep':'1/4','rho':'3/4','terminal_squared_clock_exponent':'4/5','buffered_alpha_range':'A^(9/5)<alpha_j<=A','prefix_mean_order':4,'tail_mean_order':6,'final_own_mean_order':4,'log_cutoff':'eta_j=alpha_j^2 K_log^2; K_log=C L^b fixed'},
 'outer_composition':'explicit normalized local bridge, exact reference contraction and geometric sum, finite-mode moment closure, positive terminal with delta_terminal<=alpha_T',
 'native_compilers_numerically_executed':False,
 'author_assertions':author['assertions'],'independent_assertions':checks['assertions'],
 'provenance_exceptions':author['provenance_exceptions'],
 'limitations':['every literal native/source/root/caller/gap/precision guard remains required','no numerical-constant A^4 bound','no endpoint grade above four','no terminal LAW to RAW conversion','no arbitrary-order or order-uniform cost theorem','no production arbitrary-precision scalar Gauss implementation claimed'],
 'files':files,'inputs':json.loads((P/'INPUT-PINS.json').read_text())['inputs']}
(P/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
entries=files|{'MANIFEST.json':hashlib.sha256((P/'MANIFEST.json').read_bytes()).hexdigest()}
(P/'SHA256SUMS').write_text(''.join(f'{v}  {k}\n' for k,v in sorted(entries.items())))
print(entries['MANIFEST.json'])
