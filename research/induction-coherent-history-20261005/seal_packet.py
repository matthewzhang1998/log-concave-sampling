import hashlib,json
from pathlib import Path
root=Path(__file__).parent
sources=[
'/workspace/shared/history-conditional-resummation-20261005/TRUE-PAIR-HISTORY-RESUMMATION.md',
'/workspace/shared/next-history-m4-port-20261005/NEXT-HISTORY-M4-AND-MARKED-CURRENT-PORT.md',
'/workspace/shared/dyadic-pair-secant-port-20261005/COHERENT-BLOCK-ORDER-TWO-EXTENSION.md',
'/workspace/shared/dyadic-pair-secant-port-20261005/FINITE-ORDER-TWO-C2-TAIL.md',
'/workspace/shared/dyadic-marked-history-obstruction-20261005/DYADIC-MARKED-HISTORY-NONCLOSURE.md',
'/workspace/shared/dyadic-prefix-law-join-20261005/ACTUAL-DYADIC-PREFIX-LAW-JOIN.md',
'/workspace/shared/positive-fourth-order-endpoint-20261005/FOURTH-ORDER-POSITIVE-ENDPOINT.md',
'/workspace/shared/marked-shifted-local-current-20261005/BOUNDED-SHIFTED-CURRENT.md',
'/workspace/shared/retuned-gram-native-20261005/RETUNED-NATIVE-AND-CONDITIONAL-GRAM-CONTRACT.md',
]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
(root/'INPUT-PINS.json').write_text(json.dumps([{'path':p,'sha256':sha(p)} for p in sources],indent=2)+'\n')
files=sorted(p for p in root.rglob('*') if p.is_file() and p.name not in ['MANIFEST.json','SHA256SUMS'] and '__pycache__' not in str(p))
manifest={'date':'2026-10-05','status':'constructive_partial_local_current_port','completed_full_history_induction':False,'new_results':['retuned fourteen-VALUE pair source first','nineteen-VALUE coherent marked current correction','retained-carrier clipped native current port','raw small-energy near-gradient marker invariant'],'boundaries':['true-history base/singleton fields and Delta3 source unsupplied','no unrestricted-D history comparison','no nonlinear marker-reuse or unweighted-law inference','native compiler imported with literal guards, not executed'],'files':[{'path':str(p.relative_to(root)),'sha256':sha(p)} for p in files]}
(root/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
allfiles=files+[root/'MANIFEST.json']
(root/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.relative_to(root)}\n' for p in allfiles))
print(json.dumps({'files':len(allfiles),'manifest_sha256':sha(root/'MANIFEST.json')}))
