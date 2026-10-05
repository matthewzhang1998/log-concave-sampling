from pathlib import Path
import hashlib,json
p=Path(__file__).resolve().parent
H=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
sources=[
'/workspace/shared/law-only-variance-join-20261005/LAW-ONLY-SHRINKING-BUFFER-JOIN.md',
'/workspace/shared/combined-bridge-skew-20261005/COMBINED-BRIDGE-SKEW-JOIN.md',
'/workspace/shared/recovery-20261004/cost/endpoint-stein-quadrature-20261004/positive-law-repair/higher-cumulant-gate/order-reentry/RECTANGULAR-FIRST-COEFFICIENT-MEAN-AND-GRAM-RETURN.md',
'/workspace/shared/recovery-20261004/cost/endpoint-stein-quadrature-20261004/positive-law-repair/higher-cumulant-gate/same-carrier-feedback/GAUSSIAN-COVARIANCE-MIXTURE-ONE-ENERGY-LAW.md',
'/workspace/shared/v9-curation-work/frozen/prerequisites/research-source/High Acc Ideas/ai-bucket/30_low_acc.tex',
'/workspace/shared/fourth-cumulant-return-20261005/POSITIVE-SMOOTHED-FOURTH-CUMULANT-PACKET.md',
'/workspace/shared/fourth-cumulant-return-20261005/independent-audit/AUDIT.md',
'/workspace/shared/fourth-cumulant-return-20261005/MANIFEST.json',
]
files=[f for f in p.rglob('*') if f.is_file() and '__pycache__' not in f.parts and f.name not in ('MANIFEST.json','SHA256SUMS')]
manifest={'date':'2026-10-05','status':'Independent PASS under the stated source-qualified native/path/variance guards','main':'POSITIVE-CORRECTED-GRAM-MIXTURE.md','main_sha256':H(p/'POSITIVE-CORRECTED-GRAM-MIXTURE.md'),'result':'Positive unheated corrected Gram integrated-Y W2 <= Lambda sqrt(D)[A^4/u+A^6/u^(5/2)]+floors; parent terminal debt A^5/v+A^7/v^(5/2). Native order seven for u>=c A^(3/2), generalized fixed order as stated.','nonclaims':['No zero-buffer endpoint return','No all-rank finite compiler','No improved full outer rate without its other separately proved components'],'checks':{'native_independent':839,'correction_independent':516,'author_finite_and_stress':516},'source_pins':{f:H(Path(f)) for f in sources},'artifacts':{str(f.relative_to(p)):H(f) for f in sorted(files)}}
(p/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
files.append(p/'MANIFEST.json')
(p/'SHA256SUMS').write_text(''.join(f'{H(f)}  {f.relative_to(p)}\n' for f in sorted(files)))
print('MANIFEST',H(p/'MANIFEST.json'))
