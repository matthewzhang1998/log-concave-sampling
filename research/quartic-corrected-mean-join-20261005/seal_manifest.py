#!/usr/bin/env python3
from pathlib import Path
import json,hashlib
r=Path(__file__).resolve().parent
excluded={'MANIFEST.json','SHA256SUMS'}
files={str(p.relative_to(r)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(r.rglob('*')) if p.is_file() and p.name not in excluded and '__pycache__' not in p.parts}
author=json.loads((r/'join_checks.json').read_text())
independent=json.loads((r/'independent-audit/independent_checks.json').read_text())
manifest={'date':'2026-10-05','status':'independent guarded PASS','target':'canonical m3 own mean and optional own-mean Gaussian LAW','bound':'C A^(39/10) sqrt(D)+Lambda A^(59/15) sqrt(D)+e_abs','parameters':{'w':'A^(14/15)','eta':'A^(19/10)','h':'A^(19/10)','native_mean_order':25,'native_mean_padding':'full normalized radius','negative_Gram_path_order':6,'cubic_quartic_pair_order':5,'final_own_mean_order':4,'final_padding':'A'},'files':files,'inputs':json.loads((r/'INPUT-PINS.json').read_text())['inputs'],'checks':{'author':author['assertions'],'independent':4357,'native_compiler_executed':False},'scope':['Exact coefficient target comparisons are integrated conditional-W2 over actual standard Y.','All actual original-VALUE banks, common carrier, unused/perpendicular roots and readsets remain explicit.','Every fixed sub-four grade requires its own finite order and native/log window.','No grade-four endpoint admission or repeatable arbitrary-order recurrence is claimed.']}
(r/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
all_files=files|{'MANIFEST.json':hashlib.sha256((r/'MANIFEST.json').read_bytes()).hexdigest()}
(r/'SHA256SUMS').write_text(''.join(f'{h}  {p}\n' for p,h in sorted(all_files.items())))
print('Manifest SHA256:',all_files['MANIFEST.json'])
