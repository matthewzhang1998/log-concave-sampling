#!/usr/bin/env python3
"""Exercise a full production scalar preparation and the pinned VALUE callback.
This is not a native sampler execution or a certified floating Sobolev oracle.
"""
from fractions import Fraction as F
import importlib.util,json,pathlib,hashlib
import numpy as np
from cutoff_free_sector_gauss import prepare_certified_rule
from gauss_rational_core import map_rule

root=pathlib.Path(__file__).parent
p,b,panels=prepare_certified_rule(9,7,1,F(1,2),F(1,2))
assert len(panels)==p['interval_count_per_axis']
assert len(b['x'])==p['m']
assert b['node_error']<=p['node_tolerance']
assert b['weight_l1_error']<=p['weight_tolerance']
assert sum(b['w'])==2 and sum(x*w for x,w in zip(b['x'],b['w']))==0
for cell in (panels[0],panels[-1]):
    rr=map_rule(b,*cell)
    assert all(cell[0]<x<cell[1] and w>0 for x,w in rr)

source=pathlib.Path('/workspace/shared/induction-native-rank-generator-20261005/active_probe_source.py')
spec=importlib.util.spec_from_file_location('pinned_active',source)
active=importlib.util.module_from_spec(spec);spec.loader.exec_module(active)
rng=np.random.default_rng(72031);records=[]
for k in (0,1,2,3,4):
    seen=[]
    def g(x):
        seen.append(np.array(x,copy=True))
        return np.asarray(x)+np.tanh(x)
    probes=rng.normal(size=(k,2));x=rng.normal(size=2);z=rng.normal(size=2)
    value,calls=active.active_source(g,x,probes,z,.125,2.)
    assert calls==len(seen)==2**(k+1) and np.all(np.isfinite(value))
    seen.clear();zero,calls0=active.active_source(g,np.zeros(2),probes,z,.125,2.)
    assert np.array_equal(zero,np.zeros(2)) and calls0==len(seen)==2**(k+1)
    records.append({'order':k,'raw_VALUE_calls':calls,'exact_origin_zero':True})
out={'status':'PASS','production_scalar_parameters':{'J':p['J'],'m':p['m'],'intervals':len(panels),'node_count_per_tree':p['node_count_per_tree'],'certified_bits':b['certification_bits']},'active_source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'VALUE_callback_checks':records,'scope':'Full rational scalar preparation; raw diagnostic original-VALUE callback only; no native sampler or certified Sobolev oracle executed.'}
(root/'production-callback-checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
