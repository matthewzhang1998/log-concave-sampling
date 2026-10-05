#!/usr/bin/env python3
"""Scalar author diagnostics only; the same-source separator is an analytical proof."""
import math, json, hashlib
from pathlib import Path
from scipy.integrate import quad
here=Path(__file__).parent
name='RADIAL-SHELL-REFUTES-UNSHIFTED-SINGLE-HISTORY-COVARIANCE.md'
pin=hashlib.sha256((here/name).read_bytes()).hexdigest()
assert pin=='2807f7307646cc16c79049319b574259874b09e731a13f686696bdd60a05b035'
eta=lambda s:math.exp(-1/(1-s*s/4)) if abs(s)<2 else 0.
z=quad(eta,-2,2,epsabs=1e-14)[0]
b=lambda s:eta(s)/z
phi=lambda s:math.exp(-s*s)/math.sqrt(math.pi)
c=d=.25
shift=c/2
h_mean=quad(lambda s:b(s)*(phi(s+shift)-phi(s)),-2,2,epsabs=1e-14)[0]
assert z>=2*math.exp(-4/3)
assert b(0)<.5*math.exp(1/3)<1
assert c+d<=1
assert c+d/2<=1
assert 2/3-1/16>=7/12
assert (7/12)**2-.25>=13/144-1e-15
assert h_mean<0
k_min=d*c*c/2*13/144
assert k_min>0
out={'status':'PASS','assertions':9,'scope':'Author scalar constants and witness diagnostics only; independent exact-text proof review pending.',
     'source':name,'source_sha256':pin,'bump_normalization':z,'bump_max':b(0),'radial_shift':shift,
     'limiting_shell_h_mean':h_mean,'uniform_min_abs_k':k_min,
     'asymptotic_min_abs_Ax_witness_per_A2':.5*k_min*abs(h_mean)}
(here/'radial_shell_separator_constants.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
