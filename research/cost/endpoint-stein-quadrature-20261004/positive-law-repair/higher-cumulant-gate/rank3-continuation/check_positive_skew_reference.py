#!/usr/bin/env python3
"""Scalar CDF diagnostic for the positive quadratic skew reference."""
import json,math
from pathlib import Path
import numpy as np
from scipy.special import roots_hermitenorm,ndtr
from scipy.interpolate import PchipInterpolator
from numpy.polynomial.legendre import leggauss
v,w=roots_hermitenorm(160);w/=math.sqrt(2*math.pi)
f=(v+.35*(np.logaddexp(v,-v)-math.log(2)))/1.35;f-=w@f
var=w@(f*f);third=w@(f**3);energy=math.sqrt(var)
x=np.linspace(-10,10,20001);u,wu=leggauss(3600);u=(u+1)/2;wu/=2

def quant(means,sd):
    cdf=np.zeros_like(x)
    for j in range(0,len(means),80):cdf+=ndtr((x[:,None]-means[None,j:j+80])/sd)@w[j:j+80]
    cdf=np.maximum.accumulate(cdf);keep=np.r_[True,np.diff(cdf)>1e-16]
    return PchipInterpolator(cdf[keep],x[keep])(u)
rows=[]
for a in [.4,.3,.2,.15,.1,.075,.05]:
    eta=1/math.sqrt(2);s=math.sqrt(eta*eta+a*a*var)
    b=a**3*third/(6*s*s)
    qr=quant(a*f,1);qc=quant(s*v+b*(v*v-1),eta)
    d=math.sqrt(wu@((qr-qc)**2))
    rows.append({'A':a,'W2':d,'W2_over_A4':d/a**4,'one_energy_ratio_W2_over_A3e':d/(a**4*energy),'reference_variance_extra':2*b*b,'reference_third_extra':8*b**3})
assert all(np.isfinite(z['W2']) for z in rows)
assert max(z['W2_over_A4'] for z in rows)<.1
out={'rows':rows,'scope':'Numerical scalar quantile diagnostic only. Positive reference includes its true extra covariance and higher cumulants; these are not silently set to zero.'}
Path(__file__).with_name('positive_skew_reference_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
