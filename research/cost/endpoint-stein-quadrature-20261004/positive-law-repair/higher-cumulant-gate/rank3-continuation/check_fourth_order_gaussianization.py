#!/usr/bin/env python3
"""Diagnostics, not a proof of the Hilbert-valued Riesz theorem."""
import json, math
from pathlib import Path
import numpy as np
from scipy.special import roots_hermitenorm, ndtr
from scipy.interpolate import PchipInterpolator
from numpy.polynomial.legendre import leggauss

out={"assertions":0,"scalar_laws":[],"third_hermite":[],"independent_copy_coefficients":[]}
def check(condition):
    assert bool(condition)
    out["assertions"]+=1

# Scalar globally Lipschitz asymmetric source with exact centering by quadrature.
v,w=roots_hermitenorm(160); w=w/math.sqrt(2*math.pi)
g=(v+0.35*np.logaddexp(v,-v)-0.35*math.log(2))/1.35
g-=w@g
check(abs(w@g)<1e-14)
var=w@(g*g); k3=w@(g**3); k4=w@(g**4)-3*var**2
e0=math.sqrt(var)
check(abs(k3)>0.05)
check(abs(k4)>0.005)
# Two-copy symmetry is exact on the same product quadrature.
d=(g[:,None]-g[None,:])/math.sqrt(2)
ww=w[:,None]*w[None,:]
for k in [1,3,5,7]: check(abs(np.sum(ww*d**k))<1e-11)
check(abs(np.sum(ww*d*d)-var)<2e-14)
# Odd source tanh is symmetric and Lipschitz one; e controls one energy.
s=np.tanh(v); es=math.sqrt(w@(s*s))
check(abs(w@(s**3))<1e-14)

# CDF quadrature and quantile W2 for scalar buffered sources.
# This numerical diagnostic has an absolute interpolation/tail floor.
ygrid=np.linspace(-10,10,24001)
xq,wq=leggauss(5000); uq=(xq+1)/2; wq=wq/2

def invcdf(vals):
    vals=np.maximum.accumulate(vals)
    keep=np.r_[True,np.diff(vals)>2e-16]
    return PchipInterpolator(vals[keep],ygrid[keep],extrapolate=False)(uq)

def mixture_cdf(means,weights,sd):
    ans=np.zeros_like(ygrid)
    for start in range(0,len(means),256):
        ans+=ndtr((ygrid[:,None]-means[None,start:start+256])/sd)@weights[start:start+256]
    return ans

from scipy.special import ndtri
for a in [0.4,0.3,0.2,0.15,0.1,0.075,0.05]:
    raw=invcdf(mixture_cdf(a*g,w,1.0))
    gauss=math.sqrt(1+a*a*var)*ndtri(uq)
    wr=math.sqrt(wq@((raw-gauss)**2))
    sym=invcdf(mixture_cdf(a*s,w,1.0))
    sg=math.sqrt(1+a*a*es*es)*ndtri(uq)
    ws=math.sqrt(wq@((sym-sg)**2))
    bound=math.sqrt(6)/4*a**4*es
    check(ws<bound)
    out["scalar_laws"].append({"amplitude":a,"skew_raw_W2":wr,"skew_W2_over_a3":wr/a**3,"symmetric_W2":ws,"symmetric_W2_over_a4":ws/a**4,"symmetric_bound":bound,"ratio_to_bound":ws/bound})

# Exact degree-three Hermite isometry, by independent Gaussian monomials.
# H_ijk=z_i z_j z_k - z_i delta_jk-z_j delta_ik-z_k delta_ij.
# For arbitrary symmetric T, E|T:H3|^2 = 3! ||T||HS^2.
from itertools import permutations, product
rng=np.random.default_rng(2404)
from collections import Counter

def monomial_expect(exponents):
    val=1
    for e in exponents:
        if e%2:return 0
        val*=math.prod(range(1,e,2))
    return val
for dim in [1,2,3,5]:
    T=rng.normal(size=(dim,dim,dim))
    T=sum(T.transpose(p) for p in permutations(range(3)))/6
    poly={}
    def add(exp,c):poly[tuple(exp)]=poly.get(tuple(exp),0)+c
    for i,j,k in product(range(dim),repeat=3):
        ex=[0]*dim
        for z in [i,j,k]:ex[z]+=1
        add(ex,T[i,j,k])
        for a,b,c in [(i,j,k),(j,i,k),(k,i,j)]:
            if b==c:
                ex=[0]*dim;ex[a]=1;add(ex,-T[i,j,k])
    lhs=0
    for a,ca in poly.items():
        for b,cb in poly.items():lhs+=ca*cb*monomial_expect([x+y for x,y in zip(a,b)])
    rhs=6*np.sum(T*T)
    check(abs(lhs-rhs)<1e-10*max(1,rhs))
    out["third_hermite"].append({"dimension":dim,"squared_contraction":lhs,"six_times_tensor_HS2":rhs})

for coeff in [[1/math.sqrt(2),-1/math.sqrt(2)], [1/math.sqrt(8)]*4+[-1/math.sqrt(8)]*4]:
    c=np.array(coeff)
    check(abs(np.sum(c*c)-1)<1e-14)
    check(abs(np.sum(c**3))<1e-14)
    check(np.sum(c**4)>0)
    out["independent_copy_coefficients"].append({"coefficients":coeff,"second":float(np.sum(c*c)),"third":float(np.sum(c**3)),"fourth":float(np.sum(c**4))})

out["source"]={"Lipschitz_bound":1,"energy_at_amplitude_1":e0,"third_cumulant":k3,"fourth_cumulant":k4}
out["scope"]="Scalar quadrature diagnostics and finite identities; the theorem is proved analytically in the companion note. No finite mean/covariance circuit is instantiated."
Path(__file__).with_name('fourth_order_gaussianization_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
