#!/usr/bin/env python3
"""Algebra/finite Gaussian quadrature checks, not a theorem or a path algorithm."""
import json
import numpy as np
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss
from pathlib import Path

TAU=np.array([0.25,0.75]); V=np.array([0.5,0.5]); DD=np.sqrt(1-TAU**2)
KM=np.minimum.outer(TAU,TAU)/np.maximum.outer(TAU,TAU)-np.outer(TAU,TAU)
KC=np.outer(DD,DD); DK=KM-KC; aa=.08; eps=.5

def g(x): return aa*(x+eps*np.logaddexp(x,-x)-eps*np.log(2.))
def gp(x): return aa*(1+eps*np.tanh(x))
def gpp(x): return aa*eps/np.cosh(x)**2

def normal_rule(n):
    x,w=hermgauss(n); return np.sqrt(2)*x,w/np.sqrt(np.pi)
def gauss_rule(n,lo,hi):
    x,w=leggauss(n); return lo+(hi-lo)*(x+1)/2,w*(hi-lo)/2

z,w=normal_rule(48); zz=np.array(np.meshgrid(z,z,indexing='ij')).reshape(2,-1).T
ww=np.outer(w,w).ravel(); lc=np.linalg.cholesky(KM)
lam,wl=gauss_rule(50,0,1)
rows=[]
for x in [-2.,-.5,0.,1.,2.]:
    yM=TAU*x+zz@lc.T; im=g(yM)@V
    yC=TAU*x+z[:,None]*DD; ic=g(yC)@V
    jm=ww@g(x-im); jc=w@g(x-ic)
    muM=ww@im; muC=w@ic
    cm=ww@(im-muM)**2; cc=w@(ic-muC)**2
    current=0.; cov_current=0.
    for l,weight in zip(lam,wl):
        L=np.linalg.cholesky(l*KM+(1-l)*KC)
        y=TAU*x+zz@L.T; ig=g(y)@V; j=gp(y)
        coeff=np.einsum('bi,ij,bj->b',j*V,DK,j*V)
        current+=weight*.5*(ww@(gpp(x-ig)*coeff))
        cov_current+=weight*(ww@coeff)
    rows.append(dict(x=x,mean_gap=float(muM-muC),force_gap=float(jm-jc),
                     current=float(current),force_identity_error=float(jm-jc-current),
                     covariance_gap=float(cm-cc),covariance_current=float(cov_current),
                     covariance_identity_error=float(cm-cc-cov_current)))

# The centered convolution is a homotopy of two independent ACTUAL packet laws.
z3,w3=normal_rule(24)
zzz=np.array(np.meshgrid(z3,z3,z3,indexing='ij')).reshape(3,-1).T
www=np.einsum('i,j,k->ijk',w3,w3,w3).ravel()
alph,wa=gauss_rule(60,0,np.pi/2)
homotopy=[]
for x in [-1.,0.,1.]:
    im=g(TAU*x+zzz[:,:2]@lc.T)@V
    ic=g(TAU*x+zzz[:,2,None]*DD)@V
    mum=www@im; muc=www@ic; xim=im-mum; xic=ic-muc
    gap=www@(g(x-im)-g(x-ic)); integral=0.
    for alpha,weight in zip(alph,wa):
        sn=np.sin(alpha); cs=np.cos(alpha)
        mu=sn*sn*mum+cs*cs*muc
        f=mu+sn*xim+cs*xic
        deriv=2*sn*cs*(mum-muc)+cs*xim-sn*xic
        integral-=weight*(www@(gp(x-f)*deriv))
    homotopy.append(dict(x=x,force_gap=float(gap),derivative_free_current=float(integral),
                        identity_error=float(gap-integral)))

assert np.max(np.abs(np.diag(DK)))<1e-14
assert DK[0,1]<0 and DK[1,0]<0
assert abs(V@TAU-.5)<1e-14
for row in rows:
    assert abs(row['mean_gap'])<2e-9, row
    assert abs(row['force_identity_error'])<2e-9, row
    assert abs(row['covariance_identity_error'])<2e-9,row
for row in homotopy: assert abs(row['identity_error'])<2e-10,row
out=dict(scope='Finite scalar smooth C2-class Gaussian quadrature checks only; no continuous-path execution or A4 remainder certificate.',
         g_hessian_bounds=[aa*(1-eps),aa*(1+eps)],
         same_clock=rows,centered_actual_law_homotopy=homotopy,
         status='PASS',assertions=21)
path=Path(__file__).with_name('single_history_resummed_current_checks.json')
path.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
