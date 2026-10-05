#!/usr/bin/env python3
"""Scalar positive native-pair feedback and cumulant-normalization diagnostics."""
import json, math
from pathlib import Path
import numpy as np
from scipy.special import roots_hermitenorm, ndtr, ndtri
from scipy.interpolate import PchipInterpolator
from numpy.polynomial.legendre import leggauss
p,wp=roots_hermitenorm(120);wp/=math.sqrt(2*math.pi)
g,wg=roots_hermitenorm(32);wg/=math.sqrt(2*math.pi)
A=.15+.1*np.tanh(g)
meanA=wg@A; h=math.sqrt(wg@(A*A))
c=a=b=eta=.5; beta2=c*c+a*a+b*b; nu=beta2+eta*eta
P=np.broadcast_to(p[None,:],(len(g),len(p))).ravel()
M=np.clip((A[:,None]*p[None,:]).ravel(),-2,2)
w=(wg[:,None]*wp[None,:]).ravel()
xgrid=np.linspace(-9,9,16001)
uq,wq=leggauss(2400);uq=(uq+1)/2;wq/=2

def inv(cdf):
    cdf=np.maximum.accumulate(cdf);keep=np.r_[True,np.diff(cdf)>1e-16]
    return PchipInterpolator(cdf[keep],xgrid[keep])(uq)

def mix(mean,sd,weight):
    out=np.zeros_like(xgrid)
    for st in range(0,len(mean),96):
        out+=ndtr((xgrid[:,None]-mean[None,st:st+96])/sd[None,st:st+96])@weight[st:st+96]
    return out

res={"assertions":0,"checks":[]}
def check(x):
    assert bool(x);res['assertions']+=1
check(abs(meanA-.15)<1e-14)
for q in [.12,.08,.05,.03]:
    # Integrate Y, xi, Z analytically; P,Q remain true independent Gaussian banks.
    m=b*P; var=c*c+a*a+eta*eta+2*c*a*q*M
    check(var.min()>0)
    exactmean=w@m
    exactvar=w@(m*m+var)
    exactthird=w@(m**3+3*m*var)
    targetthird=6*a*b*c*q*meanA
    check(abs(exactmean)<1e-13)
    check(abs(exactvar-nu)<1e-12)
    check(abs(exactthird-targetthird)<2e-9)
    qp=inv(mix(m,np.sqrt(var),w))
    S=math.sqrt(beta2)*p
    refmean=S+q*a*b*c*meanA*(S*S/beta2**2-1/beta2)
    qr=inv(mix(refmean,np.full_like(S,eta),wp))
    dist=math.sqrt(wq@((qp-qr)**2))
    res['checks'].append({'q':q,'pair_variance':exactvar,'pair_third_moment':exactthird,'target_third':targetthird,'W2_to_positive_reference':dist,'W2_over_q2_h':dist/(q*q*h)})
# Verify resolvent hierarchy polynomial fixture, by exact rational arithmetic.
from fractions import Fraction as F
k_a=F(16,45);k_b=F(32,135)
check(5*k_a==F(16,9));check(3*k_b-2*k_a==0)
res['hierarchy_fixture']={'g':'x^2','kappa3_x2_coefficient':str(k_a),'kappa3_constant':str(k_b)}
res['scope']='Scalar positive Gaussian-pair CDF diagnostics, including retained-center-public dependence and owned coarse-Q averaging. Not a numerical execution of the imported finite native VALUE compiler.'
Path(__file__).with_name('positive_rank3_packet_checks.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
