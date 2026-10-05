#!/usr/bin/env python3
"""Independent execution of equation (16) and the imported zero-clock VALUE action."""
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss

class Source:
    def __init__(self,delta,r,a,b,xi,B=None):
        self.delta,self.r,self.a,self.b,self.xi=delta,r,a,b,xi
        self.B=B; self.calls=0
        d=delta/2; self.sigma=np.sqrt(np.tanh(d)); alpha=1/(2*np.cosh(d))
        x,w=leggauss(24); u=np.r_[d*(x+1)/2,d+d*(x+1)/2]
        weights=np.tile(d*w/2,2)
        left=u<d; v=np.where(left,u,u-d)
        h0=np.sinh(d-v)/np.sinh(d); h1=np.sinh(v)/np.sinh(d)
        psi=np.where(left,h1,h0)
        self.ca=np.where(left,h0+alpha*h1,alpha*h0)
        self.cb=np.where(left,alpha*h1,h1+alpha*h0)
        kap2=(-np.expm1(-2*v))*(-np.expm1(-2*(d-v)))/(-np.expm1(-2*d))
        self.gamma=weights*np.exp(-u)*psi
        self.gamma*=d*np.exp(-d)/self.gamma.sum()
        self.psi=psi; self.tau=np.sqrt(kap2+self.sigma**2*psi**2*(1-r*r))
        self.z=self.ca[:,None]*a+self.cb[:,None]*b+r*self.sigma*psi[:,None]*xi
        self.ell=self.sigma*d*np.exp(-d)
        self.vectors=np.array([[1.,0.],[.6,.8]])
    def g(self,x):
        self.calls+=len(np.atleast_2d(x))
        if self.B is not None:return x@self.B.T
        return .55*x+.09*np.sin(x@self.vectors.T)@self.vectors
    def H(self,x):
        if self.B is not None:return np.tile(self.B,(len(x),1,1))
        return .55*np.eye(2)[None,:,:]+.09*np.einsum('ni,ij,ik->njk',np.cos(x@self.vectors.T),self.vectors,self.vectors)
    def f(self,N):
        shifted=self.z+self.tau[:,None]*N
        return self.sigma*np.einsum('i,ij->j',self.gamma/self.tau,self.g(shifted)-self.g(self.z))
    def dN(self,N):return self.sigma*np.einsum('i,ijk->jk',self.gamma,self.H(self.z+self.tau[:,None]*N))
    def dc(self,N,co):
        return self.sigma*np.einsum('i,ijk->jk',self.gamma*co/self.tau,self.H(self.z+self.tau[:,None]*N)-self.H(self.z))

def filt(K):
    b2=(.25+1/(4*K*K))/2+(.25-1/(4*K*K))/2*np.cos(np.arange(K)*np.pi/(K-1))
    d=np.array([np.prod([-b2[j]/(b2[i]-b2[j]) for j in range(K) if j!=i]) for i in range(K)])
    return np.sqrt(b2),d

def square(src,p,z1,z0,z2,mu=.1,K=6):
    ell=src.ell; s0=.2; s=np.sqrt(mu/100)
    f=lambda u:s0*src.f(u)/ell
    bs,ds=filt(K)
    I=sum(d*(f(np.sqrt(1-b*b)*z1+b*p)-f(np.sqrt(1-b*b)*z1-b*p))/(2*b) for b,d in zip(bs,ds))
    R=lambda w:(f(.5*w+np.sqrt(.75)*z2)-f(-.5*w+np.sqrt(.75)*z2))
    return ell*ell*(R(z0+s*I)-R(z0-s*I))/(2*s)

def main():
    rng=np.random.default_rng(5347); maxder=0.; maxenergy=0.; maxcap=0.; maxparity=0.; maxquad=0.; counts=[]
    for delta in [.0001,.01,.2,.69]:
      for r in [0.,.5,.999999]:
        a,b,xi,N=rng.normal(size=(4,2)); src=Source(delta,r,a,b,xi)
        h=2e-5; dfd=np.column_stack([(src.f(N+h*np.eye(2)[j])-src.f(N-h*np.eye(2)[j]))/(2*h) for j in range(2)])
        maxder=max(maxder,np.linalg.norm(dfd-src.dN(N))/src.ell)
        assert np.linalg.eigvalsh(src.dN(N)).min()>-1e-12
        assert np.linalg.norm(src.dN(N),2)<=src.ell*(1+1e-12)
        maxenergy=max(maxenergy,np.linalg.norm(src.f(N))/(src.ell*np.linalg.norm(N)))
        maxcap=max(maxcap,np.linalg.norm(src.dc(N,src.ca),2)/delta,np.linalg.norm(src.dc(N,src.cb),2)/delta,np.linalg.norm(src.dc(N,r*src.sigma*src.psi),2)/delta**1.5)
        p,z1,z0,z2=rng.normal(size=(4,2)); src.calls=0
        c=square(src,p,z1,z0,z2); calls=src.calls; counts.append(calls)
        # Conservative uncached count: 2 Nu times (2 Kf + 4)
        assert calls==2*48*(2*6+4)
        maxparity=max(maxparity,np.linalg.norm(c+square(src,-p,z1,z0,z2)))
        assert np.array_equal(src.f(np.zeros(2)),np.zeros(2))
        assert np.array_equal(square(src,np.zeros(2),z1,z0,z2),np.zeros(2))
        B=np.array([[.6,.19],[.19,.4]]); lin=Source(delta,r,a,b,xi,B=B)
        J=lin.ell*B; val=square(lin,p,z1,z0,z2)
        maxquad=max(maxquad,np.linalg.norm(val-.2**2*J@J@p)/max(1e-100,np.linalg.norm(J@J@p)))
    assert maxder<1e-6 and maxenergy<=1 and maxcap<2 and maxparity<1e-12 and maxquad<1e-6
    result={'scope':'Literal anchored VALUE-source and zero-clock square execution; no Monte Carlo calibration inference','max_relative_private_derivative_error':maxder,'max_energy_to_pointwise_radius_ratio':maxenergy,'max_scaled_captured_derivative':maxcap,'max_absolute_oddness_error':maxparity,'max_relative_linear_matrix_calibration_error':maxquad,'original_value_calls_per_uncached_bank':sorted(set(counts)),'inner_nodes':48,'first_filter_nodes':6}
    Path(__file__).with_name('literal_graph_checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
