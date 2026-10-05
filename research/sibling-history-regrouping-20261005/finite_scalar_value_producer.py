#!/usr/bin/env python3
"""Actual finite scalar VALUE sampler for the 64-history interior-parent packet.
No g derivative is called. Numerical clock arithmetic is finite and its bias is
not included in the exact-real Monte Carlo/cutoff certificate printed below.
This is a direct estimator, not an execution of an imported native pair program.
"""
import math, json
from pathlib import Path
import numpy as np


def root_mass(r):
    r=float(r)
    if r==0:return 0.
    z=r*r
    if z<.25:
        # Positive series avoids catastrophic subtraction of the closed form.
        total=0.;power=z**4
        for j in range(256):
            term=(j+1)*(j+2)*power/(4*(j+4))
            total+=term
            if term<max(1e-300,total*1e-17):break
            power*=z
        return total
    u=1-z
    return .75+1/(4*u*u)-1.5/u-1.5*math.log(u)+u/2


def inverse_root_cdf(U,rstar,steps=64):
    I=root_mass(rstar);lo=np.zeros_like(U);hi=np.full_like(U,rstar)
    # Fixed finite scalar geometry computation, independent of the force.
    for _ in range(steps):
        mid=(lo+hi)/2
        mass=np.fromiter((root_mass(r) for r in mid),dtype=float,count=len(mid))
        below=mass<I*U
        lo=np.where(below,mid,lo);hi=np.where(below,hi,mid)
    return (lo+hi)/2


def cutoff(x,R):return np.minimum(1.,np.maximum(0.,R+1-np.abs(x)))


def sample_packet(g,A,y,N=4096,rstar=.75,R=7.,seed=20261005):
    if not 0<rstar<1:raise ValueError('rstar must be strictly between zero and one')
    if not A>0:raise ValueError('A must be positive; A=0 is handled as the zero source')
    rng=np.random.default_rng(seed)
    U=rng.random((N,8));r=np.sqrt((1-U)*(1+U));sig=U/math.sqrt(2)
    rho=np.empty((N,7))
    rho[:,1:]=rng.random((N,6))**(1/8)
    rho[:,0]=inverse_root_cdf(rng.random(N),rstar)
    G=rng.normal(size=(N,7));W=rng.normal(size=(N,8));Z=rng.normal(size=(N,8))
    q=np.zeros((N,8));parent=np.full(N,float(y))
    for n in range(8,1,-1):
        j=n-2
        parent=rho[:,j]*parent+np.sqrt((1-rho[:,j])*(1+rho[:,j]))*G[:,j]
        if n>=3:q[:,n-1]=r[:,n-1]*parent+sig[:,n-1]*W[:,n-1]
    q[:,:2]=r[:,:2]*parent[:,None]+sig[:,:2]*W[:,:2]
    # Exactly two original g VALUE entries per force occurrence. No /sigma.
    J=math.sqrt(2)*Z*(g(q+sig*Z)-g(q))/A
    H=G[:,0]**6-15*G[:,0]**4+45*G[:,0]**2-15
    C=math.factorial(8)*A**8*root_mass(rstar)/8**6
    T=C*H*cutoff(G[:,0],R)*np.prod(J*cutoff(Z,R),axis=1)
    E=C*math.sqrt(720)
    return dict(estimate=float(np.mean(T)),empirical_standard_error=float(np.std(T,ddof=1)/math.sqrt(N)),
                exact_real_rms_bound=E/math.sqrt(N),exact_real_cutoff_bias_bound=E*math.sqrt(18)*math.exp(-R*R/4),
                source_energy_bound=E,retained_caller_L2_first_bound=16*math.sqrt(2)*E,
                original_g_value_count=16*N,normal_draw_count=23*N,scalar_clock_draw_count=15*N,
                samples=N,rstar=rstar,root_normalizer=root_mass(rstar),R=R,
                numerical_geometry_note='64-step bisection and float arithmetic executed; propagated geometry/force precision debt is separate from exact-real bounds.')


if __name__=='__main__':
    rows=[]
    rows.append({'source':'g(x)=x; true packet is exactly zero',**sample_packet(lambda x:x,1.,.3)})
    for K in (4,16,64,256):
        rows.append({'source':f'g(x)=x/2+sin({K}x)/(4*{K})','K':K,
                     **sample_packet(lambda x,K=K:.5*x+np.sin(K*x)/(4*K),.75,.3)})
    out={'scope':'Actual finite original-VALUE scalar sampler; RMS bounds are not deterministic realization error guarantees. No native theorem.', 'runs':rows}
    Path(__file__).with_name('producer-runs.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
