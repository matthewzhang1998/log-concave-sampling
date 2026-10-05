import json, math
import numpy as np
import mpmath as mp
mp.mp.dps=70
from numpy.polynomial.legendre import leggauss

def rule(epsilon):
    # Start with epsilon/3 to allow exact first-moment correction.
    e=epsilon/3
    a=e/8; T=math.log(8/e)
    n=math.ceil(math.log(16*T/e,4))
    base_nodes,base_weights=leggauss(n)
    corr=[1.0,0.0]; weights=[-math.expm1(-a),math.exp(-T)]
    s=a
    intervals=0
    while s<T:
        t=min(2*s,T)
        nodes=(s+t)/2+(t-s)*base_nodes/2
        cs=(t-s)*base_weights/2
        corr.extend(np.exp(-nodes))
        weights.extend(cs*np.exp(-nodes))
        s=t; intervals+=1
    r=np.array(corr); p=np.array(weights); p/=p.sum()
    m=float(p@r)
    if m>0.5:
        theta=1-0.5/m; p*=1-theta; p[1]+=theta
    elif m<0.5:
        theta=(0.5-m)/(1-m); p*=1-theta; p[0]+=theta
    return r,p,n,intervals

def phi(t,A):
    alpha=beta=.5
    cv=alpha*(1-1/math.sqrt(2)); tau=cv/2
    eta=A*A
    def hinge(x):
        if x<=-eta:return 0.
        if x>=eta:return x
        return (x+eta)**2/(4*eta)
    return alpha*hinge(t-.5)+beta*hinge(t-(1-A*tau))

def row_g(v,A):
    v=np.asarray(v,float); norm=np.linalg.norm(v)
    return np.zeros_like(v) if norm==0 else A*phi(norm,A)*v/norm

out=[]
for A in [1/36,1e-2,1e-3,1e-4,1e-6,1e-8]:
    eps=A**.75; r,p,n,intervals=rule(eps)
    w=math.sqrt(A); q=1-w; sigma=math.sqrt(1-q*q)
    # High-precision radial row arithmetic avoids cancellation in the A^3 witness.
    AA=mp.mpf(str(A)); ww=mp.sqrt(AA); qq=1-ww; ss=mp.sqrt(1-qq*qq)
    beta=mp.mpf(str(float(p@np.sqrt(1-r*r))))
    xx=mp.matrix([1,0,0])
    def ph(t):
        eta=AA*AA; tau=(1-1/mp.sqrt(2))/4
        def hinge(z):
            if z<=-eta:return mp.mpf(0)
            if z>=eta:return z
            return (z+eta)**2/(4*eta)
        return (hinge(t-mp.mpf('.5'))+hinge(t-(1-AA*tau)))/2
    def rg(v):
        norm=mp.sqrt(sum(a*a for a in v))
        return AA*ph(norm)/norm*v
    V=qq*AA*ph(1)*mp.matrix([qq/2,ss/2,beta])
    u=xx-V
    terminal=rg(u-ww*rg(u))
    cU=ph(1); q2=mp.sqrt(1-AA*cU+AA*AA*cU*cU/2); c2=ph(q2)/q2
    vv=mp.matrix([mp.mpf('.5'),mp.mpf('.5'),0]); zz=mp.matrix([mp.mpf('.25'),mp.mpf('.5'),mp.mpf('.25')])
    exact_terminal=rg(xx-AA*c2*vv+AA*AA*c2*cU*zz)
    witness=(terminal[0]-exact_terminal[0])/2
    k=np.unique(np.concatenate((np.arange(1001),np.geomspace(1,1e10,1000).astype(np.int64))))
    max_spectral=float(np.max(np.abs((r[None,:]**k[:,None])@p-1/(k+1))))
    theorem_constant=(A**3+(2*math.sqrt(2)/3)*A*A*w**1.5+q*w*A**3+q*A*A*eps+q*q*A**3*(.5+1/sigma))/A**2.75
    out.append(dict(A=A,nodes=len(r),gauss_order=n,dyadic_intervals=intervals,mass=float(p.sum()),first_moment=float(p@r),min_weight=float(p.min()),operator_epsilon=eps,sampled_spectral_error=max_spectral,theorem_bias_constant=theorem_constant,radial_first_chaos_witness_over_A2=float(witness/AA**2),radial_first_chaos_witness_over_A3=float(witness/AA**3)))
print(json.dumps(out,indent=2))
