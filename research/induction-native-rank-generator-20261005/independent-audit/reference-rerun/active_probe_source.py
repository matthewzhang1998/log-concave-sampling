"""All-rank original-VALUE active source. No derivative oracle is executed."""
from __future__ import annotations
import itertools
import math
import numpy as np


def radial_geometry(x, R):
    """Return y=psi(x) and J=c I+d uu^T in O(D) storage/work."""
    x=np.asarray(x,dtype=float)
    r=float(np.linalg.norm(x))
    if R<=0: raise ValueError('R must be positive')
    if r<=R: return x.copy(),1.0,0.0,np.zeros_like(x)
    if r<2*R:
        v=(r-R)/R
        h=R*(1+v-v**3+v**4/2)
        hp=1-3*v*v+2*v**3
    else:
        h=1.5*R;hp=0.0
    u=x/r
    return h*u,h/r,hp-h/r,u


def radial_pullback(x, R):
    """Dense diagnostic Jacobian; the production source uses its O(D) action."""
    y,c,d,u=radial_geometry(x,R)
    return y,c*np.eye(y.size)+d*np.outer(u,u)


def active_source(g, x, probes, z, t, A, R=None):
    """Return (vector, original_VALUE_count). All sign terms share the bank.

    Optional R executes the exact gradient-preserving radial pullback.
    Floating-point arithmetic is diagnostic, not a certified arithmetic oracle.
    """
    x=np.asarray(x,dtype=float);z=np.asarray(z,dtype=float)
    P=np.asarray(probes,dtype=float)
    if P.size==0: P=np.empty((0,x.size))
    if P.ndim!=2 or P.shape[1]!=x.size: raise ValueError('probe shape')
    if t<=0: raise ValueError('t must be positive')
    if A<0: raise ValueError('A must be nonnegative')
    if A==0:return np.zeros_like(x),0
    k=len(P); a=1/math.sqrt(k+1)
    if R is None: y=x;c=1.;d=0.;u=None
    else:y,c,d,u=radial_geometry(x,R)
    out=np.zeros_like(x);calls=0
    for eps in itertools.product((-1.,1.),repeat=k):
        shift=np.asarray(eps)@P if k else np.zeros_like(x)
        q=z+a*t*shift
        out+=math.prod(eps)*(np.asarray(g(q+a*t*y))-np.asarray(g(q)))
        calls+=2
    pulled=c*out if u is None or d==0 else c*out+d*float(u@out)*u
    return pulled/(A*t*2**k),calls


def projection_cut_constant(k,p):
    q=k+2-p
    if not 1<=p<=k+1:raise ValueError('proper cut required')
    a=1/math.sqrt(k+1)
    return a**(k+1)*2**(k/2)*math.sqrt(math.factorial(p-1)*math.factorial(q-1))


def response_floor(k,D,R):
    if R<math.sqrt(D):raise ValueError('R must be at least sqrt(D)')
    a=1/math.sqrt(k+1)
    return 2.5*a*(D*(D+2))**.25*math.exp(-(R-math.sqrt(D))**2/8)
