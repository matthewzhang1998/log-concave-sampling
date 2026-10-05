"""Finite original-gradient source and rank-two current arithmetic.

This module implements the 19-call LOCAL source once its exact finite Gaussian
rows/clocks are supplied. It does not sample true nonlinear histories, execute
an imported native compiler, or supply unresolved base/singleton coefficients.
"""
from dataclasses import dataclass
from typing import Callable
import numpy as np
Array=np.ndarray
Gradient=Callable[[Array],Array]

@dataclass
class ReplicaRows:
    x: Array
    xp: Array
    z: Array
    zp: Array
    smooth_outer: Array
    smooth_inner: Array

@dataclass
class LocalMarkRows:
    x0: Array
    x1: Array
    x2: Array
    a: float=1.0
    b: float=1.0


def project(x:Array,radius:float)->Array:
    if radius<0: raise ValueError('negative cap')
    norm=float(np.linalg.norm(x))
    return x.copy() if norm<=radius else x*(radius/norm)


def pair_secant(g:Gradient, rows:ReplicaRows, epsilon:float, delta:float, scale:float)->Array:
    """Exactly 7 original gradient VALUES, with saved coherent anchors."""
    if epsilon<=0 or delta<=0: raise ValueError('positive smoothing/secant scales required')
    anchor=g(epsilon*rows.smooth_inner)
    u=g(rows.z+epsilon*rows.smooth_inner)-anchor
    up=g(rows.zp+epsilon*rows.smooth_inner)-anchor
    outer=epsilon*rows.smooth_outer
    a=g(rows.x+outer-delta*u)
    b=g(rows.xp+outer-delta*u)
    c=g(rows.x+outer-delta*up)
    d=g(rows.xp+outer-delta*up)
    return scale*(a-b-c+d)/(2*delta)


def local_coherent_mark(g:Gradient,rows:LocalMarkRows)->Array:
    """Exactly 5 original gradient VALUES; genuine local Delta=U-V."""
    h0=g(rows.x0)
    h1=g(rows.x1-rows.a*h0)
    h1b=g(rows.x1)
    u=g(rows.x2-rows.b*h1)
    v=g(rows.x2-rows.b*h1b)
    return u-v


def h_action(f:Array,fp:Array,x:Array)->Array:
    return (f*np.dot(fp,x)+fp*np.dot(f,x))/2


def first_current_polynomial(f:Array,fp:Array,carrier:Array,v:float)->float:
    """T1(E,G), E=H/v+H²/(4v²), using O(D) work and no dense matrix."""
    if v<=0: raise ValueError('positive buffer required')
    hg=h_action(f,fp,carrier)
    trh=float(np.dot(f,fp))
    trh2=(float(np.dot(f,f))*float(np.dot(fp,fp))+trh*trh)/2
    quadratic=float(np.dot(carrier,hg))/v+float(np.dot(hg,hg))/(4*v*v)
    trace=trh/v+trh2/(4*v*v)
    return (quadratic-trace)/2


def exact_likelihood(f:Array,fp:Array,carrier:Array,v:float)->float:
    """Exact rank-two likelihood using a guarded Woodbury 2x2 solve.

    No inverse Gram matrix, eigenspace derivative or rank threshold is used.
    carrier is standard G; physical endpoint is sqrt(v)*G.
    """
    if v<=0: raise ValueError('positive buffer required')
    R=max(float(np.linalg.norm(f)),float(np.linalg.norm(fp)))
    if R*R/(2*v)>1/8+1e-14: raise ValueError('rank-two likelihood guard failed')
    U=np.column_stack((f,fp))
    V=np.column_stack((fp,f))/(4*v)
    small=np.eye(2)+V.T@U
    sign,logdet=np.linalg.slogdet(small)
    if sign<=0: raise ArithmeticError('positive determinant guard failed')
    displacement=U@np.linalg.solve(small,V.T@carrier)
    loglike=-logdet+float(np.dot(carrier,displacement))-float(np.dot(displacement,displacement))/2
    return float(np.exp(loglike))


def covariance_transport(f:Array,fp:Array,carrier:Array,v:float)->Array:
    return np.sqrt(v)*carrier+h_action(f,fp,carrier)/(2*np.sqrt(v))


def nineteen_value_source(g:Gradient,replica0:ReplicaRows,replica1:ReplicaRows,
                          mark_rows:LocalMarkRows,carrier:Array,*,
                          epsilon:float,delta:float,scale:float,
                          pair_cap:float,mark_cap:float,carrier_cap:float,v:float,
                          return_raw_marker:bool=False):
    """One bounded correction source evaluation: exactly 19 gradient calls.

    Freeze all clock labels and carrier as native caller inputs. The two
    replica rows must have the prescribed shared geometry and independent
    nuisance Gaussian rows. The mark rows retain their supplied shared
    ancestry. Native replays must rebuild all affected rows and calls.
    """
    f=project(pair_secant(g,replica0,epsilon,delta,scale),pair_cap)
    fp=project(pair_secant(g,replica1,epsilon,delta,scale),pair_cap)
    mraw=local_coherent_mark(g,mark_rows)
    m=project(mraw,mark_cap)
    gc=project(carrier,carrier_cap)
    if pair_cap*pair_cap/(2*v)>1/8: raise ValueError('buffer/cap guard failed')
    correction=m*first_current_polynomial(f,fp,gc,v)
    return mraw+correction if return_raw_marker else correction
