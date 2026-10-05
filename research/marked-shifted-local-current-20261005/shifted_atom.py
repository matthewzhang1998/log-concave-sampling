"""Seven original-VALUE calls for a bounded shifted atom; no history oracle."""
from dataclasses import dataclass
import numpy as np

@dataclass
class PairRecord:
    h0: np.ndarray
    h1: np.ndarray
    h1b: np.ndarray
    u: np.ndarray
    v: np.ndarray
    delta: np.ndarray
    query_sites: tuple

def coherent_pair(g, x0, x1, x2, a=1.0, b=1.0):
    q0 = np.asarray(x0)
    h0 = g(q0)
    q1 = np.asarray(x1) - a*h0
    h1 = g(q1)
    q1b = np.asarray(x1)
    h1b = g(q1b)
    q2 = np.asarray(x2) - b*h1
    u = g(q2)
    q2b = np.asarray(x2) - b*h1b
    v = g(q2b)
    return PairRecord(h0,h1,h1b,u,v,u-v,(q0,q1,q1b,q2,q2b))

def shifted_atom(g, x0, x1, x2, x3, a=1.0, b=1.0, c=1.0):
    record=coherent_pair(g,x0,x1,x2,a,b)
    q3=np.asarray(x3)-c*record.u
    q3b=np.asarray(x3)-c*record.v
    e=g(q3)-g(q3b)
    return e,record,record.query_sites+(q3,q3b)

def marked_current(record, keep, variance, q=1.0, theta=1.0):
    if not variance>0 or not (0<=q<=1 and 0<=theta<=1):
        raise ValueError('Require variance>0 and q, theta in [0,1].')
    z=np.sqrt(variance)*np.asarray(keep)-q*((1-theta)*record.v+theta*record.u)
    return record.delta,z,-q*record.delta
