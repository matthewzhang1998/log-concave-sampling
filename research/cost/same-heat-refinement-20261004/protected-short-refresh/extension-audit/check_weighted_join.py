#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction as F
import numpy as np
import json, math
OUT=Path(__file__).resolve().parent

def kernel(n,T):
    t=np.arange(n)*T/n
    W=np.zeros((n,n))
    for i in range(n):
        lo=t[i]-t[:i]; hi=lo-T/n
        W[i,:i]=2*np.sin((lo+hi)/2)*np.sin((lo-hi)/2)
    lo=T-t; hi=lo-T/n
    w=2*np.sin((lo+hi)/2)*np.sin((lo-hi)/2)
    return t,W,w

def build(nq,ns,a=.01,Mq=2,Ms=2):
    h=a**(19/30)
    tq,Wq,wq=kernel(nq,math.pi/2); ts,Ws,ws=kernel(ns,h)
    n=3+Mq*nq+Ms*ns+1
    A=np.zeros((n,n)); C=np.zeros((n,4)); ds=np.ones(n)
    A[:3,:3]=np.array([[0,0,0],[.012,0,0],[.008,.006,0]])
    C[:3]=np.array([[.3,0,0,.2],[.2,0,0,.1],[.4,0,0,.3]])
    oldR=a*np.array([.2,.1,.1])
    qgroups=[slice(3+k*nq,3+(k+1)*nq) for k in range(Mq)]
    sgroups=[slice(3+Mq*nq+k*ns,3+Mq*nq+(k+1)*ns) for k in range(Ms)]
    for k,g in enumerate(qgroups):
        ds[g]=1/math.sqrt(nq)
        C[g,0]=np.cos(tq)/math.sqrt(nq); C[g,1]=np.sin(tq)/math.sqrt(nq); C[g,3]=1/math.sqrt(nq)
        A[g,:3]=np.outer(np.cos(tq)/math.sqrt(nq),oldR)
        if k: A[g,qgroups[k-1]]=-a*Wq
    for k,g in enumerate(sgroups):
        ds[g]=1/math.sqrt(ns)
        C[g,1]=np.cos(ts)/math.sqrt(ns); C[g,2]=np.sin(ts)/math.sqrt(ns); C[g,3]=1/math.sqrt(ns)
        A[g,qgroups[-1]]=-a*np.outer(np.cos(ts)/math.sqrt(ns),wq*math.sqrt(nq))
        if k: A[g,sgroups[k-1]]=-a*Ws
    C[-1]=[0,math.cos(h),math.sin(h),1]
    A[-1,qgroups[-1]]=-a*math.cos(h)*wq*math.sqrt(nq)
    A[-1,sgroups[-1]]=-a*ws*math.sqrt(ns)
    return A,C,ds,tq,Wq,wq,ts,Ws,ws,h

cases=[]
for nq,ns in [(4,3),(17,7),(64,9),(127,23)]:
    A,C,d,tq,Wq,wq,ts,Ws,ws,h=build(nq,ns)
    K=2*math.sin(h/2)**2; a=.01
    assert max(Wq.sum(axis=0).max(),Wq.sum(axis=1).max())<=1+1e-13
    assert max(Ws.sum(axis=0).max(),Ws.sum(axis=1).max())<=K+1e-13
    assert abs(wq.sum()-1)<1e-13 and abs(ws.sum()-K)<1e-13
    assert max(wq)<=math.pi/(2*nq)+1e-13
    assert max(ws)<=h*h/ns+1e-13
    Anorm=float(np.linalg.norm(abs(A),2)); Cnorm=float(np.linalg.norm(C,2))
    terminal=float(np.linalg.norm(A[-1],2))
    Lnorm=float(np.linalg.norm(C[:-1,2]))
    assert Anorm<.1 and Cnorm<4 and terminal<2*a
    assert Lnorm<=math.sqrt(2)*h+1e-13
    # Similarity transform is exact, including output scaling and primitive width.
    Cphys=C/d[:,None]; Mphys=A*d[None,:]/d[:,None]
    x=np.array([.4,-.8,1.1,.2]); p=np.zeros(len(d)); y=np.zeros(len(d))
    f=lambda z:.5*z+.25*np.tanh(z)
    for i in range(len(d)):
        p[i]=d[i]*f((C[i]@x+A[i]@p)/d[i])
        y[i]=f(Cphys[i]@x+Mphys[i]@y)
    scale_error=float(np.max(abs(p/d-y)))
    assert scale_error<2e-14
    # Weighted initial forcing uses input rank, not nonlinear-node dimension.
    affine_fro=float(np.linalg.norm(C,'fro'))
    assert affine_fro<4
    cases.append({'Nq':nq,'Ns':ns,'nodes':len(d),'absolute_edge_norm':Anorm,
      'affine_operator_norm':Cnorm,'affine_frobenius_norm':affine_fro,
      'terminal_norm_over_a':terminal/a,'restricted_L_norm_over_h':Lnorm/h,
      'exact_scaled_program_error':scale_error})

# Exact gradient property of the normalized, projected, symmetric twin reference.
A,C,d,*_=build(17,7)
h=.01**(19/30); eps=.01**.9; gamma=math.sin(h)**2
C[:-1,2]=0
n=len(d); et=np.zeros((n,1)); et[-1]=1
S=np.block([[A+A.T,A.T],[A,np.zeros_like(A)]])
Ct=np.block([[C,np.zeros((n,1))],[C,eps*et]])
beta=1/math.sqrt(1/gamma+1/eps**2)
B=beta*np.array([[0,0,math.sin(h)/gamma,0,-1/eps]])
assert np.max(abs(Ct@B.T-beta*np.vstack([et,np.zeros_like(et)])))<2e-14
assert abs((B@B.T)[0,0]-1)<2e-14
D=np.concatenate([d,d]); signs=np.concatenate([np.ones(n),-np.ones(n)])
def solve(x):
    p=np.zeros(2*n)
    for _ in range(100):
        z=(Ct@x+S@p)/D
        nxt=signs*D*(.5*z+.25*np.tanh(z))
        if np.max(abs(nxt-p))<1e-15:
            p=nxt; break
        p=nxt
    else: raise AssertionError('not converged')
    return Ct.T@p
x=np.array([.3,-.8,1.1,.2,.4]); step=2e-6
Jac=np.column_stack([(solve(x+step*np.eye(5)[i])-solve(x-step*np.eye(5)[i]))/(2*step) for i in range(5)])
curl=float(np.linalg.norm(Jac-Jac.T,2))
assert curl<1e-8

v=F(19,15); J=F(29,10)
child=[]
for a in [.1,.01,.001]:
    for ratio in [.5,1/math.log(1/a)]:
        eta=a*ratio; he=eta**float(v/2)
        gamhidden=ratio*math.sin(he)**2
        expected=ratio**float(1+v)*a**float(v)
        force=math.sqrt(ratio)*eta**float(J)
        assert (25/36)*expected<=gamhidden<=expected*(1+1e-13)
        assert abs(force/(a**float(J))-ratio**float(J+F(1,2)))<1e-12
        child.append({'a':a,'eta_over_a':ratio,'gamma_over_a_v':gamhidden/a**float(v),'force_over_r_a_J':force/a**float(J)})

rec=[]
k=F(0); P=F(7)
for _ in range(12):
    R=P+F(1,2); k=max(k,R-F(3,2)); assert k==R-F(3,2)
    rec.append({'R':str(R),'k_R':str(k),'accuracy_exponent':str((R-F(1,2))/(2*(R-1))),
      'kappa_exponent':str(2+1/(2*(R-1)))})
    P=R
result={'status':'PASS bounded weighted-graph, exact-scaling, gradient-reference, child-scale, and recurrence checks',
        'weighted_cases':cases,'symmetric_reference_curl_fd':curl,'child_scaling':child,'known_center_recurrence':rec,
        'scope':'Relative admission only; existing seed and finite numerical source restoration remain explicit obligations.'}
(OUT/'weighted_join_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(result['status']); print('weighted cases:',len(cases),'reference curl:',curl)
