#!/usr/bin/env python3
"""Diagnostics of the literal finite source, not an implementation of mean compilers."""
import json
import math
from pathlib import Path
import numpy as np

rng=np.random.default_rng(20261005)
checks=0
def check(test):
    global checks
    assert bool(test)
    checks+=1

def row(s):
    r=math.exp(-s); a=2*s*r; b=2*s*(s-1)*r
    c=math.sqrt(max(0.,1-r*r-a*a-b*b))
    return np.array([r,a,b,c])

times=np.array([.2,1.2])
rr=np.exp(-times)
p=np.array([(.5-rr[1])/(rr[0]-rr[1]),(rr[0]-.5)/(rr[0]-rr[1])])
check(np.all(p>0));check(abs(p.sum()-1)<1e-15);check(abs(p@rr-.5)<1e-15)

class Source:
    def __init__(self,A,D,K=None):
        self.A=A;self.D=D;self.K=K;self.calls=0
        self.a=np.array([1.,0.]);self.b=np.array([.6,.8])
    def g(self,x):
        self.calls+=1
        if self.K is not None:return self.K@x
        z=x.reshape(-1,2)
        return (self.A*(z/2+(np.sin(z@self.a)[:,None]*self.a+np.sin(z@self.b)[:,None]*self.b)/8)).reshape(-1)
    def h(self,x):
        if self.K is not None:return self.K
        out=np.zeros((self.D,self.D))
        for j in range(0,self.D,2):
            z=x[j:j+2]
            out[j:j+2,j:j+2]=self.A*(np.eye(2)/2+(math.cos(z@self.a)*np.outer(self.a,self.a)+math.cos(z@self.b)*np.outer(self.b,self.b))/8)
        return out

def node(src,root,first=False):
    x,N,M,L=root
    v=x/2+N/2;w=x/4+N/2+M/4
    gv=src.g(v);gw=src.g(w);ggw=src.g(gw);S=x-gv+ggw
    out=src.g(S);D=src.D
    if first:
        V=src.h(v);J=src.h(gw);W=src.h(w)
        P=np.eye(D)-V/2+J@W/4;Q=-V/2+J@W/2;R=J@W/4
        Sjac=np.concatenate([P,Q,R,np.zeros((D,D))],axis=1)
        jac=src.h(S)@Sjac
    for weight,s in zip(p,times):
        r,a,b,c=row(s);U=r*x+a*N+b*M+c*L
        gu=src.g(U);delta=gv-gu;gp=src.g(S+delta);gm=src.g(S-delta)
        out+=weight*(gp-gm)/2
        if first:
            H=src.h(U)
            Djac=np.concatenate([V/2-r*H,V/2-a*H,-b*H,-c*H],axis=1)
            Hp=src.h(S+delta);Hm=src.h(S-delta)
            jac+=weight*((Hp-Hm)@Sjac+(Hp+Hm)@Djac)/2
    base=src.g(x)
    if first:
        jac[:,:D]-=src.h(x)
        return out-base,jac
    return out-base

max_fd=0.; max_first_ratio=0.;max_curl_ratio=0.;max_quadratic=0.
for s in np.r_[np.logspace(-9,3,300),rng.uniform(0,25,300)]:
    r,a,b,c=row(s)
    check(abs(r*r+a*a+b*b+c*c-1)<1e-13)
    check(c>0)
    check(abs(a/2-s*math.exp(-s))<1e-13)
    check(abs(a/2+b/4-(s*s+s)*math.exp(-s)/2)<1e-13)

for D in [2,4,8]:
  for A in [.001,.03,.1,.5]:
    for repeat in range(4):
      src=Source(A,D);roots=rng.normal(size=(4,D))
      val,jac=node(src,roots,True)
      check(src.calls==5+3*len(p))
      h=2e-6; numeric=np.zeros_like(jac)
      for j in range(4*D):
        v=np.zeros(4*D);v[j]=h
        numeric[:,j]=(node(src,roots+v.reshape(4,D))-node(src,roots-v.reshape(4,D)))/(2*h)
      err=np.max(abs(numeric-jac));max_fd=max(max_fd,float(err));check(err<2e-8)
      beta=math.sqrt(3)/2
      # x-derivative becomes beta times the private G block.
      private=jac.copy();private[:,:D]*=beta
      lift=np.zeros((4*D,4*D));lift[:D,:]=private
      first=math.sqrt(2)*np.linalg.norm(private,2)
      curl=math.sqrt(2)*np.linalg.norm(lift-lift.T,2)
      max_first_ratio=max(max_first_ratio,first/A)
      max_curl_ratio=max(max_curl_ratio,curl/(A*A))
      check(first<=4*A+1e-12);check(curl<=12*A*A+1e-12)
      # Noncommutativity really occurs in the claimed source class.
      H1=src.h(roots[0]);H2=src.h(roots[1])
      check(np.linalg.norm(H1@H2-H2@H1)>1e-16)
      # Literal coherent zero.
      check(np.array_equal(node(src,np.zeros((4,D))),np.zeros(D)))
      O,_=np.linalg.qr(rng.normal(size=(D,D)))
      K=O@np.diag(rng.uniform(0,A,size=D))@O.T
      quad=Source(A,D,K)
      out=node(quad,roots)+K@roots[0]
      Ubar=sum(weight*(row(s)@roots) for weight,s in zip(p,times))
      w=roots[0]/4+roots[1]/2+roots[2]/4
      expected=K@roots[0]-K@K@Ubar+K@K@K@w
      qerr=np.max(abs(out-expected));max_quadratic=max(max_quadratic,float(qerr));check(qerr<1e-12)
      mean=node(quad,np.array([roots[0],np.zeros(D),np.zeros(D),np.zeros(D)]))+K@roots[0]
      target=(K-K@K/2+K@K@K/4)@roots[0]
      check(np.max(abs(mean-target))<1e-12)

# Verify uniform declared envelopes without sampling a source.
for A in np.linspace(.000001,.5,1000):
    eta=A/2+A*A/4;q=A/2+A*A/2;r0=A*A/4;beta=math.sqrt(3)/2
    Cx=1.5+1.5*eta+1.5*A;Cn=1.5*q+1.5*A;Cm=1.5*r0+A
    first=math.sqrt(2)*A*math.sqrt(beta*beta*Cx*Cx+Cn*Cn+Cm*Cm+A*A)
    curl=math.sqrt(2)*(beta*(3*A*eta+3*A*A)+A*math.sqrt(Cn*Cn+Cm*Cm+A*A))
    check(first<=4*A);check(curl<=12*A*A)

report={"status":"PASS","assertions":checks,"maximum_finite_difference_error":max_fd,
        "maximum_observed_normalized_first_over_A":max_first_ratio,
        "maximum_observed_normalized_curl_over_A_squared":max_curl_ratio,
        "maximum_quadratic_identity_error":max_quadratic,
        "test_rule_nodes":times.tolist(),"test_rule_weights":p.tolist(),
        "scope":"Raw graph and analytic envelopes only; no numerical implementation of the imported mean compilers or high-accuracy bridge quadrature is claimed."}
Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
