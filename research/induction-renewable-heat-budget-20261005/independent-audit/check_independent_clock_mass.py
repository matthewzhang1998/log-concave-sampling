#!/usr/bin/env python3
"""Independent local tests. Does not import the author implementation or execute native compilers."""
from pathlib import Path
from fractions import Fraction as F
import hashlib, json, math, random
import numpy as np
from numpy.polynomial.legendre import leggauss

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
rng=random.Random(104729)
count=0

def check(x, msg):
    global count
    count+=1
    if not x: raise AssertionError(msg)

def q(k,tag): return k if tag==0 else max(k-2,0)

def rule(J,m,j=1):
    xi,wi=leggauss(m)
    nodes=[]
    for p in range(J):
        a=2.**(-p-1)
        c=1-1.5*a
        for x,w in zip(xi,wi):
            r=c+a*x/2
            v=a*w/2
            nodes.append((r,v*r**(j-1),a,p,v))
    return nodes

# Exact Bernstein ellipse geometry: with rho=3/2, horizontal semiaxis
# is 13a/24, so the maximal modulus is its rightmost real point.
for p in range(256):
    a=F(1,2**(p+1)); c=1-F(3,2)*a; A=F(13,24)*a; B=F(5,24)*a
    check(A*A-B*B==a*a/4,'ellipse focal identity')
    check(c>=0 and c+A==1-F(23,24)*a and c+A<1,'ellipse modulus')

# Positivity, original resolvent masses, common-node rule, and monomial errors.
# Monomial Gauss errors are nonnegative (f^(2m)>=0 or exact lower degree).
max_excess=0.0
for J in [2,5,10,18]:
  for m in [2,4,8,16]:
    for j in [1,2,3,8,16]:
      rr=rule(J,m,j)
      rs=np.array([x[0] for x in rr]); ws=np.array([x[1] for x in rr])
      for r,w,a,p,leb in rr:
        check(0<w<=leb*(1+1e-13),'positive density mass')
        check(leb<=a*(1+1e-13) and a<=1-r+1e-14 and 1-r<=1-r*r+1e-14,'mass certificate')
      check(ws.sum()<1+1e-14,'total rule mass')
      for ell in list(range(100))+[10**k for k in range(2,11)]:
        err=1/(ell+j)-float(ws@np.power(rs,ell))
        check(err>=-1e-12,'positive monomial remainder')
        # A safe analytic approximation constant; only checks the predicted envelope.
        bound=2.**(-J)+24*(1.5)**(-2*m)
        check(err<=bound+1e-12,'uniform multiplier envelope')
        max_excess=max(max_excess,err-bound)

# Direct product preservation: amplitudes include only their own primitive weights.
allocation_max_relative=0.0
for trial in range(1000):
    N=rng.randrange(2,11)
    ks=[rng.randrange(0,9) for _ in range(N)]
    tags=[rng.randrange(2) for _ in range(N)]
    gammas=[F(rng.randrange(1,8),28) for _ in range(N)]
    betas=[F(rng.randrange(1,4),28) for _ in range(N)]
    alpha=10**rng.uniform(-3,-.2)
    a=F(20)
    H=sum(g*q(k,tag) for g,k,tag in zip(gammas,ks,tags))
    d=a-sum(betas)-H
    root=0
    rho=[]; original=alpha**float(a)
    remaining=-.31
    for k,tag,g,b in zip(ks,tags,gammas,betas):
        r,w,_,_,_=rng.choice(rule(8,5,rng.randrange(1,9)))
        h=alpha**float(g); t=math.sqrt((1-r*r+h*h)/2)
        bk=(k+1)**(-(k+1)/2)
        z=t**(-k)/bk
        original*=w*z
        if tag==0:
            remaining*=w
            rho.append(alpha**float(b+g*k)*z)
        else:
            rho.append(alpha**float(b+g*q(k,tag))*w*z)
    rho[root]*=alpha**float(d)*remaining
    # Include the identical signed structural/history coefficient.
    original*=(-.31)
    rel=abs(math.prod(rho)/original-1)
    allocation_max_relative=max(allocation_max_relative,rel)
    check(rel<2e-12,'literal allocated product')

# Exact mixed recurrence, including repeated hits at the same occurrence.
for trial in range(1500):
    Gamma=F(1,3)
    r=rng.randrange(2,7)
    packets=[]
    for i in range(r):
        N=rng.randrange(1,7)
        tags=[rng.randrange(2) for _ in range(N)]
        ks=[rng.randrange(0,10) for _ in range(N)]
        gs=[F(rng.randrange(1,9),24) for _ in range(N)]
        bs=[F(rng.randrange(1,4),28) for _ in range(N)]
        aa=F(20)
        H=sum(g*q(k,t) for g,k,t in zip(gs,ks,tags))
        psi=aa-sum(bs)-H-2*Gamma
        packets.append([aa,ks,tags,gs,bs,psi])
    bridge_D=[]
    for child in range(1,r):
        parent=rng.randrange(child)
        delta=F(0)
        for pi in (parent,child):
            pp=packets[pi]; v=rng.randrange(len(pp[1])); k=pp[1][v]
            delta+=pp[3][v]*(q(k+1,pp[2][v])-q(k,pp[2][v]))
            pp[1][v]+=1
        bridge_D.append(delta)
    child_psi=sum(p[0] for p in packets)-sum(sum(p[4]) for p in packets)-sum(sum(g*q(k,t) for g,k,t in zip(p[3],p[1],p[2])) for p in packets)-2*Gamma
    rhs=sum(p[5] for p in packets)+sum(2*Gamma-D for D in bridge_D)
    check(child_psi==rhs,'mixed recurrence exact')
    check(child_psi>=sum(p[5] for p in packets),'mixed recurrence superadditive')

# Conditional side omission, with heterogeneous root types.
for trial in range(1000):
    Gamma=F(1,3); roots=[]; retained=[]
    for h in range(rng.randrange(2,7)):
        N=rng.randrange(1,9); aa=F(rng.randrange(10,21))
        gs=[F(rng.randrange(1,8),24) for _ in range(N)]
        bs=[F(rng.randrange(1,4),28) for _ in range(N)]
        ks=[rng.randrange(8) for _ in range(N)]; ts=[rng.randrange(2) for _ in range(N)]
        d=aa-sum(bs)-sum(g*q(k,t) for g,k,t in zip(gs,ks,ts))
        roots.append(d)
        # Every root is retained; arbitrary side occurrences may be absent.
        for v in range(N):
            if v==0 or rng.randrange(2): retained.append((bs[v],gs[v],ks[v],ts[v]))
    aa=sum(roots)+sum(b+g*q(k,t) for b,g,k,t in retained)
    dd=aa-sum(b for b,g,k,t in retained)-sum(g*q(k,t) for b,g,k,t in retained)
    check(dd==sum(roots),'side omission keeps entire amplitude')
    check(dd-2*Gamma==sum(d-2*Gamma for d in roots)+2*Gamma*(len(roots)-1),'heterogeneous conditional reserve')

# Nodewise two-power sharpness on an endpoint panel, independently of old allocation.
sharp=[]
for m in [1,3,8,16]:
  vals=[]
  for p in [4,8,12,16,20,24]:
    h=2.**(-p/2)
    rr=[x for x in rule(p+1,m,8) if x[3]==p]
    masses=[w for r,w,a,pp,leb in rr]
    weighted=[w*((1-r*r+h*h)/2)**(-3) for r,w,a,pp,leb in rr]
    panel=sum(masses)
    check(panel>0.05*h*h,'endpoint panel mass comparable to h squared')
    check(max(weighted)>=panel/m*max((1-r*r+h*h)/2 for r,w,a,pp,leb in rr)**(-3)*(1-1e-12),'pigeonhole lower bound')
    vals.append({'p':p,'mass_over_h2':panel/h**2,'two_power_scaled':max(weighted)*h**4,'false_four_power_scaled':max(weighted)*h**2})
  # A fictitious extra two powers fails exponentially as h shrinks.
  check(vals[-1]['false_four_power_scaled']>1e5*vals[0]['false_four_power_scaled'],'no second mass credit')
  sharp.append({'m':m,'values':vals})

# Exact advertised grades, cost, conservative floors and fixed-rank admission.
betas=[F(1,14),F(1,28),F(1,28)]; gammas=[F(1,8),F(1,7),F(1,5)]; Ws=[6,6,4]; Gamma=F(1,5)
d=[8-8*b-g*w for b,g,w in zip(betas,gammas,Ws)]
psi=[x-2*Gamma for x in d]; root=[b+x for b,x in zip(betas,d)]
check(psi==[F(879,140),F(226,35),F(228,35)],'reserve grades')
check(root==[F(27,4),F(193,28),F(139,20)],'root grades')
check([root[i]-gammas[i] for i in range(3)]==[F(53,8),F(27,4),F(27,4)],'caller grades')
check(psi[1]-psi[0]==F(5,28) and psi[2]-psi[1]==F(2,35),'successive reserve gains')
check(8+gammas[2]==9-4*gammas[2]==F(41,5),'heat old mixing balance')
check(F(41,5)-8+6*gammas[2]==F(7,5),'conservative floor')
check(F(7,5)+gammas[2]==F(8,5),'VALUE floor')
T=[0,1]; V=[0,1]
for n in range(2,9):
    T.append(sum(i*(n-i)*T[i]*T[n-i] for i in range(1,n)))
    V.append(sum((n-i)*(i+1)*V[i]*T[n-i]+i*(n-i+1)*T[i]*V[n-i] for i in range(1,n)))
check(T[8]==794880 and V[8]==27020800,'independent history VALUE recurrence')
check(19*V[8]==513395200,'two-shell raw bill')
for n in range(4,101):
    g=F(1,n-3); b=F(1,4*(n-1)); rr=n-F(1,4)-(n-4)*g; P=n+g
    check(rr-g==n-F(5,4),'rank caller grade')
    check(2*rr>P and 2*n-(2*n-7)*g>P,'rank admission own self')

inputs=[ROOT/'CLOCK-MASS-TWO-STEP-REENTRY.md', ROOT/'WEIGHTED-GENERATOR-CONTRACT.md',
Path('/workspace/shared/induction-heat-defect-reentry-20261005/HEAT-SHELL-REENTRY.md'),
Path('/workspace/shared/induction-mixed-queue-20261005/WEIGHTED-MIXED-QUEUE.md'),
Path('/workspace/shared/induction-mixed-queue-20261005/SIMULTANEOUS-FAMILY-JOIN.md'),
Path('/workspace/shared/induction-mixed-queue-20261005/COMMON-CARRIER-AMALGAMATION.md'),
Path('/workspace/shared/induction-native-rank-generator-20261005/ACTIVE-PROBE-GENERATOR.md'),
Path('/workspace/shared/induction-native-rank-generator-20261005/RANK8-SMOOTHED-NATIVE-RETURN.md'),
Path('/workspace/shared/induction-shared-variables-20261005/marked-spanning-tree-addendum/FROZEN-COEFFICIENT-WIDTH-ZERO-PORT.md')]
hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
result={'assertions':count,'status':'PASS','native_compiler_executed':False,'author_implementation_imported':False,'allocation_max_relative_error':allocation_max_relative,'max_multiplier_error_excess':max_excess,'sharpness':sharp,'input_sha256':hashes}
(HERE/'independent-checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('sharpness','input_sha256')},indent=2))
