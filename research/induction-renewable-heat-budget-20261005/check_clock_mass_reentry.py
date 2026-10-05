#!/usr/bin/env python3
"""Exact ledgers and finite diagnostics; no native compiler is executed."""
from fractions import Fraction as F
from pathlib import Path
import hashlib, itertools, json, math
import numpy as np
OUT=Path(__file__).resolve().parent

def compositions(total,n):
    if n==1:
        yield (total,); return
    for v in range(total+1):
        for rest in compositions(total-v,n-1): yield (v,)+rest

Gamma=F(1,5)
a=F(8)
betas=[F(1,14),F(1,28),F(1,28)]
gammas=[F(1,8),F(1,7),F(1,5)]
losses=[6,6,4]
psi=[a-8*b-g*k-2*Gamma for b,g,k in zip(betas,gammas,losses)]
roots=[a-7*b-g*k for b,g,k in zip(betas,gammas,losses)]
first=[r-g for r,g in zip(roots,gammas)]
assert psi==[F(879,140),F(226,35),F(228,35)]
assert [psi[1]-psi[0],psi[2]-psi[1]]==[F(5,28),F(2,35)]
assert roots==[F(27,4),F(193,28),F(139,20)]
assert first==[F(53,8),F(27,4),F(27,4)]
P=F(41,5); E=a-4*Gamma; J=a-5*Gamma
assert E+1==P and a+Gamma==P
assert E+J==F(71,5)>P
assert 2*min(roots)==F(27,2)>P
assert P-F(57,7)==F(2,35)
assert P-a+6*Gamma==F(7,5)
assert P-a+7*Gamma==F(8,5)

# Degree sequences of labelled rank-eight trees: all k>=0 with sum k=6
# are feasible because degree=k+1 and sum degree=14. Counting sequences
# counts degree vectors, not trees or histories.
degree_vectors=0; hits=0; maxW=0
for ks in compositions(6,8):
    W=sum(max(k-2,0) for k in ks); maxW=max(maxW,W)
    assert W<=4
    for v in range(8):
        kp=list(ks); kp[v]+=1
        assert sum(max(k-2,0) for k in kp)<=5
        assert 0<=max(kp[v]-2,0)-max(ks[v]-2,0)<=1
        hits+=1
    degree_vectors+=1
assert degree_vectors==1716 and maxW==4

# Actual positive dyadic GL rule, with common nodes for every R_j.
x,w=np.polynomial.legendre.leggauss(12)
rows=[]
for p in range(12):
    length=2.**(-p-1); lo=1-2*length; hi=1-length
    rr=(lo+hi)/2+length*x/2; ww=length*w/2
    assert abs(sum(ww)-length)<1e-14
    for r,base in zip(rr,ww):
        assert base>0 and base<=length<=1-r+1e-15
        rows.append((r,base))
max_weight_ratio=0.; max_multiplier_error=0.; weighted_tests=0
for order in range(1,18):
    for r,base in rows:
        omega=base*r**(order-1)
        assert 0<omega<=1-r*r+1e-14
        for heat in [.001,.013,.15,.7,1.]:
            t=math.sqrt((1-r*r+heat*heat)/2)
            for k in range(15):
                q=max(k-2,0)
                bound=2**(1+q/2)*heat**(-q)
                ratio=omega*t**(-k)/bound
                assert ratio<=1+1e-12
                max_weight_ratio=max(max_weight_ratio,ratio)
                weighted_tests+=1
    for ell in list(range(64))+[100,1000,10000,100000]:
        actual=sum(base*r**(ell+order-1) for r,base in rows)
        err=abs(actual-1/(ell+order))
        max_multiplier_error=max(max_multiplier_error,err)
        assert err<=2**(-12)+1e-8

# Both charge tags survive arbitrary finite repeated bank hits.
charge_tests=0
for oldtag,newtag in itertools.product([False,True],repeat=2):
    charge=lambda k,tag: max(k-2,0) if tag else k
    for ku,kv,hu,hv in itertools.product(range(12),range(12),range(1,5),range(1,5)):
        du=charge(ku+hu,oldtag)-charge(ku,oldtag)
        dv=charge(kv+hv,newtag)-charge(kv,newtag)
        assert 0<=du<=hu and 0<=dv<=hv
        charge_tests+=1

# Fixed-rank direct weighted return margins; these do not certify a
# complete endpoint or arbitrary-order recursion.
rank_ledgers=[]
for n in range(4,101):
    g=F(1,n-3); b=F(1,4*(n-1)); W=n-4
    r=F(n)-(n-1)*b-W*g
    e=F(n)-W*g; j=F(n)-(W+1)*g; p=F(n)+g
    ps=F(n)-n*b-(W+2)*g
    assert e+1==p and 2*r>p and e+j>p and r-g==F(n)-F(5,4)>1 and ps>0
    if n<=12: rank_ledgers.append(dict(n=n,gamma=str(g),P=str(p),root=str(r),Psi=str(ps)))

# Exact raw source census and no accidental extra clock dimension.
T=[0,1]; V=[0,1]
for n in range(2,9):
    T.append(sum(i*(n-i)*T[i]*T[n-i] for i in range(1,n)))
    V.append(sum((n-i)*(i+1)*V[i]*T[n-i]+i*(n-i+1)*T[i]*V[n-i] for i in range(1,n)))
assert T[8]==794880 and V[8]==27020800
assert 19*V[8]==513395200

# Two literal rank-eight scalar force-slot telescopes, sharing each
# stage's endpoints. This is exact algebra, not a native simulation.
rng=np.random.default_rng(20261005)
telescope_error=0.
for _ in range(200):
    old,middle,fine=rng.normal(size=(3,8))
    def shell(cold,warm):
        return sum((cold[m]-warm[m])*math.prod(cold[:m])*math.prod(warm[m+1:]) for m in range(8))
    actual=math.prod(old)+shell(middle,old)+shell(fine,middle)
    err=abs(actual-math.prod(fine))
    assert err<1e-12
    telescope_error=max(telescope_error,err)

result=dict(status='PASS finite algebra/numerical diagnostics; no native execution',
    Psi=list(map(str,psi)),gains=[str(psi[1]-psi[0]),str(psi[2]-psi[1])],
    roots=list(map(str,roots)),actual_root_firsts=list(map(str,first)),
    stage2=dict(heat=str(P),outside_old=str(E+1),self_bank=str(E+J),own=str(2*min(roots)),clock_floor='7/5',value_floor='8/5'),
    degree_vectors_checked=degree_vectors,one_hit_vectors_checked=hits,max_weighted_width_degree=maxW,
    weighted_source_tests=weighted_tests,max_weight_ratio=max_weight_ratio,positive_clock_nodes=len(rows),
    max_sampled_scalar_multiplier_error=max_multiplier_error,omitted_endpoint_length=2**(-12),
    mixed_charge_tests=charge_tests,rank_checks_through=100,selected_rank_ledgers=rank_ledgers,
    rank8_histories=T[8],rank8_base_values=V[8],rank8_two_shell_values=19*V[8],clocks_per_history=15,
    max_two_telescope_error=telescope_error,
    limitations=['No native compiler executed','No arbitrary-order heat restoration','Public-log constants retained','No whole-cumulant lower bound inferred from the star mass certificate'])
(OUT/'checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
