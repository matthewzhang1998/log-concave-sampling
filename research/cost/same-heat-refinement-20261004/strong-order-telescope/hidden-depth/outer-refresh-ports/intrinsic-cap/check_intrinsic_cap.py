#!/usr/bin/env python3
"""Finite algebra diagnostics for the fixed-subspace cap extension, not its proof."""
from fractions import Fraction as F
import hashlib, json
from pathlib import Path
import numpy as np

rng=np.random.default_rng(20261004)
checks=0
max_residual=0.0
witnesses=[]

def check(ok, text):
    global checks
    checks+=1
    if not bool(ok): raise AssertionError(text)

def near(a,b,text,tol=4e-10):
    global max_residual
    err=float(np.linalg.norm(a-b))
    max_residual=max(max_residual,err)
    check(err<=tol,text)

def op(a): return np.linalg.norm(a,2)
def sym(a): return (a+a.T)/2

def root(a):
    val, vec=np.linalg.eigh(sym(a))
    check(val.min()>0,"positive spectral gap")
    return (vec*np.sqrt(val))@vec.T

for n,d in [(7,1),(11,2),(17,3),(31,4)]:
    for rep in range(8):
        E,_=np.linalg.qr(rng.normal(size=(n,d)))
        P=E.T; Q=E@E.T; Qp=np.eye(n)-Q
        near(P@P.T,np.eye(d),'coisometry')
        near(Q@Q,Q,'projection')
        A=rng.normal(size=(d,n)); A*=.19/op(A)
        M=E@A
        near(Q@M,M,'fixed response output')
        check(np.linalg.norm(M,'fro')<=np.sqrt(d)*op(M)+1e-12,'one Hilbert response')
        C=M-M.T
        check(np.linalg.matrix_rank(C,tol=1e-10)<=2*d,'curl rank')
        check(np.linalg.norm(C,'fro')<=np.sqrt(2*d)*op(C)+1e-12,'curl Hilbert estimate')
        pair_cov=np.eye(n)-M@M.T
        projected_cov=Q@pair_cov@Q+Qp
        near(projected_cov,pair_cov,'conditional pair covariance fixed')
        x=rng.normal(size=n)
        near(Q@(M@x),M@x,'conditional mean fixed')
        near(Q@M,M,'incoming/output cross covariance fixed')
        mean=E@rng.normal(size=d)
        near(Q@mean,mean,'mean target fixed')
        near(Q@np.eye(n)@Q+Qp,np.eye(n),'mean covariance fixed')

        # Forward root perturbation and every matrix Taylor coefficient remain in S.
        B=M@M.T
        H=root(np.eye(n)-B)-np.eye(n)
        near(Q@H,H,'forward root perturbation in S')
        for j in range(1,5):
            Bj=np.linalg.matrix_power(B,j)
            near(Q@Bj,Bj,'root polynomial coefficient in S')
            check(np.linalg.norm(Bj,'fro')<=np.sqrt(d)*op(Bj)+1e-12,'root coefficient Hilbert')
        # Complete orientation cancellation; individual covariances have cross blocks.
        forward=np.eye(n)-sym(M@M)
        defect=M@M.T-sym(M@M)
        near(forward-defect,pair_cov,'ordered self/orientation cancellation')
        cross=float(np.linalg.norm(Qp@forward@Q,'fro'))
        check(cross>1e-6,'separate helper is not individually fixed')
        witnesses.append(cross)

        # A force-output tensor with all other legs flattened: only one d factor.
        T0=rng.normal(size=(d,3*n)); T0/=op(T0)
        T=E@T0
        check(np.linalg.norm(T,'fro')<=np.sqrt(d)+1e-12,'output-cut to full tensor norm')
        C1=rng.normal(size=(3*n,2*n)); C1/=op(C1)
        check(np.linalg.norm(T@C1,'fro')<=np.linalg.norm(T,'fro')+1e-12,'one-Hilbert contraction')
        # Root action after physical readout: no ambient norm penalty.
        R1=root(np.eye(n)-.2*sym(rng.normal(size=(n,n)))/np.sqrt(n))
        R2=np.eye(n)
        check(np.linalg.norm(P@(R1-R2),'fro')<=np.sqrt(d)*op(R1-R2)+1e-12,'projected root finite error')

# Explicit averaging warning: rank-one matrices with varying ranges average to full rank.
n=9
average=sum(np.outer(np.eye(n)[j],np.eye(n)[j]) for j in range(n))/n
check(np.linalg.matrix_rank(average)==n,'average rank is not generally preserved')
check(np.linalg.norm(average,'fro')<=1,'Hilbert averaging survives')

B=F(52,5); u=F(3,4); canonical=B*u; interior=2*canonical+1
check(canonical==F(39,5),'exact final mixed cap exponent')
check(interior==F(83,5),'exact intrinsic rank-four interior grade')
for c in [0,1,3,7,11]:
    check(2*(canonical-F(c,2))+1==interior-c,'ambient loss cannot be changed by child order')
# Simple sufficient strict cap margins, conservative e<=2 eta/2048.
eta=F(100,99); emax=2*eta/F(2048); cstar=F(199,200)
check(F(11)-9*emax>B,'force-eleven current margin')
check(F(12)-10*emax>B+cstar,'force-twelve direct margin')
check(F(16)>B+cstar,'sixteenth helper margin')
# At k>1/7 and c<=cstar, helper−(B−2+2k+c min(k,1))>0.
k=F(1,7)
check(F(8)+8*k-(B-2+2*k+cstar*k)>0,'frozen helper small-k margin')

# Public adaptation: only Exact32 is bundled. External pins are recorded,
# not counted as freshly verified by this portable diagnostic.
research_root = next(p for p in Path(__file__).resolve().parents if p.name == 'research')
sources={'Exact32':research_root/'exact-slack/exact32-proof-v2.tex'}
expected={
 'LOW30':'7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8',
 'LOW31':'3de62d349338304cce43af4e49c8862e8c260b1f5816234d63d6abfbc1d004ba',
 'Exact32':'c7487ca24ffa44ea73b52e42e4b6f715ad964171fa5bd0a61cf9323eb646cbff'}
pins={name:hashlib.sha256(path.read_bytes()).hexdigest() for name,path in sources.items()}
for name,h in pins.items(): check(h==expected[name],name+' pin')
external_pins={name:expected[name] for name in ('LOW30','LOW31')}
result={'status':'PASS','assertions':checks,'max_algebra_residual':max_residual,
        'min_individual_helper_cross_block':min(witnesses),
        'final_mixed_canonical_cap':str(canonical),'intrinsic_cap_interior':str(interior),
        'ambient_cap_interior':'83/5−c','source_sha256':pins,
        'external_source_pins_not_reverified':external_pins,
        'scope':'Finite algebra diagnostics only; not execution of the full marked family.'}
print(json.dumps(result,indent=2))
