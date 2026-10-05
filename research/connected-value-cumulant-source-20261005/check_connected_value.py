#!/usr/bin/env python3
"""Exact finite-copy cumulant checks and bounded-Hessian witness diagnostics.
The mathematical claims and limitations are stated in the companion note.
"""
import json, math
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy.integrate import quad

OUT=Path(__file__).resolve().parent

def partitions(n):
    def go(i, blocks):
        if i==n:
            yield tuple(tuple(b) for b in blocks); return
        for j in range(len(blocks)):
            blocks[j].append(i); yield from go(i+1,blocks); blocks[j].pop()
        blocks.append([i]); yield from go(i+1,blocks); blocks.pop()
    yield from go(0,[])

rows=[]
S={(0,0):1}
for n in range(1,9):
    counts={k:0 for k in range(1,n+1)}
    for p in partitions(n): counts[len(p)]+=1
    for k in range(1,n+1):
        S[n,k]=S.get((n-1,k-1),0)+k*S.get((n-1,k),0)
        assert counts[k]==S[n,k]
    cn=sum(Fraction(counts[k]*math.factorial(k-1)**2*math.factorial(n-k),math.factorial(n)) for k in counts)
    rows.append(dict(rank=n,partitions=sum(counts.values()),off_diagonal_variance_factor=str(cn),decimal=float(cn)))
assert rows[3]['off_diagonal_variance_factor']=='10/3'
assert rows[7]['off_diagonal_variance_factor']=='46212/35'

# Direct enumeration of every physical-slot-to-replica map, ranks 2 through 5.
# Monomials are mutually orthogonal for distinct physical coordinate indices.
import itertools
for n in range(2,6):
    ss=Fraction(0)
    for assignment in itertools.product(range(n),repeat=n):
        k=len(set(assignment))
        falling=math.factorial(n)//math.factorial(n-k)
        c=Fraction((-1)**(k-1)*math.factorial(k-1),falling)
        ss+=c*c
    assert str(ss)==rows[n-1]['off_diagonal_variance_factor']

# Known-centered rank-four two-copy version: 8 orthogonal monomials, coefficients +/-1/2.
coeff={tuple([0]*4):Fraction(1,2),tuple([1]*4):Fraction(1,2)}
for inds in itertools.combinations(range(4),2):
    inds=set(inds); coeff[tuple(int(i not in inds) for i in range(4))]=Fraction(-1,2)
assert len(coeff)==8 and sum(c*c for c in coeff.values())==2

# Analytic Gram formula for g(x)=a*x+b*sin(x)+c*(1-cos(x)), including antipodal and
# correlated rows. Positive definiteness is proved in the note; eigenvalues are a diagnostic.
a,b,c=.5,.125,1/64
L=np.array([[1.,0.],[-1.,0.],[0.,1.],[0.,-1.],[2**-.5,2**-.5],[2**-.5,-2**-.5],[.4,.2]])
K=L@L.T; v=np.diag(K); vi=v[:,None]; vj=v[None,:]
E=np.exp(-(vi+vj)/2)
M=a*a*K+a*b*K*(np.exp(-vi/2)+np.exp(-vj/2))+b*b*E*np.sinh(K)+c*c*(1-np.exp(-vi/2)-np.exp(-vj/2)+E*np.cosh(K))
eig=np.linalg.eigvalsh(M)
assert eig.min()>0

# Conditional OU history: Y=int exp(-t)X_t dt, X_0=0, Var(Y)=1/4;
# Cov(Y,X_t)=t*exp(-t). Cumulant derivative at epsilon=0 for g_eps.
assert abs(quad(lambda t:t*np.exp(-2*t),0,np.inf)[0]-.25)<1e-12
witness=[]
for n in [4,8]:
    integral,err=quad(lambda t:t**(n-1)*np.exp(-n*t-(1-np.exp(-2*t))/2),0,np.inf,epsabs=1e-14)
    lo=math.exp(-.5)*math.factorial(n-1)/n**n
    hi=math.factorial(n-1)/n**n
    assert lo<integral<hi
    deriv=n*a**(n-1)*(-1)**(n//2-1)*integral
    assert deriv<0
    witness.append(dict(rank=n,I_n=integral,quadrature_reported_error=err,rigorous_lower_bound=lo,rigorous_upper_bound=hi,kappa_derivative_at_zero=deriv))

# Finite tensor rank/nuclear-norm inequality checks with random outer products.
rng=np.random.default_rng(20261005)
rank_checks=0
for D in [3,5,8]:
    for Mcount in [1,2,3]:
        V=rng.normal(size=(D,Mcount))
        for n in [3,4]:
            C=rng.normal(size=(Mcount,)*n)
            letters='abcdefghijklmnop'
            left=''.join(letters[:n])
            expr=left+','+','.join(chr(65+j)+letters[j] for j in range(n))+'->'+''.join(chr(65+j) for j in range(n))
            T=np.einsum(expr,C,*([V]*n))
            mat=T.reshape(D,-1)
            sv=np.linalg.svd(mat,compute_uv=False)
            assert np.count_nonzero(sv>1e-9)<=Mcount
            assert sum(sv)**2<=Mcount*np.sum(sv*sv)+1e-8
            rank_checks+=1

result=dict(status='PASS',scope='Exact combinatorics and diagnostics; analytical source exclusion is proved in the companion note. No all-source impossibility claim.',u_statistic_rows=rows,known_centered_rank4_variance_factor=2,antithetic_dictionary_gram_eigenvalues=eig.tolist(),conditional_ou_witness=witness,rank_nuclear_checks=rank_checks)
(OUT/'checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
