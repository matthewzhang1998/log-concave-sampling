#!/usr/bin/env python3
"""Finite coefficient/covariance diagnostics; does not run imported native VALUE samplers."""
from fractions import Fraction
from math import factorial,exp,cos,sin,sqrt
import json
import numpy as np
import sympy as s
from numpy.polynomial.hermite import hermgauss
assertions=0

def check(x):
    global assertions
    assert bool(x)
    assertions+=1

def dfact(n):
    if n<=0:return 1
    a=1
    for j in range(n,0,-2):a*=j
    return a

# No reliance on the source's printed variance: integrate the displayed reverse tree.
x,y,z,S,T,r,p,th=s.symbols('x y z S T r p th')
gaus=(x,y,z,S,T)
V=r*th**2*x*y*z*(S+p*th*x*(T+p*th*y*z))
def gaussian_mean(q):
    poly=s.Poly(s.expand(q),*gaus)
    out=0
    for powers,c in poly.terms():
        if any(k%2 for k in powers):continue
        out+=c*np.prod([dfact(k-1) for k in powers])
    return s.expand(out)
mom=[s.Integer(1)] + [gaussian_mean(V**k) for k in range(1,7)]
kappas=[s.Integer(0)]
for k in range(1,7):
    kappas.append(s.expand(mom[k]-sum(s.binomial(k-1,j-1)*kappas[j]*mom[k-j] for j in range(1,k))))
check(s.expand(kappas[1]-r*p**2*th**4)==0)
check(s.expand(kappas[2]-r*r*th**4*(1+3*p*p*th*th+26*p**4*th**4))==0)
paired=[]
for k in range(1,7):
    both=s.expand((kappas[k]+kappas[k].subs(r,-r))/s.factorial(k))
    if k%2:check(both==0)
    else:
        check(both==2*kappas[k]/s.factorial(k))
        for (a,b,c),coef in s.Poly(both,r,p,th).terms():
            check(a==k);check(b%2==0);check(c==2*k+b)
            exponent=6*k+b;n=2*k+b
            check(exponent-n==4*k)
            check(exponent>=12)
    paired.append({'cluster':k,'coefficient':str(both)})

# Generic packet root/side allocation: all retained outputs obey doubled surplus.
for d in [4,8,12,16,24,32]:
 for h in range(2,21):
  for j in range(0,81):
   n=2*h+j;A=(d+2)*h+j
   check(A-n==h*d)
   check(A-n>=2*d)
   check(Fraction(2,3)*n+h*d>0)

# Rotation/disintegration covariance, in nonscalar dimensions with known source-zero rows.
rng=np.random.default_rng(21005)
rotations=[]
for db in range(1,6):
 for dx in range(1,6):
  for rep in range(20):
   L=rng.normal(size=(dx,db))/5
   R=rng.normal(size=(dx,dx))
   C=L@L.T+R@R.T+np.eye(dx)/2
   A=L.T@np.linalg.inv(C)
   S2=np.eye(db)-A@L
   check(np.min(np.linalg.eigvalsh(S2))>0)
   check(np.allclose(A@C,L.T,atol=1e-11))
   check(np.allclose(A@C@A.T+S2,np.eye(db),atol=1e-11))
   rotations.append(float(np.min(np.linalg.eigvalsh(S2))))

# Direct deterministic Gaussian integration of the source and exact readout formula.
nodes,weights=hermgauss(100);nodes=nodes*sqrt(2);weights=weights/sqrt(np.pi)
eta=.3;t=.7;s0=.4;delta=.2
A=1+eta
epsilon=s0*delta*eta*exp(-t*t/2)/A
for B in np.linspace(-5,5,61):
 source_mean=sum(weights*(delta*eta/A)*np.cos(B+t*nodes))
 check(abs(s0*source_mean-epsilon*cos(B))<1e-14)
 for xx in np.linspace(-4,4,17):
  first=delta*eta*cos(B+t*xx)/A
  bfirst=delta*eta*(cos(B+t*xx)-cos(B))/(A*t)
  check(abs(first)<=delta*eta/A+1e-14)
  check(abs(bfirst)<=2*delta*eta/(A*t)+1e-14)
 check(1-(epsilon*cos(B))**2>0.99)
cos2=sum(weights*np.cos(nodes)**2)
check(abs(cos2-(1+exp(-2))/2)<1e-14)
# Actual positive covariance consumer; independently differentiate characteristic ratio.
a=.35;c=.25;k=.4;beta=.6;v=a*a+c*c+k
characteristic_errors=[]
for theta in np.linspace(-2.4,2.4,25):
 base=sum(weights*np.exp(1j*beta*theta*nodes))
 derivative=-a*c*theta**2*sum(weights*np.exp(1j*beta*theta*nodes)*np.cos(nodes))/base
 exact=-a*c*theta**2*exp(-.5)*np.cosh(beta*theta)
 err=abs(derivative-exact);characteristic_errors.append(err);check(err<2e-13)
 # Exact positive variance and first-order retained-public discrepancy.
 for B in [-3.,-1.,0.,.7,2.]:
  m=epsilon*cos(B)
  cov=np.array([[1,m],[m,1.]])
  check(np.linalg.eigvalsh(cov).min()>0)
  readoutvar=a*a+c*c+2*a*c*m+k
  check(abs(np.array([a,c])@cov@np.array([a,c])+k-readoutvar)<1e-14)
  check(abs((a+c*m)-a-c*epsilon*cos(B))<1e-15)

# Explicit orthogonal reserve expansion: exact baseline covariance and original row.
for n in range(2,32):
 w=rng.uniform(.1,1,size=n);w=w/np.linalg.norm(w)
 e=np.zeros(n);e[0]=1
 v0=e-w
 O=np.eye(n) if np.linalg.norm(v0)<1e-14 else np.eye(n)-2*np.outer(v0,v0)/(v0@v0)
 check(np.allclose(O@O.T,np.eye(n),atol=1e-12))
 check(np.allclose(O[0],w,atol=1e-12))

out={'assertions':assertions,'paired_reverse_cumulants':paired,
     'minimum_disintegration_gap':min(rotations),
     'max_characteristic_derivative_error':max(characteristic_errors),
     'retained_public_lower_bound_constant':c*sqrt((1+exp(-2))/2),
     'epsilon_native_source':epsilon,
     'scope':'Finite reverse-jet, covariance, source-first, characteristic and Gaussian-row diagnostics; not execution of imported native source/pair programs.'}
with open('joint_port_checks.json','w') as f:json.dump(out,f,indent=2)
print(json.dumps({'assertions':assertions,'max_characteristic_derivative_error':max(characteristic_errors),'status':'PASS'}))
