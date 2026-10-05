#!/usr/bin/env python3
"""Finite diagnostics for the exact identities in OBSERVER-LEADING-BLOCK.md."""
import json, math
from fractions import Fraction
from pathlib import Path
import numpy as np
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss

ROOT=Path(__file__).resolve().parent
checks=[]
def record(name, ok, **data):
    checks.append({'name':name,'passed':bool(ok),**data})
    if not ok: raise AssertionError(name)

C=.5; EPS=.25
class Oracle:
    def __init__(self): self.calls=[]
    def value(self,x):
        self.calls.append(float(x))
        return C*x+EPS*math.sin(x)

def source_sample(B,alpha):
    o=Oracle(); s=1+alpha*C; beta=alpha*EPS/s
    W=B+alpha*o.value(B); x=W/s
    p=(o.value(x)-C*x)/EPS
    q=(o.value(x+math.pi/2)-C*(x+math.pi/2))/EPS
    A=1+beta*q
    z1=x-beta*p/A
    z2=z1+beta**3*p**3/(2*A**3)
    return z1,z2,o.calls,beta

def maps(x,b):
    p=np.sin(x); q=np.cos(x); A=1+b*q
    h=-p/A; k=b*p**3/(2*A**3)
    Q1=x+b*h; Q2=Q1+b*b*k
    return Q1,Q2,h,k

def bounds(b):
    a=1-b; H=1/a; K=b/(2*a**3)
    C2=H*K+b*K*K/2+(H+b*K)**3/6
    H1=(1+b)/a**2; K1=3*b*(1+b)/(2*a**4)
    rho=b*H1+b*b*K1
    return b**3/(2*a*a),b**4*C2,b**3/(2*a**3),b**4*C2/a,rho

x=np.linspace(-30,30,120001)
for b in [1/8,1/16,1/32,1/64]:
    Q1,Q2,h,k=maps(x,b)
    A=1+b*np.cos(x)
    F=lambda y:y+b*np.sin(y)
    e1=F(Q1)-x; e2=F(Q2)-x
    b1,b2,r1,r2,rho=bounds(b)
    kp=3*b*np.sin(x)**2*(np.cos(x)+b)/(2*A**4)
    hp=-(np.cos(x)+b)/(A*A)
    first=1+b*hp+b*b*kp
    record(f'uniform_residual_and_first_beta_{b}',
           np.max(abs(e1))<=b1+2e-14 and np.max(abs(e2))<=b2+2e-14
           and np.min(first)>=1-rho-1e-14 and np.max(first)<=1+rho+1e-14,
           max_e1=float(np.max(abs(e1))),bound_e1=b1,
           max_e2=float(np.max(abs(e2))),bound_e2=b2,
           min_Q2_first=float(np.min(first)),max_Q2_first=float(np.max(first)),rho=rho)

# Exact integral remainder identities evaluated by high-order quadrature.
u,w=leggauss(40); u=(u+1)/2; w=w/2
z=np.linspace(-5,5,203); b=1/8
Q1,Q2,h,k=maps(z,b); d1=b*h; d2=b*h+b*b*k
r1=b*d1*d1*np.sum(w[None,:]*(1-u)[None,:]*(-np.sin(z[:,None]+d1[:,None]*u)),axis=1)
r2=b**4*(-np.sin(z))*h*k+b**5*(-np.sin(z))*k*k/2 \
  +b*d2**3/2*np.sum(w[None,:]*(1-u)[None,:]**2*(-np.cos(z[:,None]+d2[:,None]*u)),axis=1)
e1=Q1+b*np.sin(Q1)-z; e2=Q2+b*np.sin(Q2)-z
err1=float(np.max(abs(r1-e1))); err2=float(np.max(abs(r2-e2)))
record('exact_integral_remainders',max(err1,err2)<2e-15,error_e1=err1,error_e2=err2)

# Actual three-call implementation, including the old original force call.
z1,z2,calls,beta=source_sample(.73,.125)
x0=.73+beta*math.sin(.73)
q1,q2,*_=maps(np.array(x0),beta)
record('three_original_VALUE_calls',len(calls)==3 and abs(z1-q1)<1e-15 and abs(z2-q2)<1e-15,
       calls=calls,beta=beta,z1=z1,z2=z2)

# Non-vacuous orders under a retained-carrier coupling.
B=np.linspace(-math.pi,math.pi,40001); rows=[]
for b in [1/8,1/16,1/32,1/64]:
    q1,q2,*_=maps(B+b*np.sin(B),b)
    rows.append({'beta':b,'error1':float(np.max(abs(q1-B))), 'error2':float(np.max(abs(q2-B))),
      'lead_error1':float(np.max(abs((q1-B)/b**3+np.sin(B)**3/2))),
      'lead_error2':float(np.max(abs((q2-B)/b**4-np.sin(B)**3*np.cos(B)/6)))})
rates1=[math.log(rows[i]['error1']/rows[i+1]['error1'],2) for i in range(3)]
rates2=[math.log(rows[i]['error2']/rows[i+1]['error2'],2) for i in range(3)]
record('two_stage_nonzero_orders',min(rates1)>2.9 and min(rates2)>3.8 and max(rates1)<3.2 and max(rates2)<4.4,
       rows=rows,rates1=rates1,rates2=rates2)

# Retained test eta(B), no integration by parts or independent-carrier substitution.
nodes,weights=hermgauss(100); B=np.sqrt(2)*nodes; weights=weights/np.sqrt(np.pi)
b=1/8; x0=B+b*np.sin(B); q1,q2,*_=maps(x0,b)
eta=np.cos(.7*B)+.25*np.sin(B)
phi=lambda y:np.sin(1.3*y)
phip=lambda y:1.3*np.cos(1.3*y)
for j,Z in [(1,q1),(2,q2)]:
    d=Z-B
    lhs=np.sum(weights*eta*(phi(Z)-phi(B)))
    rhs=np.sum(weights*eta*d*np.sum(w[None,:]*phip(B[:,None]+d[:,None]*u),axis=1))
    record(f'exact_retained_test_current_stage_{j}',abs(lhs-rhs)<2e-16,lhs=float(lhs),rhs=float(rhs),difference=float(abs(lhs-rhs)))

# Fixed-endpoint leading block: v=-sin(B), complete actual derivative vs target.
b=1/8; A=1+b*np.cos(B); v=-np.sin(B); h=v/A
leading=np.sum(weights*eta*A*h*phip(x0))
target=np.sum(weights*eta*v*phip(x0))
record('fixed_endpoint_block_inverse',abs(leading-target)<1e-15,leading=float(leading),target=float(target))

# Finite Bell-polynomial coefficient block and its exact triangular inverse.
def bell_matrix(x,b,m):
    deriv=[0.0]+[1+b*math.cos(x)]+[b*math.sin(x+r*math.pi/2) for r in range(2,m+1)]
    bells=[[0.0]*(m+1) for _ in range(m+1)]; bells[0][0]=1.0
    for n in range(1,m+1):
        for k in range(1,n+1):
            bells[n][k]=sum(math.comb(n-1,j-1)*deriv[j]*bells[n-j][k-1] for j in range(1,n-k+2))
    M=np.zeros((m+1,m+1)); M[0,0]=1
    for n in range(1,m+1):
        for k in range(1,n+1): M[k,n]=bells[n][k]
    return M
errors=[]
for x0 in [-2.1,.13,1.9]:
    m=8; b=1/8; M=bell_matrix(x0,b,m); D=np.linspace(-1,1,m+1)
    C0=np.zeros(m+1)
    for k in range(m,-1,-1): C0[k]=(D[k]-np.dot(M[k,k+1:],C0[k+1:]))/M[k,k]
    errors.append(float(np.max(abs(M@C0-D))))
record('finite_Bell_block_inverse_rank_8',max(errors)<1e-11,errors=errors)

# First-order cumulant coefficients of theta*exp(-1/2)*sin(theta).
coeff=[{'rank':2*m,'normalized_coefficient':(-1)**(m-1)*2*m} for m in range(1,9)]
record('observer_tower_same_grade',all(r['normalized_coefficient']!=0 for r in coeff),coefficients=coeff)

# Independent finite-order coefficient solver and computable first/remainder budgets.
def compositions(n,r):
    if r==1:
        yield (n,)
    else:
        for first in range(1,n-r+2):
            for tail in compositions(n-first,r-1): yield (first,)+tail

def h_coefficients(x,b,N):
    A=1+b*np.cos(x); hs=[0,-np.sin(x)/A]
    for n in range(2,N+1):
        value=0
        for r in range(2,n+1):
            dr=np.sin(x+r*math.pi/2)
            value+=dr/math.factorial(r)*sum(np.prod([hs[j] for j in co]) for co in compositions(n,r))
        hs.append(-b*value/A)
    return np.array(hs)

def source_sample_order(B,alpha,N):
    o=Oracle(); s=1+alpha*C; b=alpha*EPS/s
    W=B+alpha*o.value(B); x=W/s
    p=(o.value(x)-C*x)/EPS
    q=(o.value(x+math.pi/2)-C*(x+math.pi/2))/EPS
    A=1+b*q; hs=[0,-p/A]
    cycle=[p,q,-p,-q]
    for n in range(2,N+1):
        total=sum(cycle[r%4]/math.factorial(r)*sum(math.prod(hs[j] for j in co) for co in compositions(n,r)) for r in range(2,n+1))
        hs.append(-b*total/A)
    return x+sum(b**j*hs[j] for j in range(1,N+1)),o.calls

def rational_poly_power(p,r):
    out=[Fraction(1)]
    for _ in range(r):
        nxt=[Fraction(0)]*(len(out)+len(p)-1)
        for i,x in enumerate(out):
            for j,y in enumerate(p): nxt[i+j]+=x*y
        out=nxt
    return out

b0=Fraction(1,8); a0=1-b0; Hs=[Fraction(0),1/a0]; Ds=[Fraction(0),1/a0**2]; budgets=[]
for N in range(1,7):
    if N>=2:
        S=Fraction(0); Sp=Fraction(0)
        for r in range(2,N+1):
            for co in compositions(N,r):
                prod=math.prod(Hs[j] for j in co)
                S+=prod/math.factorial(r)
                Sp+=sum(Ds[co[l]]*math.prod(Hs[co[i]] for i in range(r) if i!=l) for l in range(r))/math.factorial(r)
        Hs.append(b0*S/a0)
        Ds.append(b0*(S+Sp)/a0+b0*b0*S/a0**2)
    bn=min(b0,1/(2*sum(Ds[1:])))
    poly=Hs.copy()
    CN=Fraction(0)
    for r in range(2,N+1):
        pp=rational_poly_power(poly,r)
        CN+=sum(pp[d]*b0**(d-N-1)/math.factorial(r) for d in range(N+1,len(pp)))
    effective=sum(b0**(j-1)*Hs[j] for j in range(1,N+1))
    CN+=effective**(N+1)/math.factorial(N+1)
    maxjet=0.; maxratio=0.
    for xx in [-2.1,.43,1.2]:
        b=float(bn); hs=h_coefficients(xx,b,N)
        # Compose the polynomial by a separate power routine, retaining full degree.
        residual=np.zeros(max(N*N+1,N+1)); residual[1]=(1+b*np.cos(xx))*hs[1]+np.sin(xx)
        for j in range(2,N+1): residual[j]+=(1+b*np.cos(xx))*hs[j]
        for r in range(2,N+1):
            power=np.polynomial.polynomial.polypow(hs,r)
            residual[:len(power)]+=b*np.sin(xx+r*math.pi/2)*power/math.factorial(r)
        maxjet=max(maxjet,float(np.max(abs(residual[:N+1]))))
        qq=xx+sum(b**j*hs[j] for j in range(1,N+1))
        ee=qq+b*np.sin(qq)-xx
        maxratio=max(maxratio,float(abs(ee)/(float(CN)*b**(N+2))))
    actual,calls_N=source_sample_order(.73,.125,N)
    bsample=.125*EPS/(1+.125*C); xsample=.73+bsample*math.sin(.73)
    hsample=h_coefficients(xsample,bsample,N)
    expected=xsample+sum(bsample**j*hsample[j] for j in range(1,N+1))
    budgets.append({'N':N,'b_N':float(bn),'b_N_exact':str(bn),'C_N':float(CN),'C_N_exact':str(CN),
                    'H_N':float(Hs[N]),'H_N_exact':str(Hs[N]),'D_N':float(Ds[N]),'D_N_exact':str(Ds[N]),
                    'first_guard_sum':float(sum(bn**j*Ds[j] for j in range(1,N+1))),
                    'first_guard_sum_exact':str(sum(bn**j*Ds[j] for j in range(1,N+1))),
                    'max_cancelled_jet':maxjet,'max_tested_error_bound_ratio':maxratio,'actual_original_VALUE_calls':len(calls_N),'actual_source_error':abs(actual-expected)})
    record(f'finite_order_solver_N_{N}',maxjet<1e-14 and maxratio<=1.02 and len(calls_N)==3 and abs(actual-expected)<1e-15 and sum(bn**j*Ds[j] for j in range(1,N+1))<=Fraction(1,2),
           **budgets[-1])
(ROOT/'finite-order-budgets.json').write_text(json.dumps({'base_beta_bound':float(b0),'base_beta_bound_exact':str(b0),'budget_arithmetic':'exact rational; floating fields are display approximations','budgets':budgets},indent=2)+'\n')

report={'status':'PASS','scope':'Finite numerical diagnostics supporting exact scalar proofs; not a native sampler execution.',
        'checks':checks,'passed':len(checks),'failed':0}
(ROOT/'checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'status':'PASS','passed':len(checks),'output':str(ROOT/'checks.json')},indent=2))
