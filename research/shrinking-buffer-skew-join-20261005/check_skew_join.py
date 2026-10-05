#!/usr/bin/env python3
"""Finite diagnostics for the proof's new identities; no native compiler execution."""
import itertools, json, math
from fractions import Fraction
from pathlib import Path
import numpy as np
from numpy.polynomial.hermite import hermgauss

rng=np.random.default_rng(67320261005)
counts={}
def check(ok, family):
    assert bool(ok), family
    counts[family]=counts.get(family,0)+1

def sym3(T):
    return sum(T.transpose(p) for p in itertools.permutations(range(3)))/6

def root(C):
    d,U=np.linalg.eigh(C)
    return (U*np.sqrt(d))@U.T

# Exact conditional bulk geometry and strictly positive full variance shares.
for A in np.geomspace(1e-8,0.1,24):
    w=A;q=1-w;sig2=1-q*q
    for t in np.linspace(0,q,57):
        c2=1-t*t;s=1-q*q*t*t;v=c2*sig2/s
        check(v>=w*(1-1e-8) and v<=2*w*(1+1e-8),'bulk_geometry')
        shares=np.array([v/8]*4+[v/2])
        check(abs(shares.sum()-v)<1e-13 and shares.min()>0,'variance_shares')
        beta=q*c2/math.sqrt(s)
        check(abs(beta*beta+v-c2)<1e-12,'known_carrier_row')

# Cross-covariance Hilbert inequality for arbitrary finite joint laws.
for D in range(1,9):
    for _ in range(20):
        n=113;R=rng.normal(size=(n,D));Q=rng.normal(size=(n,D*D))
        Q+=R@rng.normal(size=(D,D*D))*0.4
        p=rng.uniform(size=n);p/=p.sum()
        R-=p@R;Q-=p@Q
        C=R.T@(p[:,None]*Q);CQ=Q.T@(p[:,None]*Q)
        lhs=np.sum(C*C);rhs=np.linalg.eigvalsh(CQ)[-1]*np.sum(p[:,None]*R*R)
        check(lhs<=rhs*(1+1e-11),'one_energy_covariance_inequality')

# Smooth Lipschitz scalar coherent fixture: exact high-degree GH integration.
x,wx=hermgauss(100);x=x*math.sqrt(2);wx=wx/math.sqrt(math.pi)
h=np.tanh(x)+0.2*np.logaddexp(x,-x)-0.2*math.log(2)
f=np.tanh(0.7*x)
h-=wx@h;f-=wx@f
for A in np.geomspace(1e-5,0.1,40):
    V=A*h;U=V+A*A*f;R=U-V
    L=1.2*A+0.7*A*A
    eR=math.sqrt(wx@(R*R))
    kapdiff=abs(wx@(U**3-V**3))
    check(kapdiff<=6*L*L*eR*(1+1e-12),'third_tensor_stability_fixture')
    for a,b in [(U,U),(V,U),(V,V)]:
        Q=a*b;var=wx@((Q-wx@Q)**2)
        check(var<=4*L**4*(1+1e-11),'quadratic_covariance_bound')
    check(kapdiff/A**4<10,'third_tensor_A4_scaling')

# Exact native three-channel target under common radius/output scaling.
for u in np.geomspace(1e-7,0.2,17):
    for N in [1,3,11,37]:
        cbar=abar=bbar=0.5/math.sqrt(N)
        for _ in range(12):
            A=u*0.0001;q=1-A;alpha=q*A/math.sqrt(u)
            sig=rng.uniform(0.03,0.7);weight=rng.uniform(1e-5,0.1)*sig*sig
            rho1=-4*weight*alpha/(cbar*abar*bbar*sig)
            c,a,b=math.sqrt(u)*np.array([cbar,abar,bbar])
            native_log3=c*a*b*rho1*alpha*alpha
            target_log3=-4*q**3*weight*A**3/sig
            check(abs(native_log3/target_log3-1)<2e-12,'cubic_exact_rescaling')
            check(abs(math.sqrt(u)*alpha-q*A)<1e-16,'actual_first_rescaling')
            check(abs(math.sqrt(u)*alpha**6/(q**6*A**6*u**(-2.5))-1)<2e-12,'feedback_v_power')
            check(abs(math.sqrt(u)*alpha**5/(q**5*A**5*u**(-2))-1)<2e-12,'prior_v_power')

# Exact anisotropic Gaussian conditional quadratic reference regression.
for D in range(2,7):
    for _ in range(30):
        n=4;bs=rng.uniform(0.01,0.1,n)
        M=rng.normal(size=(D,D));J=M@M.T/D+np.eye(D)*0.7
        C=J+sum(bs)*np.eye(D);Ci=np.linalg.inv(C)
        S=rng.normal(size=D)
        for beta2 in bs:
            N=rng.normal(size=(D,D,D));N=(N+N.swapaxes(1,2))/2
            mu=beta2*Ci@S;V=beta2*np.eye(D)-beta2**2*Ci
            raw=np.einsum('iab,ab->i',N,(np.outer(mu,mu)+V)/beta2**2-np.eye(D)/beta2)
            target=np.einsum('iab,ab->i',N,np.outer(Ci@S,Ci@S)-Ci)
            check(np.linalg.norm(raw-target)<2e-10,'anisotropic_reference_regression')

# Full-output symmetrization has zero first-order characteristic response.
def char_quad(C,eta2,N,theta):
    Rt=root(C);Ci=np.linalg.inv(C)
    B=np.einsum('i,iab,ac,bd->cd',theta,N,Ci,Ci)
    M=Rt@B@Rt;h=Rt@theta
    K=np.eye(len(theta),dtype=complex)-2j*M
    const=np.einsum('i,iab,ab',theta,N,Ci)
    return np.exp(-0.5*eta2*(theta@theta)-1j*const-0.5*np.log(np.linalg.det(K))-0.5*h@np.linalg.solve(K,h))
ratios=[]
for dim in [2,3,4]:
    for _ in range(20):
        M=rng.normal(size=(dim,dim));C=np.eye(dim)+M@M.T/dim
        T=rng.normal(size=(dim,dim,dim))*0.1;N=(T+T.swapaxes(1,2))/2
        K=sym3(N);D=N-K;theta=rng.normal(size=dim)*0.6
        check(np.linalg.norm(sym3(D))<1e-12,'zero_full_symmetry')
        check(abs(np.einsum('iab,i,a,b',D,theta,theta,theta))<1e-12,'leading_rank3_cancellation')
        errs=[]
        for lam in [0.02,0.01,0.005]:
            err=abs(char_quad(C,0.7,lam*N,theta)-char_quad(C,0.7,lam*K,theta))
            errs.append(err)
            check(err<=10*lam*lam,'positive_symmetry_feedback_quadratic')
        if errs[1]>1e-12:
            r=errs[2]/errs[1];ratios.append(r)
            check(0.20<r<0.30,'positive_symmetry_feedback_order')

# Grade ledger at w=A; all public-log and absolute floors stay separate.
terms={'prefix':(2,Fraction(3,2)),'prefix_commutation':(3,1),'endpoint':(3,Fraction(1,2)),
       'mean_intrinsic':(4,0),'skew_and_mean_gram':(5,Fraction(-3,2)),
       'kappa_mismatch':(5,-1),'mixed_K':(6,Fraction(-7,6)),
       'native_prior_b5':(6,-2),'native_and_reference_feedback':(7,Fraction(-5,2)),
       'K_mixture':(7,Fraction(-3,2))}
exponents={}
for k,(a,b) in terms.items():
    exponents[k]=str(Fraction(a)+b)
    check(Fraction(a)+b>=Fraction(7,2),'grade_ledger')
check(Fraction(2)+Fraction(3,2)==Fraction(5)-Fraction(3,2),'prefix_skew_balance')

out={'assertions':sum(counts.values()),'families':counts,'grade_exponents':exponents,
     'symmetry_feedback_ratio_range':[min(ratios),max(ratios)],
     'scope':'Algebra, Gaussian covariance inequalities, exact anisotropic regression, positive quadratic characteristic laws, buffer/radius powers and exponent diagnostics. Does not numerically execute imported native VALUE compilers or certify their guards for an unspecified A.'}
path=Path(__file__).with_name('skew_join_checks.json');path.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
