#!/usr/bin/env python3
"""Independent chord audit; never imports or executes the author's checker.

Separate scipy quadrature plus analytic coefficient extraction, literal VALUES,
and a translated original-gradient replay test. Numerical assertions are
regression diagnostics, not substitutes for the analytic remainder estimates.
"""
from pathlib import Path
import json, math
from fractions import Fraction
import numpy as np
from scipy.special import roots_legendre, roots_hermitenorm

OUT=Path(__file__).resolve().parent
nchecks=0
def ck(value, label):
    global nchecks
    if not value: raise AssertionError(label)
    nchecks += 1

def rule(kind):
    if kind=='dyadic_midpoint':
        x,v=roots_legendre(14)
        panels=[]
        for j in range(24):
            lo=1-2.**(-j);hi=1-2.**(-j-1)
            panels.append((lo+(hi-lo)*(x+1)/2,(hi-lo)*v/2))
        panels.append((np.array([1-2.**(-25)]),np.array([2.**(-24)])))
        r=np.concatenate([p[0] for p in panels]);w=np.concatenate([p[1] for p in panels])
    elif kind=='dyadic_exact_tail':
        # Same dyadic panels; tail is a 3-point Gaussian rule, exact through r^5.
        x,v=roots_legendre(14);panels=[]
        for j in range(24):
            lo=1-2.**(-j);hi=1-2.**(-j-1)
            panels.append((lo+(hi-lo)*(x+1)/2,(hi-lo)*v/2))
        xt,vt=roots_legendre(3);h=2.**(-24)
        panels.append((1-h+h*(xt+1)/2,h*vt/2))
        r=np.concatenate([p[0] for p in panels]);w=np.concatenate([p[1] for p in panels])
    else:
        x,v=roots_legendre(32);r=(x+1)/2;w=v/2
    return r,w,np.sqrt(1-r*r)

q=.25;eps=.1;k=1.;tau=math.exp(-k*k/2)
u=lambda x:q*x*x+eps*(np.sin(k*x)-k*x)
f=lambda x:2*q*x+eps*k*(np.cos(k*x)-1)
fp=lambda x:2*q-eps*k*k*np.sin(k*x)
ck(.4 <= 2*q-eps*k*k and 2*q+eps*k*k <= .6,'global potential Hessian bounds')
ck(f(0.)==0,'anchored force')

# A mathematically nonzero finite-tail floor that IEEE double cannot see.
h=Fraction(1,2**24)
m2_exact=Fraction(1,3)-h**3/Fraction(12)
ck(m2_exact != Fraction(1,3),'fixed midpoint rule is not exactly quadratic-exact')
first_floor=3*eps*k**3*tau*float(m2_exact-Fraction(1,3))
ck(first_floor<0,'tiny finite-tail linear cumulant coefficient is negative')

z,prob=roots_hermitenorm(96);prob/=math.sqrt(2*math.pi)
Z=z[:,None];G=z[None,:];P=prob[:,None]*prob[None,:]
E=lambda a:float(np.sum(P*a))
results={}
for kind in ('dyadic_midpoint','dyadic_exact_tail','global_gauss'):
    r,w,s=rule(kind);beta=float(w@s);lam=1.5-2*beta*beta
    m1=float(w@r);m2=float(w@r**2);m4=float(w@r**4)
    m11=float(w@(r*s));m31=float(w@(r**3*s))
    K=m2+2*beta*m11;J=m4+2*beta*m31
    target1=eps*k**3*tau;target2=q*eps*tau*(k**5-6*k**3)
    pred1=3*eps*k**3*tau*m2
    pred2=3*q*eps*tau*(lam*(k**5*J-3*k**3*m2)-2*k**3*K)
    delta=pred2-target2
    ck(abs(m1-.5)<1e-14 and abs(sum(w)-1)<1e-14,'mass and first moment')
    ck(0<=lam<=1 and min(w)>0,'positive rule and chord range')
    H=np.zeros((len(z),len(z)));Hz=np.zeros_like(H)
    for a,b,c in zip(r,s,w):
        H+=c*f(a*Z+b*G);Hz+=c*a*fp(a*Z+b*G)
    # Obtain coefficients directly from third centered moments of the map
    # Z+A*(-H)+A²*(lambda Hz H), without using the simplified identities.
    L=-H;R=lam*Hz*H
    direct1=3*E((Z*Z-1)*L)
    direct2=3*E((Z*Z-1)*R)+3*E(Z*L*L)-6*E(L)*E(Z*L)
    ck(abs(direct1-pred1)<2e-13,'first cumulant from direct moment expansion')
    ck(abs(direct2-pred2)<2e-13,'second cumulant from direct moment expansion')
    # Tilt partition derivatives at A=0, independently from author formulas.
    U=u(z);eU=float(prob@U)
    raw=[]
    for power in range(4):
        a=z**power
        zer=float(prob@a)
        one=zer*eU-float(prob@(a*U))
        two=float(prob@(a*U*U))/2-float(prob@(a*U))*eU+zer*(eU*eU-float(prob@(U*U))/2)
        raw.append((zer,one,two))
    t1=raw[3][1]-3*raw[1][1]
    t2=raw[3][2]-3*raw[1][2]-3*raw[1][1]*raw[2][1]
    ck(abs(t1-target1)<3e-13 and abs(t2-target2)<3e-13,'independent partition coefficient extraction')
    rows=[]
    for A in (1/16,1/32,1/64,1/128,1/256):
        H1=np.zeros_like(H)
        for a,b,c in zip(r,s,w): H1+=c*f(a*(Z-A*H)+b*G)
        Y=Z-A*H-lam*A*(H1-H)
        ym=E(Y);yv=E((Y-ym)**2)
        pp=prob*np.exp(-A*U);pp/=pp.sum()
        tm=float(pp@z);tv=float(pp@(z-tm)**2);tt=float(pp@(z-tm)**3)
        repaired=tm+math.sqrt(tv/yv)*(Y-ym)
        rt=E((repaired-tm)**3)
        fourth=max(E((repaired-tm)**4),float(pp@(z-tm)**4))
        coeff=(rt-tt)/(A*A)
        ck(abs(E(repaired)-tm)<1e-13,'exact mean repair diagnostic')
        ck(abs(E((repaired-tm)**2)-tv)<1e-13,'exact variance repair diagnostic')
        ck(fourth<3.1,'fourth-moment numerical regression')
        ck(abs(coeff-delta)<.02*A,'literal VALUE map A² defect convergence')
        # Moment-based W2 floor is valid for any coupling with these moments.
        w2_floor=abs(rt-tt)/(3*math.sqrt(fourth))
        rows.append(dict(A=A,defect_over_A2=coeff,repair_factor=math.sqrt(tv/yv),fourth_max=fourth,W2_moment_floor=w2_floor))
        t=.5
        # Independent centered Gaussian N kills cubic cross terms exactly.
        by3=t**3*rt;bt3=t**3*tt
        by4=t**4*E((repaired-tm)**4)+6*t*t*(1-t*t)*tv+3*(1-t*t)**2
        bt4=t**4*float(pp@(z-tm)**4)+6*t*t*(1-t*t)*tv+3*(1-t*t)**2
        ck(abs((by3-bt3)/(A*A)-coeff/8)<1e-13,'bridge preserves scaled third defect')
        ck(max(by4,bt4)<3.1,'bridge fourth-moment numerical regression')
    results[kind]=dict(nodes=len(r),beta=beta,lambda_chord=lam,first_coefficient=pred1,second_coefficient=pred2,delta=delta,rows=rows)

# Matrix calibration retains one shared G at every node.
r,w,s=rule('global_gauss');beta=float(w@s);lam=1.5-2*beta*beta
rng=np.random.default_rng(941323)
quadratic=[]
for d in (1,3,7):
    Q,_=np.linalg.qr(rng.normal(size=(d,d)))
    B=Q@np.diag(rng.uniform(.01,.25,d))@Q.T;I=np.eye(d)
    # Polynomial coefficient audit, not merely an O(A³) tolerance.
    C2=(.25+beta*beta+lam/2)*(B@B)
    ck(np.linalg.norm(C2-B@B)<1e-14,'full matrix A² covariance coefficient')
    Z0=rng.normal(size=d);G0=rng.normal(size=d)
    hh=sum(c*(B@(a*Z0+b*G0)) for a,b,c in zip(r,s,w))
    hh1=sum(c*(B@(a*(Z0-hh)+b*G0)) for a,b,c in zip(r,s,w))
    literal=Z0-hh-lam*(hh1-hh)
    algebra=(I-B/2+lam*(B@B)/4)@Z0+(-beta*B+lam*beta*(B@B)/2)@G0
    ck(np.linalg.norm(literal-algebra)<1e-13,'literal matrix VALUE graph')
    quadratic.append(dict(d=d,pointwise_error=float(np.linalg.norm(literal-algebra))))

# Actual original-gradient VALUE/JVP/VJP ports with a moving captured anchor.
# Only first matrices at recorded VALUE points enter the replay.
vectors=np.array([[1.,0.],[.6,.8]])
def F(x): return .55*x+.08*np.sum((np.cos(vectors@x)-1)[:,None]*vectors,axis=0)
def DF(x): return .55*np.eye(2)-.08*sum(math.sin(float(v@x))*np.outer(v,v) for v in vectors)

def value_first(A,b,z0,g0):
    anchor=F(b);ra=math.sqrt(A);r,w,s=rule('global_gauss')
    def packet(x):
        hv=np.zeros(2);Jx=np.zeros((2,2));Jg=np.zeros((2,2));Jb=np.zeros((2,2));sites=[]
        for ri,wi,si in zip(r,w,s):
            site=b+ra*(ri*x+si*g0);fv=F(site);D=DF(site)
            hv+=wi*ra*(fv-anchor)
            Jx+=wi*A*ri*D;Jg+=wi*A*si*D
            Jb+=wi*ra*(D-DF(b));sites.append(site)
        return hv,Jx,Jg,Jb,sites
    H,Hz,Hg,Hb,sites0=packet(z0)
    H1,Hz1,Hg1,Hb1,sites1=packet(z0-H)
    ez=Hz1@(np.eye(2)-Hz)-Hz
    eg=Hg1-Hz1@Hg-Hg
    eb=Hb1-Hz1@Hb-Hb
    Y=z0-H-lam*(H1-H)
    D=np.concatenate((-Hb-lam*eb,np.eye(2)-Hz-lam*ez,-Hg-lam*eg),axis=1)
    return Y,D,len(sites0)+len(sites1)+1
ports=[]
for A in (.25,.1,.025):
    for j in range(7):
        b=rng.normal(size=2);z0=rng.normal(size=2);g0=rng.normal(size=2)
        Y,D,count=value_first(A,b,z0,g0);theta=np.r_[b,z0,g0]
        numeric=np.empty_like(D);step=1e-5
        for idx in range(6):
            e=np.zeros(6);e[idx]=step
            yp=value_first(A,*(theta+e).reshape(3,2))[0]
            ym=value_first(A,*(theta-e).reshape(3,2))[0]
            numeric[:,idx]=(yp-ym)/(2*step)
        err=float(np.linalg.norm(D-numeric))
        ck(err<2e-9,'recorded-site first equals literal VALUE finite difference')
        aa=rng.normal(size=6);vv=rng.normal(size=2)
        ck(abs(float(vv@(D@aa))-float(aa@(D.T@vv)))<1e-13,'JVP and adjoint duality')
        zero,Dzero,_=value_first(A,b,np.zeros(2),np.zeros(2))
        ck(np.array_equal(zero,np.zeros(2)),'same-site anchor gives literal private zero')
        ck(count==65,'two complete node banks and one captured anchor')
        phys=np.eye(2)+math.sqrt(A)*D[:,:2]
        ck(np.linalg.norm(phys-np.eye(2),2)<4*A,'physical caller first I plus O(A) regression')
        ports.append(dict(A=A,replay_error=err,original_value_sites=count))

# Whole-root gradient lift and the separate fixed-u projected curl.
lift_rows=[]
r,w,s=rule('global_gauss');LQ=float(w@(1/r));P0=np.concatenate((np.eye(2),np.zeros((2,2))),axis=1);Pi=P0.T@P0
for A in (.25,.08,.02):
    for trial in range(8):
        W=rng.normal(size=4)
        def lift(V):
            J=np.zeros(4);B=np.zeros((4,4));K=np.zeros(2);DK=np.zeros((2,4))
            for ri,wi,si in zip(r,w,s):
                R=np.concatenate((ri*np.eye(2),si*np.eye(2)),axis=1)
                fv=A*F(R@V);first=A*DF(R@V)
                J+=(wi/ri)*(R.T@fv);B+=(wi/ri)*(R.T@first@R)
                K+=wi*fv;DK+=wi*first@R
            return J,B,K,DK
        J0,B0,K0,DK0=lift(W);J1,B1,K1,DK1=lift(W-P0.T@K0)
        EJ=J1-J0;DEJ=B1-B0-B1@Pi@B0
        expected_skew=-B1@Pi@B0+B0@Pi@B1
        ck(np.linalg.norm(P0@J0-K0)<1e-13,'full lift projects to actual packet')
        ck(np.linalg.norm(P0@EJ-(K1-K0))<1e-13,'full chord lift projects to actual chord')
        ck(np.linalg.norm(DEJ-DEJ.T-expected_skew)<1e-13,'full lift exact noncommuting curl identity')
        ck(np.linalg.norm(DEJ-DEJ.T,2)<=2*A*A+1e-13,'whole-root lift quadratic curl bound')
        # Freezing u does not lose small curl: K_v differences are symmetric.
        Kv0=DK0[:,2:];Kv1=DK1[:,2:];Ku1=DK1[:,:2]
        Ev=Kv1-Kv0-Ku1@Kv0
        evskew=-Ku1@Kv0+Kv0@Ku1
        ck(np.linalg.norm(Ev-Ev.T-evskew)<1e-13,'fixed-u projected exact small-curl identity')
        ck(np.linalg.norm(Ev-Ev.T,2)<=A*A*beta+1e-13,'fixed-u projected quadratic curl bound')
        lift_rows.append(dict(A=A,LQ=LQ,full_curl=float(np.linalg.norm(DEJ-DEJ.T,2)),fixed_u_curl=float(np.linalg.norm(Ev-Ev.T,2))))

# Exact mean/variance Gaussianization misses a first-order bounded Lipschitz test.
gaussianization=[]
for t in (.25,.5):
    expected=-eps*math.exp(-1)*(math.sinh(t)-t)
    A=1/1024;pp=prob*np.exp(-A*u(z));pp/=pp.sum()
    tm=float(pp@z);tv=float(pp@(z-tm)**2)
    target=math.exp(-(1-t*t)/2)*float(pp@np.sin(t*z))
    reference=math.exp(-((1-t*t)+t*t*tv)/2)*math.sin(t*tm)
    observed=(target-reference)/A
    ck(abs(observed-expected)<.003*abs(expected),'bounded sine-test first defect')
    gaussianization.append(dict(t=t,predicted=expected,observed=observed))

report=dict(status='PASS',assertions=nchecks,
    midpoint_tail=dict(width=str(h),m2_exact=str(m2_exact),linear_cumulant_floor=first_floor,
        repair='Use an A-linked tail, retain a quadrature floor, or replace tail midpoint by a positive 3-point Gauss rule.'),
    scalar=results,matrix=quadratic,replay=ports,lift=lift_rows,gaussianization=gaussianization,
    scope='Independent diagnostics of a finite VALUE chord obstruction only; no all-order impossibility or general near-gradient-consumer admission.')
(OUT/'chord_independent_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(dict(status=report['status'],assertions=nchecks,midpoint_linear_floor=first_floor,
    defects={key:val['delta'] for key,val in results.items()},max_replay_error=max(x['replay_error'] for x in ports)),indent=2))
