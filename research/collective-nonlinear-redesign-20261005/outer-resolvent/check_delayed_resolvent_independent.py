"""Independent algebra/numeric diagnostics; not a substitute for the analytic proof."""
from fractions import Fraction as Q
import json, math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from pathlib import Path
OUT=Path(__file__).parent

def rule(delta):
    K=math.ceil(math.log2(4/delta)); n=max(1,math.ceil(math.log(16/delta)/math.log(4)))
    z,p=leggauss(n); ts=[]; ws=[]
    for j in range(K):
        a=1-2.**(-j); b=1-2.**(-j-1)
        ts.extend((a+b)/2+(b-a)*z/2);ws.extend(p*(b-a)/2)
    h=2.**(-K);ts.append(1-h/2);ws.append(h)
    return np.array(ts),np.array(ws),K,n

qchecks=[]
for delta in (0.2,0.01,1e-4,1e-8):
    t,p,K,n=rule(delta)
    ks=np.unique(np.r_[np.arange(200),np.geomspace(200,1e11,1000).astype('int64')])
    errs=np.abs(np.exp(np.outer(ks,np.log(t)))@p-1/(ks+1.))
    proof=2*2.**(-K)+8*4.**(-n)
    assert p.min()>0 and t.min()>0 and t.max()<1
    assert abs(p.sum()-1)<1e-14 and abs(p@t-.5)<1e-14
    assert errs.max()<=proof<=delta
    qchecks.append(dict(delta=delta,nodes=len(t),exact_real_certificate=proof,
        sampled_multiplier_max=float(errs.max()),mass_error=float(abs(p.sum()-1)),first_moment_error=float(abs(p@t-.5))))

# Original smooth radial obstruction in the deterministic Gaussian-row limit.
def hinge(v,eta):
    if v<=-eta:return 0.
    if v>=eta:return v
    return (v+eta)**2/(4*eta)
def phi(r,A):
    cv=.5*(1-1/math.sqrt(2));tau=cv/2
    return .5*hinge(r-.5,A*A)+.5*hinge(r-(1-A*tau),A*A)
def radial_g(row,A):
    norm=np.linalg.norm(row)
    return np.zeros_like(row) if norm==0 else A*phi(norm,A)/norm*row
radial=[]
for A in (1e-2,1e-3,1e-4,1e-5,1e-6):
    w=math.sqrt(A);q=1-w;sig=math.sqrt(1-q*q)
    t,p,_,_=rule(A**.75)
    x=np.array([1.,0.,0.]); y=np.array([q,sig,0.]); ell=np.array([0.,0.,1.])
    rows=t[:,None]*y+np.sqrt(1-t*t)[:,None]*ell
    assert np.max(np.abs(np.sum(rows*rows,axis=1)-1))<1e-14
    B=q*sum((pj*radial_g(row,A) for pj,row in zip(p,rows)),np.zeros(3))
    u=x-B;terminal=radial_g(u-w*radial_g(u,A),A)
    # Genuine nested OU reduction, with actual retained linear histories v,w.
    vh=np.array([.5,.5,0.]);wh=np.array([.25,.5,.25])
    cU=phi(1,A)
    rr=math.sqrt(1-A*cU+.5*A*A*cU*cU)
    c2=phi(rr,A)/rr
    target=radial_g(x-A*c2*vh+A*A*c2*cU*wh,A)
    firstchaos=.5*(terminal[0]-target[0])
    radial.append(dict(A=A,nodes=len(t),outer_first_chaos_defect=float(firstchaos),
        over_A2=float(firstchaos/A**2),over_A3=float(firstchaos/A**3),
        over_A_11_4=float(abs(firstchaos)/A**2.75)))
assert all(abs(radial[i+1]['over_A2'])<abs(radial[i]['over_A2']) for i in range(len(radial)-1))

# Exact quadratic coefficient ledger, including optional cheap correction.
quadratic=[]
for A in (.08,.01,1e-4):
    w=math.sqrt(A);q=1-w
    d2=-(w+q*q/2)+.5
    d3=w*q*q/2-.25
    assert abs(d2+w*w/2)<1e-15
    corr2=w*w/2;corr3=.25-w*q*q/2
    assert abs(d2+corr2)<1e-15 and abs(d3+corr3)<1e-15
    assert corr2>0 and corr3>0
    quadratic.append(dict(A=A,inner_K2_defect=d2,inner_K3_defect=d3,correction_K2=corr2,correction_K3=corr3))

# Noncommuting C-infinity convex gradients: analytic firsts of literal finite graph.
a=np.array([1.,0.]);b=np.array([.6,.8]); eye=np.eye(2)
def source(z,A):
    return A*(z/2+(math.sin(a@z)*a+math.sin(b@z)*b)/8)
def H(z,A):
    return A*(eye/2+(math.cos(a@z)*np.outer(a,a)+math.cos(b@z)*np.outer(b,b))/8)
def inner(x,N,L,A,t,p):
    w=math.sqrt(A);q=1-w;sig=math.sqrt(1-q*q);Y=q*x+sig*N
    B=np.zeros(2);Bx=np.zeros((2,2));BN=Bx.copy();BL=Bx.copy()
    for tj,pj in zip(t,p):
        cj=math.sqrt(1-tj*tj);U=tj*Y+cj*L;HH=H(U,A)
        B+=q*pj*source(U,A)
        Bx+=q*q*pj*tj*HH;BN+=q*sig*pj*tj*HH;BL+=q*pj*cj*HH
    u=x-B;v=u-w*source(u,A);T=source(v,A)
    HT=H(v,A)@(eye-w*H(u,A))
    return T-source(x,A), HT@(eye-Bx)-H(x,A), -HT@BN, -HT@BL
rng=np.random.default_rng(41417); ports=[]
for A in (.08,.01):
    it,ip,_,_=rule(A**.75);ot,op,_,_=rule(A*A)
    maxfirst=maxcurl=maxcaller=0.
    for ktest in range(8):
        z,G,N,L=rng.normal(size=(4,2))
        Dpriv=np.zeros((2,6));Dz=np.zeros((2,2))
        for t,p in zip(ot,op):
            c=math.sqrt(1-t*t);e,dx,dn,dl=inner(t*z+c*G,N,L,A,it,ip)
            Dpriv+=p*np.hstack((c*dx,dn,dl));Dz+=p*t*dx
        lift=np.zeros((6,6));lift[:2,:]=Dpriv
        maxfirst=max(maxfirst,math.sqrt(2)*np.linalg.norm(Dpriv,2)/A)
        maxcurl=max(maxcurl,math.sqrt(2)*np.linalg.norm(lift-lift.T,2)/(A*A))
        maxcaller=max(maxcaller,math.sqrt(2)*np.linalg.norm(Dz,2)/A)
    assert maxfirst<2 and maxcurl<3 and maxcaller<1
    ports.append(dict(A=A,inner_nodes=len(it),outer_nodes=len(ot),sampled_normalized_first_over_A=maxfirst,
        sampled_normalized_curl_over_A2=maxcurl,sampled_uncaptured_caller_over_A=maxcaller))

# Uniform grade constant from the full displayed error, before outer/completion floors.
AA=np.geomspace(1e-18,1/12,1000)
ratios=[]
for A in AA:
    w=math.sqrt(A);q=1-w;sig=math.sqrt(1-q*q) if q<1 else math.sqrt(2*w-w*w)
    err=A**3+(2*math.sqrt(2)/3)*A*A*w**1.5+q*w*A**3+q*A*A*A**.75+q*q*A**3*(.5+1/sig)
    ratios.append(err/A**2.75)
assert max(ratios)<5

backup_const=Q(3,2)*(Q(11,7)*Q(101,96)**2*Q(7,4)+Q(13,24)*Q(101,96)*Q(21,20))+Q(13,24)*Q(101,96)+Q(13,8)
result=dict(status='passed',positive_ou_quadrature=qchecks,radial_C2_row_limit=radial,
    exact_quadratic_ledger=quadratic,noncommuting_source_ports=ports,
    delayed_grade_max_sampled_constant=max(ratios),grade_three_backup_rational_constant=str(backup_const),
    limitations=['Numerical multiplier checks do not replace the uniform ellipse proof.',
    'Radial rows test the exact source formula in the Gaussian-row limit, not an infinite-dimensional sampled gradient.',
    'Ports are analytic in the proof; sampled source checks are diagnostics only.',
    'No completed native compiler or finite-precision production implementation was executed.'])
(OUT/'delayed_resolvent_independent_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
