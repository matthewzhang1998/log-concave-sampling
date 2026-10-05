#!/usr/bin/env python3
"""Independent diagnostics; imports no author checker or completed mean compiler."""
from __future__ import annotations
import hashlib, json, math
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.signal import fftconvolve
from scipy.special import gammaln, gammaincc


# Publication-only source locator and original/public pin mapping.
import sys
_ARCHIVE=next(p for p in Path(__file__).resolve().parents if (p/'publication_sources.py').is_file())
sys.path.insert(0,str(_ARCHIVE))
from publication_sources import source_path, source_label, verify_source_pin, pin_report

HERE=Path(__file__).resolve().parent
COUNTS={}; DETAILS={}
def check(x, group, message):
    COUNTS[group]=COUNTS.get(group,0)+1
    if not bool(x): raise AssertionError(f'{group}: {message}')
def close(x,y,tol,group,message): check(np.max(np.abs(np.asarray(x)-np.asarray(y)))<=tol,group,message)
def a2(x): return x/2-x*x/4
def step(u):
    u=np.asarray(u); v=np.clip(u,0,1)
    return 10*v**3-15*v**4+6*v**5
def h(u):
    u=np.asarray(u); v=np.clip(u,0,1)
    return v**6-3*v**5+2.5*v**4+np.maximum(u-1,0)
def chi_pdf(x,D):
    if x<=0: return 0.
    return math.exp((D-1)*math.log(x)-x*x/2-(D/2-1)*math.log(2)-gammaln(D/2))
def tail_moment(D,s,R,which='hinge'):
    f=(lambda u: u*u) if which=='hinge' else (lambda u: float(h(u))**2)
    val,err=quad(lambda u:f(u)*chi_pdf((R+u)/s,D)/s,0,np.inf,epsabs=1e-30,epsrel=3e-10,limit=300)
    return val,err

def symbolic():
    l,A=sp.symbols('lambda A', real=True)
    V=sp.Matrix([[1,sp.Rational(1,2),sp.Rational(1,4)],
                 [sp.Rational(1,2),sp.Rational(1,2),sp.Rational(3,8)],
                 [sp.Rational(1,4),sp.Rational(3,8),sp.Rational(3,8)]])
    row=sp.Matrix([[1,-l,0]])
    check(sp.expand((row*V*row.T)[0])==1-l+l*l/2,'symbolic','genuine future variance')
    check(V.det()>0,'symbolic','joint endpoint/history covariance positive')
    C=V[1:3,1:3]-V[1:3,0:1]*V[0:1,1:3]
    check(C==sp.Matrix([[sp.Rational(1,4),sp.Rational(1,4)],[sp.Rational(1,4),sp.Rational(5,16)]]),'symbolic','conditional backbone covariance')
    row=sp.Matrix([[l,-l*l]])
    check(sp.expand((row*C*row.T)[0])==l*l/4-l**3/2+5*l**4/16,'symbolic','same resummed b2')
    x=sp.symbols('x',real=True)
    p=10*x**3-15*x**4+6*x**5
    check(sp.factor(sp.diff(p,x))==30*x*x*(x-1)**2,'smooth_step','monotone polynomial step')
    for k in range(3):
        for endpoint,want in [(0,0),(1,1 if k==0 else 0)]:
            check(sp.diff(p,x,k).subs(x,endpoint)==want,'smooth_step','matching derivatives')
    H=sp.integrate(p,(x,0,x))
    check(H==x**6-3*x**5+sp.Rational(5,2)*x**4,'smooth_step','exact integrated transition')
    check(sp.integrate(1-p,(x,0,1))==sp.Rational(1,2),'smooth_step','kappa=one half')
    alpha,beta=A/3,2*A/3
    p0,p1=alpha/2-alpha**2/4,beta/2-beta**2/4
    q=(alpha**2*p0+beta**2*p1)/(alpha**2+beta**2)
    min_per_D=sp.simplify((alpha**2*(q-p0)**2+beta**2*(q-p1)**2)/8)
    check(sp.simplify(min_per_D-A**4*(2-A)**2/(36**2*10))==0,'scalar_obstruction','exact least square minimum')
    DETAILS['symbolic']={'sigma1_squared':'1-lambda+lambda^2/2','conditional_covariance':[['1/4','1/4'],['1/4','5/16']],
        'linear_obstruction_minimum':'A^2 (2-A) sqrt(D)/(36 sqrt(10))'}

def linear_checks():
    rng=np.random.default_rng(83021)
    max_ratio=0.; max_bias=0.
    for A in [1e-4,.003,.02,.1,.5]:
        for D in [1,2,8,101]:
            k=rng.uniform(0,A,size=D)
            for lam in np.linspace(0,A,9):
                b=k-lam; e0=np.linalg.norm(b); sigma=math.sqrt(1-lam+lam*lam/2)
                L=max(lam,A-lam); e1=sigma*e0; bound=e1+(lam+L)*e0
                c=2*lam*b+b*b
                var=b*b/2-3*b*c/4+3*c*c/8
                check(np.min(var)>=-1e-15,'linear_RMS','exact force difference has nonnegative variance')
                err=math.sqrt(max(0.,sum(var)))
                check(err<=bound+1e-14,'linear_RMS','full correlated R2 bound')
                bias=np.linalg.norm(k*(a2(lam)-a2(k)))/2
                check(bias<=A*bound+1e-14,'linear_RMS','canonical m3 RMS target bound')
                check(np.max(np.abs(b))<=L+1e-15,'linear_RMS','global residual Lipschitz constant')
                if bound: max_ratio=max(max_ratio,err/bound); max_bias=max(max_bias,bias/(A*bound))
        alpha,beta=A/3,2*A/3
        q=(alpha**2*a2(alpha)+beta**2*a2(beta))/(alpha**2+beta**2)
        lam=1-math.sqrt(1-4*q)
        check(alpha-1e-14<=lam<=beta+1e-14,'scalar_obstruction','optimizing scalar lies in interval')
        got=math.sqrt((alpha**2*(a2(lam)-a2(alpha))**2+beta**2*(a2(lam)-a2(beta))**2)/8)
        expected=A*A*(2-A)/(36*math.sqrt(10))
        close(got,expected,1e-14,'scalar_obstruction','minimum value numeric check')
    # A Gaussian variance contraction need not lower a general residual RMS.
    A=.5; lam=A/2; eta=A/2; width=.1; s=math.sqrt(1-lam+lam*lam/2)
    e=lambda v:eta*v/(1+2*v*v/(width*width))**.75
    check(e(s)>e(1),'scope','epsilon1 cannot generically be replaced by epsilon0')
    DETAILS['linear_checks']={'max_R2_to_bound_ratio':max_ratio,'max_m3_bias_to_bound_ratio':max_bias,
        'contracted_Gaussian_residual_RMS_counterexample_ratio':e(s)/e(1),
        'counterexample':'r(x)=eta x exp(-x^2/(2 width^2)), lambda=eta=A/2, width=0.1'}

def tails():
    rows=[]
    for A in [.5,.2,.05,.01]:
        lam=d=A/2; s=math.sqrt(1-lam+lam*lam/2); T=math.sqrt(8*math.log(1/A))
        for D in [1,5,100,10000]:
            R=math.sqrt(D)+T
            e=[]
            for scale in [1.,s]:
                hinge,err=tail_moment(D,scale,R)
                exact,_=tail_moment(D,scale,R,'h')
                check(exact<=hinge*(1+1e-8)+1e-25,'tail_certificate','h is bounded by hinge')
                check(hinge<=2*math.exp(-T*T/2)*(1+1e-8)+1e-25,'tail_certificate','dimension uniform concentration certificate')
                coarse=scale*scale*D*gammaincc(D/2+1,R*R/(2*scale*scale))
                check(hinge<=coarse*(1+1e-8)+1e-25,'tail_certificate','incomplete gamma second moment certificate')
                e.append(d*math.sqrt(exact))
            check(e[1]<=e[0]*(1+1e-8)+1e-20,'tail_certificate','monotone radial fixture epsilon1<=epsilon0')
            cert=A**3/math.sqrt(2)
            check(max(e)<=cert*(1+1e-8),'tail_certificate','epsilon bound A^3/sqrt2')
            bias=A*(e[1]+A*e[0]); grade=A**4*math.sqrt(D)
            check(bias<=((1+A)/math.sqrt(2))*grade*(1+1e-8),'tail_certificate','order-four grade for every D>=1')
            rows.append({'A':A,'D':D,'R':R,'epsilon0_numeric':e[0],'epsilon1_numeric':e[1],
                'epsilon_each_certificate':cert,'bias_to_A4sqrtD':bias/grade,
                'old_bounded_certificate_to_grade':A*d*(R+.5)*(1+A)/grade})
            for rho in [0,R/2,R,R+.01,R+.5,R+1,R+10,1e6]:
                u=max(0,rho-R); radial=lam+d*float(step(u))
                tangential=lam+(d*float(h(u))/rho if rho else 0)
                check(lam-1e-14<=radial<=A+1e-14,'tail_Hessian','radial eigenvalue within interval')
                check(lam-1e-14<=tangential<=A+1e-14,'tail_Hessian','tangential eigenvalue within interval')
            test_rho=R+10
            old_epsilon=d*(test_rho-float(h(test_rho-R)))
            close(old_epsilon,d*(R+.5),1e-11,'tail_Hessian','exact best bounded-remainder amplitude at slope A')
    DETAILS['tail_fixture']=rows

def ancestry_diagnostics():
    lam=.25
    errors=[]
    for dt in [.2,.1,.05,.025,.0125]:
        q=math.exp(-dt)
        covariance=1/(1+q)
        variance=(1+q*q)/(1+q)**2
        sigma2=1-2*lam*covariance+lam*lam*variance
        errors.append(abs(sigma2-(1-lam+lam*lam/2)))
    check(all(errors[i+1]<errors[i] for i in range(len(errors)-1)),'ancestry','genuine geometric future covariance converges')
    check(abs((1+lam*lam/2)-(1-lam+lam*lam/2)-lam)<1e-15,'ancestry','independent ancestry has a different marginal')
    rng=np.random.default_rng(10052026)
    samples=7000; n=200; dt=.04; A=.5; d=.25; R=.3
    q=math.exp(-dt); w=q**np.arange(n); w=w/w.sum()
    X=np.empty((samples,2*n-1)); X[:,0]=rng.normal(size=samples)
    for j in range(1,2*n-1): X[:,j]=q*X[:,j-1]+math.sqrt(1-q*q)*rng.normal(size=samples)
    residual=lambda x:d*np.sign(x)*h(np.abs(x)-R)
    g=lambda x:lam*x+residual(x)
    F1=fftconvolve(g(X),w[None,::-1],mode='valid',axes=1)
    H1=fftconvolve(X,w[None,::-1],mode='valid',axes=1)
    F2=np.sum(w*g(X[:,:n]-F1),axis=1)
    backbone=lam*(X[:,:n]@w)-lam*lam*(X@np.convolve(w,w))
    explicit=-lam*((F1-lam*H1)@w)+residual(X[:,:n]-F1)@w
    algebra=np.max(np.abs(F2-backbone-explicit))
    check(algebra<1e-13,'ancestry','exact finite genuine-ancestry residual decomposition')
    Cov=q**np.abs(np.arange(n)[:,None]-np.arange(n)[None,:])
    sigma=math.sqrt(1-2*lam*(w@(q**np.arange(n)))+lam*lam*(w@Cov@w))
    e0=d*math.sqrt(tail_moment(1,1,R,'h')[0]); e1=d*math.sqrt(tail_moment(1,sigma,R,'h')[0])
    L=A/2; bound=e1+(lam+L)*e0
    rms=lambda y:float(np.sqrt(np.mean(y*y)))
    delta=rms(F1[:,0]-lam*H1[:,0]); err=rms(F2-backbone)
    # Monte Carlo observations are diagnostics, not exact expectation proofs.
    check(delta<e0,'ancestry','sampled first residual stays comfortably below analytical RHS')
    check(err<bound,'ancestry','sampled second residual stays comfortably below analytical RHS')
    DETAILS['finite_ancestry']={'samples':samples,'nodes_per_future_window':n,'dt':dt,
        'max_decomposition_error':float(algebra),'sigma_discrete':sigma,'sigma_continuum':math.sqrt(1-lam+lam*lam/2),
        'epsilon0_exact_quadrature':e0,'epsilon1_discrete_exact_quadrature':e1,
        'sampled_Delta1_RMS':delta,'sampled_R2_RMS':err,'R2_bound':bound,
        'geometric_marginal_errors':errors,'qualification':'Finite-window Monte Carlo diagnostic; continuum proof is in report.'}

def raw_graph_checks():
    from numpy.polynomial.legendre import leggauss
    u,w=leggauss(7); t=(u+1)/2; w=w/2; c=np.sqrt(1-t*t); beta=w@c
    rng=np.random.default_rng(97201); max_fd=0.; max_comm=0.
    for A in [.01,.1,.4]:
        D=3; lam=d=A/2; R=math.sqrt(D)+math.sqrt(8*math.log(1/A))
        aa=a2(lam); bb=math.sqrt(lam*lam/4-lam**3/2+5*lam**4/16); alpha=1-aa
        def value(x):
            rho=np.linalg.norm(x)
            return lam*x+(d*float(h(rho-R))*x/rho if rho else 0)
        def Hess(x):
            rho=np.linalg.norm(x)
            if rho==0: return lam*np.eye(D)
            tang=lam+d*float(h(rho-R))/rho
            rad=lam+d*float(step(rho-R)); v=x/rho
            return tang*np.eye(D)+(rad-tang)*np.outer(v,v)
        sites=[]
        def query(x): sites.append(np.array(x)); return value(x)
        def raw(z,g,n,kind):
            out=np.zeros(D)
            for ti,ci,wi in zip(t,c,w):
                x=ti*z+ci*g
                if kind=='B': out+=wi*query(x)
                elif kind=='F': out+=wi*query(alpha*x-bb*n)
                else: out+=wi*(query(alpha*x-bb*n)-query(x))
            return out
        for scale in [.1,1,3,10]:
            Z,G,N=rng.normal(size=(3,D))*R*scale
            sites.clear(); F=raw(Z,G,N,'F')
            check(len(sites)==len(t),'raw_counts','F uses one original leaf per node')
            sites.clear(); B=raw(Z,G,N,'B')
            check(len(sites)==len(t),'raw_counts','B uses one original leaf per node')
            sites.clear(); E=raw(Z,G,N,'E')
            check(len(sites)==2*len(t),'raw_counts','E uses two original leaves per node')
            close(E,F-B,1e-13,'raw_graph','literal same-record split')
            JG=np.zeros((D,D)); JN=JG.copy(); JZ=JG.copy()
            xs=t[:,None]*Z+c[:,None]*G
            for ti,ci,wi,x in zip(t,c,w,xs):
                H1,H0=Hess(alpha*x-bb*N),Hess(x)
                max_comm=max(max_comm,np.linalg.norm(H1@H0-H0@H1))
                diff=alpha*H1-H0
                JG+=wi*ci*diff; JN-=wi*bb*H1; JZ+=wi*ti*diff
            close(JG,JG.T,1e-14,'raw_graph','G block is symmetric without Hessian commutation')
            check(np.linalg.norm(JG,2)<=A*beta+1e-14,'raw_graph','G first bound')
            check(np.linalg.norm(JN,2)<=A*bb+1e-14,'raw_graph','N first bound')
            check(np.linalg.norm(JZ,2)<=A/2+1e-14,'raw_graph','captured first bound')
            J=np.hstack((JG,JN)); lift=np.vstack((J,np.zeros_like(J)))
            close(np.linalg.norm(lift-lift.T,2),np.linalg.norm(JN,2),1e-14,'raw_graph','full lift curl equals N block')
            check(np.linalg.norm(J,2)<=A*math.sqrt(beta*beta+bb*bb)+1e-14,'raw_graph','complete private first')
            amp=A*(aa*(w@np.linalg.norm(xs,axis=1))+bb*np.linalg.norm(N))
            check(np.linalg.norm(E)<=amp+1e-13,'raw_graph','chord energy envelope')
            v=rng.normal(size=(3,D)); v/=np.linalg.norm(v); dz,dg,dn=v
            eps=1e-5
            fd=(raw(Z+eps*dz,G+eps*dg,N+eps*dn,'E')-raw(Z-eps*dz,G-eps*dg,N-eps*dn,'E'))/(2*eps)
            err=np.linalg.norm(fd-JZ@dz-JG@dg-JN@dn); max_fd=max(max_fd,err)
            check(err<1e-8,'raw_graph','actual directional first sweep')
            zero=np.zeros(D); sites.clear(); E0=raw(Z,zero,zero,'E')
            check(len(sites)==2*len(t),'raw_counts','caller-only origin fully billed')
            check(np.linalg.norm(E0)<=A*aa*np.linalg.norm(Z)/2+1e-13,'raw_graph','caller-only origin amplitude')
            close(raw(zero,zero,zero,'E'),zero,0,'raw_graph','literal complete zero')
    check(max_comm>1e-12,'raw_graph','tail fixture has genuinely noncommuting local Hessians')
    DETAILS['raw_finite_graph']={'max_directional_first_error':max_fd,'max_Hessian_commutator':max_comm,
        'nodes':len(t),'qualification':'Original radial tail fixture, fixed scalar parameters; complete mean compilers imported.'}

def main():
    pins=[
        (HERE.parent/'GAUSSIAN-RMS-RESIDUAL-RESUMMED-M3.md','18640f9c6c567acba4fe830c7b0d64394924887ee4335a7a49f57f817435dccf'),
        (HERE.parent.parent/'resummed-linear-backbone/BOUNDED-REMAINDER-RESUMMED-M3-VALUE-CONSUMER.md','798d20b87124e7069e465ce4e7bcf92c5af3734d86c8a62aeb288b5fec8eaef3'),
        (HERE.parent.parent/'resummed-linear-backbone/independent-audit/INDEPENDENT-BOUNDED-RESUMMED-M3-AUDIT.md','49d135ba6692c54012acd9d0d94f603c0015fca98b5d1be1164457ff45d1213a')]
    for p,digest in pins: check(hashlib.sha256(p.read_bytes()).hexdigest()==digest,'provenance','frozen source/import pin: '+p.name)
    symbolic(); linear_checks(); tails(); ancestry_diagnostics(); raw_graph_checks()
    sources=[HERE.parent/'GAUSSIAN-RMS-RESIDUAL-RESUMMED-M3.md']
    result={'status':'PASS','assertions':sum(COUNTS.values()),'assertions_by_group':COUNTS,'details':DETAILS,
        'author_files':[{'path':source_label(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sources],
        'audited_imports':[{'path':source_label(p),'sha256':digest} for p,digest in pins],
        'limits':['No author checker imported.','No completed gradient/near-gradient compiler is executed.',
        'Numerical checks supplement the written proof; sampled estimates are not theorem certificates.']}
    (HERE/'gaussian_rms_independent_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'assertions':result['assertions'],'groups':COUNTS,'ancestry':DETAILS['finite_ancestry']},indent=2))
if __name__=='__main__': main()
