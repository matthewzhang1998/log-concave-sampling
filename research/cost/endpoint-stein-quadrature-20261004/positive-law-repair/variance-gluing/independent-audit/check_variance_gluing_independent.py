#!/usr/bin/env python3
"""Independent finite diagnostics; not a time-order or compiler proof.

No author/upstream checker is imported. Source claims are checked against
direct Gaussian block conditioning, independently integrated nonlinear laws,
and literal small synthetic finite graphs. NumPy/SciPy only.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.hermite_e import hermegauss
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad_vec
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "EXACT-REVERSE-OU-BRIDGE-AND-BUFFER-LEDGER.md"
EXPECTED = "6846e7e565181a0d8acdc25d8d011b1155a4d350a2b2ca516f6d7cab24eea01a"
SEED = 610042114
rng = np.random.default_rng(SEED)
counts = {}
metrics = {}


def ck(group, condition, detail=""):
    counts[group] = counts.get(group, 0) + 1
    if not bool(condition):
        raise AssertionError(f"{group}: {detail}")


def close(group, actual, expected, tol=2e-10, detail=""):
    error = float(np.max(np.abs(np.asarray(actual) - np.asarray(expected))))
    ck(group, error <= tol, f"{detail}; error={error}, tol={tol}")
    metrics[group + "_max_absolute_comparison_error"] = max(
        metrics.get(group + "_max_absolute_comparison_error", 0), error)


def psd(group, mat, floor=-1e-10):
    ck(group, np.linalg.eigvalsh((mat + mat.T) / 2).min() >= floor)


def bridge(r, t):
    s2 = 1-r*r
    delta = t*t-r*r
    return (r*(1-t*t)/(t*s2), delta/(t*s2),
            (1-t*t)*delta/(t*t*s2), r/t, delta/t, delta/(t*t))


def sqrt_psd(mat):
    e, v = np.linalg.eigh(mat)
    return (v*np.sqrt(np.maximum(e, 0))) @ v.T


def gaussian_w2(m1, c1, m2, c2):
    s2 = sqrt_psd(c2)
    q = np.linalg.norm(m1-m2)**2 + np.trace(c1+c2-2*sqrt_psd(s2@c1@s2))
    return math.sqrt(max(0, q))


def quadratic_checks():
    group = "quadratic_block_conditioning"
    parameters = [(0., .1), (0., 1.), (.01, .04), (.2, .7),
                  (.5, math.sqrt(.25+.25)), (.8, .99), (.999, 1.)]
    for dim in [1, 2, 5, 9]:
        for A in [.02, .2, .5]:
            q, _ = np.linalg.qr(rng.normal(size=(dim, dim)))
            B = (q*rng.uniform(0, A, dim))@q.T
            I = np.eye(dim)
            sigma = np.linalg.inv(I+B)
            for r, t in parameters:
                z = rng.normal(size=dim)
                a,b,v,c,d,v0 = bridge(r,t)
                s2 = 1-r*r
                posterior = np.linalg.inv(I/s2+B)
                mux = posterior@(r*z/s2)
                m = B@mux
                H = B-B@posterior@B
                Ctt = t*t*sigma+(1-t*t)*I
                Crr = r*r*sigma+s2*I
                Ctr = (r/t)*Ctt
                direct_mean = Ctr@np.linalg.solve(Crr,z)
                direct_cov = Ctt-Ctr@np.linalg.solve(Crr,Ctr.T)
                close(group, a*z+b*mux, direct_mean, 5e-10)
                close(group, b*b*posterior+v*I, direct_cov, 5e-10)
                close(group, c*z-d*m, direct_mean, 5e-10)
                close(group, v0*I-d*d*H, direct_cov, 5e-10)
                close(group, H, B@np.linalg.inv(I+s2*B))
                close(group, posterior, s2*I-s2*s2*H)
                psd(group,H)
                psd(group,A*I-H)
                psd(group,direct_cov-v0*(1-A*(t*t-r*r))*I)
                ck(group, v >= -1e-14 and v0-d*d >= -1e-14)
                # The exact mean-only defect is matrix-valued.
                close(group,v0*I-direct_cov,d*d*H,5e-10)
                if t*t-r*r <= .25+1e-14:
                    reserve=v0*I-d*d*(I+H)
                    psd(group,reserve-.625*v0*I)
                    close(group,d*d*I+reserve,direct_cov,5e-10)
    metrics["quadratic_cases"] = 4*3*len(parameters)


# This potential is genuinely C2 but not C3 at x=shift. Its Hessian is
# A min(sqrt(abs(x-shift)),1). All definitions enforce g(0)=0 exactly.
def primitive1(u):
    u=np.asarray(u); a=np.abs(u)
    return np.sign(u)*np.where(a<=1, (2/3)*a**1.5, a-1/3)


def primitive2(u):
    a=np.abs(np.asarray(u))
    return np.where(a<=1, (4/15)*a**2.5, .5*a*a-a/3+.1)


def scalar_fields(x,A,shift=.37):
    g=A*(primitive1(x-shift)-primitive1(-shift))
    U=A*(primitive2(x-shift)-primitive2(-shift)-primitive1(-shift)*x)
    H=A*np.minimum(np.sqrt(np.abs(x-shift)),1)
    return U,g,H


def scalar_posterior(A,r,z,shift=.37,discrete=False):
    s2=1-r*r; s=math.sqrt(s2); a=r*z
    mode=brentq(lambda x:x+s2*scalar_fields(x,A,shift)[1]-a,
                a-8, a+8,xtol=1e-14)
    Um,gm,_=scalar_fields(mode,A,shift)
    def vector(y):
        x=mode+s*y
        U,g,H=scalar_fields(x,A,shift)
        w=np.exp(-y*y/2-(U-Um-gm*s*y))
        return w*np.array([1.,x,x*x,x**3,x**4,g,g*g,H,g*x])
    points=[(p-mode)/s for p in [shift-1,shift,shift+1]]
    points=sorted(p for p in points if -14<p<14)
    val,err=quad_vec(vector,-14,14,epsabs=2e-11,epsrel=2e-11,points=points)
    val/=val[0]
    ex,ex2,ex3,ex4,eg,eg2,eh,egx=val[1:]
    variance=ex2-ex*ex
    k3=ex3-3*ex*ex2+2*ex**3
    k4=ex4-4*ex*ex3+6*ex*ex*ex2-3*ex**4-3*variance**2
    result=dict(mean=ex,var=variance,k3=k3,k4=k4,m=eg,
                H=eh-(eg2-eg*eg),covgx=egx-eg*ex,error=float(err))
    if discrete:
        nodes,weights=leggauss(150)
        yy=[]; ww=[]
        borders=[-14]+points+[14]
        for lo,hi in zip(borders[:-1],borders[1:]):
            y=(lo+hi)/2+(hi-lo)*nodes/2
            x=mode+s*y; U,_,_=scalar_fields(x,A,shift)
            w=weights*(hi-lo)/2*np.exp(-y*y/2-(U-Um-gm*s*y))
            yy.extend(x); ww.extend(w)
        result["nodes"]=np.array(yy)
        result["weights"]=np.array(ww)/np.sum(ww)
    return result


def nonlinear_scalar_checks():
    group="nonlinear_C2_conditional_identities"
    max_k3=0.; max_k4=0.; max_covg=0.
    gn,gw=hermegauss(16); gw/=math.sqrt(2*math.pi)
    for A in [.03,.2,.5]:
        for r in [0.,.25,.8,.98]:
            for z in [-1.7,.6]:
                p=scalar_posterior(A,r,z,discrete=True); s2=1-r*r
                close(group,p["mean"],r*z-s2*p["m"],2e-9)
                close(group,p["var"],s2-s2*s2*p["H"],2e-9)
                close(group,p["covgx"],s2*p["H"],2e-9)
                ck(group,-1e-9<=p["H"]<=A+1e-9)
                step=5e-5
                derivative=(scalar_posterior(A,r,z+step)["m"]-
                            scalar_posterior(A,r,z-step)["m"])/(2*step)
                close(group,derivative,r*p["H"],2e-7)
                t=(r+1)/2
                a,b,v,c,d,v0=bridge(r,t)
                y=a*z+b*p["nodes"][:,None]+math.sqrt(v)*gn[None,:]
                w=p["weights"][:,None]*gw[None,:]
                mean=float(np.sum(w*y)); centered=y-mean
                variance=float(np.sum(w*centered**2))
                k3=float(np.sum(w*centered**3))
                k4=float(np.sum(w*centered**4)-3*variance**2)
                close(group,mean,c*z-d*p["m"],2e-8)
                close(group,variance,v0-d*d*p["H"],2e-8)
                close(group,k3,b**3*p["k3"],2e-8)
                close(group,k4,b**4*p["k4"],2e-8)
                max_k3=max(max_k3,abs(k3)); max_k4=max(max_k4,abs(k4))
                # E Dg alone omits a real same-law covariance term.
                U,g,H=scalar_fields(p["nodes"],A)
                covg=float(p["weights"]@(g*g)-(p["weights"]@g)**2)
                max_covg=max(max_covg,covg)
    ck(group,max_k3>1e-5 and max_k4>1e-5)
    ck(group,max_covg>1e-3)
    metrics.update(nonlinear_scalar_cases=24,max_nonlinear_bridge_third_cumulant=max_k3,
                   max_nonlinear_bridge_fourth_cumulant=max_k4,
                   max_retained_scalar_covariance_of_force=max_covg)


V=np.array([[1.,.1],[.55,.9],[-.25,.8]])
BIAS=np.array([.2,-.5,.4])
GRAM=np.linalg.norm(V,2)**2


def matrix_fields(x,A):
    u=x@V.T+BIAS
    logcosh=lambda a:np.logaddexp(a,-a)-math.log(2)
    U=A/GRAM*np.sum(logcosh(u)-logcosh(BIAS)-(x@V.T)*np.tanh(BIAS),axis=-1)
    g=A/GRAM*((np.tanh(u)-np.tanh(BIAS))@V)
    hs=A/GRAM*np.einsum("...k,ki,kj->...ij",1-np.tanh(u)**2,V,V)
    return U,g,hs


def matrix_posterior(A,r,z,n=70):
    nodes,weights=hermegauss(n); weights/=math.sqrt(2*math.pi)
    grid=np.stack(np.meshgrid(nodes,nodes,indexing="ij"),axis=-1).reshape(-1,2)
    weights=np.outer(weights,weights).reshape(-1)
    x=r*z+math.sqrt(1-r*r)*grid
    U,g,hs=matrix_fields(x,A)
    weights*=np.exp(-U); weights/=sum(weights)
    mx=weights@x; mg=weights@g
    covx=np.einsum("n,ni,nj->ij",weights,x-mx,x-mx)
    covg=np.einsum("n,ni,nj->ij",weights,g-mg,g-mg)
    H=np.einsum("n,nij->ij",weights,hs)-covg
    return mx,mg,covx,H,covg


def nonlinear_matrix_checks():
    group="noncommuting_nonlinear_matrix"
    max_comm=0.
    for A in [.1,.5]:
        for r in [0.,.45,.9]:
            z=np.array([1.1,-.7])
            mx,m,covx,H,covg=matrix_posterior(A,r,z)
            fine=matrix_posterior(A,r,z,100)
            for lhs,rhs in zip((mx,m,covx,H,covg),fine):
                close(group,lhs,rhs,2e-8)
            s2=1-r*r; I=np.eye(2)
            close(group,mx,r*z-s2*m,2e-8)
            close(group,covx,s2*I-s2*s2*H,2e-8)
            psd(group,H); psd(group,A*I-H)
            derivative=[]
            for j in range(2):
                dz=np.eye(2)[j]*1e-4
                derivative.append((matrix_posterior(A,r,z+dz)[1]-
                                   matrix_posterior(A,r,z-dz)[1])/(2e-4))
            close(group,np.array(derivative).T,r*H,2e-8)
            _,_,h1=matrix_fields(np.array([1.,-.5]),A)
            _,_,h2=matrix_fields(np.array([-.7,1.3]),A)
            max_comm=max(max_comm,np.linalg.norm(h1@h2-h2@h1,2))
    ck(group,max_comm>1e-4)
    metrics["max_noncommuting_Hessian_commutator_norm"]=max_comm


def density_unnormalized(A,r,z):
    s=math.sqrt(1-r*r)
    def f(y):
        U,g,_=scalar_fields(r*z+s*y,A)
        w=math.exp(-y*y/2-float(U))
        return np.array([w,w*g])
    points=[(x-r*z)/s for x in [-.63,.37,1.37]]
    val,_=quad_vec(f,-14,14,epsabs=2e-11,epsrel=2e-11,
                   points=[p for p in points if -14<p<14])
    rho=math.exp(-z*z/2)*val[0]/(2*math.pi)
    return rho,val[1]/val[0]


def density_checks():
    group="density_score_continuity_and_FP"
    for A in [.1,.5]:
        for r in [.2,.6,.95]:
            for z in [-1.2,.4]:
                h=1e-4
                rho,m=density_unnormalized(A,r,z)
                rp,mp=density_unnormalized(A,r,z+h)
                rm,mm=density_unnormalized(A,r,z-h)
                dt=(density_unnormalized(A,r+h,z)[0]-
                    density_unnormalized(A,r-h,z)[0])/(2*h)
                current=(rp*mp-rm*mm)/(2*h)
                dx=(rp-rm)/(2*h)
                dxx=(rp-2*rho+rm)/(h*h)
                bp=-(z+h)-(1+r)*mp; bm=-(z-h)-(1+r)*mm
                fp=-(bp*rp-bm*rm)/(2*h)+dxx
                close(group,dx/rho,-z-r*m,2e-7)
                close(group,dt,current,3e-7)
                close(group,dt,fp,3e-7)


def buffer_and_graph_checks():
    group="buffer_positivity_and_actual_graph"
    beta=.25; dim=3
    # A concrete private tape of 4D coordinates. It is deliberately not two
    # roots. This tests affine chain rules, not the imported mean compiler.
    J=np.zeros((dim,4*dim)); J[:,:dim]=np.eye(dim)
    L=np.zeros_like(J); L[:,2*dim:3*dim]=np.eye(dim)
    max_gap_ratio=0.; max_first_ratio=0.
    for r in np.linspace(0,1,17):
        for h in [1e-6,.001,.03,.1,.25,math.log(5/3)]:
            c=math.exp(-h); d=(1+r)*(1-c); q2=1-c*c-d*d
            ck(group,q2>=-1e-14)
            ck(group,q2+1e-14 >= (1-c)*(5*c-3))
            if h<=.25:
                ck(group,q2+1e-14>=9*h/16 and q2<=2*h+1e-14)
                x=rng.normal(size=dim); w=rng.normal(size=4*dim)
                n=rng.normal(size=dim); q=math.sqrt(max(q2,0))
                # M=-beta*r*x+Jw+beta*tanh(Lw), a complete explicit graph.
                M=lambda xx,ww:-beta*r*xx+J@ww+beta*np.tanh(L@ww)
                stage=lambda xx,ww,nn:c*xx-d*M(xx,ww)+q*nn
                Dcaller=(c+d*beta*r)*np.eye(dim)
                Dprivate=np.hstack([-d*(J+beta*np.diag(1-np.tanh(L@w)**2)@L),q*np.eye(dim)])
                close(group,np.linalg.norm(Dcaller,2),c+d*beta*r)
                ck(group,np.linalg.norm(Dcaller,2)<=(1+c)/2+1e-14)
                ck(group,np.linalg.norm(Dprivate,2)<=math.sqrt(d*d*(1+beta)**2+q2)+1e-14)
                vx=rng.normal(size=dim); vw=rng.normal(size=4*dim); vn=rng.normal(size=dim)
                eps=1e-5
                fd=(stage(x+eps*vx,w+eps*vw,n+eps*vn)-
                    stage(x-eps*vx,w-eps*vw,n-eps*vn))/(2*eps)
                close(group,fd,Dcaller@vx+Dprivate@np.r_[vw,vn],3e-9)
                max_gap_ratio=max(max_gap_ratio,q2/h)
                max_first_ratio=max(max_first_ratio,np.linalg.norm(Dprivate,2)/math.sqrt(h))
    # Threshold is a uniform-in-r sufficient condition, even outside r+h<=1.
    # A value slightly above it fails at r=1; this is not a claim of optimality
    # when the admissible r+h<=1 restriction is also imposed.
    h=math.log(5/3)+1e-5; c=math.exp(-h)
    ck(group,1-c*c-4*(1-c)**2<0)
    metrics.update(max_numerical_reserve_variance_over_h=max_gap_ratio,
                   max_actual_fresh_first_over_sqrt_h=max_first_ratio)


def law_gluing_checks():
    group="conditional_law_gluing_and_negative_controls"
    dim=4; I=np.eye(dim); e=.04; bias=np.linspace(-.02,.02,dim)
    for r,t in [(0.,.4),(.2,.4),(.6,.7),(.8,.9)]:
        a,b,v,c,d,v0=bridge(r,t)
        B=np.diag(np.linspace(.02,.4,dim))
        H=B@np.linalg.inv(I+(1-r*r)*B)
        reserve=v0*I-d*d*(I+H)
        psd(group,reserve)
        S=sqrt_psd(reserve)
        m=rng.normal(size=dim); z=rng.normal(size=dim)
        target_mean=c*z-d*m; target_cov=v0*I-d*d*H
        # Both independently supplied services have nonzero errors.
        actual_mean=target_mean-d*bias
        actual_cov=d*d*(1+e)**2*I+(S+e*I)@(S+e*I)
        epsM=math.sqrt(np.linalg.norm(bias)**2+dim*e*e)
        epsR=e*math.sqrt(dim)
        ck(group,gaussian_w2(actual_mean,actual_cov,target_mean,target_cov)
           <=d*epsM+epsR+2e-12)
        ck(group,np.linalg.norm(actual_cov-target_cov)>1e-4)
        ck(group,np.linalg.norm(actual_mean-target_mean)>0)
        # A correct marginal law does not allow retaining the service root.
        actual_joint=np.block([[I,I],[I,I]])
        forbidden_joint=np.block([[I,np.zeros_like(I)],[np.zeros_like(I),I]])
        ck(group,np.linalg.norm(actual_joint-forbidden_joint)>1)
    # Heterogeneous, random entering callers: the reference remains Gaussian
    # conditional on the caller although the rotation depends on that caller.
    n=30000; caller=rng.normal(size=n); w=rng.normal(size=(n,4)); fresh=rng.normal(size=(n,2))
    theta=caller
    G=np.c_[np.cos(theta)*w[:,0]-np.sin(theta)*w[:,1],
            np.sin(theta)*w[:,0]+np.cos(theta)*w[:,1]]
    residual=e*np.tanh(w[:,2:])/math.sqrt(2)
    mean=np.c_[.07*caller,.03*caller]
    M=mean+G+residual
    r=.4; h=.13; c=math.exp(-h); d=(1+r)*(1-c); q=math.sqrt(1-c*c-d*d)
    exposed=np.c_[caller,caller*caller]
    actual=c*exposed-d*M+q*fresh
    reference=c*exposed-d*mean-d*G+q*fresh
    close(group,actual-reference,-d*residual,2e-14)
    actual_squared_error=np.mean(np.sum((actual-reference)**2,axis=1))
    ck(group,actual_squared_error<=d*d*e*e)
    close(group,d*d+q*q,1-c*c)
    normalized=(-d*G+q*fresh)/math.sqrt(1-c*c)
    close(group,np.cov(normalized,rowvar=False),np.eye(2),.025)
    # Fully expanded synthetic operation census: every fresh semantic bank is
    # charged, even when the captured numerical caller is identical.
    value_cost=37; hvp_sites=37; private_dimension=4*dim
    records=[dict(caller="same",bank=j) for j in range(7)]
    ck(group,len({(q["caller"],q["bank"]) for q in records})==7)
    ck(group,sum(value_cost for q in records)==259)
    ck(group,sum(hvp_sites for q in records)==259)
    ck(group,sum(private_dimension+dim for q in records)==140)
    metrics["random_caller_observed_squared_coupling_error"]=float(actual_squared_error)
    metrics["random_caller_certified_squared_coupling_bound"]=d*d*e*e


def main():
    ck("source_integrity",hashlib.sha256(SOURCE.read_bytes()).hexdigest()==EXPECTED)
    quadratic_checks()
    nonlinear_scalar_checks()
    nonlinear_matrix_checks()
    density_checks()
    buffer_and_graph_checks()
    law_gluing_checks()
    result={"status":"PASS","seed":SEED,"source_sha256":EXPECTED,
            "assertions":sum(counts.values()),"assertions_by_group":counts,
            "metrics":metrics,
            "scope":"Finite diagnostics only. No author checkers imported. No numerical time-order, all-order nonlinear closure, original compiler execution, or c(P)/P claim.",
            "program_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    output=HERE/"variance_gluing_independent_checks.json"
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
