#!/usr/bin/env python3
"""Bounded diagnostic checks for the same-heat short refresh; not a theorem oracle."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import numpy as np

OUT = Path(__file__).resolve().parent
v, J = F(19, 15), F(29, 10)
assert 1 + F(3, 2)*v == J
assert F(3, 2) + F(3, 2)*v == F(17, 5)
assert max(v/2, J-2) == F(9, 10)
assert 1-max(v/2,J-2) == F(1,10)
checks = {"exponents": {"v": str(v), "projection_mark": str(J), "short_quadrature_base": "17/5", "beta_power": "9/10", "full_radius_at_r_a": "1/10"}}

quadratic = []
for a in [.2, .05, .01, .001]:
    h = a**float(v/2)
    gamma = math.sin(h)**2
    eps = a**.9
    beta = 1/math.sqrt(1/gamma+1/eps**2)
    assert (25/36)*a**float(v) <= gamma <= a**float(v)
    assert 5/math.sqrt(61)*eps <= beta <= eps
    for lam in [.2, .7, 1.0]:
        om = math.sqrt(1+a*lam)
        c, s = math.cos(om*h), math.sin(om*h)/om
        target_var = a/(1+a*lam)
        actual_var = c*c*target_var+a*s*s
        invariant_residual = abs(actual_var-target_var)
        assert invariant_residual < 2e-15
        K = 2*math.sin(h/2)**2
        coupling_coeff = math.sqrt(a)*abs(s-math.sin(h))
        bound = math.sqrt(a)*a*K*math.sin(h)/(1-a*K)
        assert coupling_coeff <= bound*(1+1e-10)
        protected_var = c*c*target_var+a*math.sin(h)**2
        assert protected_var >= target_var-1e-15
        quadratic.append({"a":a,"lambda":lam,"invariance_residual":invariant_residual,
                          "projection_coupling_bound_ratio":coupling_coeff/bound,
                          "projection_changes_law_variance":protected_var-target_var,
                          "gamma":gamma,"beta":beta})
checks["quadratic"] = quadratic

# Nonlinear bounded-Hessian primitive: f'=1/2+(1/4)sech^2 in [1/2,3/4].
def force(x): return .5*x+.25*np.tanh(x)
def hessian(x): return .5+.25/np.cosh(x)**2

def finite(a,N,M,x0,L,protected=False,y=0.,dx0=1.):
    h=a**float(v/2)
    t=np.linspace(0,h,N+1)
    w=np.zeros((N+1,N+1))
    # Stable exact-kernel weights: cos(t_i-t_{j+1})-cos(t_i-t_j).
    for i in range(1,N+1):
        lo=t[i]-t[:i]; hi=t[i]-t[1:i+1]
        w[i,:i]=2*np.sin((lo+hi)/2)*np.sin((lo-hi)/2)
    c,s=np.cos(t),np.sin(t)
    base=y+c*(x0-y)+math.sqrt(a)*s*(0 if protected else L)
    # Differentiation in y, with dx0 the actual entering caller first.
    dbase=1+c*(dx0-1)
    q=base.copy(); dq=dbase.copy()
    for _ in range(M):
        dq=dbase-a*w@(hessian(q)*dq)
        q=base-a*w@force(q)
    val=q[-1]+(math.sqrt(a)*s[-1]*L if protected else 0)
    return val,dq[-1],q,w

nonlinear=[]
for a in [.2,.05,.01]:
    h=a**float(v/2); K=2*math.sin(h/2)**2
    for N in [4,16,64]:
        for L in [-2.,.5,1.7]:
            x0=.3*math.sqrt(a)
            full,dy,_,w=finite(a,N,4,x0,L,False,dx0=1+.1*a)
            prot,_,_,_=finite(a,N,4,x0,L,True,dx0=1+.1*a)
            bound=math.sqrt(a)*a*K*math.sin(h)*abs(L)/(1-a*K)
            assert abs(full-prot)<=bound*(1+1e-10)
            assert np.all(w>=0)
            assert abs(w[-1].sum()-K)<1e-14
            # Explicit actual finite caller first, independent of VALUE convergence.
            assert abs(dy-1)<=2*a
            zero1=finite(a,N,4,x0,0,False)[0]
            zero2=finite(a,N,4,x0,0,True)[0]
            assert zero1==zero2
            nonlinear.append({"a":a,"N":N,"L":L,"projection_ratio":abs(full-prot)/bound,"caller_first":dy,"origin_equal":True})
checks["nonlinear_finite"] = nonlinear

# Exact protected row / signed-twin identities including a full center port U.
a=.05; h=a**float(v/2); eps=a**.9; gam=math.sin(h)**2
beta=1/math.sqrt(1/gam+1/eps**2)
# Input columns: old Gaussian W, fresh L, external normalized center U.
C=np.array([[1.,0.,1.],[math.cos(h/2),0.,1.],[math.cos(h),math.sin(h),1.]])
P=np.array([[math.cos(h),math.sin(h),0.]])
Pi=np.diag([0.,1.,0.])
et=np.array([[0.],[0.],[1.]])
Ctw=np.block([[C,np.zeros((3,1))],[C,eps*et]])
R=np.concatenate([(P@Pi)/gam,np.array([[-1/eps]])],axis=1)
B=beta*R
expect=np.array([[0.,0.,1.,0.,0.,0.]])
assert np.max(np.abs(R@Ctw.T-expect))<1e-13
assert np.max(np.abs(Ctw@B.T-beta*expect.T))<1e-13
assert abs((B@B.T)[0,0]-1)<1e-13
assert abs((P@Pi@P.T)[0,0]-gam)<1e-15
checks["row_identities"]={"R_Ctw_transpose_error":float(np.max(np.abs(R@Ctw.T-expect))),"B_coisometry_error":float(abs((B@B.T)[0,0]-1)),"external_center_column_annihilated":bool(R[0,2]==0)}

# Finite-column lemma with deliberately varying Hessians at successive iterates.
A=np.array([[0.,0.,0.],[.008,0.,0.],[.011,.009,0.]])
S=np.block([[A+A.T,A.T],[A,np.zeros_like(A)]])
q=float(np.linalg.norm(S,2)); cn=float(np.linalg.norm(Ctw,2))
Jac=np.zeros((6,4)); firsts=[]
for k in range(1,21):
    hs=np.array([.6+.2*math.sin(k+j) for j in range(3)]+[-(.6+.2*math.cos(k+j)) for j in range(3)])
    Jac=np.diag(hs)@(Ctw+S@Jac)
    DG=(a/beta)*Ctw.T@Jac
    right=float(np.linalg.norm(DG@B.T,2))
    right_bound=a*cn/(1-q)
    assert right<=right_bound*(1+1e-12)
    assert float(np.max(np.abs(B@DG-a*Jac[2:3,:])))<1e-13
    firsts.append({"iteration":k,"physical_column":right,"lemma_bound":right_bound,"physical_row":float(np.linalg.norm(B@DG,2))})
checks["finite_firsts"]={"q":q,"iterations":firsts}

# Naive stack caveat: this is an actual norm growth, not a source-width theorem.
stack=[]
for N in [16,64,256,1024]:
    t=np.linspace(0,h,N,endpoint=False)
    norm=float(np.linalg.norm(np.sin(t)))
    stack.append({"N":N,"unweighted_restricted_norm_over_h":norm/h,"scaled_by_sqrtN":norm/(h*math.sqrt(N))})
checks["unweighted_stack_caveat"] = stack
checks["status"]="PASS: bounded algebraic and numerical diagnostics; weighted graph normalization remains a separate premise."
(OUT/'protected_short_refresh_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
print(checks['status'])
print('quadratic cases',len(quadratic),'nonlinear finite cases',len(nonlinear),'finite first iterates',len(firsts))
