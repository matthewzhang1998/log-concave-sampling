#!/usr/bin/env python3
"""Independent diagnostics, not a reimplementation of imported mean compilers.

Checks use exact symbolic covariance algebra and independent smooth source
fixtures. The analytic audit supplies proofs; numerical checks cannot certify
an arbitrary black-box source's Gaussian residual norms.
"""
from pathlib import Path
import hashlib
import json
import math
import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite import hermgauss

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
rng = np.random.default_rng(610052026)
checks = 0
metrics = {}

def check(condition, description):
    global checks
    checks += 1
    if not bool(condition):
        raise AssertionError(description)

def close(a, b, atol=1e-11, rtol=1e-10):
    return np.allclose(a, b, atol=atol, rtol=rtol)

def op(a):
    return float(np.linalg.norm(a, 2))

# Independently derive the conditional covariance from OU kernel integration.
a, b = sp.symbols('a b', positive=True)
L = (1/(a+b)) * (1/(a+1) + 1/(b+1))
sub = {a: 1, b: 1}
unconditional = sp.Matrix([
    [L.subs(sub), (-sp.diff(L,b)).subs(sub)],
    [(-sp.diff(L,a)).subs(sub), sp.diff(L,a,b).subs(sub)]])
reg = sp.Matrix([sp.Rational(1,2), sp.Rational(1,4)])
conditional = unconditional - reg*reg.T
root = sp.Matrix([[sp.Rational(1,2), 0],
                  [sp.Rational(1,2), sp.Rational(1,4)]])
check(unconditional == sp.Matrix([[sp.Rational(1,2), sp.Rational(3,8)],
                                 [sp.Rational(3,8), sp.Rational(3,8)]]), 'OU unconditional block')
check(conditional == root*root.T, 'two-root exact conditional block')
check(conditional.det() == sp.Rational(1,64), 'nondegenerate two-root covariance')
metrics['conditional_covariance'] = [[str(x) for x in row] for row in conditional.tolist()]

# A positive interior rule. This tests finite-source identities, not the
# infinite-Hermite multiplier bound of the separately imported rule.
t, weights = leggauss(7)
t, weights = (t+1)/2, weights/2
c = np.sqrt(1-t*t)
beta = float(weights @ c)
check(np.all(weights > 0) and np.all((t > 0) & (t < 1)), 'positive interior nodes')
check(close(weights.sum(),1) and close(weights @ t,.5), 'mass and first moment')
check(beta<=math.sqrt(3)/2+1e-14,'Jensen beta bound')
for aa in np.linspace(0,.5,101):
    eta=aa/2+aa*aa/4; qq=aa/2+aa*aa/2; rr=aa*aa/4
    bb=math.sqrt(3)/2
    normalized_first=math.sqrt(2)*aa*math.sqrt(bb*bb*(1+eta)**2+qq*qq+rr*rr)
    normalized_curl=math.sqrt(2)*(2*bb*aa*eta+aa*math.sqrt(qq*qq+rr*rr))
    check(normalized_first<=2*aa+1e-14,'safe declared ell_E=2A')
    check(normalized_curl<=3*aa*aa+1e-14,'safe normalized curl <=3A squared')

class SmoothSource:
    """Anchored smooth convex gradient with noncommuting local Hessians."""
    def __init__(self, dim, A):
        self.dim, self.A, self.count = dim, A, 0
        dirs = rng.normal(size=(dim+3,dim))
        dirs /= np.linalg.norm(dirs,axis=1)[:,None]
        self.dirs = dirs
        self.scale = .8 / op(dirs.T @ dirs)
    def value(self,x):
        self.count += 1
        return self.A*(.1*x + self.scale*(np.tanh(self.dirs @ x) @ self.dirs))
    def hess(self,x):
        tt = np.tanh(self.dirs @ x)
        return self.A*(.1*np.eye(self.dim) + self.scale*self.dirs.T @ ((1-tt*tt)[:,None]*self.dirs))

class QuadraticSource:
    def __init__(self,K):
        self.K, self.dim, self.count = K,K.shape[0],0
    def value(self,x):
        self.count += 1
        return self.K @ x
    def hess(self,x):
        return self.K

def evaluate(src,z,G,N,M, which='E', derivatives=False):
    d = len(z)
    ans = np.zeros(d)
    Jz,JG,JN,JM = [np.zeros((d,d)) for _ in range(4)]
    tape=[]
    for ti,ci,wi in zip(t,c,weights):
        x=ti*z+ci*G
        if which == 'B':
            ans += wi*src.value(x)
            if derivatives:
                hx=src.hess(x)
                Jz += wi*ti*hx
                JG += wi*ci*hx
            continue
        v=x/2+N/2
        w=x/4+N/2+M/4
        gv=src.value(v)
        gw=src.value(w)
        ggw=src.value(gw)
        T=x-gv+ggw
        ans += wi*src.value(T)
        hx=np.zeros((d,d))
        if which == 'E':
            ans -= wi*src.value(x)
            hx=src.hess(x)
        if derivatives:
            H,V,J,W=map(src.hess,(T,v,gw,w))
            dx=H@(np.eye(d)-V/2+J@W/4)-hx
            dn=H@(-V/2+J@W/2)
            dm=H@J@W/4
            Jz += wi*ti*dx
            JG += wi*ci*dx
            JN += wi*dn
            JM += wi*dm
            tape.append((ti,ci,wi,H,V,J,W,hx))
    return (ans,(Jz,JG,JN,JM),tape) if derivatives else ans

def reverse(tape,u):
    # Actual reverse sweep uses Hessian actions only at original VALUE sites.
    bz,bG,bN,bM=[np.zeros_like(u) for _ in range(4)]
    for ti,ci,wi,H,V,J,W,hx in tape:
        q=H@(wi*u)
        bv=-V@q
        bw=W@(J@q)
        bx=q+bv/2+bw/4-hx@(wi*u)
        bz += ti*bx
        bG += ci*bx
        bN += bv/2+bw/2
        bM += bw/4
    return np.concatenate((bz,bG,bN,bM))

max_fd=0.
max_noncommuting_skew=0.
max_private_ratio=max_curl_ratio=0.
for A in (0.,1e-4,.03,.17,.5):
    for dim in (1,2,5,9):
        src=SmoothSource(dim,A)
        for rep in range(4):
            z,G,N,M=rng.normal(size=(4,dim))*(1+rep)
            val,js,tape=evaluate(src,z,G,N,M,derivatives=True)
            Jz,JG,JN,JM=js
            Htest=src.hess(rng.normal(size=dim))
            ev=np.linalg.eigvalsh(Htest)
            check(ev.min() >= -1e-13 and ev.max() <= A+1e-13, 'convex source Hessian interval')
            q=1+A/2+A*A/4
            n=A*A*(1+A)/2
            m=A**3/4
            private=np.hstack((JG,JN,JM))
            radius=math.sqrt((A*beta*q)**2+n*n+m*m)
            lifted=np.zeros((3*dim,3*dim)); lifted[:dim,:]=private
            skew=lifted-lifted.T
            curl_bound=beta*(A*A+A**3/2)+math.sqrt(n*n+m*m)
            check(op(JG)<=A*beta*q+1e-12, 'G first')
            check(op(JN)<=n+1e-12, 'N first')
            check(op(JM)<=m+1e-12, 'M first')
            check(op(Jz)<=A*q/2+1e-12, 'captured caller first')
            check(op(private)<=radius+1e-12, 'complete private first')
            check(op(skew)<=curl_bound+1e-12, 'complete lifted curl')
            if A:
                max_private_ratio=max(max_private_ratio,op(private)/radius)
                max_curl_ratio=max(max_curl_ratio,op(skew)/curl_bound)
            max_noncommuting_skew=max(max_noncommuting_skew,op(JG-JG.T))
            amp=A*A*((.25+A/8)*np.linalg.norm(z)
                +(beta/2+A*beta/4)*np.linalg.norm(G)
                +(.5+A/2)*np.linalg.norm(N)+A*np.linalg.norm(M)/4)
            check(np.linalg.norm(val)<=amp+1e-12,'pointwise amplitude')
            origin,origin_js,_=evaluate(src,z,np.zeros(dim),np.zeros(dim),np.zeros(dim),derivatives=True)
            check(np.linalg.norm(origin)<=A*A*(.25+A/8)*np.linalg.norm(z)+1e-12,'actual caller origin amplitude')
            check(op(Jz-origin_js[0])<=A*q+1e-12,'anchored caller first')
            zero=evaluate(src,*[np.zeros(dim) for _ in range(4)])
            check(np.array_equal(zero,np.zeros(dim)),'literal total zero')
            base,bjs,_=evaluate(src,z,G,N,M,which='B',derivatives=True)
            check(close(bjs[1],bjs[1].T),'baseline genuine gradient')
            check(op(bjs[1])<=A*beta+1e-12,'baseline first')
            direction=rng.normal(size=4*dim); direction/=np.linalg.norm(direction)
            fields=np.concatenate((z,G,N,M)); hh=2e-6
            plus=evaluate(src,*np.split(fields+hh*direction,4))
            minus=evaluate(src,*np.split(fields-hh*direction,4))
            exact=np.hstack(js)@direction
            fd=(plus-minus)/(2*hh)
            diff=np.linalg.norm(fd-exact)
            max_fd=max(max_fd,diff)
            check(diff<=4e-9,'full directional chain rule')
            u=rng.normal(size=dim)
            check(close(reverse(tape,u),np.hstack(js).T@u),'complete adjoint replay')

check(max_noncommuting_skew>1e-7,'fixture detects nonsymmetric G derivative')
metrics.update(max_directional_chain_rule_error=max_fd,
               max_observed_G_skew=max_noncommuting_skew,
               max_private_bound_ratio=max_private_ratio,
               max_curl_bound_ratio=max_curl_ratio)

# Every anisotropic quadratic is handled by four VALUES per outer node.
max_quad_error=0.
for dim in (1,2,4,8):
    for A in (1e-4,.05,.25,.5):
        Q,_=np.linalg.qr(rng.normal(size=(dim,dim)))
        K=Q@np.diag(np.linspace(0,A,dim) if dim>1 else [A*.7])@Q.T
        src=QuadraticSource(K)
        z,G,N,M=rng.normal(size=(4,dim))
        I=np.eye(dim); K2=K@K
        P=I-K/2+K2/4
        C=-K/2+K2/2
        D=K2/4
        expected=K@(P@(z/2+beta*G)+C@N+D@M)
        actual=evaluate(src,z,G,N,M,which='F')
        max_quad_error=max(max_quad_error,float(np.linalg.norm(actual-expected)))
        check(close(actual,expected),'exact anisotropic quadratic values')
        check(src.count==4*len(t),'four original VALUES per terminal occurrence')
        src.count=0; evaluate(src,z,G,N,M,which='B')
        check(src.count==len(t),'baseline leaf count')
        src.count=0; evaluate(src,z,G,N,M)
        check(src.count==5*len(t),'five original VALUES per correction occurrence')
        # Matrix covariance calculation uses analytical K; no K is a producer leaf.
        conditional_cov=K2/4-K@K2/2+5*(K2@K2)/16
        check(close(C@C.T+D@D.T,conditional_cov),'all quadratic covariance directions')
        own_mean=K@P@z/2
        check(close(own_mean,.5*(K-.5*K2+.25*K2@K)@z),'exact canonical quadratic mean')
        S=I-K+K2/2
        check(np.linalg.eigvalsh(S).min()>0 and op(S)<=1+1e-12,'explicit anisotropic history covariance')
        # Independent roots at every node preserve mean but change raw covariance.
        shared=beta*beta*K@P@P@K+K@(C@C+D@D)@K
        iid=sum(wi*wi*(ci*ci*K@P@P@K+K@(C@C+D@D)@K) for wi,ci in zip(weights,c))
        if A>.01:
            check(op(shared-iid)>1e-10,'shared-level ancestry has material raw covariance')
metrics['max_quadratic_value_error']=max_quad_error

# Original VALUE absolute perturbations propagate through all four calls.
src=SmoothSource(4,.5)
z,G,N,M=rng.normal(size=(4,4)); nu=3e-7
true=evaluate(src,z,G,N,M,which='F')
class Perturbed:
    def __init__(self,base): self.base=base; self.dim=base.dim
    def value(self,x):
        # A bounded deterministic error at each actual, possibly perturbed site.
        v=np.cos(np.arange(self.dim)+float(x.sum())); v/=np.linalg.norm(v)
        return self.base.value(x)+nu*v
approx=evaluate(Perturbed(src),z,G,N,M,which='F')
check(np.linalg.norm(approx-true)<=(1+src.A)**2*nu+1e-12,'nested absolute VALUE floor')

# Analytic fixed-D obstruction: h(x)=(x+sin x)/2. Evaluate its exact
# first-Hermite coefficient independently by Gaussian and t quadrature.
tl,wl=leggauss(160); tl=(tl+1)/2; wl=wl/2
zh,wh=hermgauss(120); zh=zh*math.sqrt(2); wh=wh/math.sqrt(math.pi)
inner=np.sin(zh[:,None]*tl[None,:])@(wl*np.exp(-(1-tl*tl)/2))
qfun=.25*(1+np.cos(zh))*(inner-math.exp(-.125)*np.sin(zh/2))
hermite_coefficient=float(wh@(zh*qfun))
exact_coefficient=.25*(.5*math.exp(-.5)+2*math.exp(-1)-.5
    -1.5*math.exp(-2)-.25*math.exp(-.25)-.75*math.exp(-1.25))
check(abs(hermite_coefficient-exact_coefficient)<1e-12,'explicit obstruction Hermite coefficient')
check(exact_coefficient<-.0183,'nonzero first Hermite obstruction')
# Exact rational Taylor enclosures establish the sign and decimal lower bound.
def exp_interval_negative(x):
    x=sp.Rational(x)
    lower=sum((-x)**j/sp.factorial(j) for j in range(32))
    upper=sum((-x)**j/sp.factorial(j) for j in range(33))
    return lower,upper
terms=[(sp.Rational(1,8),sp.Rational(1,2)),(sp.Rational(1,2),1),
       (-sp.Rational(3,8),2),(-sp.Rational(1,16),sp.Rational(1,4)),
       (-sp.Rational(3,16),sp.Rational(5,4))]
coef_lo=coef_hi=-sp.Rational(1,8)
for mult,arg in terms:
    lo,hi=exp_interval_negative(arg)
    coef_lo+=mult*(lo if mult>0 else hi)
    coef_hi+=mult*(hi if mult>0 else lo)
check(coef_hi<0,'exact rational sign of obstruction coefficient')
check(-coef_hi/2>sp.Rational('0.0091946882585889'),'rigorous decimal lower bound')
check(coef_lo<sp.Rational(str(exact_coefficient))<coef_hi or abs(float(coef_lo)-exact_coefficient)<1e-15,'numerical coefficient matches certified enclosure')
A=.5
Cmax=1+math.sqrt(3)/4*(1+A)**2+math.sqrt(3/8)+math.sqrt(3)/4*(1/math.sqrt(2)+A*math.sqrt(3/8))**2
check(Cmax<3.04,'uniform Taylor remainder constant')
metrics['obstruction']={'q_first_hermite':exact_coefficient,
 'R1_q_L2_lower_bound':abs(exact_coefficient)/2,
 'Taylor_remainder_constant_at_A_half':Cmax,
 'certified_lower_bound':'(0.009194688258588945) A^2 - 3.04 A^3'}

pins={}
for path in sorted(BASE.glob('*.md')):
    pins[str(path.relative_to(BASE))]=hashlib.sha256(path.read_bytes()).hexdigest()
for rel in ('../gaussian-rms-backbone/GAUSSIAN-RMS-RESIDUAL-RESUMMED-M3.md',
            '../resummed-linear-backbone/BOUNDED-REMAINDER-RESUMMED-M3-VALUE-CONSUMER.md',
            '../resummed-linear-backbone/independent-audit/INDEPENDENT-BOUNDED-RESUMMED-M3-AUDIT.md'):
    p=BASE/rel
    pins[rel]=hashlib.sha256(p.read_bytes()).hexdigest()
result={'status':'PASS','assertions':checks,'seed':610052026,
 'scope':'Independent raw-source/covariance/chain-rule/port/cost/obstruction diagnostics; not execution of imported complete mean compilers.',
 'metrics':metrics,'source_sha256':pins}
(HERE/'independent_matrix_free_m3_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
