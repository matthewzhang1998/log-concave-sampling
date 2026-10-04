#!/usr/bin/env python3
"""Independent finite diagnostics; imports no author diagnostic code.

Analytical proofs are in INDEPENDENT-CLOCK-GRADIENT-AUDIT.md. These finite
checks are not a proof of the LOW30 consumer or an all-order law compiler.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite import hermgauss
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent / 'POSITIVE-CLOCK-CALIBRATION-AND-DIAGONAL-GATE.md'
LOW30_PIN = '7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8'
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--low30-source', type=Path, help='Optional original LOW30 file for archival hash verification; not bundled.')
args = parser.parse_args()
# Publication adaptation: this was only a report-time hash read, not a
# mathematical assertion. All 1060 original mathematical checks are retained.
LOW30 = args.low30_source
LOW30_ACTUAL = hashlib.sha256(LOW30.read_bytes()).hexdigest() if LOW30 is not None else None
if LOW30_ACTUAL is not None and LOW30_ACTUAL != LOW30_PIN:
    raise ValueError('LOW30 pin mismatch')
rng = np.random.default_rng(620041004)
checks = 0
maxima = {'factor_error': 0., 'readout_error': 0., 'potential_fd_error': 0.,
          'private_fd_error': 0., 'caller_fd_error': 0., 'centered_caller_fd_error': 0.}

def check(ok, text):
    global checks
    assert bool(ok), text
    checks += 1

def small(x, tol, text):
    check(np.max(np.abs(x)) <= tol, text)

def op(x):
    return float(np.linalg.norm(x, 2))

def positive_rule(panels, order):
    nodes, weights = leggauss(order)
    out = []
    for j in range(panels):
        left, right = 2. ** (-j-1), 2. ** (-j)
        for x, w in zip(nodes, weights):
            t = (left + right)/2 + (right-left)*x/2
            out.append((1-t, (right-left)*w/2))
    out.append((1-2. ** (-panels-1), 2. ** (-panels)))
    out.sort(reverse=True)
    return tuple(np.array(z) for z in zip(*out))

def brownian_factor(r):
    q = np.expm1(-2*np.log(r))/2
    prev = np.r_[0., q[:-1]]
    dq = q-prev
    check(np.all(dq > 0), 'strictly ordered positive Brownian increments')
    L = np.sqrt(2)*r[:,None]*np.sqrt(dq)[None,:]*np.tri(len(r))
    # Algebraically identical to (sqrt(q)-sqrt(prev))/sqrt(dq), but stable.
    t = np.sqrt(dq)/(np.sqrt(q)+np.sqrt(prev))
    c = np.sqrt((1-r)*(1+r))
    R = np.minimum.outer(r,r)/np.maximum.outer(r,r)
    e = np.max(abs(L@L.T-(R-np.outer(r,r))))
    maxima['factor_error'] = max(maxima['factor_error'], float(e))
    small(e, 3e-13, 'Brownian factor gives exact OU conditional covariance')
    small(L@t-c, 3e-13, 'telescoping L t equals conditional row norms')
    small(np.linalg.norm(L,axis=1)-c, 3e-13, 'conditional row norm c')
    bound = 1 + .25*np.log(q[-1]/q[0])
    check(t@t <= bound+1e-13, 'logarithmic readout bound')
    if len(r)>1:
        small(t[1:]**2-np.tanh(.25*np.diff(np.log(q))), 4e-14, 'exact tanh increment formula')
    check(abs(t[0]-1)<1e-14, 'first readout increment equals one')
    return L,t,R,q,bound

def calibrated_factor(r,w):
    L,t,R,q,bound = brownian_factor(r)
    c = np.sqrt((1-r)*(1+r))
    independent = np.diag(c)
    common = c[:,None]
    S = float(w@R@w)
    s0 = float((w@r)**2 + np.sum(w*w*c*c))
    sm = float((w@r)**2 + (w@c)**2)
    check(s0 < .5 < sm, 'strict calibration bracket')
    aux, branch = (independent,'independent') if S>=.5 else (common,'common')
    sa = s0 if S>=.5 else sm
    eta = (.5-S)/(sa-S) if S!=.5 else 0.
    check(0 <= eta <= 1, 'positive covariance interpolation')
    C = (1-eta)*R+eta*(np.outer(r,r)+aux@aux.T)
    Lcal = np.column_stack((math.sqrt(1-eta)*L, math.sqrt(eta)*aux))
    tcal = np.r_[t/math.sqrt(1-eta), np.zeros(aux.shape[1])]
    small(C-(np.outer(r,r)+Lcal@Lcal.T), 4e-13, 'literal concatenated factor, no ridge or covariance oracle')
    small(np.diag(C)-1, 4e-13, 'standard calibrated marginals')
    small(Lcal@tcal-c, 4e-13, 'calibrated readout exactly recovers c')
    small(w@C@w-.5, 4e-13, 'exact calibrated quadratic variance')
    check(np.linalg.eigvalsh(C-np.outer(r,r)).min()>=-4e-13, 'conditional covariance PSD')
    small(tcal@tcal-(t@t)/(1-eta), 5e-13, 'exact calibration norm multiplier')
    if eta<=.5:
        check(tcal@tcal <= 2*bound+2e-13, 'factor-two guard when eta <= 1/2')
    return Lcal,tcal,C,{'raw_S':S,'S0':s0,'Smax':sm,'eta':eta,'branch':branch,
                        'roots':Lcal.shape[1],'t_norm2':float(tcal@tcal),'log_bound':bound/(1-eta)}

class C2Potential:
    """Uniformly convex C2, generally non-C3 nonorthogonal ridge fixture."""
    def __init__(self,d,A):
        self.d,self.A = d,A
        self.B = rng.normal(size=(d+2,d))/math.sqrt(d+2)
        self.scale = A/(1+op(self.B)**2)
    def U(self,x):
        s = abs(np.asarray(x)@self.B.T)
        h = np.where(s<=1, 4*s**2.5/15, .5*s*s-s/3+.1)
        return self.scale*(np.sum(np.asarray(x)**2,axis=-1)/2+h.sum(axis=-1))
    def g(self,x):
        z = np.asarray(x)@self.B.T
        s = abs(z)
        hp = np.sign(z)*np.where(s<=1, 2*s**1.5/3,s-1/3)
        return self.scale*(np.asarray(x)+hp@self.B)
    def hess(self,x):
        z = np.asarray(x)@self.B.T
        hpp = np.minimum(np.sqrt(abs(z)),1.)
        return self.scale*(np.eye(self.d)+np.einsum('...i,ij,ik->...jk',hpp,self.B,self.B))

def gradient_values(pot,z,G,r,w,L):
    c = np.sqrt((1-r)*(1+r))
    x = r[:,None]*z+L@G
    f = pot.g(x)
    j = L.T@((w/c)[:,None]*f)
    j0 = L.T@((w/c)[:,None]*pot.g(r[:,None]*z))
    h = pot.hess(x)
    h0 = pot.hess(r[:,None]*z)
    priv = np.einsum('ia,ib,ide,i->adbe',L,L,h,w/c).reshape(G.size,G.size)
    dz = np.einsum('ia,ide,i->ade',L,h,w*r/c).reshape(G.size,z.size)
    dzc = np.einsum('ia,ide,i->ade',L,h-h0,w*r/c).reshape(G.size,z.size)
    return j,j0,x,f,priv,dz,dzc

source_cases=[]
for K,m in [(2,2),(3,3),(5,5),(8,4)]:
    r,w=positive_rule(K,m)
    check(np.all(w>0),'positive weights')
    small([w.sum()-1,w@r-.5],5e-14,'exact constant/linear clocks')
    raw= brownian_factor(r)
    cal= calibrated_factor(r,w)
    for name,L,t in [('OU',raw[0],raw[1]),('calibrated',cal[0],cal[1])]:
        c=np.sqrt((1-r)*(1+r)); beta=np.linalg.norm(t)
        small(np.sum((t/beta)**2)-1,4e-14,'normalized physical readout is a coisometry')
        for d in [1,3]:
            for A in [0.,.03,.3]:
                pot=C2Potential(d,A)
                z=rng.normal(size=d); G=rng.normal(size=(L.shape[1],d))
                j,j0,x,f,priv,dz,dzc=gradient_values(pot,z,G,r,w,L)
                e=np.max(abs(t@j-w@f))
                maxima['readout_error']=max(maxima['readout_error'],float(e))
                small(e,4e-13,'original VALUE full-gradient readout H')
                small(t@j0-w@pot.g(r[:,None]*z),4e-13,'same caller-only origin readout')
                ev=np.linalg.eigvalsh(priv)
                check(ev.min()>=-5e-13,'genuine convex private gradient Hessian')
                check(ev.max()<=A*(w@c)+5e-13,'actual private first <= A sum w c')
                check(op(dz)<=A/2+5e-13,'actual endpoint caller first <= A/2')
                check(op(dzc)<=A/2+5e-13,'centered endpoint caller first also <= A/2 by PSD sandwich')
                check(np.linalg.norm(j0)<=A*np.linalg.norm(z)/2+5e-13,'nonzero origin bound')
                check(np.linalg.norm(j-j0)<=A*np.sum(w*np.linalg.norm(L@G,axis=1))+5e-13,
                      'pointwise one-energy triangle before Gaussian moments')
                zero=gradient_values(pot,np.zeros(d),np.zeros_like(G),r,w,L)
                small(zero[0],0.,'literal joint zero')
                small(zero[1],0.,'literal origin zero')
                eps=2e-6
                v=rng.normal(size=G.shape); v/=np.linalg.norm(v)
                vj=rng.normal(size=d); vj/=np.linalg.norm(vj)
                def psi(H):
                    return np.sum((w/c)*pot.U(r[:,None]*z+L@H))
                fd=(psi(G+eps*v)-psi(G-eps*v))/(2*eps)
                er=abs(fd-np.sum(j*v)); maxima['potential_fd_error']=max(maxima['potential_fd_error'],float(er))
                small(er,2e-8,'J is full gradient of the displayed analytical potential')
                jp=gradient_values(pot,z,G+eps*v,r,w,L)[0]
                jm=gradient_values(pot,z,G-eps*v,r,w,L)[0]
                er=np.max(abs((jp-jm).ravel()/(2*eps)-priv@v.ravel()))
                maxima['private_fd_error']=max(maxima['private_fd_error'],float(er))
                small(er,2e-7,'actual private HVP chain rule against VALUE directional difference')
                pp=gradient_values(pot,z+eps*vj,G,r,w,L)
                pm=gradient_values(pot,z-eps*vj,G,r,w,L)
                er=np.max(abs((pp[0]-pm[0]).ravel()/(2*eps)-dz@vj))
                maxima['caller_fd_error']=max(maxima['caller_fd_error'],float(er))
                small(er,2e-7,'actual endpoint caller derivative including every VALUE site')
                er=np.max(abs(((pp[0]-pp[1])-(pm[0]-pm[1])).ravel()/(2*eps)-dzc@vj))
                maxima['centered_caller_fd_error']=max(maxima['centered_caller_fd_error'],float(er))
                small(er,2e-7,'centered caller derivative includes all origin sites')
        # Exact L2 energy diagnostic for g(x)=A x, with arbitrary D factor canceled.
        M=L.T@((w/c)[:,None]*L)
        check(np.linalg.norm(M,'fro')<=w@c+5e-13,'quadratic one-energy bound has no sqrt(number of roots)')
        source_cases.append({'panels':K,'order':m,'nodes':len(r),'factor':name,'private_scalar_roots':L.shape[1],
                             'readout_norm2':float(t@t),'quadratic_energy_ratio':float(np.linalg.norm(M,'fro')/(w@c))})

# Common-root auxiliary branch, plus a valid calibration that fails the optional half-mixing guard.
r=np.array([.99,.5,.01]);w=np.array([.4,.2,.4])
cm=calibrated_factor(r,w)
check(cm[3]['branch']=='common' and cm[3]['eta']<.5,'explicit maximal-root auxiliary branch')
r_bad=np.array([.53,.51,.49,.47]);w_bad=np.ones(4)/4
bad=calibrated_factor(r_bad,w_bad)
check(bad[3]['eta']>.5,'bracketing alone does not imply eta <= 1/2')

# Broad range of clocks checks normalization and logarithmic behavior independently of quadrature.
for n in [1,2,5,20,80]:
    qq=np.geomspace(1e-10,1e8,n)
    r=1/np.sqrt(1+2*qq)
    L,t,R,q,bnd=brownian_factor(r)
    check(np.max(abs((L/math.sqrt(2))@t-np.sqrt(1-r*r)))>1e-7,
          'negative control catches dropping sqrt(2) from Brownian normalization')

# Finite calibrated rules, Markov triples, full matrix quadratic law, and diagonal separator.
tau,v=positive_rule(24,24)
def Cfun(u,R):
    return .5*(np.exp(-u*(1-R))+np.exp(-u*(1+R))-2*math.exp(-u))
def Ffun(u,t):
    return math.exp(-u/2)-.5*((1-t)*np.exp(-u*(1-t))+(1+t)*np.exp(-u*(1+t)))
def Iexact(u):
    return -math.expm1(-2*u)/(2*u)-math.exp(-u)
def Fexact(u):
    return math.exp(-u/2)+math.exp(-2*u)/u+math.expm1(-2*u)/(2*u*u)
cosine=[]
for K,m in [(2,2),(4,4),(8,8),(12,12),(20,20)]:
    r,w=positive_rule(K,m); L,t,C,info=calibrated_factor(r,w)
    check(info['eta']<=.5,'tested fine rules satisfy optional calibration guard')
    rows=np.column_stack((r,L))
    for i in [0,len(r)//2,len(r)-1]:
        for tauj in [.05,.5,.95]:
            joint=np.zeros((3,rows.shape[1]+1));joint[0,0]=1
            joint[1,:-1]=rows[i];joint[2,:-1]=tauj*rows[i];joint[2,-1]=math.sqrt(1-tauj*tauj)
            target=np.array([[1,r[i],r[i]*tauj],[r[i],1,tauj],[r[i]*tauj,tauj,1]])
            small(joint@joint.T-target,5e-13,'actual retained endpoint/outer/inner Markov triple')
    B0=rng.normal(size=(3,3));B=.07*(B0@B0.T)/(op(B0)**2)
    h=w@rows
    vi=float(np.sum(w*w)*np.sum(v*v*(1-tau*tau)))
    ell=np.r_[float(v@tau)*h,math.sqrt(vi)]
    h=np.r_[h,0.];z=np.zeros_like(h);z[0]=1.
    blocks=[zj*np.eye(3)-hj*B+lj*(B@B) for zj,hj,lj in zip(z,h,ell)]
    cov=sum(b@b.T for b in blocks)
    formula=np.eye(3)-B+B@B-2*(h@ell)*np.linalg.matrix_power(B,3)+(ell@ell)*np.linalg.matrix_power(B,4)
    small(cov-formula,5e-13,'full non-diagonal matrix covariance through fourth degree')
    s2=float(w@w);u=4/s2
    kq=float(w@Cfun(u,C)@w); iq=Iexact(u)
    fq=float(v@Ffun(u,tau)); fb=Fexact(u)
    dcoef=((kq-iq)/u+fq-fb)/16
    check(kq>=s2*(1-math.exp(-u))**2/2-2e-13,'positive cosine covariance diagonal survives every PSD calibration')
    check((kq-iq)/u>=s2*s2/12-2e-13,'rank-two diagonal lower bound')
    check(abs(fq-fb)<=s2*s2/48,'independent high-resolution inner rule fits nonlinear coefficient budget')
    check(dcoef>=s2*s2/256-2e-13,'nested predictor second-order variance separator')
    cosine.append({'nodes':len(r),'s2':s2,'u':u,'coefficient_debt':dcoef,'lower_bound':s2*s2/256,
                   'inner_error':fq-fb,'eta':info['eta'],'readout_norm2':float(t@t)})

# Independent Gaussian integration confirms the nonlinear current and target centering.
gh,gw=hermgauss(100);gh=math.sqrt(2)*gh;gw=gw/math.sqrt(math.pi)
a,b=.5,.25
current=[]
for u in [4.,9.,20.]:
    k=math.sqrt(u)
    W=a*gh*gh/2+b/u*np.sin(k*gh)-b/k*gh
    EW=float(gw@W)
    coef=float(.5*(gw@((gh*gh-1)*W*W))-(gw@((gh*gh-1)*W))*EW-(gw@(gh*W))**2)
    T=math.exp(-u/2)+(math.exp(-2*u)-math.exp(-u))/u
    small(coef-(a*a+b*b*T),2e-12,'true target centered second-order variance coefficient')
    for rr in [.1,.7,.99]:
        for tt in [0.,.2,.8,.99]:
            X=gh[:,None];Y=tt*X+math.sqrt(1-tt*tt)*gh[None,:]
            fp=a-b*np.sin(k*X);fy=a*Y+b/k*(np.cos(k*Y)-1)
            actual=float(rr*np.sum(gw[:,None]*gw[None,:]*X*fp*fy))
            exact=rr*(a*a*tt+b*b*Ffun(u,tt))
            small(actual-exact,3e-11,'exact Gaussian Markov rank-one current including exp(-u/2)')
            current.append(abs(actual-exact))
    integ=quad(lambda tt:Ffun(u,tt),0,1,epsabs=1e-13)[0]
    small(integ-Fexact(u),2e-13,'exact integrated current formula')
    small(Iexact(u)/u+Fexact(u)-T,2e-14,'continuous predictor and true target coefficient coincide')

report={'status':'PASS','assertions':checks,'seed':620041004,
        'scope':'Independent finite tests of the protected-endpoint gradient lift, calibration, actual firsts/origins, and cosine separator. No execution or proof of the imported high-order LOW30 consumer.',
        'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        'low30_sha256':LOW30_PIN,
        'low30_pin_status':'verified_original_bytes' if LOW30 is not None else 'external_not_reverified',
        'mathematical_assertions_removed':0,
        'maxima':maxima,'source_cases':source_cases,'cosine_cases':cosine,
        'common_auxiliary_case':cm[3],'half_mixing_guard_counterexample':bad[3],
        'max_rank_one_gaussian_integration_error':max(current),
        'consumer_dimension_warning':'Safe square-gradient consumer has n=M*D, readout beta=||t||, and normalized radius at most beta*A*sum(w*c)/sqrt(v). Its printed theorem has sqrt(M*D), not a proved sqrt(D) replacement.'}
(ROOT/'clock_gradient_independent_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
