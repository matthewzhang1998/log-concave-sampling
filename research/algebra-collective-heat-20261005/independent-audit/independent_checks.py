#!/usr/bin/env python3
"""Independent finite mathematical diagnostics; no native execution or sampling."""
import hashlib, json, math, itertools
from pathlib import Path
import sympy as S
import numpy as np
import mpmath as mp
HERE=Path(__file__).resolve().parent
BASE=HERE.parent
checks=[]
def ck(name, condition, detail=None):
    ok=bool(condition)
    checks.append(dict(name=name,passed=ok,**({'detail':detail} if detail is not None else {})))
    if not ok: raise AssertionError(name)
def ez(e): return S.expand(e)==0
# Verify each imported source independently without modifying it.
pins=json.loads((BASE/'INPUT-PINS.json').read_text())
for ent in pins['sources']:
    data=Path(ent['path']).read_bytes()
    ck('pin_'+Path(ent['path']).name,len(data)==ent['bytes'] and hashlib.sha256(data).hexdigest()==ent['sha256'])
# Vector-dimensional exact differential operators, including unequal/zero rates.
z=S.symbols('z0:6'); vv=[z[0:2],z[2:4],z[4:6]]
f=(1+z[0]*z[1]+z[0]**2)*(2+z[2]**2+z[3])*(1+z[4]+z[4]*z[5])
s,theta=S.symbols('s theta')
def ops(rates):
    def ind(f): return S.expand(sum(rates[v]**2*S.diff(f,vv[v][i],2) for v in range(3) for i in range(2))/2)
    def cross(f): return S.expand(sum(rates[u]*rates[v]*S.diff(f,vv[u][i],vv[v][i]) for u in range(3) for v in range(u+1,3) for i in range(2)))
    def common(f): return S.expand(sum(rates[u]*rates[v]*S.diff(f,vv[u][i],vv[v][i]) for u in range(3) for v in range(3) for i in range(2))/2)
    return ind,cross,common
def heat(f,op):
    out=f; d=f
    for j in range(1,10):
        d=op(d)
        if d==0: return S.expand(out)
        out+=s**j*d/S.factorial(j)
    raise AssertionError('Fixture not truncated')
for rates in [(1,1,1),(0,2,3),(S.Rational(1,2),1,3)]:
    hi,bc,hc=ops(rates)
    label='_'.join(map(str,rates))
    ck('2d_operator_'+label,ez(hi(f)-hc(f)+bc(f)))
    ck('2d_commute_'+label,ez(hi(bc(f))-bc(hi(f))))
    ck('2d_complete_finite_jet_'+label,ez(heat(f,hi)-heat(heat(f,bc).subs(s,-s),hc)))
    p=heat(f,lambda h:hi(h)+theta*bc(h))
    ck('2d_covariance_derivative_'+label,ez(S.diff(p,theta)-s*heat(bc(f),lambda h:hi(h)+theta*bc(h))))
    ck('2d_positive_homotopy_'+label,ez(p.subs(theta,0)-p.subs(theta,1)+s*S.integrate(heat(bc(f),lambda h:hi(h)+theta*bc(h)),(theta,0,1))))
# Exact vector endpoint adjoint with overlapping non-unit Gaussian directions.
b,c,k,l=S.symbols('b c k l'); xx,yy=S.symbols('x y')
directions=[(1,2),(S.Rational(1,2),-1)]
def dd(f,a):return a[0]*S.diff(f,b)+a[1]*S.diff(f,c)
def ge(f):
    ans=0
    for powers,coeff in S.Poly(S.expand(f),b,c,k,l).terms():
        if any(p%2 for p in powers):continue
        ans+=coeff*math.prod(int(S.factorial2(p-1)) if p else 1 for p in powers)
    return S.simplify(ans)
X=S.Matrix([b*b+c+S.Rational(2,3)*k,b*c+c*c+S.Rational(3,4)*l])
C=S.Matrix([1+b*c,2+b*b-c]); Q=1+b+c*c; chi=1+b*b+b*c
phi=xx**3*yy**2+xx*yy**4+xx**2
Dphi=S.Matrix([S.diff(phi,t) for t in (xx,yy)])
def sub(f):return f.subs({xx:X[0],yy:X[1]},simultaneous=True)
# h = C : grad phi. Keep the coefficient C out of adjoint endpoint differentiation.
grad=sub(Dphi)
beta=[a[0]*b+a[1]*c for a in directions]
V=[X.applyfunc(lambda e:dd(e,a)) for a in directions]
W=[v.applyfunc(lambda e:dd(e,a)) for v,a in zip(V,directions)]
score=sum(be*be-sum(q*q for q in a) for be,a in zip(beta,directions))/2
HC=C.applyfunc(lambda e:sum(dd(dd(e,a),a) for a in directions)/2)
lhs=(HC.dot(grad))
rhs=score*C.dot(grad)
for a,be,v,w in zip(directions,beta,V,W):
    for i in range(2):
        for j in range(2):
            rhs+=C[i]*(w[j]/2-be*v[j])*sub(S.diff(phi,(xx,yy)[i],(xx,yy)[j]))
            for h in range(2):rhs+=C[i]*v[j]*v[h]*sub(S.diff(phi,(xx,yy)[i],(xx,yy)[j],(xx,yy)[h]))/2
ck('vector_endpoint_full_adjoint',ge(lhs-rhs)==0)
defect=0
for a,be,v in zip(directions,beta,V):
    dchi=dd(chi,a)
    defect+=(dd(dchi,a)/2-be*dchi)*C.dot(grad)
    defect+=dchi*sum(C[i]*v[j]*sub(S.diff(phi,(xx,yy)[i],(xx,yy)[j])) for i in range(2) for j in range(2))
ck('vector_retained_test_adjoint',ge(chi*(lhs-rhs)-defect)==0)
ck('vector_naive_retention_really_fails',ge(chi*(lhs-rhs))!=0,{'missing_current_pairing':str(ge(chi*(lhs-rhs)))})
HA=lambda e:sum(dd(dd(e,a),a) for a in directions)/2
for i in range(2):ck('vector_scalar_product_'+str(i),ez(Q*HA(C[i])-HA(Q*C[i])+sum(dd(Q,a)*dd(C[i],a) for a in directions)+HA(Q)*C[i]))
# Live-width dependence must be restored, even at a frozen-center lift.
y,t=S.symbols('y t', real=True)
F=S.exp(-t*t/2)*S.cos(y); live_t=1+b*b/4
live=F.subs({y:b,t:live_t})
frozen=S.diff(F,y,2).subs({y:b,t:live_t})/2
chain=(S.diff(F,y,t)*S.diff(live_t,b)+S.diff(F,t,2)*S.diff(live_t,b)**2/2+S.diff(F,t)*S.diff(live_t,b,2)/2).subs({y:b,t:live_t})
ck('live_width_second_chain_rule',S.simplify(S.diff(live,b,2)/2-frozen-chain)==0)
ck('live_width_defect_nonzero',chain.subs(b,0)!=0)
# Exact integer factorial estimate far beyond author fixture, including high K.
for K in [0,1,2,7,12,24,50]:
 for j in [0,1,2,7,25,50,70]:
    ck(f'factorial_K{K}_j{j}',math.factorial(K+2*j+2)<=2**(K+2+4*j)*math.factorial(K+2)*math.factorial(j)**2)
# Deterministic random tensors: all proper cuts and one-Hilbert bound in multigraphs.
rng=np.random.default_rng(81123)
def cutnorm(T,indices):
    rest=tuple(i for i in range(T.ndim) if i not in indices)
    return float(np.linalg.norm(T.transpose(tuple(indices)+rest).reshape(2**len(indices),-1),2))
for fixture in range(12):
    N=3+fixture%2
    edges=[(i,i+1) for i in range(N-1)]+[(0,N-1),(0,1)]
    if fixture%3==0:edges.append((1,N-1))
    M=len(edges); tensors=[]; inds=[]; cs=[]
    for v in range(N):
        ids=[i for i,(u,w) in enumerate(edges) if v in (u,w)]+[M+v]
        T=rng.normal(size=(2,)*len(ids))
        cmax=max(cutnorm(T,I) for size in range(1,T.ndim) for I in itertools.combinations(range(T.ndim),size))
        T/=cmax; tensors.append(T); inds.append(ids);cs.append(1.)
    args=[]
    for T,I in zip(tensors,inds):args.extend([T,I])
    out=np.einsum(*args,list(range(M,M+N)))
    norm=max(cutnorm(out,I) for size in range(1,N) for I in itertools.combinations(range(N),size))
    ck('cyclic_parallel_all_cuts_'+str(fixture),norm<=1+1e-12,{'max_global_cut':norm})
    ck('cyclic_parallel_one_Hilbert_'+str(fixture),np.linalg.norm(out)<=np.linalg.norm(tensors[0])+1e-12)
# Analytic Fourier homotopy at positive private heat, including cancellation and degenerate endpoint.
mp.mp.dps=70
for xi in [[1,-1,0],[1,2,-3],[mp.mpf('0.1'),mp.mpf('0.2'),mp.mpf('0.3')],[4,-2,3]]:
    for rate in [[1,1,1],[0,1,2],[mp.mpf('.5'),2,1]]:
        for heat_s in [mp.mpf('.1'),mp.mpf('1')]:
            a=sum((r*x)**2 for r,x in zip(rate,xi))/2
            cross=sum(rate[i]*rate[j]*xi[i]*xi[j] for i in range(3) for j in range(i+1,3))
            base=mp.exp(-sum(x*x for x in xi)/2)
            P=lambda th:base*mp.exp(-heat_s*(a+th*cross))
            residual=P(0)-P(1)+heat_s*mp.quad(lambda th:-cross*P(th),[0,1])
            ck('fourier_homotopy_'+str((xi,rate,heat_s)),abs(residual)<mp.mpf('1e-60'))
# Nonpolynomial positive Gauss certificate on a product of heated C0 Hessians.
from numpy.polynomial.legendre import leggauss
for tau in [.25,1.]:
 for ratio in [.01,.125,1.,2.]:
    ss=ratio*tau*tau; Q=1.; K=0; N=3; p=3
    M=Q*2**(K+2)*math.sqrt(math.factorial(K+2))*p*tau**(-K-2)
    R0=tau*tau/(8*p*ss); J=max(1,math.ceil(1/R0)); eps=1e-8
    m=max(1,math.ceil(math.log(max(1,2*ss*M/eps),4)))
    gn,gw=leggauss(m)
    coeffs=[]
    # f_v(z)=1/2+1/4 exp(-tau^2 xi_v^2/2)cos(xi_v z), each Hessian in [1/4,3/4].
    freq=np.array([1.,-2.,3.])/tau
    for signs in itertools.product([-1,0,1],repeat=N):
        q=freq*np.array(signs); amplitude=math.prod(.5 if v==0 else .125*math.exp(-tau*tau*freq[i]**2/2) for i,v in enumerate(signs))
        independent=sum(q*q)/2; cross=sum(q[i]*q[j] for i in range(N) for j in range(i+1,N))
        coeffs.append((amplitude,independent,cross))
    def f(th):return sum(-cr*a*math.exp(-ss*(ind+th*cr)) for a,ind,cr in coeffs)
    qu=0.
    for j in range(J):
      for x,w in zip(gn,gw):qu+=w/(2*J)*f((j+.5+x/2)/J)
    p0=sum(a*math.exp(-ss*ind) for a,ind,cr in coeffs)
    p1=sum(a*math.exp(-ss*(ind+cr)) for a,ind,cr in coeffs)
    err=abs(p0-p1+ss*qu);cert=2*ss*M*4**(-m)
    ck(f'nonpoly_three_occurrence_tau{tau}_ratio{ratio}',err<=cert and cert<=eps,{'error':err,'certificate':cert,'panels':J,'m':m})
# Infinite fixed source lower bound: high-precision finite partial sum is already a lower bound.
eps=mp.mpf(1)/16;kap=mp.mpf(1)/4;u=mp.mpf('.7')
for n in [1,2,5,10,30,60]:
    f=mp.mpf(0);fw=mp.mpf(0);h2=mp.mpf(0)
    for j in range(1,n+1):
        ratio=mp.power(2,j*j-n*n);a=mp.mpf(1)/(j*j)
        summ=a*ratio**2*mp.exp(-ratio**2/2)
        f+=summ;fw+=summ*mp.exp(-u*ratio**2/2);h2+=summ*ratio**2
    # Positive partial sums of A-B and A+B separately bound the infinite product.
    heat=(kap*eps)**2*(f-fw)*(f+fw)
    lower=(kap*eps)**2/mp.mpf(n)**4*(mp.exp(-mp.mpf('.5'))-mp.exp(-(1+u)/2))*mp.exp(-mp.mpf('.5'))
    ck(f'lacunary_positive_partial_sum_n{n}',heat>=lower,{'heat_to_lower':str(heat/lower)})
    ck(f'lacunary_trace_n{n}',f*h2>=(1-mp.mpf('1e-65'))*mp.exp(-1)/mp.mpf(n)**4)
# Distinct fixed powers cannot dominate the lacunary sequence; record log-domain witnesses.
for delta in [.01,.1,1.]:
    n=1000
    ck('power_failure_delta'+str(delta),delta*n*n-4*math.log2(n)>9000)
report={'status':'PASS','test_count':len(checks),'checks':checks,'scope':'Independent algebra, multivariate Gaussian moments, graph-norm diagnostics, real positive quadrature and lacunary arithmetic. Native-source admission was audited against source contracts; no native sampler or numerical target program was executed.'}
(HERE/'independent-checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'status':report['status'],'test_count':len(checks)},indent=2))
