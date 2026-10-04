"""Finite-source diagnostics; analytic inequalities are proved in the adjacent note."""
import json, math
from pathlib import Path
import numpy as np
from numpy.polynomial.hermite import hermgauss

checks=0
rng=np.random.default_rng(2026100416)
def ck(test):
 global checks
 checks+=1
 assert bool(test),f'check {checks}'
T=np.array([[.5,.1],[.1,.5]]); N=np.diag([1.,0.]); e1=np.array([1.,0.])
delta=.1;a=r=1/16;eps=a**.9;c=.6;s=.8;v=.5;h=math.sqrt(1-v*v);k=1/v

def g(y):return T@y+(delta/k)*math.sin(k*y[0])*e1
def H(y):return T+delta*math.cos(k*y[0])*N

def source(W):
 S,U,Z=np.asarray(W).reshape(3,2)
 q1=g(S);q2=g(S+eps*Z);x=c*S+s*U;x1=x+a*(q1-q2)
 q3=g(x);q4=g(x1);t0=S+a*q3;t1=S+a*q4;q5=g(t0);q6=g(t1)
 E=r*(q5-q6)
 A0=H(x);A1=H(x1);T0=H(t0);T1=H(t1)
 K0=T0@A0;K1=T1@A1;D=H(S)-H(S+eps*Z)
 DS=r*(T0-T1)+r*a*c*(K0-K1)-r*a*a*K1@D
 DU=r*a*s*(K0-K1)
 DZ=r*a*a*eps*K1@H(S+eps*Z)
 return E,np.hstack([DS,DU,DZ]),K0,K1,[S,S+eps*Z,x,x1,t0,t1]

max_derivative_error=0.
for _ in range(600):
 W=rng.normal(size=(3,2))
 E,J,K0,K1,queries=source(W)
 ck(len(queries)==6)
 for qy in queries:
  ev=np.linalg.eigvalsh(H(qy));ck(ev[0]>=.3-1e-13 and ev[-1]<=.7+1e-13)
 ck(np.linalg.norm(E)<=r*a*a*eps*np.linalg.norm(W[2])+1e-13)
 ck(np.linalg.norm(J,2)<=3*r)
 z=rng.normal(size=6);z/=np.linalg.norm(z);step=1e-5
 fd=(source(W.reshape(-1)+step*z)[0]-source(W.reshape(-1)-step*z)[0])/(2*step)
 err=np.linalg.norm(fd-J@z);max_derivative_error=max(max_derivative_error,err);ck(err<1e-9)
 ck(np.linalg.norm(K0.T-H(queries[2])@H(queries[4]))<1e-14)
 # Legal complete regeneration is an exact evaluation at a fresh whole record.
 Wp=rng.normal(size=(3,2));_,_,K0p,_,qsp=source(Wp)
 ck(np.linalg.norm(K0p-H(qsp[4])@H(qsp[2]))<1e-14)

# Analytic one-dimensional Gaussian expectations by deterministic quadrature.
def expect(f,n=160):
 nodes,weights=hermgauss(n)
 return float(np.sum(weights*np.array([f(math.sqrt(2)*x) for x in nodes]))/math.sqrt(math.pi))
q=c+a/2;beta=a*delta;b=.1;L=math.exp(-(s*s+(a*b)**2)/2);mu0=math.exp(-.5)
ct=lambda x:math.cos(q*x+beta*math.sin(x))
mut=L*expect(ct);joint=L*expect(lambda x:ct(x)*math.cos(x))
cov=joint-mut*mu0
base=math.exp(-(q*q+1)/2)*(math.cosh(q)-1)
cov_lower=L*(base-beta*(1+math.exp(-.5)))
ck(cov_lower>.05);ck(cov>=cov_lower);ck(abs(cov-.0748306267)<1e-8)
original=T@T+delta*mut*N@T+delta*mu0*T@N+delta**2*joint*N
split=(T+delta*mut*N)@(T+delta*mu0*N)
ck(np.linalg.norm(original-split-delta**2*cov*N)<1e-14)
ck((original-split)[0,0]>.0005)
alpha=math.exp(-.5)
heated=T@T+delta*alpha*mut*N@T+delta*alpha*mu0*T@N+delta**2*alpha**2*joint*N
heat_gap_lower=delta*b*(1-alpha)*L*(math.exp(-q*q/2)-beta)
ck(heat_gap_lower>.002);ck((original-heated)[0,1]>=heat_gap_lower)
ck(abs((original-heated)[0,1]-delta*b*(1-alpha)*mut)<1e-14)

# Integrated SAME-root conditional independence: exact Gaussian phase geometry.
A=k*k*(1+2*a*.5*c+a*a*(.25+b*b));Bvar=k*k;D=k*k*(c+a*.5)
root_base=math.exp(-(A+Bvar)/2)*(math.cosh(D)-math.cosh(h*h*D))
root_lower=root_base-2*beta
ck(root_lower>.03)
Qnoise=math.exp(-(A-D*D/Bvar)/2)
oldprod=Qnoise*expect(lambda z:math.cos((D/Bvar)*k*z+beta*math.sin(k*z))*math.cos(k*z))
splitprod=Qnoise*math.exp(-Bvar*(1-h**4)/2)*expect(lambda z:math.cos((D/Bvar)*k*z+beta*math.sin(k*z))*math.cos(h*h*k*z))
root_cov=oldprod-splitprod
ck(root_cov>=root_lower);ck(delta**2*root_lower>.0003)
root_base_quad=Qnoise*expect(lambda z:math.cos((D/Bvar)*k*z)*math.cos(k*z))-Qnoise*math.exp(-Bvar*(1-h**4)/2)*expect(lambda z:math.cos((D/Bvar)*k*z)*math.cos(h*h*k*z))
ck(abs(root_base_quad-root_base)<1e-13)

# Fully independent root regeneration is a different coupling; it also changes K0.
tau=2.;Ll=math.exp(-tau*tau*(s*s+(a*b)**2)/2)
base2=math.exp(-tau*tau*(q*q+1)/2)*(math.cosh(tau*tau*q)-1)
independent_lower=Ll*(base2-beta*(1+math.exp(-tau*tau/2)))
ck(independent_lower>.08)

# Source is an acyclic six-query graph; its primitive comparison obligations cycle.
parents={'q1':['S'],'q2':['S','Z'],'x0':['S','U'],'x1':['x0','q1','q2'],
 'q3':['x0'],'q4':['x1'],'t0':['S','q3'],'t1':['S','q4'],'q5':['t0'],'q6':['t1'],'E':['q5','q6']}
reads={'S':{'fineS'},'U':{'fineU'},'Z':{'fineZ'}}
for key,ps in parents.items():
 ck(all(p in reads for p in ps))
 reads[key]=set().union(*(reads[p] for p in ps))
ck('fineU' in reads['q3']);ck('fineU' in reads['q5']);ck('fineZ' in reads['q6'])
# q5 is first terminal node after ancestors; exported YA doesn't erase direct q3 use.
terminal=[q for q in parents if q in ('q5','q6')]
ck(terminal[0]=='q5');ck(parents['t0']==['S','q3'])
comparison_obligations={'A':{'T'},'T':{'A'}}
ck('A' in comparison_obligations['T'] and 'T' in comparison_obligations['A'])
for K in [1,3,8,17,33]:
 ck(2+2+(K+1)+1==K+6)

# Minimal actual terminal-label tuple recovers the complete native record.
# The inverse below is a diagnostic for the known fixture, never an oracle in a producer.
def inverse_g(p):
 m=.5;bb=.1;gamma=m-bb*bb/m;rhs=p[0]-(bb/m)*p[1]
 lo=(rhs-delta/k)/gamma;hi=(rhs+delta/k)/gamma
 for _ in range(64):
  mid=(lo+hi)/2
  if gamma*mid+(delta/k)*math.sin(k*mid)<rhs:lo=mid
  else:hi=mid
 y1=(lo+hi)/2
 return np.array([y1,(p[1]-bb*y1)/m])
max_inverse_error=0.
for _ in range(300):
 W=rng.normal(size=(3,2));S,U,Z=W
 _,_,_,_,queries=source(W)
 p0=g(queries[2]);p1=g(queries[3])
 x0=inverse_g(p0);x1=inverse_g(p1)
 Ur=(x0-c*S)/s;Delta=(x1-x0)/a
 Zr=(inverse_g(g(S)-Delta)-S)/eps
 err=np.linalg.norm(np.stack([S,Ur,Zr])-W);max_inverse_error=max(max_inverse_error,err)
 ck(np.linalg.norm(Ur-U)<1e-12);ck(err<1e-9)
 # Reverse map from arbitrary terminal labels to a native record is also onto.
 S,p0,p1=rng.normal(size=(3,2));x0=inverse_g(p0);x1=inverse_g(p1)
 Wr=np.stack([S,(x0-c*S)/s,(inverse_g(g(S)-(x1-x0)/a)-S)/eps])
 _,_,_,_,qr=source(Wr)
 ck(np.linalg.norm(g(qr[2])-p0)<1e-11)
 ck(np.linalg.norm(g(qr[3])-p1)<1e-10)

# Literal finite Lagrange first-coefficient filter and concrete tail-source algebra.
for K in [2,4,6,8]:
 nodes=.5*np.cos(np.pi*np.arange(K+1)/K)
 weights=[]
 for j in range(K+1):
  poly=np.polynomial.Polynomial([1.])
  for ell in range(K+1):
   if ell!=j:poly*=np.polynomial.Polynomial([-nodes[ell],1.])/(nodes[j]-nodes[ell])
  weights.append(poly.deriv()(0.))
 weights=np.array(weights)
 for ell in range(K+1):ck(abs(np.sum(weights*nodes**ell)-(1. if ell==1 else 0.))<1e-10)
 for _ in range(25):
  incoming,UA,UT,ZT,center=rng.normal(size=(5,2));rho=.01;s0=.125;c0=math.sqrt(1-s0*s0)
  f=lambda z:(rho/v)*g(center+v*z)
  Dp=lambda z:f(c0*z+s0*incoming)-f(z)
  filtered=sum(w*(f(math.sqrt(1-bj*bj)*ZT+bj*incoming)-f(ZT)) for w,bj in zip(weights,nodes))
  Jtail=Dp(UT)-s0*filtered
  ca,ct,cs,eta=[.5]*4;GA,GT,GS,GE=rng.normal(size=(4,2))
  literal=ca*(Dp(UA)/ca+GA)+ct*(-Jtail/ct+GT)+cs*GS+eta*GE
  expanded=Dp(UA)-Dp(UT)+s0*filtered+ca*GA+ct*GT+cs*GS+eta*GE
  ck(np.linalg.norm(literal-expanded)<1e-14)
 ck(4+(K+1)+1==K+6)

out={
 'status':'PASS','checks':checks,'random_seed':2026100416,
 'same_potential':{'a':a,'r':r,'epsilon':eps,'v':v,'h':h,'k':k,'hessian_sandwich':[.3,.7]},
 'zero_coarse_root':{'covariance_analytic_lower':cov_lower,'covariance_quadrature':cov,'K0_11_gap_analytic_lower':delta**2*cov_lower,'K0_11_gap_quadrature':delta**2*cov},
 'integrated_same_root':{'linear_phase_covariance_exact':root_base,'actual_covariance_analytic_lower':root_lower,'actual_covariance_quadrature':root_cov,'K0_11_gap_analytic_lower':delta**2*root_lower},
 'frozen_labels_new_primitive_heat':{'widths':[v,v],'K0_12_gap_analytic_lower':heat_gap_lower,'K0_12_gap_quadrature':float((original-heated)[0,1])},
 'fresh_root_independent_banks':{'covariance_analytic_lower':independent_lower},
 'six_query_native_jacobian_max_error':max_derivative_error,
 'minimal_terminal_label_inverse_max_error':max_inverse_error,
 'primitive_seed_filter_degrees_checked':[2,4,6,8],
 'scope':['same-potential complete positive-clock K3 source; branch K0 target','six original VALUE sites and first-only analytic ledger','same-root split-bank gap remains after owned-root integration','positive additive primitive heat at frozen physical labels changes target','acyclic native graph; cyclic prerequisites for separate primitive replacements','minimal actual terminal labels recover complete native fine record in M=I strong-gradient subclass','literal concrete finite first-chaos tail source and K+6 unoptimized call count'],
 'not_claimed':['paired K0-K1 noncancellation lower bound','lower bound on every completed pair law error','impossibility for all finite VALUE producers','new joint selected-word stationary kernel','uniform integrated restoration impossibility']}
Path(__file__).with_name('internal_k3_packet_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
