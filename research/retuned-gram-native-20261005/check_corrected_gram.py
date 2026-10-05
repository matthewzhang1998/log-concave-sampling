#!/usr/bin/env python3
"""Finite identities and stress checks; does not replace native compiler proofs."""
import json, math, hashlib
from pathlib import Path
import numpy as np
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss
ROOT=Path(__file__).resolve().parent
nchecks=0
def ck(x,msg):
 global nchecks
 nchecks+=1
 if not x: raise AssertionError(msg)
def gh(n):
 x,w=hermgauss(n); return math.sqrt(2)*x,w/math.sqrt(math.pi)
x,w=gh(70)
r=np.array([.25,.7]); weights=np.array([.4,.6]); sig=math.sqrt(.75)
def B(z,A=1):
 z=np.asarray(z)
 return sum((A/2)*wi*ri*(1+.5*np.exp(-(1-ri*ri)/2)*np.cos(ri*sig*z)) for wi,ri in zip(weights,r))
def dB(z,A=1):
 return sum(-(A/4)*wi*ri*ri*sig*np.exp(-(1-ri*ri)/2)*np.sin(ri*sig*z) for wi,ri in zip(weights,r))
c=B(x)**2; dc=2*B(x)*dB(x)
mean=float(w@c); variance=float(w@(c-mean)**2)
ck(variance>0,'positive actual coarse variance')
rr,rw=leggauss(60); rr=(rr+1)/2;rw=rw/2
riesz=0.
for rho,wr in zip(rr,rw):
 xp=rho*x[:,None]+math.sqrt(1-rho*rho)*x[None,:]
 cp=2*B(xp)*dB(xp)
 riesz+=wr*np.sum(w[:,None]*w[None,:]*dc[:,None]*cp)
ck(abs(riesz-variance)<2e-13,'matrix-space covariance Riesz coefficient')
# Per-node four-path product-rule coefficients for two distinct primitive clocks.
rng=np.random.default_rng(417)
for k in range(400):
 G,H=rng.normal(size=2);rho=rng.uniform(0,1);Gp=rho*G+math.sqrt(1-rho*rho)*H
 b=lambda ri,z:(ri/2)*(1+.5*np.exp(-(1-ri*ri)/2)*np.cos(ri*sig*z))
 db=lambda ri,z:-(ri*ri*sig/4)*np.exp(-(1-ri*ri)/2)*np.sin(ri*sig*z)
 expanded=0.
 for i in range(2):
  for j in range(2):
   for k0 in range(2):
    for l in range(2):
     wt=weights[i]*weights[j]*weights[k0]*weights[l]
     expanded+=wt*(db(r[i],G)*b(r[j],G)+b(r[i],G)*db(r[j],G))*(db(r[k0],Gp)*b(r[l],Gp)+b(r[k0],Gp)*db(r[l],Gp))
 ck(abs(expanded-(2*B(G)*dB(G))*(2*B(Gp)*dB(Gp)))<1e-15,'four hit-tree terms and literal weights')
# Exact scalar Gaussian-Hermite feedback moments.
H3=x**3-3*x
moments={str(j):float(w@((x**i)*(H3**j))) for j,i in [(1,3),(2,2),(3,1),(4,0)]}
ck(abs(moments['1']-6)<1e-10,'P^3 H3=6')
ck(abs(moments['2']-42)<1e-10,'P^2 H3^2=42')
ck(abs(moments['3']-324)<1e-9,'P H3^3=324')
ck(abs(moments['4']-3348)<1e-8,'H3^4=3348')
feedback=[]
for A in [.03,.06,.12,.24]:
 u=.8;up=u/4;ug=u/2;uk=u/4;s=math.sqrt(up)
 C=A*A*c;cm=float(w@C);cv=float(w@(C-cm)**2);Q=-cv/(8*s**3)
 V=ug+uk+C; L=s*x+Q*H3
 EW2=float(w@V+w@(L*L)); targetvar=u+cm
 EW4=float(3*w@(V*V)+6*(w@V)*(w@(L*L))+w@(L**4))
 linear=3*cv+24*s**3*Q
 ck(abs(linear)<1e-16,'same-law fourth cancellation scalar coefficient -1/8')
 ck(abs(EW2-targetvar-6*Q*Q)<1e-14,'actual cubic covariance feedback charged')
 expected=6*(ug+uk+cm)*6*Q*Q+252*s*s*Q*Q+1296*s*Q**3+3348*Q**4
 ck(abs(EW4-3*targetvar**2-expected)<2e-12,'quartic feedback exact polynomial')
 lb=3*cv/(4*math.sqrt(15)*(u+float(C.max()))**1.5)
 feedback.append({'A':A,'u':u,'mixture_variance_noise':cv,'unchanged_mixture_W2_lower_bound':lb,'quartic_feedback_stable_polynomial':expected,'quadrature_feedback_cancellation_error':EW4-3*targetvar**2-expected,'quartic_feedback_over_A8':expected/A**8})
# Actual j_Q (one inner finite source) converges to baseline coefficient, not an artificial Hessian oracle.
X,WX=gh(28);Z,WZ=gh(32);HH,WH=gh(32); ri=.55;coarse_sig=.7;ti=.4
jrows=[]
for A in [.01,.02,.04,.08]:
 g=lambda z:A/2*(z+.5*np.sin(z));dg=lambda z:A/2*(1+.5*np.cos(z))
 bx=[]
 for xx in X:
  site=ri*coarse_sig*xx+math.sqrt(1-ri*ri)*Z[:,None]
  inner=ti*site+math.sqrt(1-ti*ti)*HH[None,:]
  nested=g(inner); derivative=dg(site-nested)*(1-ti*dg(inner))
  bx.append(ri*np.sum(WZ[:,None]*WH[None,:]*derivative))
 cv=np.asarray(bx)**2; vm=float(WX@(cv-float(WX@cv))**2)
 ck(vm>0,'actual j covariance noise positive')
 jrows.append({'A':A,'actual_var_C_over_A4':vm/A**4})
ck(max(z['actual_var_C_over_A4'] for z in jrows)/min(z['actual_var_C_over_A4'] for z in jrows)<1.5,'actual j noise retains A4 scale')
# All exponent comparisons are separate substantive rows.
for beta in [0,.8,1.4,1.5]:
 target=4-beta
 rows=[7-3*beta,4-beta,4-beta,4-.5*beta,4.5-1.25*beta,5-1.5*beta]
 for exponent in rows:ck(exponent+1e-12>=target,'retuned seven-order domination')
 for A in [1e-8,1e-5,.003]:
  u=A**beta;alpha=A/math.sqrt(u)
  ck(A**6/u**2.5<=1.00001*A**5/u**2,'rank6 below packet alpha5')
  ck(A**8/u**3.5<=1.00001*A**5/u**2,'cubic feedback below packet alpha5')
for beta,b in [(1.5,7),(1.8,13),(1.9,23)]:
 ck(b-beta*(b-1)/2>=4-beta-1e-12,'fixed order formula')
# Exact dyadic path-only inverse-shield sums and root/profile bounds.
path_panels=[]
for panels in [4,8,16,32,64,100]:
 delta=2.0**(-np.arange(1,panels+1,dtype=float))
 ww=delta.copy();sigma=np.sqrt(delta)
 s1=float(np.sum(ww/sigma));s2=float(np.sum(ww/sigma**2))
 q2=float(np.sum(ww**2/sigma**2));q3=float(np.sum(ww**2/sigma**3))
 ck(s1<2.5,'unsquared leaf inverse width integrates')
 ck(abs(s2-panels)<1e-12,'unsquared hit inverse width squared is logarithmic')
 ck(q2<1.01 and q3<2.5,'squared path feedback inverse widths sum')
 # Two hit slots, two leaves; feedback is fully factored over four clocks.
 z0=float(np.sum(ww**2));z1=float(np.sum(ww**2/sigma))
 feedback_bound=2*q3*q2*z0*z0+2*q2*q2*z1*z0
 ck(feedback_bound<4,'complete squared-weight old-root feedback remains bounded')
 ck(s2*s1*float(np.sum(ww))**2<2.5*panels,'complete hit caller only logarithmic')
 path_panels.append({'panels':panels,'leaf_sum':s1,'hit_sum':s2,'squared_feedback':feedback_bound})
# Positive finite Riesz rule: operator from Hilbert-valued input to Gaussian-root HS.
# Its chaos-n multiplier is sqrt(n)[1/n-sum w*rho^(n-1)].
riesz_rules=[]
for panels in [8,12,16,20]:
 gl,gw=leggauss(12);nodes=[];mass=[]
 for j in range(panels):
  lo=1-2.0**(-j);hi=1-2.0**(-(j+1))
  nodes.extend(((lo+hi)/2+(hi-lo)*gl/2).tolist());mass.extend(((hi-lo)*gw/2).tolist())
 tail=2.0**(-panels);nodes.append(1-tail/2);mass.append(tail)
 nodes=np.array(nodes);mass=np.array(mass)
 ns=np.unique(np.concatenate((np.arange(1,500),np.rint(np.geomspace(500,1e9,2000))))).astype(float)
 err=np.sqrt(ns)*np.abs(1/ns-np.exp((ns[:,None]-1)*np.log(nodes)[None,:])@mass)
 mx=float(err.max())
 ck(abs(float(mass.sum())-1)<1e-14,'Riesz rule mass exact positive')
 ck(np.all(mass>0) and np.all((nodes>0)&(nodes<1)),'Riesz interior positive nodes')
 ck(mx<math.sqrt(tail),'Riesz Hilbert operator error controlled by square-root cutoff')
 riesz_rules.append({'panels':panels,'nodes':len(nodes),'operator_error_sampled':mx,'sqrt_tail':math.sqrt(tail)})
res={'assertions':nchecks,'status':'PASS','scope':'Finite identities and numerical stress checks only; analytical packet audits are separate.','coarse_variance':variance,'covariance_Riesz_integral':float(riesz),'four_path_error':float(abs(riesz-variance)),'hermite_feedback_moments':moments,'corrected_mixture_scalar_rows':feedback,'actual_nested_j_rows':jrows,'unsmoothed_path_panel_checks':path_panels,'positive_Riesz_operator_checks':riesz_rules}
(ROOT/'corrected_gram_checks.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
