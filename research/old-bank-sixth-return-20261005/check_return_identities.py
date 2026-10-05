import numpy as np, math, json
from numpy.polynomial.hermite_e import hermegauss, hermevander
from numpy.polynomial.legendre import leggauss

def gauss(n):
 x,w=hermegauss(n); return x,w/math.sqrt(2*math.pi)
q,wq=gauss(180); p,wp=gauss(130)
r,wr=leggauss(90); r=(r+1)/2;wr=wr/2
Rsin=np.sum(wr[None,:]*np.exp(-(1-r*r)[None,:]/2)*np.cos(q[:,None]*r[None,:]),axis=1)
v=(1-math.exp(-2))/2
h=Rsin*np.cos(q)-v
# L2 Gaussian Hermite expansion. R H_n/sqrt(n!)=H_(n-1)/sqrt((n-1)! n).
V=hermevander(q,40)/np.sqrt(np.array([math.factorial(i) for i in range(41)],dtype=float))[None,:]
coef=V.T@(wq*h)
Rh=sum(coef[n]*V[:,n-1]/math.sqrt(n) for n in range(1,41))
checks=[]
for a in [.03,.11,.24]:
 for t in [.0,.2,.6,1.0]:
  for om in [.3,.8,1.4]:
   H=p*p-1; base=.45*p +.07*(p*p-1); kap=.65
   F=a*np.sin(q[:,None])*H[None,:]
   C=a*a*H*H*v
   phase=np.exp(1j*om*(base[None,:]+t*F)-.5*om*om*(kap+(1-t*t)*C[None,:]))
   def ev(x):return np.sum(wq[:,None]*wp[None,:]*x)
   direct=ev((1j*om*F+om*om*t*C[None,:])*phase)
   one=ev(t*a*a*h[:,None]*H[None,:]**2*(1j*om)**2*phase)
   two=ev(t*t*a**3*Rh[:,None]*np.cos(q[:,None])*H[None,:]**3*(1j*om)**3*phase)
   err=max(abs(direct-one),abs(one-two))
   checks.append({'a':a,'t':t,'omega':om,'error':float(err)})
assert max(c['error'] for c in checks)<2e-10
# H2^2=H4+4H2+2; Gaussian current then has 1/2,2,1 coefficients.
x=np.linspace(-5,5,101);H2=x*x-1;H4=x**4-6*x*x+3
assert np.max(abs(H2*H2-H4-4*H2-2))<1e-12
# Native adapter scalar normalization for all 36 physical permutation pairs.
s0=.1;b0=s0/math.sqrt(2);b2=s0/2;d=32/(9*s0**6)
assert abs(36*d*b0**4*b2**2-8)<1e-12
# Positive dyadic original center and new bridge mass bounds, retained exact weights.
clock=[]
for J in [8,16,32,64,100]:
 Delta=2.**(-np.arange(1,J+1)); w=Delta.copy();sig2=Delta
 for m in [1.,2.**-10,2.**-30]:
  h=np.sqrt(Delta);u=Delta.copy()
  S0=S1=0.
  for hh,uu in zip(h,u):
   tc2=sig2+hh*hh*m
   S0+=uu*np.sum(w*w/tc2**2)
   S1+=uu*np.sum(w*w/tc2**2.5)*math.sqrt(m)
  clock.append({'J':J,'m':m,'mass':float(sum(u)),'u_over_h':float(sum(u/h)),'weighted_beta':float(S0),'weighted_derivative_times_sqrt_m':float(S1)})
assert max(c['u_over_h'] for c in clock)<2.42
assert max(c['weighted_derivative_times_sqrt_m'] for c in clock)<8
out={'two_riesz_identity_max_error':max(c['error'] for c in checks),'identity_cases':checks,'wick_covariance_trace_coefficients':[.5,2.,1.],'native_normalization_sum':float(36*d*b0**4*b2**2),'clock_cases':clock,'scope':'Finite scalar diagnostics for the two-Riesz positive interpolation, Wick trace cancellation, exact adapter scalar and clock ledgers. Imported native programs are not numerically executed.'}
open('/workspace/shared/old-bank-sixth-return-20261005/return_identity_checks.json','w').write(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['two_riesz_identity_max_error','wick_covariance_trace_coefficients','native_normalization_sum']},indent=2))
