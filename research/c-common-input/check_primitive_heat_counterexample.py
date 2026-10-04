# Publication copy; original/public SHA-256 mapping is in INVENTORY.json.
# New diagnostic; not a byte-identical copy of an earlier checker.
import numpy as np, math, json
from pathlib import Path
from scipy.integrate import quad
from scipy.special import gammaln
EPS=1e-8;C=100.
dh=(math.exp(-.2)-1)/math.sqrt(5)
lead=math.exp(-1)*abs(dh)/(2*math.sqrt(2))
rem=4*C*C*EPS
assert lead>.0105 and rem<.0005
assert lead/2-rem>.004
assert 2+16*math.sqrt(2)+68<100
rows=[]
for m in [100,1000,10000,100000,1000000]:
 sm=math.sqrt(m);s=m**(-.25)
 def density(t):
  x=m+sm*t
  if x<=0:return 0.
  return math.exp((m/2-1)*math.log(x)-x/2-(m/2)*math.log(2)-gammaln(m/2)+math.log(sm))
 h0=quad(lambda t:math.exp(-t*t)*density(t),max(-12,-sm),12,epsabs=1e-10)[0]
 hs=quad(lambda t:math.exp(-((1+s*s)*t+s*s*sm)**2)*density(t),max(-12,-sm),12,epsabs=1e-10)[0]
 a=.001;ba=math.exp(-((1+a)**2+1)/2)
 coef=ba*abs(math.exp(-s*s)*hs-h0)/(2*math.sqrt(2))
 assert coef-rem>.004
 rows.append({'m':m,'s':s,'h0':h0,'hs':hs,'main_coefficient_div_eps_kappa_squared':coef,'lower_after_full_Riesz_remainder':coef-rem})
assert abs(rows[-1]['h0']-1/math.sqrt(5))<1e-4
assert abs(rows[-1]['hs']-math.exp(-.2)/math.sqrt(5))<1e-3
rng=np.random.default_rng(73405);worst=0
for _ in range(1000):
 a=rng.uniform(.0001,.2);h0=rng.uniform(.1,.6);hshift=rng.uniform(.1,.6);s=rng.uniform(.001,.1)
 b0=math.exp(-1);ba=math.exp(-((1+a)**2+1)/2)
 D=np.diag([1.,0.]);K=np.array([[0.,-1.],[1.,0.]])
 # The exact constant Riesz component is .5D.
 delta=ba*(math.exp(-s*s)*hshift-h0)*K
 main=-(D@delta+(D@delta).T)/4
 target=ba*abs(math.exp(-s*s)*hshift-h0)/(2*math.sqrt(2))
 er=abs(np.linalg.norm(main,'fro')-target);worst=max(worst,er);assert er<1e-14
out={'status':'PASS','provenance':'New diagnostic for exact preserved source a025c546; no historical test count is reused.','analytic_limit_coefficient_div_eps_kappa_squared':lead,'uniform_full_Riesz_remainder_div_eps_kappa_squared':rem,'certified_eventual_lower_coefficient':.004,'radial_quadrature_rows':rows,'matrix_identity_cases':1000,'max_matrix_identity_error':worst,'scope':'The proof uses CLT and a uniform full-mark remainder bound. These numerical radial integrals and matrix checks do not replace that proof or imply an actual CW7 batch admission.'}
Path(__file__).with_name('primitive_heat_counterexample_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
