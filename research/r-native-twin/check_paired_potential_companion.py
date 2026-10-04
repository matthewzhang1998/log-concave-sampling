from pathlib import Path
import numpy as np,json,math
rng=np.random.default_rng(812625);checks=0;records=[]
m=.5;d=.1;L=.6;c=.6;s=.8
# Deterministic envelope proof and exact Gaussian leading norm.
for a in np.geomspace(1e-6,1/16,60):
 for eps in [1.,a**.9,.01]:
  for fac in [1.,2.,8.]:
   k=fac/(a*eps);t=m*k*a*eps
   leading2=d*d*m*m*(1-math.exp(-t*t/2)+(1-2*math.exp(-t*t/2)+math.exp(-2*t*t))/2*(1-4*c*c*k*k)*math.exp(-2*k*k))
   error=.09*a+.1296*a*a*eps
   assert math.sqrt(leading2)-error>.01;checks+=1
# Same-potential actual native companion, sampled only as diagnostics.
S,U,Z=rng.normal(size=(3,160000));x=c*S+s*U
for a in [1/16,.03,.01]:
 eps=a**.9;r=a
 for fac in [1.,3.,11.]:
  k=fac/(a*eps)
  def g(y):return m*y+d*np.sin(k*y)/k
  def H(y):return m+d*np.cos(k*y)
  def Hp(y):return -d*k*np.sin(k*y)
  delta=g(S)-g(S+eps*Z);x1=x+a*delta;t0=S+a*g(x);t1=S+a*g(x1)
  C=H(x)*g(t0)-H(x1)*g(t1)
  E=(g(t0)-g(t1))/(a*a*eps)
  phi=k*x;t=m*k*a*eps
  lead=d*m*S*(np.cos(phi)-np.cos(phi-t*Z))
  err=float(np.sqrt(np.mean((C-lead)**2)))
  measured=float(np.sqrt(np.mean(C*C)))
  mark=float(np.std(E));assert measured>.01;assert mark<L**3+.001;checks+=2
  # Uniform first witness, with the exact same phase perturbation retained.
  SS=1.;xx=np.pi/(2*k);ZZ=np.pi/(m*k*a*eps)
  dd=g(SS)-g(SS+eps*ZZ);xx1=xx+a*dd;tt0=SS+a*g(xx);tt1=SS+a*g(xx1)
  derivative=Hp(xx)*g(tt0)-Hp(xx1)*g(tt1)+a*H(xx)**2*H(tt0)-a*H(xx1)**2*H(tt1)
  ratio=abs(float(derivative))/k;assert ratio>.09;checks+=1
  records.append({'a':a,'epsilon':eps,'frequency_multiple':fac,'companion_over_ra':measured,'actual_e_over_ra2epsilon':mark,'measured_L2_approx_error':err,'analytic_error_allowance':.09*a+.1296*a*a*eps,'first_over_ra_k':ratio})
out={'status':'PASS','checks':checks,'exact_proof_lower_companion_over_ra':.01,'records':records,'scope':'Deterministic Gaussian phase/envelope and exact first witnesses; Monte Carlo entries are diagnostics only.'}
Path(__file__).with_name('paired_potential_companion_checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':'PASS','checks':checks,'min_measured_companion_over_ra':min(v['companion_over_ra'] for v in records),'min_first_over_ra_k':min(v['first_over_ra_k'] for v in records)},indent=2))
