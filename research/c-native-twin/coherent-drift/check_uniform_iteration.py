# Publication copy; original/public SHA-256 mapping is in INVENTORY.json.
from pathlib import Path
import numpy as np,json
rng=np.random.default_rng(41238);lam=.5;de=be=.025;c=s=2**-.5;checks=0;rows=[]
for N in [200,1000,5000]:
 a=N**(-5/9);ep=N**(-.5);dim=int(round((2*np.pi*N/(a*de))**2));df=float(dim-1);n=12000
 S,U,Z=rng.normal(size=(3,n));x=c*S+s*U
 l00=np.sqrt(rng.chisquare(df,n));l10=rng.normal(size=n);l20=rng.normal(size=n);l11=np.sqrt(rng.chisquare(df-1,n));l21=rng.normal(size=n);l22=np.sqrt(rng.chisquare(df-2,n))
 SS=l00*l00;UU=l10*l10+l11*l11;ZZ=l20*l20+l21*l21+l22*l22;SU=l00*l10;SZ=l00*l20
 xx=c*c*SS+s*s*UU+2*c*s*SU
 def feedback(q):
  R=np.sqrt(1+SS+q*q);Rp=np.sqrt(1+SS+2*ep*SZ+ep*ep*ZZ+(q+ep*Z)**2)
  rd=(2*ep*(q*Z+SZ)+ep*ep*(Z*Z+ZZ))/(Rp+R)
  return a*(-lam*ep*Z-de*rd+de*(q*q/R-(q+ep*Z)**2/Rp)+be*(np.sin(q)-np.sin(q+ep*Z)))
 def L(h):
  y=x+h;R=np.sqrt(1+xx+y*y)
  g=lam*y+de*(R-1+y*y/R)+be*np.sin(y)
  return feedback(S+a*g)
 zero=np.zeros(n);h={1:zero,2:feedback(S)};L0=L(zero);q=.72*a*a;worst=0.
 for j in range(3,42):
  h[j]=L(h[j-2]);err=np.abs(h[j]-L0)
  assert np.max(err-q*np.abs(h[j-2]))<5e-11
  assert np.max(np.abs(h[j])-.6*a*ep*np.sqrt(Z*Z+ZZ))<1e-10
  worst=max(worst,float(np.sqrt(np.mean((h[j]+np.pi)**2))));checks+=2
 assert worst<.012;checks+=1
 rows.append({'N':N,'maximum_h_j_plus_pi_L2_for_3_to_41':worst,'h2_plus_pi_L2':float(np.sqrt(np.mean((h[2]+np.pi)**2))),'supremum_Lh_minus_L0_L2':max(float(np.sqrt(np.mean((h[j]-L0)**2))) for j in range(3,42))})
out={'checks':checks,'rows':rows,'scope':'Literal two-step finite twin recurrence on exact radial Gaussian sufficient statistics; no iteration is replaced by a law reference.'}
print(json.dumps(out,indent=2));open(str(Path(__file__).with_name('uniform_iteration_checks.json')),'w').write(json.dumps(out,indent=2)+'\n')
