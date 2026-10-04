# Publication copy; original/public SHA-256 mapping is in INVENTORY.json.
from pathlib import Path
import numpy as np,json
rng=np.random.default_rng(41237)
lam=.5;delta=beta=.025;c=s=2**-.5
cU=2*beta*s*(lam*np.exp(-.5)+beta*np.exp(-1)*np.cosh(c))
rows=[]; count=180000
for N in [200,1000,5000,20000,100000]:
 a=N**(-5/9);eps=N**(-.5); D=int(round((2*np.pi*N/(a*delta))**2)); df=float(D-1)
 S,U,Z=rng.normal(size=(3,count));x=c*S+s*U
 # Exact 3x3 Wishart Gram of the remaining S,U,Z coordinates via Bartlett.
 l00=np.sqrt(rng.chisquare(df,size=count));l10=rng.normal(size=count);l20=rng.normal(size=count)
 l11=np.sqrt(rng.chisquare(df-1,size=count));l21=rng.normal(size=count);l22=np.sqrt(rng.chisquare(df-2,size=count))
 SS=l00*l00;UU=l10*l10+l11*l11;ZZ=l20*l20+l21*l21+l22*l22
 SU=l00*l10;SZ=l00*l20
 xx=c*c*SS+s*s*UU+2*c*s*SU
 RS=np.sqrt(1+S*S+SS);RSp=np.sqrt(1+(S+eps*Z)**2+SS+2*eps*SZ+eps*eps*ZZ)
 rad_diff=(2*eps*(S*Z+SZ)+eps*eps*(Z*Z+ZZ))/(RSp+RS)
 h=a*(-lam*eps*Z-delta*rad_diff+delta*(S*S/RS-(S+eps*Z)**2/RSp)+beta*(np.sin(S)-np.sin(S+eps*Z)))
 RX=np.sqrt(1+x*x+xx);RXh=np.sqrt(1+(x+h)**2+xx)
 rdiff=-h*(2*x+h)/(RX+RXh)
 d0=-lam*h+delta*(rdiff-h*(2*x+h)/RX-(x+h)**2*rdiff/(RX*RXh))+beta*(np.sin(x)-np.sin(x+h))
 K=a*delta*np.sqrt(float(D))
 phase=(K-2*np.pi*N)+S+a*(lam*x+delta*(RX-np.sqrt(float(D))-1+x*x/RX)+beta*np.sin(x))
 t0=2*np.pi*N+phase;t1=t0-a*d0
 R0=np.sqrt(1+SS+t0*t0);R1=np.sqrt(1+SS+t1*t1)
 hrad=(t0+t1)/(R0+R1)+(t0+t1)/R0-t1*t1*(t0+t1)/((R0+R1)*R0*R1)
 mid=phase-a*d0/2
 E0=d0*(lam+delta*hrad+beta*np.cos(mid)*np.sinc((a*d0/2)/np.pi))
 L=1/R0-t1*(t0+t1)/((R0+R1)*R0*R1)
 bulk2=delta**2*d0*d0*L*L*SS
 Flim=(lam*np.pi+2*beta*np.sin(x))*(lam+beta*np.cos(S))
 row={'N':N,'a':float(a),'epsilon':float(eps),'dimension':str(D),'mean_shift':float(h.mean()),'shift_L2_error':float(np.sqrt(np.mean((h+np.pi)**2))),'coordinate_L2_limit_error':float(np.sqrt(np.mean((E0-Flim)**2))),'hidden_first_chaos':float(np.mean(U*(E0-E0.mean()))),'limit_first_chaos':float(cU),'energy_over_ra':float(np.sqrt(E0.var()+bulk2.mean())),'samples':count}
 assert abs(row['mean_shift']+np.pi)<.01
 assert row['energy_over_ra']<.2 and row['energy_over_ra']>.01
 assert abs(row['hidden_first_chaos']-cU)<.0006
 rows.append(row)
out={'rows':rows,'analytic_positive_first_chaos':float(cU),'orientation_limit_lower_bound':float(cU*cU),'scope':'Exact low-dimensional sampling of native radial sufficient statistics; proof uses Gaussian concentration and heat-gradient isometry, not this Monte Carlo diagnostic.'}
print(json.dumps(out,indent=2));open(str(Path(__file__).with_name('coherent_drift_checks.json')),'w').write(json.dumps(out,indent=2)+'\n')
