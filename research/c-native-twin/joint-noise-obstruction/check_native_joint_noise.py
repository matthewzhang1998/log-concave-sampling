# Publication copy; original/public SHA-256 mapping is in INVENTORY.json.
from pathlib import Path
import numpy as np,json
from numpy.polynomial.hermite import hermgauss
from scipy.optimize import least_squares
z,w=hermgauss(120);z=z*np.sqrt(2);w=w/np.sqrt(np.pi)
m=.5;d=.1;lo=m-d;hi=m+d
vo=(1-np.exp(-2))/2-np.exp(-1)
ve=(1+np.exp(-2))/2-np.exp(-1)
ck=d*np.sqrt(vo); checks=0; maxerr=0.; worst=np.inf
for theta in np.linspace(0,2*np.pi,401):
 k=-m*z-d*(np.sin(theta+z)-np.sin(theta))
 mu=d*np.sin(theta)*(1-np.exp(-.5));sigma=m+d*np.cos(theta)*np.exp(-.5)
 val=np.dot(w,(k-(mu-sigma*z))**2)
 exact=d*d*(np.cos(theta)**2*vo+np.sin(theta)**2*ve)
 maxerr=max(maxerr,abs(val-exact));assert abs(val-exact)<1e-13
 assert val>=ck*ck-1e-13; checks+=2
for eps in [1.,.1,.003]:
 for a in [.01,.04]:
  r=.03;rho=a*eps
  g=lambda x:m*x+d*eps*np.sin(x/eps)
  H=lambda x:m+d*np.cos(x/eps)
  for theta in [0.,.5,1.57,2.4]:
   S=eps*theta
   for U in [-.7,.4]:
    x=(S+U)/np.sqrt(2)
    k=-rho*(m*z+d*(np.sin(theta+z)-np.sin(theta)))
    def action(t):return r*(g(S+a*g(x))-g(S+a*g(x+t)))
    E=action(k); mu=rho*d*np.sin(theta)*(1-np.exp(-.5));sig=rho*(m+d*np.cos(theta)*np.exp(-.5))
    # Optimize the Gaussian inserted BEFORE the same nonlinear paired observable.
    def residual(p):return np.sqrt(w)*(action(rho*(p[0]-np.exp(p[1])*z))-E)/(r*a*rho)
    fit=least_squares(residual,[mu/rho,np.log(sig/rho)],xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=200)
    err=np.linalg.norm(residual(fit.x));worst=min(worst,err)
    assert err>=lo*lo*ck-2e-7
    assert np.max(np.abs(E)-hi**3*r*a*rho*np.abs(z))<1e-13
    deriv=r*a*rho*H(S+a*g(x+k))*H(x+k)*H(S+eps*z)
    assert np.min(deriv)>=lo**3*r*a*rho*(1-1e-12)
    checks+=3
out={'checks':checks,'odd_residual_variance':vo,'even_residual_variance':ve,'uniform_shift_gaussian_W2_constant':ck,'uniform_joint_pair_constant':ck/np.sqrt(2),'paired_E_constant_in_units_r_a_rho':lo*lo*ck,'paired_E_constant_relative_to_actual_e_upper':lo*lo*ck/hi**3,'max_projection_identity_error':maxerr,'smallest_optimized_paired_E_ratio':worst,'scope':'One original strongly convex primitive; exact conditional Gaussian insertion only. Gaussian mixtures or weaker same-E coefficient comparisons are not ruled out.'}
print(json.dumps(out,indent=2));open(str(Path(__file__).with_name('native_joint_noise_checks.json')),'w').write(json.dumps(out,indent=2)+'\n')
