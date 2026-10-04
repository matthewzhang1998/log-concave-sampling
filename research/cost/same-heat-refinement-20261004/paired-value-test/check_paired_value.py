"""Numerical/algebraic checks of the scoped paired VALUE identities, not proofs."""
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite_e import hermegauss

out = {}

# Exact affine and formal-quadratic cancellation.
d1, d2 = .17, -.31
q, k0, k1, k2 = .23, -.2, .7, .4
b = lambda x: k0 + k1*x + k2*x*x/2
actual = 2*b(q+(d1+d2)/2)-(b(q+d1)+b(q+d2))/2
expected = b(q)+(k1+k2*q)*(d1+d2)/2+k2*d1*d2/2
out['quadratic_packet_identity_abs_error'] = float(abs(actual-expected))
assert abs(actual-expected)<1e-14

# Smooth strict-Hessian-sandwich near kink.
x, w = leggauss(220)
p = w/2
d1, d2 = x[:,None], x[None,:]
wp = p[:,None]*p[None,:]
a, c = .3, .2
limits=[]
for eta in [.1, .01, .001, .0001]:
    b=lambda z: a*(z/2+c*(np.sqrt(z*z+eta*eta)-eta))
    packet=2*b((d1+d2)/2)-(b(d1)+b(d2))/2
    val=float(np.sum(wp*packet)/(a*c))
    limits.append({'eta_over_epsilon':eta,'scaled_pair_mean':val})
out['near_kink_limit_target']=1/6
out['near_kink_values']=limits
assert abs(limits[-1]['scaled_pair_mean']-1/6)<2e-4

# Exact all-layer paired law-current identity for a non-polynomial force.
x,w=leggauss(60); p=w/2
d1=.2*x[:,None]; d2=.2*x[None,:]; A=(d1+d2)/2
wp=p[:,None]*p[None,:]
q=.15; a=.3; c=.2; eta=.08; weight=.7; B=.4
b=lambda z:a*(z/2+c*(np.sqrt(z*z+eta*eta)-eta))
bp=lambda z:a*(.5+c*z/np.sqrt(z*z+eta*eta))
phi=lambda y:np.exp(-y*y/2)
phip=lambda y:-y*np.exp(-y*y/2)
phipp=lambda y:(y*y-1)*np.exp(-y*y/2)
Y0=B-weight*b(q)
Y1=B-weight*(2*b(q+A)-(b(q+d1)+b(q+d2))/2)
lhs=float(np.sum(wp*phi(Y1))-phi(Y0))
J0=-weight*bp(q)*A
t,wt=leggauss(100);t=(t+1)/2;wt=wt/2
rank1=0.;rank2=0.
for u,wu in zip(t,wt):
    Y=B-weight*(2*b(q+u*A)-(b(q+u*d1)+b(q+u*d2))/2)
    J=-weight*(2*bp(q+u*A)*A-(bp(q+u*d1)*d1+bp(q+u*d2)*d2)/2)
    rank1+=wu*float(np.sum(wp*(J-J0)*phip(Y)))
    rank2+=wu*(1-u)*float(np.sum(wp*J0*J*phipp(Y)))
out['paired_current']={'direct':lhs,'rank_one':rank1,'rank_two':rank2,
                       'abs_residual':abs(lhs-rank1-rank2)}
assert abs(lhs-rank1-rank2)<1e-13

# Bounded odd C1 cap, identity on [-Bcap,Bcap], then C1 saturation.
def cap(z,Bcap=2.):
    az=np.abs(z);u=np.maximum(az-Bcap,0.)
    return np.sign(z)*np.where(az<=Bcap,az,
                              Bcap+np.where(u<=1,u-u*u/2,.5))
k,wk=hermegauss(200);pk=wk/np.sqrt(2*np.pi)
ck=cap(k)
r=float(pk@ck**2);s=float(pk@(k*ck))
sigma=.9;v=sigma/np.sqrt(2)
H=np.array([-.17,.02,.21,.31]);pH=np.array([.2,.3,.4,.1])
mu=float(pH@H);U=H-mu
S=(H[:,None]-H[None,:])**2/2
prob=pH[:,None]*pH[None,:]
gap=v*v*s*s/(2*r)
assert np.max(S)<gap
d=(-v*s+np.sqrt(v*v*s*s-r*S))/r
algebra=2*v*s*d+r*d*d+S
fillvar=2*v*v+2*v*s*d+r*d*d
mean=mu
var=float(pH@U**2+np.sum(prob*fillvar))
third_original=float(pH@U**3)
fourth_original=float(pH@U**4)-3*float(pH@U**2)**2
third_gauss=float(np.sum(prob*(U[:,None]**3+3*U[:,None]*(sigma*sigma-S))))
fourth_gauss=float(np.sum(prob*(U[:,None]**4+6*U[:,None]**2*(sigma*sigma-S)
                               +3*(sigma*sigma-S)**2)))-3*sigma**4
out['capped_fill']={'r':r,'s':s,'max_S':float(S.max()),'gap_limit':gap,
                    'quadratic_equation_abs_error':float(np.abs(algebra).max()),
                    'output_mean':mean,'input_mean':mu,'output_variance':var,
                    'target_variance':sigma*sigma,
                    'source_zero_d':float((-v*s+np.sqrt(v*v*s*s))/r)}
out['gaussian_fill_higher_moments']={
    'third':third_gauss,'minus_half_original_third':-.5*third_original,
    'fourth_cumulant':fourth_gauss,'minus_half_original_fourth_cumulant':-.5*fourth_original}
assert np.abs(algebra).max()<1e-15
assert abs(var-sigma*sigma)<1e-14
assert abs(third_gauss+.5*third_original)<1e-14
assert abs(fourth_gauss+.5*fourth_original)<1e-13

out['R4_ledger']={
    'sigma_exponent':1.25,'physical_protected_variance_exponent':3.5,
    'nearest_query_exponent':1.,'old_final_query_exponent':11/12,
    'optional_scalar_calibrated_final_query_exponent':2/3,
    'remote_M3_query_exponent':1/3,
    'calibrated_covariance_branch_normalized_energy_exponent':2.75,
    'calibrated_covariance_branch_physical_energy_exponent':3.25}
path=Path(__file__).with_name('paired_value_checks.json')
path.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
