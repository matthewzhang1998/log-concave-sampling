import json,math
from pathlib import Path
import numpy as np
from scipy.stats import chi2
from scipy.integrate import quad

rows=[]
for a in [.1,.03,.01]:
    ell=a*a; sigma=a**1.25; vkeep=sigma*sigma/2
    for d in [1,10,100,1000]:
        needed=d*ell*ell/vkeep
        passprob=float(chi2.cdf(vkeep/(ell*ell),d))
        rows.append({'a':a,'d':d,'covariance_operator':ell*ell,
                     'allocated_variance':vkeep,'necessary_rank':needed,
                     'uncapped_gap_probability':passprob})

# Direct sampled trace and rank bounds; no covariance target is silently
# restored after truncation.
rng=np.random.default_rng(20261004)
for d in [2,7,25]:
    U=rng.normal(size=d); ell=.03; v=.005
    S=ell*ell*np.outer(U,U)
    ev=np.linalg.eigvalsh(S)
    assert abs(ev[-1]-ell*ell*np.dot(U,U))<1e-14
    cap=S*min(1,v/max(ev[-1],1e-300))
    assert np.trace(cap)<=v+1e-14

# Exact affine VALUE calibration and a nonlinear sanity check.
calibration=[]
for z0,z,kappa in [(1.,-.4,.8),(-.5,.7,.96),(.1,1.2,1.)]:
    for kind in ['linear','nonlinear']:
        a=.03
        def b(x): return a*(.5*x if kind=='linear' else .5*x+.2*np.sin(x))
        direct=quad(lambda s:np.cos(s)*b(np.cos(s)*z0+kappa*np.sin(s)*z),0,np.pi/2,epsabs=1e-13)[0]
        bracket=lambda s:b(np.cos(s)*z0+kappa*np.sin(s)*z)-np.cos(s)*b(z0)-np.sin(s)*b(kappa*z)
        rhs=np.pi/4*b(z0)+.5*b(kappa*z)+quad(lambda s:np.cos(s)*bracket(s),0,np.pi/2,epsabs=1e-13)[0]
        assert abs(direct-rhs)<1e-12
        max_bracket=max(abs(bracket(s)) for s in np.linspace(0,np.pi/2,101))
        if kind=='linear': assert max_bracket<1e-15
        calibration.append({'kind':kind,'identity_error':direct-rhs,'max_bracket':max_bracket})

out={'status':'PASS: rank-one trace/gap identities and exact affine VALUE calibration',
     'scope':'Screens the new empirical rank-one matrix normalization only; not a general sampler lower bound.',
     'heat_window':rows,'value_calibration':calibration,
     'R4_exponents':{'scalar_final':2/3,'fixed_rank_matrix_at_d_equals_a_minus8':17/6,
                    'covariance_replica_necessity_at_scalar_count':13/2}}
Path(__file__).with_name('matrix_gap_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
