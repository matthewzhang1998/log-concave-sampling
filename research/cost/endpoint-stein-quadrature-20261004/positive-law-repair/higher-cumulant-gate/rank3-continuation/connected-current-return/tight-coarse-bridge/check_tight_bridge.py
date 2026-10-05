"""Finite diagnostic checks; not a proof of the open grouped-current port."""
import itertools, json, math
import numpy as np
from numpy.polynomial.hermite import hermgauss

checks=[]
def check(name, value):
    assert bool(value), name
    checks.append(name)

# Noncommuting matrix channel: no covariance oracle is executed by this check.
rng=np.random.default_rng(17)
for D in (1,2,5,9):
    Mp=rng.normal(size=(D,D)); Mm=rng.normal(size=(D,D))
    Mp*=.23/np.linalg.norm(Mp,2); Mm*=.19/np.linalg.norm(Mm,2)
    Diff=(Mp-Mm)/2; Sum=(Mp+Mm)/2
    cov=Diff@Diff.T+Sum@Sum.T+(2*np.eye(D)-Mp@Mp.T-Mm@Mm.T)/2
    check(f'difference_covariance_{D}',np.max(np.abs(cov-np.eye(D)))<1e-12)
    check(f'cross_from_original_input_{D}',np.max(np.abs(Diff-(Mp-Mm)/2))<1e-12)

# Exact bounded witness and conditional covariance at a fixed old input caller.
x,w=hermgauss(120); x=x*math.sqrt(2); w=w/math.sqrt(math.pi)
for delta in (.5,.2,.05,.01):
    for eps in (.1,.25):
        m=eps*np.sin(delta*x)
        v=float(w@(m*m)); exact=eps**2*(-math.expm1(-2*delta*delta))/2
        check(f'sine_variance_{delta}_{eps}',abs(v-exact)<1e-13)
        check(f'sine_odd_mean_{delta}_{eps}',abs(float(w@m))<1e-13)
        for r in (0.,.5,1.5,2.):
            direct=float(w@((r*m)**2+1-m*m))
            check(f'reroot_covariance_{delta}_{eps}_{r}',abs(direct-(1+(r*r-1)*exact))<1e-12)
            check(f'reroot_not_standard_{delta}_{eps}_{r}',abs(direct-1)>1e-8)

# All choices of the five possible inverse-delta nonroot rescalings.
# Three nonroots on a four-vertex spine; two side leaves are omitted.
p=.5
for flags in itertools.product((0,1),repeat=5):
    root=1+2*p+p*sum(flags)
    nonroots=[1-p*f for f in flags]
    main=root+sum(nonroots)
    spine=root+sum(nonroots[:3])
    check(f'main_preserved_{flags}',abs(main-7)<1e-12)
    check(f'spine_only_offspine_gain_{flags}',abs(spine-(5+p*sum(flags[3:])))<1e-12)
    check(f'positive_sources_{flags}',min([root]+nonroots)>0)

# Leading zero-shift fourth current of four independent Gaussian probes.
# E[P Q S T]=0; E[(P Q S T)^2]=1 exactly by factorization.
e0=float(w@x); e2=float(w@(x*x))
check('four_probe_mean_zero',abs(e0**4)<1e-14)
check('four_probe_variance_one',abs(e2**4-1)<1e-12)
check('two_spines_cycle_rank_three',10-8+1==3)

# Endpoint ledger: beta~1 for Delta=sigma^2; two derivative hits cost sigma^-2.
for k in (4,8,12,16):
    sigma=2.**(-k/2); Delta=sigma*sigma
    beta=Delta**2/sigma**4
    check(f'tight_beta_{k}',abs(beta-1)<1e-12)
    check(f'own_covariance_loss_{k}',abs((Delta**4/sigma**10)/(sigma**-2)-1)<1e-12)

out={"status":"PASS", "checks":len(checks), "scope":"finite matrix, scalar Gaussian, amplitude and endpoint-ledger diagnostics only", "open":"positive whole-coarse compiler, higher joint physical tilts and full source-qualified closure"}
with open(__file__.replace('check_tight_bridge.py','tight_bridge_checks.json'),'w') as f: json.dump(out,f,indent=2); f.write('\n')
print(json.dumps(out,indent=2))
