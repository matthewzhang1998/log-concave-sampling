import json, math
import numpy as np
rng=np.random.default_rng(905)
min_ratio=math.inf
max_error=0.
for _ in range(20000):
    q,s,r1,r2,r3=1-np.exp(rng.uniform(-20,0,5))
    sig=np.sqrt((1-np.array([r1,r2,r3])**2)/2)
    vy=1-q*q; vx=1-s*s*q*q
    G=np.array([[r1*r1*vy+sig[0]**2,r1*r2*s*vy,r1*r3*s*vy],
      [r1*r2*s*vy,r2*r2*vx+sig[1]**2,r2*r3*vx],
      [r1*r3*s*vy,r2*r3*vx,r3*r3*vx+sig[2]**2]])
    m=min(sig[0]**2,sig[2]**2,1-s*s)
    min_ratio=min(min_ratio,float(np.linalg.eigvalsh(G)[0]/m))
    r=rng.uniform(); h2=1-r*r
    t2=sig**2+h2*m/32
    R=h2*(G-m/32*np.eye(3))
    max_error=max(max_error,float(np.max(abs(R+np.diag(t2)-(h2*G+np.diag(sig**2))))))
# In normalized Hermite basis, L=sum c_k h_k.
# Cov P_a L=sum c_k^2 a^(2k); derivative semigroup density is
# 2r sum k*c_k^2*r^(2k-2).
cs=rng.normal(size=50); a=.93
cov=float(np.sum(cs**2)); cov_smoothed=float(np.sum(cs**2*a**(2*np.arange(1,51))))
tail=float(np.sum(cs**2*(1-a**(2*np.arange(1,51)))))
# Endpoint sums with literal proxy mass w=Delta and sigma^2=Delta.
ratios=[]
for h in [2.**(-j) for j in range(1,15)]:
    d=2.**(-np.arange(0,80))
    ratios.append(dict(h=h, first=float(h*np.sum(d*d/(d+h*h)**2.5)), covariance=float(h*h*np.sum(d**4/(d+h*h)**5))))
result=dict(covariance_redistribution_max_error=max_error, minimum_observed_spectral_ratio=min_ratio, proved_spectral_ratio=1/16, ou_identity_error=abs(cov-cov_smoothed-tail), dyadic_scaled_ratios=ratios, status='Numerical sanity checks only; the accompanying note supplies proofs and lists readout/producer gates.')
with open('/workspace/shared/old-bank-sixth-return-20261005/exact_ou_grouping_checks.json','w') as f: json.dump(result,f,indent=2)
print(json.dumps({k:v for k,v in result.items() if k!='dyadic_scaled_ratios'},indent=2))
