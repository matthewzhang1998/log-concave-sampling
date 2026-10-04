#!/usr/bin/env python3
"""Deterministic constant/census checks, not a replacement for the proof."""
import hashlib,json,math
from pathlib import Path

ROOT=Path(__file__).resolve().parent
checks=0
def check(ok,msg):
    global checks
    if not ok: raise AssertionError(msg)
    checks+=1

fixtures=[]
for rho in [.001,.01,.1,1.,10.,100.,1e4,1e6]:
    for log_inverse_delta in [2.01,3.,5.,10.,50.,100.,1000.]:
        L=math.log(math.e+rho)+log_inverse_delta
        R=math.sqrt(64*L)
        eps=1/(16*R)
        K=math.ceil(8*L);B=math.ceil(8*L)
        M=math.ceil(2**20*(1+rho*rho)*R*R*L)
        check(eps*R+eps*eps/2 < math.log(5/4),'uniform positivity guard')
        pilot_exponent=M/(32*rho*rho*256*R*R)
        check(pilot_exponent>=8*L,'Hoeffding exponent')
        check(K*math.log(16)>=8*L,'rank-uniform remainder')
        check(B*math.log(3)/2>=4*L,'finite fallback W2 exponent')
        check(math.log(1+rho*rho)-8*L <= -2*log_inverse_delta,'bad-event moment exponent')
        check(math.log(4*rho)-8*L <= -log_inverse_delta,'good clipping bias exponent')
        count=(K+1)*M
        bound=10*2**20*64*(1+rho*rho)*L**3
        check(count<=bound,'literal provider census polynomial bound')
        fixtures.append({'rho':rho,'log_inverse_delta':log_inverse_delta,
                         'K':K,'B':B,'M':M,'complete_provider_count':count})

# Finite Bernoulli bounded-source clipping example, independent of asymptotic M.
def binomial_distribution(M,p):
    return [(2*k/M-1,math.comb(M,k)*p**k*(1-p)**(M-k)) for k in range(M+1)]
for M in [2,4,8,16]:
    for p in [.1,.3,.5,.8]:
        mu=2*p-1;dist=binomial_distribution(M,p)
        check(abs(sum(pr for _,pr in dist)-1)<1e-13,'batch probability')
        check(abs(sum(x*pr for x,pr in dist)-mu)<1e-13,'batch mean')
        for b,_ in dist:
            a=.5
            clipped=sum(min(b+a,max(b-a,x))*pr for x,pr in dist)
            check(-1-1e-14<=clipped<=1+1e-14,'clipped mean original interval')
            if abs(b-mu)<=a/4:
                tail=sum(pr for x,pr in dist if abs(x-mu)>3*a/4)
                check(abs(clipped-mu)<=2*tail+1e-13,'exact good-pilot clipping bound')

note=ROOT/'POLYLOG-SCALAR-MEAN-LAW-REFINEMENT.md'
out={'status':'PASS','assertions':checks,
     'source_sha256':hashlib.sha256(note.read_bytes()).hexdigest(),
     'parameter_fixtures':fixtures,
     'scope':'Checks constants, sample census, and finite clipping identities only.'}
(ROOT/'polylog_scalar_refinement_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'assertions':checks,'fixtures':len(fixtures)},indent=2))
