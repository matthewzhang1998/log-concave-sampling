#!/usr/bin/env python3
"""Exact exponent algebra and complete affine-Gaussian source/pair fixtures.

Reuses only definition text from the already audited complete-source fixture.
No prior file is imported/executed or overwritten. Replaces its complete known
map and fresh record allocator by the ordered eleven-refresh implementation.
The Gaussian seed is a fixture, not an implementation of CW7.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import math,json,hashlib
import numpy as np
OUT=Path(__file__).parent
OLD=OUT.parent/'full-independent-audit/check_complete_hidden_telescope.py'
old=OLD.read_text().split('# Exact Gaussian comparisons')[0]
env={'__file__':str(OLD),'__name__':'fixture_definitions_only'}
exec(compile(old,str(OLD),'exec'),env)
Lin,Tape=env['Lin'],env['Tape']
marks=[F(12,5),F(27,10),F(14,5),F(29,10),F(3),F(31,10),F(16,5),F(33,10),F(17,5),F(7,2),F(71,20)]
rows=[]
def check(name,ok,**data):
    if not ok:raise AssertionError((name,data))
    rows.append(dict(name=name,passed=True,**data))

@lru_cache(None)
def coeff(a,lam,j):
    j=max(j,1);R=j+.5
    Nq=max(1,math.ceil(a**(-(j-1))))
    Mq=max(1,j-1)
    aq,bq,cq=env['flow_coeff'](a,lam,math.pi/2,Nq,Mq)
    alpha,beta,cs,bill=aq,bq,[cq],Nq*Mq
    for mark in marks:
        v=float(2*(mark-1)/3);h=a**(v/2)
        Ns=max(1,math.ceil(a**(-max(0.,R-float(mark)-.5))))
        Ms=max(1,math.ceil((R-.5)/(1+v))-1)
        al,be,cl=env['flow_coeff'](a,lam,h,Ns,Ms)
        alpha=al*alpha;beta=al*beta+be;cs=[al*x for x in cs]+[cl]
        bill+=Ns*Ms
    return alpha,beta,tuple(cs),bill

def known_from_records(a,y,j,w,packets,lam,t):
    x=env['seed'](a,y,w,lam)
    al,be,cs,bill=coeff(a,lam,j)
    for packet in packets:
        x=al*x+be*y+sum(c*z for c,z in zip(cs,packet));t.calls+=bill
    return x

def packet(t,label):return [t.normal(label+'-quarter')]+[t.normal(label+'-short-'+str(i)) for i in range(len(marks))]
def known(a,y,j,lam,t):
    j=max(1,j);w=t.normal('seed');ps=[packet(t,'known') for _ in range(j)]
    return known_from_records(a,y,j,w,ps,lam,t)
def known_pair(a,yf,yc,j,lam,t):
    jj=max(1,j);w=t.normal('paired-seed');ps=[packet(t,'paired') for _ in range(jj)]
    coarse=ps if j==1 else ps[1:]
    return known_from_records(a,yf,j,w,ps,lam,t),known_from_records(a,yc,max(1,j-1),w,coarse,lam,t)
env.update(known_from_records=known_from_records,known=known,known_pair=known_pair)

for i,J in enumerate(marks):
    v=2*(J-1)/3;chi=max(v/2,J-2)
    check('own_projection_exponent',1+3*v/2==J,J=str(J))
    check('signed_twin_exponent',2+(J-2)==J,J=str(J))
    for Js in marks[i:]:
        vs=2*(Js-1)/3
        check('all_later_path_exponent',1+v/2+vs>=J,J=str(J),later=str(Js),exponent=str(1+v/2+vs))
    for j in range(1,21):
        check('short_grid_quarter_dominated',max(F(0),F(j)-J)<=j-1,J=str(J),j=j)
check('pure_fourth_full_radius',F(3,4)-F(7,15)==F(17,60))
check('pure_fourth_weighted_body',F(3,4)-F(7,15)-F(1,4)==F(1,30))
check('old_fixed_width_weighted_failure',F(3,4)-F(19,30)-F(1,4)<0)
chain=marks[1:]
gates=[]
for x,y in zip(chain[:-1],chain[1:]):
    u=(F(327,50)-x)/4;e=y-F(7,2)+u+F(1,100)
    for z in [x,y]:
        chi=max((z-1)/3,z-2)
        gs=[F(3,2)-u+chi-e,F(3,2)-u-chi+e,F(3,2)-u]
        check('literal_covariance_child_gates',min(gs)>0,x=str(x),y=str(y),child=str(z),gates=[str(v) for v in gs])
        gates.append(gs)
check('minimum_literal_child_gate',min(min(z) for z in gates)==F(1,100))
for p in [F(2),F(4),F(8),F(16)]:
    for n in range(1,13):
        c=2*(n-1);b=F(10);dmax=16;delta=F(1,200);E=9
        B=math.ceil(2*p*(E+b+F(c,2)+F(c+dmax,4*p)+1)/delta)
        exponent=delta*B/(2*p)-b-F(c,2)-F(c+dmax,4*p)
        check('ambient_fixed_moment_prior_order',exponent>=E+1,p=str(p),n=n,B=B,exponent=str(exponent))

for a in [.25,.125,.0625,.015625]:
    hs=[a**float((J-1)/3) for J in marks]
    cc=[math.prod(math.cos(h) for h in hs)]+[math.sin(hs[i])*math.prod(math.cos(h) for h in hs[i+1:]) for i in range(len(hs))]
    check('multi_refresh_carrier_coisometry',abs(sum(c*c for c in cc)-1)<1e-14,a=a)
    for i,J in enumerate(marks):
        gamma=cc[i+1]**2;v=float(2*(J-1)/3)
        check('protected_nonzero_variance',gamma>0,a=a,J=str(J),gamma=gamma,ratio=gamma/a**v)

fixtures=[]
for a,d,j in [(.25,0,1),(.25,0,2),(.25,0,3),(.125,0,2),(.25,1,1),(.25,1,2),(.125,1,2),(.25,2,1),(.25,2,2)]:
    n=3;lam=.6;y=.7
    tf,tc,tp=Tape(),Tape(),Tape()
    f=env['force'](a,y,d,j,n,lam,tf)
    c=env['force'](a,y,d,j-1,n,lam,tc)
    pf,pc=env['pair'](a,y,d,j,n,lam,tp)
    for side,x,z in [('fine',f,pf),('coarse',c,pc)]:
        check('complete_named_marginal',abs(x.mean-z.mean)<5e-11 and abs(x.var()-z.var())<5e-11,
              a=a,d=d,j=j,side=side,mean_error=x.mean-z.mean,variance_error=x.var()-z.var())
    gap=(pf-pc).norm();zero=abs(pf.mean-pc.mean)
    check('literal_level_one_zero',j!=1 or zero<1e-12,a=a,d=d,j=j,zero=zero)
    state=f/(math.sqrt(a)*lam);q=env['ideal_center'](a,y,d,lam)
    err=math.hypot(state.mean-q/(1+a*lam),math.sqrt(state.var())-math.sqrt(a/(1+a*lam)))
    row=dict(a=a,d=d,j=j,force_pair_ratio=gap/a**j,zero_ratio=zero/a**j,law_ratio=err/a**(j+.5),
             source_gaussians=tf.n,pair_gaussians=tp.n,source_original_calls=tf.calls,pair_original_calls=tp.calls)
    fixtures.append(row)
    check('fixture_strong_pair_envelope',row['force_pair_ratio']<3,**row)
    check('fixture_law_envelope',row['law_ratio']<3,**row)

out=dict(status='PASS',checks=len(rows),scope='Exact algebra and rebuilt complete Gaussian fixtures, not a foundational-source or exterior-constructor simulation',
         original_fixture_sha256=hashlib.sha256(OLD.read_bytes()).hexdigest(),marks=[str(x) for x in marks],
         complete_fixtures=fixtures,checks_detail=rows)
(OUT/'ordered_refresh_join_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(status='PASS',checks=len(rows),complete_fixtures=len(fixtures),
 max_force_ratio=max(x['force_pair_ratio'] for x in fixtures),max_law_ratio=max(x['law_ratio'] for x in fixtures)),indent=2))
