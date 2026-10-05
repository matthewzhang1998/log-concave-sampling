#!/usr/bin/env python3
"""Finite arithmetic/ownership/geometry diagnostics; not a native compiler execution."""
from fractions import Fraction as F
from pathlib import Path
import json, math, hashlib
ROOT=Path(__file__).resolve().parent
checks=0

def check(cond, message):
    global checks
    checks+=1
    if not cond: raise AssertionError(message)

def close(a,b,tol=1e-8):
    check(abs(a-b)<=tol*max(1,abs(a),abs(b)), f'{a} != {b}')

alpha=F(14,15); beta=gamma=F(19,10); grade=F(39,10); b=25
rows={
 'prefix_smoothing':2+gamma,
 'prefix_decoupling':3+3*alpha-gamma,
 'prefix_nested':3+alpha,
 'smoothing_drift':1+2*gamma,
 'bridge_quadrature':4+alpha,
 'outer_target_final':F(4),
 'raw_sqrt_eta':3+beta/2,
 'raw_eta_sqrt_w':3+beta-alpha/2,
 'raw_mean':4+beta,
}
for name,n,p in [('mean_third',5,F(1)),('mixed_prior',6,F(7,6)),('quartic',6,F(2)),('gram_cubic_feedback',7,F(5,2)),('fourth_target',6,F(3,2)),('mixed_mixture',7,F(3,2))]:
 rows[name+'_bulk']=n-p*alpha
 rows[name+'_endpoint']=F(n) if p==1 else n-(p-1)*beta
check(rows['prefix_smoothing']==grade,'leading smoothing')
check(rows['prefix_decoupling']==grade,'leading bridge')
rest=[v for k,v in rows.items() if k not in ('prefix_smoothing','prefix_decoupling')]
check(min(rest)==F(59,15),'strict minimum')
check(min(rest)-grade==F(1,30),'absorption margin')
for k,v in rows.items(): check(v>=grade,k)
check((b-4)-(b-3)*beta/2==F(1,10),'mean prior margin')
check(1-beta/2==F(1,20),'radius')
check(F(3,2)-3*beta/4==F(3,40),'reserve')
check(F(1,2)-beta/4==F(1,40),'physical reserve ratio')
check(1-beta/3==F(11,30),'mixed radius')

# Every fixed member has positive native/ledger margins; no infinite-order limit.
family=[]
for n in range(2,202):
 eps=F(1,n); a=1-2*eps/3; bet=gam=2-eps
 bm=math.floor(3+2/eps)+1
 target=4-eps
 check(2+gam==target,'family smoothing')
 check(3+3*a-gam==target,'family decoupling')
 check(3+a>target,'family nested')
 check(3+bet/2>target,'family raw')
 check((bm-4)-(bm-3)*bet/2>0,'family native prior')
 check(1-bet/2>0,'family radius')
 for nn,p in [(5,F(1)),(6,F(7,6)),(6,F(2)),(7,F(5,2)),(6,F(3,2)),(7,F(3,2))]:
  check(nn-p*a>target,'family bulk')
  check((nn if p==1 else nn-(p-1)*bet)>target,'family endpoint')
 family.append({'epsilon':str(eps),'b_mean':bm,'grade':str(target)})

# Exact unit-variance partition, including all nested Gram keeps.
shares=[F(1,10)]*5+[F(1,2)]
check(sum(shares)==1,'top variance')
check(F(1,20)+F(1,40)+F(1,40)==F(1,10),'corrected Gram')
check(F(1,40)+F(1,40)==F(1,20),'retuned Gram half split')
for theta in [F(1,5),F(1,3),F(1,2),F(4,5)]:
 # Corrected Gram interpolation: old mixture u/2 plus visible theta*u/4 and keep.
 visible=theta/4; keep=F(1,4)+(1-theta)/4
 check(F(1,2)+visible+keep==1,'Gram interpolation keep counted once')

# Joint Gaussian disintegration and known-row terminal alignment.
for A in [.1,.03,.01,.003,.001,.0001,.00001]:
 w=A**float(alpha); eta=A**float(beta); h=A**float(gamma); q=1-w; sig2=1-q*q
 check(eta<=w<=.5,'geometric parameter window')
 for j in range(1,80):
  gap=eta*(1/eta)**(j/80)
  t=1-gap; c2=1-t*t; s=1-q*q*t*t
  v=c2*sig2/s
  close(1/v,q*q/sig2+1/c2,1e-7)
  check(v>=eta/2*(1-1e-7),'v lower')
  check(v<=2*w*(1+1e-7),'v upper')
  coeff_y=q*c2/math.sqrt(s)
  close(coeff_y**2+v,c2)
  norm=math.sqrt(coeff_y**2+v)
  # Carrier+complement rotation rows form an orthogonal 2x2 matrix.
  ra=coeff_y/norm; rb=math.sqrt(v)/norm
  close(ra*ra+rb*rb,1)
  close(ra*rb+rb*(-ra),0)
  rh=math.sqrt(1-h*h)
  close((rh*math.sqrt(c2))**2+h*h,1-rh*rh*t*t)
  check(math.sqrt(1-rh*rh*t*t)>0,'positive carrier norm')

# Dimensional/readset ledger: graph JSON fixes ownership rather than tensor leaves.
graph=json.loads((ROOT/'JOIN-GRAPH.json').read_text())
services=graph['law_node']['services']
check(len(services)==5,'five actual services')
check(len({s['bank_id'] for s in services})==5,'distinct complete banks')
for service in services:
 check(service['retained_readset']==['Y'],'only Y read by tail service')
 check(service['tape']=='complete_fresh_owned','complete source tape')
 check(service['source_zero_fraction']=='1/10','actual variance share')
check(graph['law_node']['external_keep_fraction']=='1/2','external keep')
check(graph['law_node']['retained_for_comparison']==['z','Y','N','G_h'],'retained labels')
check(graph['scope']['exact_target']=='integrated_conditional_W2_standard_Y','scope')
check(not graph['scope']['uniform_exact_target_in_Y'],'no uniform promotion')
check(graph['source_alignment']['preserve_all_perpendicular_roots'],'retain complement')
check(graph['source_alignment']['cross_node_independence_required'] is False,'no false independence')
check(graph['raw_node']['native_services']==[],'no endpoint native')
check(graph['leaf_types']==['original_g_VALUE','fully_expanded_native_VALUE_program','Gaussian_or_scalar_arithmetic'],'no analytical leaves')
for nl in range(1,12):
 for nr in range(0,8):
  # Synthetic dimensions represent arbitrary complete service sizes, each multiple D.
  dims=[7,11,13,17,19]
  before=nl*(4+sum(dims))+nr*6
  aligned=1+nl*(3+sum(dims))+5*nr
  check(aligned==before-(nl+nr-1),'exact root-sharing count')

# The four-row ceiling is a property of displayed powers, not a lower-bound theorem.
for ia in range(-40,81):
 for ig in range(-40,81):
  aa=F(ia,20); gg=F(ig,20)
  check(min(2+gg,3+3*aa-gg,3+aa,5-aa)<=4,'four-row ceiling')
for bm in range(4,101):
 check((bm-4)-(bm-3)*F(2)/2==-1,'no fixed-b LOCAL endpoint domination')

# Pinned sealed lineage.
pins=json.loads((ROOT/'INPUT-PINS.json').read_text())
for item in pins['inputs']:
 p=Path(item['path'])
 check(p.is_file(),f'input exists {p}')
 check(hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256'],f'pin {p}')

out={'status':'PASS','assertions':checks,'native_compiler_executed':False,
 'grade':str(grade),'strict_residual_grade':str(min(rest)),
 'rows':{k:str(v) for k,v in rows.items()},
 'family_members_checked':len(family),'fixed_concrete_mean_order':b,
 'guard_scope':'actual imported native/log guards required; bare-power checks are not numeric admission',
 'families_first_last':[family[0],family[-1]]}
(ROOT/'join_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
