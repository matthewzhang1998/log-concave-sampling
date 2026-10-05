#!/usr/bin/env python3
"""Checks scalar geometry and bounds only. No OU path/grid/history is executed."""
from pathlib import Path
import json, hashlib, math
import mpmath as mp
mp.mp.dps=70
count=0

def check(cond, label):
    global count
    assert cond, label
    count += 1

def close(x,y,label):
    check(abs(x-y) <= mp.mpf('1e-55')*max(1,abs(x),abs(y)),label)

A=mp.mpf('0.01'); l=A/4; u=3*A/4
mmin=mp.mpf(3)/4/mp.cosh(2*(1+u))**2
check(mmin>mp.mpf(1)/20,'conditional sech lower bound')
close((l/20-u*u)/A,mp.mpf(11)/1600,'uniform negative bracket')
close((A/4)*(l/20-u*u)/A**2,mp.mpf(11)/6400,'mixed derivative constant')

geometries=[]
for d in [mp.mpf('0.01'),mp.mpf('0.1'),mp.mpf('0.25'),mp.mpf('0.75'),mp.mpf('2')]:
  for a in [mp.mpf('0'),mp.mpf('0.25'),mp.mpf('1.5')]:
    for gap in [mp.mpf('0'),mp.mpf('0.5')]:
      b=a+d+gap
      def hat(t,base):
        if t<base or t>base+d: return mp.mpf(0)
        return mp.sinh(min(t-base,base+d-t))/mp.sinh(d/2)
      hi=mp.quad(lambda t:hat(t,a),[a,a+d/2,a+d])
      he=mp.quad(lambda t:mp.exp(-t)*hat(t,a),[a,a+d/2,a+d])
      ke=mp.quad(lambda t:mp.exp(-t)*hat(t,b),[b,b+d/2,b+d])
      close(hi,2*mp.tanh(d/4),'hat mass')
      close(he,(d/2)*mp.exp(-a-d/2),'weighted earlier hat')
      close(ke,(d/2)*mp.exp(-b-d/2),'weighted later hat')
      J=hi*ke; K=(hi-he)*ke
      check(J>0,'strict positive overlap integral')
      check(0<=K<J,'K bounded by J')
      Kdirect=ke*mp.quad(lambda t:(1-mp.exp(-t))*hat(t,a),[a,a+d/2,a+d])
      close(K,Kdirect,'ordered Tonelli contraction')
      for av in [mp.mpf('0.0001'),mp.mpf('0.001'),mp.mpf('0.01')]:
        lv=av/4; uv=3*av/4
        bound=(av/4)*(-lv*J/20+uv*uv*K)
        certified=-mp.mpf(11)/6400*av*av*J
        check(bound<=certified,'negative actual-history expectation envelope')
        # Exact zero-path variation for arbitrary g'(0)=av/2, g''(0)=av/4.
        direct=(av/4)*(-(av/2)*J+(av/2)**2*K)
        check(direct<0,'exact zero-path mixed differential')
      geometries.append({'a':str(a),'b':str(b),'d':str(d),'J':str(J),'K_over_J':str(K/J)})

max_beta_sq=mp.mpf(0); max_tail_ratio=mp.mpf(0)
for j in range(1,101):
  av=mp.mpf(j)/200
  coeff=[mp.mpf(3),3*av,av*av/2]
  # Integral beta/A squared is a finite scalar gamma moment sum.
  beta_sq=sum(coeff[r]*coeff[s]*mp.factorial(r+s)/2**(r+s+1) for r in range(3) for s in range(3))
  max_beta_sq=max(max_beta_sq,beta_sq)
  check(beta_sq<8,'uniform one-HS sensitivity energy')
  first=sum(coeff[r]*mp.factorial(r) for r in range(3))
  check(first<=mp.mpf(19)/4,'complete first kernel')
  for T in [mp.mpf('0'),mp.mpf('0.25'),mp.mpf('0.5'),mp.mpf('1'),mp.mpf('2'),mp.mpf('10'),mp.mpf('100')]:
    tailpoly=3+3*av*(T+1)+(av*av/2)*(T*T+2*T+2)
    ratio=tailpoly**2/(1+T**4)
    max_tail_ratio=max(max_tail_ratio,ratio)
    check(ratio<80,'analytic coherent-future tail envelope')

sources=[
 ('dyadic covariance','/workspace/shared/dyadic-bridge-square-20261005/MANIFEST.json','2365fc791230cd2d73be22964abf45c071825771ca9d6e860b483a7843a67a0a'),
 ('prefix LAW join','/workspace/shared/dyadic-prefix-law-join-20261005/MANIFEST.json','130131317183c6f00ea3566c187296912db7850081904322215a3ec9d5eb37f2'),
 ('next m4 port','/workspace/shared/next-history-m4-port-20261005/MANIFEST.json','ac9845847d49f85da1ce7cb9728e18c2f0a6ed8f7da4de1befe685878ca7250f'),
 ('local shifted source','/workspace/shared/marked-shifted-local-current-20261005/MANIFEST.json','787ddcc7f03de97085e2fcb93d20b631c705950c3d393b43eadd1c6551eabb34')]
pins=[]
for name,path,sha in sources:
 actual=hashlib.sha256(Path(path).read_bytes()).hexdigest()
 check(actual==sha,'input manifest '+name)
 pins.append({'name':name,'path':path,'sha256':actual})
out=Path(__file__).parent
(out/'INPUT-PINS.json').write_text(json.dumps(pins,indent=2)+'\n')
result={'status':'PASS','assertions':count,'precision_decimal_digits':mp.mp.dps,'conditional_sech_minimum_bound':str(mmin),'max_sensitivity_squared_ratio':str(max_beta_sq),'max_checked_tail_ratio':str(max_tail_ratio),'geometries':geometries,'scope':'Scalar geometry and deterministic analytic envelopes only; no path-grid, force-history program, native action, or block-LAW compiler was executed.'}
(out/'scalar_obstruction_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='geometries'},indent=2))
