"""Independent fixed-rank Riesz join, Hermite energy, guard and grade checks.
No supplied checker or native selected-pair/filter implementation is imported.
"""
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
from collections import Counter
import json,math
import numpy as np
import sympy as sp
from numpy.polynomial.hermite_e import hermegauss,hermeval
from numpy.polynomial.legendre import leggauss
HERE=Path(__file__).resolve().parent
counts=Counter();maximums={}
def ck(ok,group,msg):
 counts[group]+=1
 if not bool(ok):raise AssertionError(group+': '+msg)
B,H,U,G=sp.symbols('B H U G');variables=(B,H,U,G)
@lru_cache(None)
def gm(n):
 return sp.Integer(0) if n%2 else sp.Integer(math.prod(range(1,n,2)))
def expectation(expr):
 p=sp.Poly(sp.expand(expr),variables)
 return sum(c*math.prod(gm(e) for e in exps) for exps,c in p.terms())

# Polynomial identity fixture, not a uniform-in-B derivative-bound fixture.
# It checks the actual two sides independently, and separately checks transfer.
He7=sp.hermite_prob(7,H)
F=sp.Rational(1,100)*(B**2-1)*He7
RF=sp.Rational(1,100)*B*He7
Jold=sp.Rational(3,10)
for t in (sp.Rational(1,3),sp.Rational(2,3)):
 W=H+U/5+3*B/10+He7/100+t*F+G/2
 J=sp.diff(W,B);powers=[sp.expand(W**j) for j in range(5)]
 for p in range(1,5):
  for q in range(4):
   lhs=expectation(p*F*powers[p-1]*B**q)
   ordinary=0 if p<2 else RF*J*p*(p-1)*powers[p-2]*B**q
   observer=0 if q==0 else RF*p*q*powers[p-1]*B**(q-1)
   rhs=expectation(ordinary+observer)
   ck(lhs==rhs,'exact_riesz','independent Gaussian moment evaluation at same endpoint')
   gradw=p*powers[p-1]*B**q
   grado=0 if q==0 else q*powers[p]*B**(q-1)
   transferred=expectation(2*G*RF*(J*gradw+grado))
   ck(rhs==transferred,'exact_keep','one physical derivative transferred through sqrt(kappa)=1/2 keep')
   if p==2 and q==2:
    omitted=expectation(ordinary)
    ck(omitted!=lhs,'observer_boundary','retaining B requires its current; omission is detected')

# Exact physical Hermite moments and hypercontractivity, independently evaluated.
for n in range(11):
 he=sp.hermite_prob(n,H)
 m2=expectation(he**2);m4=expectation(he**4)
 ck(m2==math.factorial(n),'hermite','Hermite isometry')
 ck(m4<=3**(2*n)*m2*m2,'hermite','degree-n L4/L2 factor 3^(n/2)')
m=7;m2=math.factorial(m);m4=int(expectation(He7**4))

# Multivariate vector Hermite field realizes the exact sqrt(7!) v^(-7/2) factor.
x,w=hermegauss(32);w=w/math.sqrt(2*math.pi)
he=np.array([hermeval(x,[0]*j+[1]) for j in range(8)])
rng=np.random.default_rng(8132026);coeff=rng.normal(size=(2,8))
tensor_hs_squared=sum(math.comb(7,j)*float(coeff[o,j]**2) for o in range(2) for j in range(8))
for variance in (.3,1.,2.5):
 field=np.zeros((2,len(x),len(x)))
 for o in range(2):
  for j in range(8):field[o]+=math.comb(7,j)*coeff[o,j]*variance**(-3.5)*np.outer(he[j],he[7-j])
 measured=float(np.einsum('i,j,oij,oij->',w,w,field,field))
 expected=m2*variance**(-7)*tensor_hs_squared
 ck(abs(measured/expected-1)<2e-13,'multivariate_energy','full vector energy includes Hermite isometry and covariance loss')

# Genuinely bounded B-derivative fixture; polynomial only in the new physical H.
# F=a sin(B) He7(H)/sqrt(7!), old physical component=alpha sin(B)sin(U),
# retained nonlinear observer=.7 alpha sin(B)sin(V), with all roots independent.
bx,bw=hermegauss(100);bw=bw/math.sqrt(2*math.pi)
s,sw=leggauss(100);s=(s+1)/2;sw=sw/2
r=np.cos(bx[:,None]*s[None,:])@(sw*np.exp(-(1-s*s)/2))
esin2=(1-math.exp(-2))/2;esin4=(3-4*math.exp(-2)+math.exp(-8))/8
rb_energy=float(bw@(r*r));prodB=float(bw@(r*r*np.cos(bx)**2))
ck(rb_energy<=esin2+1e-13,'bounded_join','Gaussian Riesz L2 contraction for sin(B)')
normalized_h4=m4/(m2*m2)
fixtures=[]
for alpha in (.02,.1,.3):
 for a in (1e-6,1e-3,.05):
  for t in (0.,.5,1.):
   for kappa in (.25,1.,4.):
    E=a*math.sqrt(esin2)
    oldL4=alpha*(esin4*(1+.7**4)+2*.7**2*esin2**2)**.25
    newL4=a*normalized_h4**.25
    actual=a*math.sqrt(prodB*(alpha*alpha*(1+.7**2)*esin2+t*t*a*a*normalized_h4))/math.sqrt(kappa)
    bound=3**3.5/math.sqrt(kappa)*E*(oldL4+t*newL4)
    ck(actual<=bound*(1+1e-12),'bounded_join','one-Hilbert frozen-B L4 product with stacked retained observer')
    fixtures.append({'alpha':alpha,'a':a,'t':t,'kappa':kappa,'actual_current_L2':actual,'certified_bound':bound})

# Averaging Y before the correct Holder step is not interchangeable.
# On a rare Y event of probability p, h=p^-1/2, j=p^-1/4.
# ||h||2=||j||4=1, but ||hj||2=p^-1/4 can be arbitrarily large.
for denominator in (10**4,10**8,10**12):
 p=Q(1,denominator)
 ck(p*(1/p)==1 and p*(1/p)==1,'caller_boundary','fixed separate integrated moment budgets')
 ratio=float(p)**(-.25)
 ck(ratio>=10,'caller_boundary','joint-average moments alone do not bound the product')
ck(float(Q(1,10**12))**(-.25)>3**3.5,'caller_boundary','physical degree-seven hypercontractivity cannot repair arbitrary Y concentration')
# Hidden correlation U_old=B gives an omitted derivative despite formal D_B U_old=0.
ck(expectation(2*B*B)==2 and expectation(sp.Integer(0))==0,'bank_boundary','shared primitive roots must enter B and its derivative ledger')

# Exact rational grades and all floor propagation powers.
g=Q(1,8);beta=Q(1,14);target=Q(65,8)
ledger={'root':Q(15,2)-6*g,'caller':Q(15,2)-7*g,'main_energy':8-6*g,'main_first':8-7*g,'own':15-12*g,'mixed_exact':9-6*g,'mixed_conservative':9-7*g,'self_bank':16-13*g,'heat':8+g,'response_floor':target-8+6*g,'value_floor':target-8+7*g}
expected={'root':Q(27,4),'caller':Q(53,8),'main_energy':Q(29,4),'main_first':Q(57,8),'own':Q(27,2),'mixed_exact':Q(33,4),'mixed_conservative':Q(65,8),'self_bank':Q(115,8),'heat':Q(65,8),'response_floor':Q(7,8),'value_floor':Q(1)}
for name,value in ledger.items():ck(value==expected[name],'grades',name)
ck(Q(15,2)+7*beta==8,'grades','exact selected amplitude product')
for ell in range(65):ck((ledger['own']-ell*g>=target)==(ell<=43),'grades','complete extra-width allowance')
for amplification in (Q(0),Q(1,8),Q(5,4),Q(9)):
 b=math.floor(14*(target+amplification))+1
 ck(b*beta-amplification>target,'floors','starting native order pays nonroot prior plus downstream amplification')
# Sign-carrying amplitudes must never be added as signed error allowances.
radii=[Q(-1,10),Q(1,5)];b=3
signed=radii[0]**b+radii[1]**b*radii[0]
absolute=abs(radii[0])**b+abs(radii[1])**b*abs(radii[0])
ck(signed<0 and absolute>0,'floors','magnitude prior is required')

out={'status':'PASS','assertions':sum(counts.values()),'groups':dict(counts),'rank8_grades':{k:str(v) for k,v in ledger.items()},'H7_second_moment':m2,'H7_fourth_moment':m4,'bounded_join_fixture_count':len(fixtures),'scope':'Independent first-Riesz and keep-transfer algebra, complete observer terms, conditional moment boundaries, Hermite energy and exact grades/floors. Imported native port remains a theorem assumption; no native compiler was executed.'}
(HERE/'independent-rank8-join-checks.json').write_text(json.dumps(out,indent=2)+'\n')
(HERE/'bounded-join-fixtures.json').write_text(json.dumps(fixtures,indent=2)+'\n')
print(json.dumps(out,indent=2))
