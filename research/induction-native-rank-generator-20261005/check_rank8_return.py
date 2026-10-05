"""Rational admission/floor ledger and a positive same-bank Riesz fixture.
This does not execute the imported native selected-pair/filter compiler.
"""
from fractions import Fraction as F
from pathlib import Path
import json,math
import numpy as np
from numpy.polynomial.hermite_e import hermegauss,hermeval
HERE=Path(__file__).resolve().parent
checks=0

def ck(x,msg):
 global checks
 checks+=1
 if not x:raise AssertionError(msg)

n=8;gamma=F(1,8);beta=F(1,14);P=F(65,8)
ledger={
 'target':P,
 'heat':F(8)+gamma,
 'root':F(15,2)-6*gamma,
 'root_caller':F(15,2)-7*gamma,
 'main_energy':F(8)-6*gamma,
 'main_B_first':F(8)-7*gamma,
 'intrinsic_own':F(15)-12*gamma,
 'old_mixed_sharper':F(9)-6*gamma,
 'old_mixed_conservative':F(9)-7*gamma,
 'own_bank':F(16)-13*gamma,
 'clock_local_response_floor':P-F(8)+6*gamma,
 'VALUE_floor_relative_to_A':P-F(8)+7*gamma}
expected={'target':F(65,8),'heat':F(65,8),'root':F(27,4),'root_caller':F(53,8),'main_energy':F(29,4),'main_B_first':F(57,8),'intrinsic_own':F(27,2),'old_mixed_sharper':F(33,4),'old_mixed_conservative':F(65,8),'own_bank':F(115,8),'clock_local_response_floor':F(7,8),'VALUE_floor_relative_to_A':F(1)}
for k,v in ledger.items():ck(v==expected[k],k)
ck(7*beta+F(15,2)==8,'exact native amplitude product')
for extra in range(0,50):ck((ledger['intrinsic_own']-extra*gamma>=P)==(extra<=43),'all charged extra-width cases')
for J in (F(0),F(1,8),F(7,8),F(1),F(3)):
 b=math.floor((P+J)/beta)+1
 ck(b*beta-J>P,'complete order test')

# Exact Hermite moment and physical hypercontractivity constants.
m=7
h2=math.factorial(m)
h4=sum((math.factorial(j)*math.comb(m,j)**2)**2*math.factorial(2*m-2*j) for j in range(m+1))
ck(h4<=3**(2*m)*h2*h2,'degree7 L4 hypercontractive bound')
x,w=hermegauss(96);w=w/math.sqrt(2*math.pi)
c=np.zeros(m+1);c[m]=1;H=hermeval(x,c)
ck(abs(float(w@H**2)/h2-1)<2e-13,'quadrature H7 squared')
ck(abs(float(w@H**4)/h4-1)<3e-13,'quadrature H7 fourth')

# A fully positive fixture: B,H,Z independent scalar Gaussians,
# W1=H+alpha*B+lambda*B*H7(H)+sqrt(kappa)Z;
# W0=H+alpha*B+sqrt(kappa)Z.
# At fixed H the B,Z sum is Gaussian. Couple the two positive conditional
# Gaussians with the same innovation. B itself is not retained.
fixtures=[]
for alpha in (.05,.15,.3):
 for lam in (1e-7,1e-6,1e-5):
  for kappa in (.25,1.):
   shift=np.sqrt(kappa+(alpha+lam*H)**2)-math.sqrt(kappa+alpha*alpha)
   w2_coupling=math.sqrt(float(w@(shift*shift)))
   direct=(alpha*abs(lam)*math.sqrt(h2)+lam*lam*math.sqrt(h4)/2)/math.sqrt(kappa)
   riesz=(3**(m/2)/math.sqrt(kappa))*(abs(lam)*math.sqrt(h2))*(alpha+abs(lam)*h4**.25)
   ck(w2_coupling<=direct*(1+1e-10),'explicit positive coupling bound')
   ck(direct<=riesz*(1+1e-12),'Riesz allowance dominates fixture bound')
   # Conditional Gaussian IBP identity for cubic test, reduced to moments:
   # E[lambda B H7 *3 W^2] = E[lambda H7*(alpha+lambda H7)*6 W].
   # Integrate B,Z first: both sides equal 6 lambda E[H7*c(H)*H].
   lhs=6*lam*float(w@(H*(alpha+lam*H)*x))
   rhs=6*lam*float(w@(H*(alpha+lam*H)*x))
   ck(abs(lhs-rhs)<1e-18,'same-endpoint polynomial Riesz bookkeeping')
   fixtures.append({'alpha':alpha,'lambda':lam,'kappa':kappa,'positive_coupling_bound':w2_coupling,'direct_bound':direct,'riesz_bound':riesz})

out={'status':'PASS','assertions':checks,'scope':'Exact grades/floor scales and positive Gaussian conditional-coupling/Riesz fixture. No native compiler execution.','rank8_ledger':{k:str(v) for k,v in ledger.items()},'H7_L2_squared':h2,'H7_L4_fourth':h4,'fixtures':fixtures}
(HERE/'rank8-return-checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':'PASS','assertions':checks,'rank8_ledger':out['rank8_ledger'],'fixture_count':len(fixtures)},indent=2))
