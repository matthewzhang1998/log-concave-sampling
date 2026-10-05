"""Independent common-H current matching and fixed-rank admission algebra.
This is not a native program implementation.
"""
from pathlib import Path
from fractions import Fraction as Q
import math,json
import sympy as sp
HERE=Path(__file__).resolve().parent
count=0
def ck(ok,msg):
 global count
 count+=1
 if not bool(ok):raise AssertionError(msg)
P1,P2=sp.symbols('P1 P2')
def gm(n):return 0 if n%2 else math.prod(range(1,n,2))
def expectation(p):
 return sum(c*gm(a)*gm(b) for (a,b),c in sp.Poly(sp.expand(p),P1,P2).terms())
sigma1=sp.Rational(3,5);sigma2=sp.Rational(4,5);H=sigma1*P1+sigma2*P2
ck(sigma1**2+sigma2**2==1,'common Gaussian has fixed unit covariance')
for a,b in ((sp.Rational(2,7),sp.Rational(-3,11)),(sigma2**7,-sigma1**7)):
 qnative=a*sp.hermite_prob(7,P1)+b*sp.hermite_prob(7,P2)
 K=a*sigma1**7+b*sigma2**7
 qcommon=K*sp.hermite_prob(7,H)
 for power in range(1,14):
  lhs=power*expectation(qnative*H**(power-1))
  rhs=power*expectation(qcommon*H**(power-1))
  ck(lhs==rhs,'every checked r-linear test current matches at common Gaussian endpoint')
 variance_difference=expectation(qnative**2)-expectation(qcommon**2)
 expected=math.factorial(7)*(a*a+b*b-K*K)
 ck(sp.simplify(variance_difference-expected)==0,'exact quadratic comparison term')
 if K==0:
  ck(expectation(qnative**2)>0 and expectation(qcommon**2)==0,'zero total tensor does not mean independent local Hermite fields cancel')
# The common reference variance contains a mixed r1*r2 coefficient.
ck(-2*math.factorial(7)*sigma1**7*sigma2**7!=0,'cross-group quadratic term is genuinely present')
for values in ((Q(-1,5),Q(2,7)),(Q(1,4),Q(-3,8),Q(2,3)),tuple(Q(i,13) for i in range(-4,5))):
 ck(sum(abs(v) for v in values)**2<=len(values)*sum(v*v for v in values),'N times sum of squares safely charges every cross group')

rows=[]
for n in range(4,41):
 g=Q(1,n);beta=Q(1,2*(n-1));target=n+g;root=n-Q(1,2)-(n-2)*g
 row={'n':n,'tau':g,'beta':beta,'root':root,'root_caller':root-g,'main_energy':n-(n-2)*g,'main_first':n-(n-1)*g,'own':2*root,'old_mixed_exact':n+1-(n-2)*g,'old_mixed_conservative':n+1-(n-1)*g,'self_bank':2*n-(2*n-3)*g,'heat':n+g,'response_floor':(n-1)*g,'value_floor':Q(1)}
 ck(root+(n-2)*g+(n-1)*beta==n,'exact source amplitude product before width estimate')
 ck(row['own']==2*n-3+Q(4,n),'fixed-rank intrinsic own grade')
 ck(row['self_bank']==2*n-2+Q(3,n),'fixed-rank own bank grade')
 ck(row['root_caller']==n-Q(3,2)+g,'fixed-rank caller grade')
 ck(row['old_mixed_conservative']==target and row['heat']==target,'balanced old/heat target')
 ck(row['own']>target and row['self_bank']>target,'strict fixed-rank return surplus')
 ck(row['root_caller']>1 and root+beta-g>1,'actual root and nonroot caller paths')
 ck(row['main_energy']+row['response_floor']==target,'response/quadrature precision amplification')
 ck(row['response_floor']+g==1,'minimum-width VALUE floor')
 threshold=n*n-3*n+3
 ck(2*root-threshold*g==target,'complete inverse-width allowance threshold')
 for J in (Q(0),Q(1,3),Q(5,2)):
  b=math.floor(2*(n-1)*(target+J))+1
  ck(b*beta-J>target,'fixed-rank starting native order')
 rows.append({k:(v if isinstance(v,int) else str(v)) for k,v in row.items()})
out={'status':'PASS','assertions':count,'fixed_ranks_checked':[4,40],'common_reference_fixture':'Two positive polynomial Gaussian groups, with all linear moment currents checked through degree 13, a nonzero cross-group quadratic term, and a cancelling-total-tensor example.','scope':'Common-reference current algebra and fixed-rank exponent/floor consequences only. Native comparison remains imported; no arbitrary-order queue or complexity theorem.','grades':rows}
(HERE/'common-reference-fixed-rank-checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='grades'},indent=2))
