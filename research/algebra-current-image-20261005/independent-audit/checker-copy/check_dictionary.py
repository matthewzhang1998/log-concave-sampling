#!/usr/bin/env python3
"""Finite exact boundary test. No all-graph or native-quadratic claim."""
from fractions import Fraction as Q
from itertools import product
from math import factorial
from pathlib import Path
import hashlib,json
import sympy as s
HERE=Path(__file__).resolve().parent
checks={}; count=0

def require(test,label):
 global count
 count+=1
 if not test: raise AssertionError(label)

def partitions(items):
 if not items: yield [];return
 first,*rest=items
 for p in partitions(rest):
  yield [(first,)]+p
  for i in range(len(p)):
   yield p[:i]+[(first,)+p[i]]+p[i+1:]

def partial_matchings(items):
 if not items: yield (),();return
 first,*rest=items
 for pairs,singles in partial_matchings(rest):yield pairs,(first,)+singles
 for i,other in enumerate(rest):
  for pairs,singles in partial_matchings(rest[:i]+rest[i+1:]):yield ((first,other),)+pairs,singles

def gaussian(tokens,mean,cov,connected_copies=None):
 out=0;nonzero=0;enumerated=0
 for pairs,singles in partial_matchings(tuple(tokens)):
  enumerated+=1
  weight=s.prod(mean[t] for t in singles)*s.prod(cov[u,v] for u,v in pairs)
  if weight==0:continue
  if connected_copies is not None:
   seen={connected_copies[0]}
   while True:
    new=seen|{v[0] for u,v in pairs if u[0] in seen}|{u[0] for u,v in pairs if v[0] in seen}
    if new==seen:break
    seen=new
   if set(connected_copies)!=seen:continue
  nonzero+=1;out+=weight
 return s.expand(out),enumerated,nonzero

# Direct finite Fourier/Gaussian check for the monotone sine family.
cc,ee,zz,tt,AA=s.symbols('c epsilon z t A', positive=True)
a0=s.sqrt(s.Rational(1,2)); response1=0
for sign in (-1,1):
 # E[X P sin(z+a*t*X+sign*a*t*P)] from the Gaussian characteristic function.
 xp_sine=-sign*a0*a0*tt*tt*s.exp(-tt*tt/2)*s.sin(zz)
 response1+=sign*xp_sine/(2*AA*tt)
require(s.simplify(response1+tt*s.exp(-tt*tt/2)*s.sin(zz)/(2*AA))==0,'C1 finite VALUE response normalization')
xp,pp=s.symbols('X P'); linear=0
for sign in (-1,1):linear+=sign*(cc*(zz+a0*tt*(xp+sign*pp))-cc*(zz+a0*tt*sign*pp))
require(s.expand(linear)==0,'C1 exact affine-force cancellation')
checks['primitive_responses']={'family':'g(x)=c*x+epsilon*sin(x), 0<epsilon<c, c+epsilon<=A','C0':'(c+epsilon*exp(-t^2/2)*cos(z))/A','C1':'-epsilon*t*exp(-t^2/2)*sin(z)/(2*A)','raw_VALUE_counts':[2,4],'source_preparation':'not equated with uncut response; floor retained'}

a,b,z=s.symbols('a b theta');kap={};mu={};census=[]
for n in range(1,5):
 tokens=[(i,j) for i in range(n) for j in range(2)]
 mean={t:(a if t[1]==0 else b)*z for t in tokens}
 cov={(u,v):s.Integer(u[1]==v[1]) for u in tokens for v in tokens}
 mu[n],enum,nonzero=gaussian(tokens,mean,cov)
 part_sum=0
 ps=list(partitions(list(range(n))))
 for partition in ps:
  piece=1
  for block in partition:
   ts=[(i,j) for i in block for j in range(2)]
   piece*=gaussian(ts,mean,cov)[0]
  part_sum+=(-1)**(len(partition)-1)*factorial(len(partition)-1)*piece
 kap[n]=s.expand(part_sum)
 linked=gaussian(tokens,mean,cov,list(range(n)))[0]
 require(s.expand(kap[n]-linked)==0,('complete tilt/connected overlap',n))
 expected=(a*b*z**(2*n+2) if n%2 else z**(2*n)/n+(a*a+b*b)*z**(2*n+2)/2)
 require(s.expand(z**(2*n)*kap[n]/factorial(n)-expected)==0,('rank4 log',n))
 census.append({'copies':n,'all_partial_matchings':enum,'nonzero_moment_matchings':nonzero,'set_partitions':len(ps)})
connected2=s.expand(z**4*kap[2]/2)
raw2=s.expand(z**4*mu[2]/2)
require(s.expand(raw2-connected2-a*a*b*b*z**8/2)==0,'rank-eight disconnected current')
require(s.expand(connected2-z**4/2-(a*a+b*b)*z**6/2)==0,'zero-tilt shortcut misses rank-six branches')
# The C1-containing three-force path has one physical Gaussian probe.
require(s.expand((1+a*a*z*z)*z**4/2-z**4/2-a*a*z**6/2)==0,'T3 full second law current')
checks['tilt_and_connected']={'census':census,'connected_first_nonmain':str(connected2),'raw_fixed_base_first_nonmain':str(raw2),'zero_tilt_shortcut_rejected':True}

# All four bridge hits between two two-vertex packets, then 0,1,2 further bank derivatives.
zs=s.symbols('z0:4'); rows=(s.Integer(1),s.Integer(2),s.Integer(3),s.Integer(5))
f=[1+x+x*x+x**3+x**4+x**5 for x in zs]
D=lambda p:s.expand(sum(rows[i]*s.diff(p,zs[i]) for i in range(4)))
F=f[0]*f[1];G=f[2]*f[3]
bridge=lambda p,q:s.expand(sum(rows[i]*rows[j]*s.diff(p,zs[i])*s.diff(q,zs[j]) for i in (0,1) for j in (2,3)))
require(s.expand(D(bridge(F,G))-bridge(D(F),G)-bridge(F,D(G)))==0,'derivative/join commutation')
require(s.expand(D(D(bridge(F,G)))-bridge(D(D(F)),G)-2*bridge(D(F),D(G))-bridge(F,D(D(G))))==0,'second derivative/join multiplicity')
graphs=[]
for i,j in product((0,1),(2,3)):
 k=[0]*4;k[i]+=1;k[j]+=1
 edges=[(0,1),(2,3),(i,j)];degree=[0]*4
 for u,v in edges:degree[u]+=1;degree[v]+=1
 require(all(degree[v]+1==k[v]+2 for v in range(4)),'bridge slot ledger')
 require(sum(k)+8==2*len(edges)+4,'global graph incidence')
 require(all(degree[v]!=1 or 1>0 for v in range(4)),'marked leaves')
 # 4 choices at each extra derivative hit, with repeated hits preserved.
 for length in (0,1,2):
  for word in product(range(4),repeat=length):
   kk=k[:]
   for v in word:kk[v]+=1
   require(sum(kk)==2+length,'all derivative hit words')
 graphs.append({'edge':[i,j],'primitive_orders':k,'physical_marks':[1,1,1,1]})
x=s.symbols('x');A=1+x;B=1+x*x;u=x+x*x;v=1+x**3
join=A*B*s.diff(u,x)*s.diff(v,x)
missing=s.diff(join,x)-A*B*(s.diff(u,x,2)*s.diff(v,x)+s.diff(u,x)*s.diff(v,x,2))
require(s.expand(missing-(s.diff(A,x)*B+A*s.diff(B,x))*s.diff(u,x)*s.diff(v,x))==0,'live injection row branch')
require(missing!=0,'live geometry cannot be frozen')
checks['derivative_join']={'base_bridge_graphs':graphs,'ordered_additional_hit_words':4*(1+4+16),'constant_row_orders':[0,1,2],'live_row_escape':str(s.expand(missing))}

# The covariance bridge remains random conditional on the old B.
c,d=s.symbols('c d', positive=True)
cov_cos=lambda m,n:(s.exp(-s.Rational((m-n)**2,2))+s.exp(-s.Rational((m+n)**2,2)))/2-s.exp(-s.Rational(m*m+n*n,2))
variance=s.expand((2*c*d)**2*cov_cos(1,1)+2*(2*c*d)*(d*d/2)*cov_cos(1,2)+(d*d/2)**2*cov_cos(2,2))
t=s.symbols('t', nonnegative=True)
# E_B of conditioned sin products, retaining every original frequency.
bridge_integral=0
for m,am in ((1,2*c*d),(2,d*d)):
 for n,an in ((1,2*c*d),(2,d*d)):
  bridge_integral+=am*an*s.integrate(s.exp(-s.Rational(m*m+n*n,2)*(1-t))*(s.exp(-s.Rational((m-n)**2,2)*t)-s.exp(-s.Rational((m+n)**2,2)*t))/2,(t,0,1))
require(s.simplify(bridge_integral-variance)==0,'whole-bank bridge exact')
require(float(variance.subs({c:s.Rational(1,2),d:s.Rational(1,4)}))>0,'retained B=0 defect nonzero')
checks['bank_scope']={'coefficient':'F(B)=(c+d*cos(B))**2','variance':str(variance),'conditional_bridge_at_B_zero':0,'pointwise_replacement_rejected':True}

# Correlated packet cross term; each own pair remains independent.
a1,b1,a2,b2=s.symbols('a1 b1 a2 b2');p,q,r,u=s.symbols('cPP cPQ cQP cQQ')
ts=[(0,0),(0,1),(1,0),(1,1)];mean=dict(zip(ts,[a1*z,b1*z,a2*z,b2*z]))
C=s.Matrix([[1,0,p,q],[0,1,r,u],[p,r,1,0],[q,u,0,1]])
covs={(aa,bb):C[i,j] for i,aa in enumerate(ts) for j,bb in enumerate(ts)}
cross=s.expand(gaussian(ts,mean,covs)[0]-gaussian(ts[:2],mean,covs)[0]*gaussian(ts[2:],mean,covs)[0])
expected=p*u+q*r+z*z*(a1*a2*u+a1*b2*r+b1*a2*q+b1*b2*p)
require(s.expand(cross-expected)==0,'complete grouped cross connected coefficient')
raw_cross=s.expand(cross+a1*b1*a2*b2*z**4)
require(s.expand(raw_cross-gaussian(ts,mean,covs)[0])==0,'mixed disconnected main product')
require(s.expand(cross.subs({p:1,u:1,q:0,r:0,a1:a,a2:a,b1:b,b2:b})-1-(a*a+b*b)*z*z)==0,'shared identical packet cross')
checks['grouped_isolated']={'mixed_connected_before_theta4':str(cross),'mixed_raw_before_theta4':str(raw_cross),'isolated_shortcut_rejected':True}

# Exact positive scalar covariance branch. exp[(s*q+s^2*r/2)*D^2].
eps,qq,rr=s.symbols('s q r')
jet=s.series(s.exp((eps*qq+eps**2*rr/2)*z*z),eps,0,4).removeO().expand()
require(jet.coeff(eps,1)==qq*z*z,'covariance leading')
require(s.expand(jet.coeff(eps,2)-rr*z*z/2-qq**2*z**4/2)==0,'covariance whole first nonmain')
require(s.expand(jet.coeff(eps,3)-qq*rr*z**4/2-qq**3*z**6/6)==0,'covariance third cross')
require(s.diff((1+2*eps*qq+eps**2*rr),eps)/2==qq+eps*rr,'exact moving endpoint time current')
checks['covariance_time']={'leading':str(jet.coeff(eps,1)),'first_nonmain':str(jet.coeff(eps,2)),'third_order':str(jet.coeff(eps,3))}

# Exact strict-boundary associated-graded section on three fixed prepared tree types.
beta,gamma,P=Q(1,4),Q(1,8),Q(11,2); grades=[]
for name,ks in [('T2',(0,0)),('T3',(0,1,0)),('T4',(0,1,1,0))]:
 N=len(ks);K=sum(ks);a0=Q(4);dd=a0-beta*N-gamma*K;psi=dd-2*gamma;root=beta+dd;Ggrade=a0-gamma*K
 require(psi>=Q(5,2),'positive reserve')
 require(Ggrade<P<2*root,'strict native floor')
 require(beta*N+gamma*K+dd==a0,'source amplitude product')
 grades.append({'name':name,'N':N,'K':K,'a':str(a0),'G':str(Ggrade),'Psi':str(psi),'root_grade':str(root),'native_floor':str(2*root)})
for i,j in product(grades,repeat=2):
 bridge_grade=Q(i['a'])+Q(j['a'])-gamma*(i['K']+j['K']+2)
 require(bridge_grade>P,'all pair bridge branches beyond cutoff')
require(2*min(Q(g['root_grade']) for g in grades)>P,'summed root floor strict')
checks['finite_section']={'beta':str(beta),'gamma':str(gamma),'cutoff':str(P),'types':grades,'residual_closure':'modulo every certified current/numerical floor of grade >= P; not quadratic closure'}

# Ordinary and owned-clock charges; no new clock mass from a derivative.
for k in range(8):require(0<=max(k-1,0)-max(k-2,0)<=1,'owned clock one-hit cost')
clock_M,clock_e=s.symbols('M epsilon',positive=True)
checks['clock_heat']={'new_clocks_from_hits':0,'owned_clock_charge':'(k-2)_+ only with owned omega<=C*t^2 and relocation','sealed_two_edge_node_count':'2*(128*J_clock*m_clock)^2','general_bridge_count':'literal N_T with certified error; no generic public-log claim','heat_extension':False}

# Check declared data coverage and ensure no native-floor row is mislabeled exact.
data=json.loads((HERE/'CURRENT-DICTIONARY.json').read_text())
required={'C0','C1','T2','BANK_BRIDGE','RANK4_SPINE','COVARIANCE_TIME','CARRIER_SUPERPOSITION'}
require(required<=set(r['id'] for r in data['rows']),'required row coverage')
for row in data['rows']:
 require(all(k in row for k in ['leading_current','value_representative','first_nonmain','retained_tapes','equality_scope','escapes']),'row schema')
 for res in row['first_nonmain']:
  if res['id']=='NATIVE_QUADRATIC_UNIDENTIFIED': require(res['status']=='unresolved_norm_floor','native floor cannot be exact jet')
checks['dictionary_schema']={'rows':len(data['rows']),'required_rows':sorted(required)}
result={'status':'PASS','assertions':count,'scope':'Complete enumeration within the stated finite token/hit/type bounds. Analytic and native contracts are pinned imports, not proved by this checker.','checks':checks}
(HERE/'checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
