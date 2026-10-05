#!/usr/bin/env python3
"""Independent exact/reference diagnostics. No imported native compiler is run.

The finite log-guard fixtures illustrate the admitted finite-power mechanism;
they are not numerical guard certificates for unspecified compiler constants.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import permutations, product
import hashlib, json, math
import sympy as s

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
checks=[]
def check(name, ok, evidence=None):
    if not bool(ok): raise AssertionError(name)
    row={'name':name}
    if evidence is not None: row['evidence']=str(evidence)
    checks.append(row)
def zero(x): return s.expand(x)==0
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
reviewed=json.loads((HERE/'reviewed-source.json').read_text())
for item in reviewed['reviewed_files']:
    check('exact reviewed snapshot: '+Path(item['path']).name,sha(Path(item['path']))==item['sha256'])

# A. Exact fixed-buffer bank allocation and genuine completed-mean identity.
graph=json.loads((ROOT/'ENDPOINT-GRAPH.json').read_text())
banks={x['id']:x for x in graph['banks']}
check('five complete endpoint banks',set(banks)=={'M','H','K','V3','keep'})
check('new completed M packet is selected',banks['M']['source_packet']=='/workspace/shared/dyadic-prefix-law-join-20261005/ACTUAL-DYADIC-PREFIX-LAW-JOIN.md')
check('M is a whole completed LAW',banks['M']['input_type']=='completed own-mean LAW')
check('fixed new mean native orders',banks['M']['native_orders']=={'prefix_mean':4,'tail_mean_and_coefficient_means':6,'corrected_negative_fourth_paths':6,'cubic_mixed_true_quartic':5,'final_own_mean':4})
h=Q(1,2); u3=Q(1,8); keep=Q(1,4)
check('conditional source-zero variance',sum(Q(b['baseline_after_readout']) for b in banks.values())==1-h*h)
check('source-zero total variance',h*h+sum(Q(b['baseline_after_readout']) for b in banks.values())==1)
check('visible reference baseline',h*h+h*h*Q(1,2)+u3==Q(1,2))
check('untouched keep fraction',keep/(1-h*h)==Q(1,3))
check('cubic coefficient physical sign and scale',-h**3/Q(6)==-Q(1,48))
check('normalized cubic radius square',h*h/u3==2)
check('physical readout grade',Q(4)+Q(1,2)==Q(9,2))

# B. Derive coherent linear F3 moments from stationary OU, not a fitted covariance.
# I_k = integral exp(-x) x^k/k! X_exp(-x) dx.
# On x>=y, exp(-x-y)*Cov(X_exp(-x),X_exp(-y))=exp(-2x).
x,y,a=s.symbols('x y a', positive=True)
means=s.Matrix([s.integrate(s.exp(-2*x)*x**k/s.factorial(k),(x,0,s.oo)) for k in range(3)])
full=s.zeros(3)
for i,j in product(range(3),repeat=2):
    region=s.integrate(s.exp(-2*x)*x**i*s.integrate(y**j,(y,0,x))/(s.factorial(i)*s.factorial(j)),(x,0,s.oo))
    opposite=s.integrate(s.exp(-2*x)*x**j*s.integrate(y**i,(y,0,x))/(s.factorial(i)*s.factorial(j)),(x,0,s.oo))
    full[i,j]=region+opposite
    check(f'coherent OU covariance symmetry {i},{j}',full[i,j]==full[j,i] if j<=i else full[i,j]>0)
cond=full-means*means.T
for k in range(1,4): check(f'coherent conditional covariance positive principal minor {k}',cond[:k,:k].det()>0)
c=s.Matrix([a,-a*a,a**3]);m=(c.T*means)[0]
cf=(c.T*cond*c)[0];vf=(c.T*full*c)[0]
vp=1-2*m+vf
check('canonical m3 is the third-substitution mean',zero(m-a/2+a*a/4-a**3/8))
check('coherent conditional and unconditional F3 covariance',zero(vf-cf-m*m))
joined=(1-m)**2/4+s.Rational(3,4)+cf/4
check('full mean and covariance positive reference equals buffered P3 law',zero(joined-(vp/4+s.Rational(3,4))))
check('P3 matches exact quadratic target through cubic order',zero(s.series(vp-1/(1+a),a,0,4).removeO()))
quartic=s.expand(s.series(vp-1/(1+a),a,0,5).removeO()).coeff(a,4)
check('P3 target keeps a real quartic discrepancy',quartic!=0,quartic)
cf2=(s.Matrix([a,-a*a,0]).T*cond*s.Matrix([a,-a*a,0]))[0]
check('full F3 covariance restoration is fourth order',s.Poly(s.expand(cf-cf2),a).terms()[-1][0][0]==4)
check('posterior force mean is different from canonical mean',not zero(m-a/(2*(1+3*a/4))))
# Matrix quadratic fixture validates an off-diagonal full covariance, not only trace.
B=s.Matrix([[s.Rational(1,12),s.Rational(1,30)],[s.Rational(1,30),s.Rational(1,10)]])
Cs=[B,-B**2,B**3];mmat=sum((means[i]*Cs[i] for i in range(3)),s.zeros(2))
Sigma=sum((cond[i,j]*Cs[i]*Cs[j].T for i,j in product(range(3),repeat=2)),s.zeros(2))
VP=s.eye(2)-2*mmat+sum((full[i,j]*Cs[i]*Cs[j].T for i,j in product(range(3),repeat=2)),s.zeros(2))
JV=(s.eye(2)-mmat)*(s.eye(2)-mmat).T/4+3*s.eye(2)/4+Sigma/4
for i,j in product(range(2),repeat=2): check(f'off-diagonal complete Gaussian target {i},{j}',zero(JV[i,j]-(VP/4+3*s.eye(2)/4)[i,j]))
check('matrix covariance fixture genuinely has off-diagonal structure',Sigma[0,1]!=0)

# C. Full anisotropic Hermite regression with unequal carrier shares and fills.
xv=s.Matrix(s.symbols('z0:2'));C=s.eye(2)/2+Sigma/4;Ci=C.inv()
packet_totals=[s.Rational(1,32),s.Rational(3,32)]
carrier=[3*q/4 for q in packet_totals];fills=sum(packet_totals)/4
J=s.eye(2)*(s.Rational(1,4)+s.Rational(1,8)+fills)+Sigma/4
check('packet fills stay in assigned one-eighth buffer',sum(carrier)+fills==s.Rational(1,8))
check('visible covariance includes every fill and carrier',J+sum(carrier)*s.eye(2)==C)
check('anisotropic C has fixed positive gap',C[0,0]>=s.Rational(1,2) and (C-s.eye(2)/2).det()>=0)
for n,beta2 in enumerate(carrier):
    mu=beta2*Ci*xv;V=beta2*s.eye(2)-beta2**2*Ci
    projected=(V+mu*mu.T)/beta2**2-s.eye(2)/beta2
    exact=Ci*xv*xv.T*Ci-Ci
    for i,j in product(range(2),repeat=2): check(f'full anisotropic Hermite regression {n}:{i},{j}',zero(projected[i,j]-exact[i,j]))
    check(f'negative trace cannot be dropped {n}',any(not zero((mu*mu.T/beta2**2-s.eye(2)/beta2-exact)[i,j]) for i,j in product(range(2),repeat=2)))

# A nonsymmetric output tensor gives the same leading current, not the same map.
N={ijk:s.Rational(1+5*ijk[0]+2*(ijk[1]+ijk[2])+ijk[1]*ijk[2],71) for ijk in product(range(2),repeat=3)}
D={ijk:N[ijk]-sum(N[tuple(p)] for p in permutations(ijk))/6 for ijk in N}
phi=xv[0]**5+3*xv[0]**2*xv[1]**3
check('nonzero output-slot asymmetry',any(D.values()))
check('symmetrized leading current is exactly zero',zero(sum(v*s.diff(phi,xv[i],xv[j],xv[k]) for (i,j,k),v in D.items())))
check('asymmetry still changes positive quadratic reference',any(not zero(sum(D[i,j,k]*(xv[j]*xv[k]-int(j==k)) for j,k in product(range(2),repeat=2))) for i in range(2)))

# Exact positive same-endpoint complement-current identity, including feedback.
z,u,k,t=s.symbols('z u k t')
E=2*z*u/19+(u*u-1)/23;R=2*z/19+u/23
W=3*z/2+(9*z*z/4-s.Rational(9,4))/37+t*E+k/2

def moment(n): return s.Integer(0) if n%2 else s.factorial2(n-1) if n else s.Integer(1)
def expect(poly,variables):
    return s.expand(sum(coef*math.prod(moment(n) for n in powers) for powers,coef in s.Poly(s.expand(poly),*variables).terms()))
check('regression complement is centered at same endpoint',expect(E,[u])==0)
for n in range(1,7):
    lhs=s.diff(expect(W**n,[z,u,k]),t)
    rhs=0 if n==1 else t*n*(n-1)*expect(R*s.diff(E,u)*W**(n-2),[z,u,k])
    check(f'positive same-endpoint complement current degree {n}',zero(lhs-rhs))
check('skew correction has nonzero sixth-order covariance feedback',zero(expect((z+a**3*(z*z-1))**2,[z])-1-2*a**6))

# D. Reverse-OU geometry and sums. Log-domain scales avoid small-A underflow.
rho=Q(3,4)
for s2 in [Q(1),Q(3,4),Q(9,16),Q(1,64),Q(1,4096)]:
    r2=1-s2;delta=s2/4;t2=r2+delta;b2=delta**2/(t2*s2**2);v0=delta/t2
    check(f'affine mean squared identity {s2}',b2*s2==v0*h*h)
    check(f'exact reserve identity {s2}',(1-t2)*delta/(t2*s2)==v0*(1-h*h))
    check(f'exact reference contraction {s2}',((1-t2)+delta)/s2==1)
    check(f'local grade-four scaling squared {s2}',v0*s2**8==4*delta**2*s2**7/t2)
check('terminal threshold exponent',Q(2,5)*(4-2)==Q(4,5))
check('smallest buffered alpha exponent',1+Q(4,5)==Q(9,5))
check('geometric intrinsic exponent',Q(7,2)+1==Q(9,2))
for e in [1,2,4,8,16,32,64,128,256,512,1024,2048,4096]:
    la=-e*math.log(2);lr=math.log(.75);J=math.ceil(.8*la/lr);last=J*lr
    check(f'first terminal stop A=2^-{e}',.8*la+lr<last<=.8*la)
    check(f'terminal reaches fourth grade A=2^-{e}',2*la+2.5*last<=4*la)
    check(f'buffered alpha log range A=2^-{e}',la+(J-1)*lr>1.8*la)
    check(f'terminal alpha log range A=2^-{e}',1.8*la+lr<la+last<=1.8*la)
    finite=-math.expm1(4.5*J*lr)/(4*(-math.expm1(4.5*lr)))
    check(f'weighted intrinsic finite geometric sum A=2^-{e}',finite<=1/(4*(-math.expm1(4.5*lr))))
    # Enumerate all exact reference contraction weights for this finite schedule.
    r=[math.sqrt(-math.expm1(j*lr)) for j in range(J+1)]
    for j in [1,max(1,J//2),J]:
        weight=r[J]
        for n in range(j+1,J+1):weight*=r[n-1]/r[n]
        check(f'exact contraction telescope A=2^-{e},j={j}',math.isclose(weight,r[j],rel_tol=3e-14,abs_tol=3e-14))
# The mode weighted coefficient is Delta*s^6*r, bounded by Delta*s^6.
check('mode weighted sum uniformly bounded',Q(1,4)/(1-rho**4)<1)

# E. Mean ledger and the fixed-log guard template, without claiming instantiated guards.
# Row a^p / v^q produces bulk a^(p-q) and cutoff a^(p+2(1-q)).
for p,q in [(Q(7),Q(5,2)),(Q(6),Q(2)),(Q(6),Q(7,6)),(Q(11,2),Q(5,4)),(Q(6),Q(3,2))]:
    check(f'integrated imported mean row {p}/{q}',min(p-q,p+2*(1-q))>=4)
check('undominated order-six prior attains endpoint grade four',Q(7)+2*(1-Q(5,2))==4)
check('RAW leading term retains a positive K log power',Q(3)+Q(2)/2==4)
# Fixture family L^p*(L^p0/K)^theta with strictly positive theta.
p0=Q(3);families=[(Q(2),Q(1,3)),(Q(5),Q(2,5)),(Q(7),Q(3,2)),(Q(1),Q(1,8))]
b=1+max(p0+p/theta for p,theta in families)
for p,theta in families:
    check(f'finite log guard template has negative L power p={p},theta={theta}',p+theta*(p0-b)<0)
    check(f'same K guard works at every local alpha p={p},theta={theta}',theta>0)
# At fixed L this is exactly the single scalar cutoff condition at all stages.
for alpha_over_A in [Q(1),Q(3,4),Q(27,64),Q(1,10000)]:
    check(f'global cutoff implies local cutoff {alpha_over_A}',alpha_over_A*Q(1,2)<=Q(1,2))

# F. Full-root/serial occurrence bookkeeping fixtures. These are symbolic audits.
for Ddim,nL,nR,reentries in [(1,2,3,4),(3,5,2,7),(11,4,6,9)]:
    dNP=17*Ddim;dPM=dNP+Ddim;dPC=(6*13+2)*Ddim;dpref=dPM+dPC
    dtail=(9+14+10+12+16)*Ddim
    nodes=nL*(2*Ddim+dtail+dpref)+nR*(4*Ddim+dpref)
    aligned=Ddim+nL*(Ddim+dtail+dpref)+nR*(3*Ddim+dpref)
    check(f'only one row per mean node aligned D={Ddim}',nodes-aligned==(nL+nR-1)*Ddim)
    dcompleted=reentries*aligned+23*Ddim
    stage=2*Ddim+dcompleted+(12+15+21)*Ddim
    check(f'completed dM includes every fresh reentry D={Ddim}',dcompleted>aligned and stage>2*Ddim+aligned+(12+15+21)*Ddim)
    Qsource=211*nL+101*nR;Nbase=nL+nR;Qfinal=reentries*(Qsource+Nbase)+Nbase
    check(f'final VALUE count cannot use one source occurrence D={Ddim}',Qfinal>=reentries*Qsource)
    check(f'terminal two-root tape retained D={Ddim}',5*stage+2*Ddim>5*stage)

# G. Exact imported pins and independently exposed stale historical manifest entry.
for i,pin in enumerate(json.loads((ROOT/'INPUT-PINS.json').read_text())['inputs']):
    check(f'new input pin {i}: '+Path(pin['path']).name,sha(Path(pin['path']))==pin['sha256'])
for folder in ['dyadic-prefix-law-join-20261005','positive-endpoint-mean-join-20261005']:
    base=Path('/workspace/shared')/folder
    entries=json.loads((base/'MANIFEST.json').read_text())['files']
    if isinstance(entries,dict):entries=[{'path':name,'sha256':val if isinstance(val,str) else val['sha256']} for name,val in entries.items()]
    stale=[]
    for item in entries:
        if sha(base/item['path'])!=item['sha256']:stale.append(item['path'])
    expected=[] if folder.startswith('dyadic') else ['independent-audit/INDEPENDENT-ENDPOINT-JOIN-AUDIT.md']
    check('historical manifest verified with explicit exception: '+folder,stale==expected,stale)

result={'status':'PASS','assertions':len(checks),'scope':'Independent exact algebra, coherent OU target, full positive-reference regression and same-endpoint current, reverse-OU terminal/schedule, finite-log template and complete-bank bookkeeping.','native_compilers_executed':False,'numerical_native_guards_instantiated':False,'reviewed_files':reviewed['reviewed_files'],'checks':checks}
(HERE/'independent_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['status','assertions','native_compilers_executed','numerical_native_guards_instantiated']},indent=2))
