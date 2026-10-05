import itertools,json,collections
from pathlib import Path
P=Path(__file__).parent
base=((0,1),(1,2),(3,4),(4,5))

def path(edges,a,b):
    adj={i:[] for i in range(6)}
    for u,v in edges:adj[u].append(v);adj[v].append(u)
    stack=[(a,[a])]
    while stack:
        x,p=stack.pop()
        if x==b:return p
        for y in adj[x]:
            if y not in p:stack.append((y,p+[y]))
    return None

def record(kind,cross,marks):
    choices=[]
    cuts=[None] if kind=='M' else list(cross)
    for cut in cuts:
        edges=base+tuple(e for e in cross if e!=cut)
        for a in marks:
            for b in marks:
                if a>=b:continue
                p=path(edges,a,b)
                if (a<3)!=(b<3):choices.append((len(p),p,cut,edges))
    n,p,cut,edges=max(choices,key=lambda x:(x[0],x[1],x[2] or (-1,-1)))
    assert 4<=n<=6
    assert 1 in p and 4 in p
    sides=set(range(6))-set(p)
    assert len(sides)==6-n<=2 and sides.isdisjoint({1,4})
    assert all(sum(x in e for e in edges)==1 for x in sides)
    assert set(marks)&set(range(3)) and set(marks)&set(range(3,6))
    if cut:
        assert cut[0]<3 and cut[1]>=3
        assert sum(v in cut for v in range(3))==1
        assert sum(v in cut for v in range(3,6))==1
    return dict(kind=kind,cross=cross,marks=marks,cut=cut,spine=p,spine_forces=n,
                side_forces=sorted(sides),first_conditional_force=2*n,
                auxiliary_covariance_force=12 if cut else None,
                coarse_covariance_force=12)
rows=[]
for a in range(3):
 for b in range(3,6):
    rows.append(record('M',((a,b),),[x for x in range(6) if x not in (a,b)]))
for left in itertools.combinations(range(3),2):
 for right in itertools.combinations(range(3,6),2):
  for perm in itertools.permutations(right):
    cross=tuple(zip(left,perm))
    rows.append(record('N',cross,[x for x in range(6) if x not in (*left,*right)]))
assert len(rows)==27
counts={k:dict(collections.Counter(r['spine_forces'] for r in rows if r['kind']==k)) for k in ['M','N']}
assert counts['M']=={6:4,5:4,4:1}
assert counts['N']=={6:4,5:12,4:2}
# Local slot conservation. Total potential-derivative slots are 2,3,2 in each copy.
for r in rows:
 for v in range(6):
    d=sum(v in e for e in base)+sum(v in e for e in r['cross'])
    assert d+int(v in r['marks'])==[2,3,2,2,3,2][v]
# A scalar exact native conditional diagnostic for the worst four-force spine:
# X,Y iid N(0,1), aux G iid; side outputs S=alpha*G+sqrt(1-alpha²)*X,
# T=alpha*G+sqrt(1-alpha²)*Y. V=lambda theta²*S*T, lambda=s*alpha⁴.
# E_ST V=s alpha⁶ G² theta²; Var_ST(V)/2=
# alpha⁸ theta⁴/2*((1-alpha²)²+2alpha²(1-alpha²)G²).
# After G average the complete order-lambda² log coefficient is
# alpha⁸ theta⁴/2*(1+alpha⁴). First unwanted grade8, and grade12 auxiliary term.
import sympy as s
alpha,G,theta,lam=s.symbols('alpha G theta lam')
mu=alpha*G;cov=1-alpha**2
variance=cov**2+2*cov*mu**2
avg_var=s.expand(variance).subs(G**2,1)
aux_var=2*alpha**4
assert s.expand(avg_var+aux_var)==1+alpha**4
assert s.expand((alpha**8/2)*(avg_var+aux_var))==alpha**8/2+alpha**12/2
# Same-endpoint cancellation amplitude identities.
t=s.symbols('t');a=1-t**3
assert s.expand(s.diff(-a*a,t)-6*t*t*a)==0
assert s.expand(-6*t*t*a+ s.diff(-a*a,t))==0
assert s.expand(-12*t*t*a+2*s.diff(-a*a,t))==0
# Scalar Hermite-map check, exactly through skew-square order.
x,k=s.symbols('x k');H1=x;H2=x*x-1;H3=x**3-3*x
f=x+k*H2-k*k*H1-2*k*k*H3

def gauss(poly):
 out=0
 for (n,),c in s.Poly(s.expand(poly),x).terms():
  if n%2==0:out+=c*(s.factorial2(n-1) if n else 1)
 return s.expand(out)
mean=gauss(f);var=gauss(f*f)-mean*mean
m3=gauss((f-mean)**3);k4=gauss((f-mean)**4)-3*var**2
assert mean==0
assert s.expand(var).coeff(k,2)==0
assert s.expand(k4).coeff(k,2)==0
assert s.expand(m3).coeff(k,1)==6
out={'status':'PASS','graph_count':len(rows),'spine_counts':counts,'graphs':rows,
     'scalar_worst_spine_order_lambda_squared':'alpha^8/2 + alpha^12/2',
     'scalar_corrected_map_variance':str(var),
     'scope':'Finite graph/jet diagnostics; native compilers are imported, not executed.'}
(P/'six_force_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='graphs'},indent=2))
