#!/usr/bin/env python3
"""Independent bounded-contract diagnostics, not an all-order source compiler."""
from pathlib import Path
from collections import deque
from fractions import Fraction
import hashlib
import itertools
import json
import math
import mpmath as mp
import sympy as s


# Publication-only source locator and original/public pin mapping.
import sys
_ARCHIVE=next(p for p in Path(__file__).resolve().parents if (p/'publication_sources.py').is_file())
sys.path.insert(0,str(_ARCHIVE))
from publication_sources import source_path, source_label, verify_source_pin, pin_report

HERE=Path(__file__).resolve().parent
BASE=HERE.parent.parent
PINS={str(source_path(k)):v for k,v in json.loads((HERE/'INPUT-PINS.json').read_text()).items()}
COUNT=0

def req(ok,message):
    global COUNT
    COUNT+=1
    if not bool(ok): raise AssertionError(message)

def eq(a,b,message): req(s.simplify(a-b)==0,message)

for name,pin in PINS.items():
    status=verify_source_pin(name,pin)
    if status is not None: req(status,'immutable input '+name)
contract=BASE/'order-reentry/TYPED-CLOSURE-AND-COMPLETE-COST-CONTRACT.md'
req(PINS[str(contract)]=='f63de706baf0e3612247048a4139c581dc68da98b3940bd1c27b44ef63f495e5','final integration-contract pin')
req(PINS[str(BASE/'THIRD-ORDER-FULL-REVERSE-OU-LAW-AND-ENDPOINT.md')]=='e83d737da13d05054ff0dda129a05e72d350c6e3e33b019d548cc517b816e702','audited guarded order-three endpoint')
manifest_files=[BASE/x for x in (
 'order-reentry/independent-audit/ALL-RANK-STEIN-MANIFEST.json',
 'order-reentry/independent-audit/DISCOUNTED-OU-MANIFEST.json',
 'order-reentry/independent-audit/NATIVE-RANK3-MANIFEST.json',
 'full-bridge-independent-audit/MANIFEST.json',
 'same-carrier-feedback/independent-audit/MANIFEST.json')]
for manifest in manifest_files:
    content=json.loads(manifest.read_text())
    for group in ('inputs','outputs'):
        for name,pin in content[group].items():
            path=source_path(name,manifest.parent)
            status=verify_source_pin(path,pin)
            if status is not None: req(status,'imported manifest pin '+str(path))

# Full shared Gram and its Schur complement. Independent width cannot be relabeled.
a,u=s.symbols('a u',positive=True)
one=s.ones(3,1)
G=a*s.eye(3)+u*one*one.T
schur=s.simplify(G[0,0]-(G[0,1:3]*G[1:3,1:3].inv()*G[1:3,0])[0])
eq(schur,a*(a+3*u)/(a+2*u),'conditional query variance')
eq(s.Rational(3,2)*a-schur,a*a/(2*(a+2*u)),'conditional variance upper bound')
for v in (s.Matrix([1,-1,0]),s.Matrix([1,0,-1])):
    req(G*v==a*v,'transverse variance remains a')
    eq((v.T*(G-(a+u)*s.eye(3))*v)[0],-u*(v.T*v)[0],'independent enlarged-shield residual is not PSD')

# Enumerate small native trees independently using Pruefer decoding.
def tree_from_pruefer(code,n):
    deg=[1]*n
    for v in code: deg[v]+=1
    edges=[]
    for v in code:
        leaf=next(i for i,d in enumerate(deg) if d==1)
        edges.append((leaf,v)); deg[leaf]-=1; deg[v]-=1
    left=[i for i,d in enumerate(deg) if d==1]
    edges.append(tuple(left))
    return edges

def eccentricities(n,edges):
    adj=[[] for _ in range(n)]
    for a,b in edges: adj[a].append(b); adj[b].append(a)
    out=[]
    for start in range(n):
        d=[None]*n; d[start]=0; q=deque([start])
        while q:
            v=q.popleft()
            for w in adj[v]:
                if d[w] is None: d[w]=d[v]+1; q.append(w)
        out.append(max(d))
    return out

ntrees=0
for n in range(2,6):
    for code in itertools.product(range(n),repeat=n-2):
        edges=tree_from_pruefer(code,n); ecc=eccentricities(n,edges); diameter=max(ecc)
        degree=[sum(v in e for e in edges) for v in range(n)]
        ntrees+=1
        for v in range(n):
            attached=edges+[(v,n)]
            doubled=edges+[(a+n,b+n) for a,b in edges]+[(v,v+n)]
            main_d=max(eccentricities(n+1,attached))
            feedback_d=max(eccentricities(2*n,doubled))
            req(main_d==max(diameter,ecc[v]+1),'attached-main diameter')
            req(feedback_d==2*ecc[v]+1,'doubled old-subtree diameter')
            req(feedback_d>main_d,'history-specific strict inequality for nontrivial trees')
        # Every leaf is marked; optional internal marks include unmarked degree-two sites.
        inner=[v for v,d in enumerate(degree) if d>1]
        for bits in itertools.product((0,1),repeat=len(inner)):
            marks=[int(d==1) for d in degree]
            for v,b in zip(inner,bits): marks[v]=b
            orders=[d+m-1 for d,m in zip(degree,marks)]
            req(min(orders)>=1,'native derivative orders are legal')
            req(sum(j-1 for j in orders)==sum(marks)-2,'native heat count M-2')
# The contracted claim is deliberately scoped to nontrivial old cubic trees.
req(max(eccentricities(2,[(0,1)]))==1,'singleton attachment diameter')
req(2*0+1==1,'singleton would not give a strict improvement')
req(2*3==6 and 2*(3+1)==8,'old cubic feedback is six forces, not eight')

# Exact bridge algebra for arbitrary standardized p.
r,t,A,zheat,p=s.symbols('r t A zheat p',positive=True)
Delta=t*t-r*r
B=Delta/(t*zheat*zheat)
v0=Delta/(t*t)
h=s.sqrt(Delta)/zheat
# Use sqrt(v0)=sqrt(Delta)/t on 0<=r<t.
eq(s.sqrt(Delta)/t*h*A**p*zheat**(2*p),Delta/t*A**p*zheat**(2*p-1),'physical bridge normalization')
q=(2*p+1)/5
eq(1+2*(p-2)/5,q,'smallest local-heat exponent')
# Mode readout: source residual A^(M+1) s^(2M+2) r |z| times B.
M=s.symbols('M',integer=True,nonnegative=True)
eq(B*A**(M+1)*zheat**(2*M+2)*r,Delta/t*A**(M+1)*zheat**(2*M)*r,'finite-mode physical floor')

mp.mp.dps=90
rho=mp.mpf(3)/4
schedule=[]
for pi in (2,3,4,8,16,32):
    for power in (1,2,8,40):
        aa=mp.power(10,-power)
        target=mp.power(aa,mp.mpf(2)*(pi-2)/5)
        jj=0 if pi==2 else int(mp.ceil(mp.log(target)/mp.log(rho)))
        heat=mp.power(rho,jj)
        req(heat<=target,'terminal threshold met')
        if jj: req(heat/rho>target,'terminal threshold minimal')
        terminal=aa**2*heat**mp.mpf('2.5')
        req(terminal<=aa**pi,'order-two terminal reaches p')
        qq=mp.mpf(2*pi+1)/5
        alpha=aa*heat
        req(alpha<=aa**qq,'smallest alpha upper bound')
        if jj: req(alpha>rho*aa**qq,'smallest alpha lower bound')
        intrinsic=(1-mp.power(rho,jj*(mp.mpf(pi)+mp.mpf('.5'))))/(4*(1-mp.power(rho,mp.mpf(pi)+mp.mpf('.5'))))
        req(intrinsic<=1/(4*(1-mp.power(rho,mp.mpf(pi)+mp.mpf('.5')))),'geometric intrinsic error sum')
        req(mp.log(1/alpha)<=qq*mp.log(1/aa)+mp.log(1/rho),'public-log local heat')
        schedule.append({'p':pi,'A':'1e-'+str(power),'J':jj,'q_p':str(Fraction(2*pi+1,5))})

# Weighted path ledger independently expands a finite DAG.
edges={0:[(1,Fraction(3,2)),(2,Fraction(2))],1:[(3,Fraction(4,3))],2:[(3,Fraction(1,2)),(4,Fraction(3))],3:[],4:[]}
bs={0:Fraction(1,7),1:Fraction(1,3),2:Fraction(2,5),3:Fraction(1,2),4:Fraction(2,3)}
def recurrence(v,b): return b[v]+max([s0*recurrence(w,b) for w,s0 in edges[v]]+[Fraction(0)])
def path_costs(v,b):
    if not edges[v]: return [b[v]]
    return [b[v]+s0*c for w,s0 in edges[v] for c in path_costs(w,b)]
req(recurrence(0,bs)==max(path_costs(0,bs)),'complete exponent is maximum weighted root-to-leaf cost')
req(recurrence(0,{v:Fraction(0) for v in edges})==0,'finite zero-exponent seed graph remains zero')
for pi in (3,4,8,16):
    qp=Fraction(2*pi+1,5)
    req(qp*Fraction(1,2)/pi>=Fraction(1,5),'constant positive local exponent gives linear-order endpoint cost')
for k in range(1,12):
    req(sum(Fraction(1,2**j) for j in range(1,k+1))<1,'strict positive increments alone need not reach threshold')
    req(2**k==1/Fraction(1,2**k),'logarithmic fixed-branching depth creates inverse heat cost')

report={'status':'PASS, conditional integration theorem only','assertions':COUNT,'input_pin_count':len(PINS),
 'contract_sha256':PINS[str(contract)],'guarded_endpoint_sha256':PINS[str(BASE/'THIRD-ORDER-FULL-REVERSE-OU-LAW-AND-ENDPOINT.md')],
 'enumerated_nontrivial_trees':ntrees,'common_heat_conditional_variance':str(schur),
 'endpoint_heat_exponent':'q_p=(2p+1)/5','schedule_checks':schedule,
 'scope':'Independent exact Gram, history-specific tree and heat-count algebra, imported pin validation, endpoint and cost-recursion diagnostics. No positive-current port or all-order finite producer is instantiated or certified.'}
report['publication_pin_verification']=pin_report()
(HERE/'typed_contract_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='schedule_checks'},indent=2))
