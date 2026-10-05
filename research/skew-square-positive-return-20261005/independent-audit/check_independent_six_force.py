"""Independent finite graph, exact side-centering and anisotropic tensor checks.
These diagnostics do not execute the imported native selected-pair compilers.
"""
from itertools import combinations, permutations, product
from collections import deque, Counter
from pathlib import Path
import json
import sympy as s

out={}; base={(0,1),(1,2),(3,4),(4,5)}; verts=set(range(6))
def shortest(edges,a,b):
    q=deque([[a]])
    while q:
        p=q.popleft()
        if p[-1]==b:return p
        for edge in edges:
            if p[-1] in edge:
                nxt=edge[0] if edge[1]==p[-1] else edge[1]
                if nxt not in p:q.append(p+[nxt])
    raise AssertionError('Disconnected')
rows=[]
for r in (1,2):
    for left in combinations(range(3),r):
        for right in combinations(range(3,6),r):
            for rr in permutations(right):
                cross=set(zip(left,rr)); marks=verts-set(left)-set(right)
                candidates=[]
                for cut in ([None] if r==1 else cross):
                    edges=base|cross
                    if cut is not None:edges=edges-{cut}
                    for a,b in combinations(sorted(marks),2):
                        if a//3==b//3:continue
                        p=shortest(edges,a,b)
                        candidates.append((len(p),p,cut,edges))
                n,p,cut,edges=max(candidates,key=lambda z:(z[0],z[1]))
                side=verts-set(p)
                assert n>=4 and {1,4}.issubset(p)
                assert len(side)==6-n and side.isdisjoint({1,4})
                assert all(sum(v in e for e in edges)==1 for v in side)
                for copy in (set(range(3)),set(range(3,6))):
                    assert marks&copy
                    if cut is not None:assert len(set(cut)&copy)==1
                for v in verts:
                    # Original potential derivative order equals incident source
                    # edges + free physical marks + exposed auxiliary marks.
                    slots=sum(v in e for e in edges)+(v in marks)+(cut is not None and v in cut)
                    assert slots==(3 if v in (1,4) else 2)
                rows.append((r,n,6-n))
assert len(rows)==27
out['graph_count']=len(rows)
out['spine_histograms']={str(r):dict(Counter(n for rr,n,_ in rows if rr==r)) for r in (1,2)}
out['minimum_own_conditional_force']=min(2*n for _,n,_ in rows)

# A fixed contraction of two six-term physical permutation averages induces
# these topology multiplicities; the census must not be summed a second time.
Mcount=Counter();Ncount=Counter()
for lp in permutations(range(3)):
    for rp0 in permutations(range(3)):
        rp=tuple(x+3 for x in rp0)
        Mcount[(lp[1],rp[1])]+=1
        Ncount[tuple(sorted(((lp[1],rp[1]),(lp[2],rp[2]))))]+=1
assert len(Mcount)==9 and set(Mcount.values())=={4}
assert len(Ncount)==18 and set(Ncount.values())=={2}
s0=s.symbols('s0',positive=True)
normalization=s.simplify((s.Rational(16,36))/((s0/s.sqrt(2))**4*(s0/s.sqrt(3))**2))
assert normalization==16/(3*s0**6)
out['physical_permutation_pairs']=36
out['induced_topology_multiplicities']={'M':4,'N':2}
out['root_scalar_d']=str(normalization)

# Exact two-side centered innovation with arbitrary noncommuting matrices.
z=s.symbols('z0:4'); z1=s.Matrix(z[:2]);z2=s.Matrix(z[2:])
L=s.Matrix([[s.Rational(1,2),s.Rational(1,5)],[0,s.Rational(2,3)]])
R=s.Matrix([[s.Rational(1,3),0],[s.Rational(1,7),s.Rational(1,2)]])
B=s.Matrix([[2,-1],[3,4]])
m=s.Matrix([s.Rational(2,7),s.Rational(-3,5)])
n=s.Matrix([s.Rational(-1,3),s.Rational(4,9)])
def normal(poly,variables,variances=None):
    if variances is None:variances=[1]*len(variables)
    ans=0
    for powers,c in s.Poly(s.expand(poly),*variables).terms():
        value=c
        for p,v in zip(powers,variances):
            if p%2:value=0;break
            if p:value*=s.factorial2(p-1)*v**(p//2)
        ans+=value
    return s.simplify(ans)
F=((m+L*z1).T*B*(n+R*z2))[0]
E=s.expand(F-(m.T*B*n)[0])
assert normal(E,z)==0
CL=L*L.T;CR=R*R.T
var=s.trace(B*CR*B.T*CL)+(m.T*B*CR*B.T*m)[0]+(n.T*B.T*CL*B*n)[0]
assert normal(E**2,z)==var
out['exact_side_innovation_mean']=0
out['noncommuting_side_variance_identity']=True

# Genuinely coupled 2D skew tensor, anisotropic active carrier, independent buffer.
x=s.symbols('x0:2'); reserve=s.symbols('r0:2'); allvars=x+reserve
C=[s.Rational(2),s.Rational(3)]; q=[1/v for v in C]
varsigma=C+[s.Rational(1,2),s.Rational(4,3)]
K={idx:[s.Rational(1,7),s.Rational(-2,11),s.Rational(3,13),s.Rational(1,5)][sum(idx)] for idx in product(range(2),repeat=3)}
y=[q[i]*x[i] for i in range(2)]
def H2(a,b):return y[a]*y[b]-(q[a] if a==b else 0)
def H3(a,b,c):return y[a]*y[b]*y[c]-(q[a]*y[c] if a==b else 0)-(q[a]*y[b] if a==c else 0)-(q[b]*y[a] if b==c else 0)
p=[sum(K[i,a,b]*H2(a,b) for a,b in product(range(2),repeat=2)) for i in range(2)]
N={(i,j):sum(K[i,a,b]*q[a]*q[b]*K[j,a,b] for a,b in product(range(2),repeat=2)) for i,j in product(range(2),repeat=2)}
M={(i,j,b,d):sum(K[i,a,b]*q[a]*K[j,a,d] for a in range(2)) for i,j,b,d in product(range(2),repeat=4)}
correction=[-sum(N[i,j]*y[j] for j in range(2))-2*sum(M[i,j,b,d]*H3(j,b,d) for j,b,d in product(range(2),repeat=3)) for i in range(2)]
X=[x[i]+reserve[i] for i in range(2)]
def coeff_moment(indices,degree):
    ans=0
    for levels in product(range(3),repeat=len(indices)):
        if sum(levels)!=degree:continue
        term=1
        for i,l in zip(indices,levels):term*=([X,p,correction][l])[i]
        ans+=normal(term,allvars,varsigma)
    return s.simplify(ans)
for ij in product(range(2),repeat=2):assert coeff_moment(ij,2)==0
for ijk in product(range(2),repeat=3):assert coeff_moment(ijk,1)==6*K[ijk]
for ijkl in product(range(2),repeat=4):assert coeff_moment(ijkl,2)==0
out['anisotropic_coupled_2D_skew_square_cancellation']=True
out['third_cumulant_first_coefficient']='6 K'
out['scope']='Finite algebra/graph diagnostics; conditional and native analytic bounds require the written certificates.'
out['status']='PASS'
path=Path(__file__).with_name('independent_checks.json');path.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
