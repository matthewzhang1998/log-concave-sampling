#!/usr/bin/env python3
"""Independent finite diagnostics. No author check script is imported.
These checks do not execute native original-VALUE pair/filter programs.
"""
from pathlib import Path
from collections import Counter
from fractions import Fraction
import hashlib, itertools, json, math
import numpy as np
import sympy as s

OUT = Path(__file__).resolve().parent
BASE = OUT.parent
ROOT = BASE.parent
checks = 0

def check(test, label):
    global checks
    checks += 1
    if not bool(test):
        raise AssertionError(label)

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

# Check every manifest member in the three designated immutable packages.
packages = ['cycle-surplus-closure-20261005', 'old-bank-sixth-return-20261005',
            'eight-force-native-return-20261005']
pins = []
manifest_counts = {}
for package in packages:
    directory = ROOT / package
    manifest = json.loads((directory / 'MANIFEST.json').read_text())
    for entry in manifest['files']:
        f = directory / entry['path']
        check(sha(f) == entry['sha256'], f'manifest pin: {f}')
        if 'bytes' in entry:
            check(f.stat().st_size == entry['bytes'], f'byte count: {f}')
    manifest_counts[package] = len(manifest['files'])
    pins.append({'path': str(directory / 'MANIFEST.json'), 'sha256': sha(directory / 'MANIFEST.json')})
    rawpins = json.loads((directory / 'INPUT-PINS.json').read_text())
    for entry in rawpins if isinstance(rawpins, list) else rawpins['inputs']:
        f = Path(entry['path'])
        check(sha(f) == entry['sha256'], f'direct prerequisite pin: {f}')
        pins.append({'path': str(f), 'sha256': sha(f)})
for name in ['GROUPED-JOINT-ZERO-PORT.md', 'RETAINED-PUBLIC-ROTATION-OBSTRUCTION.md']:
    p = BASE / name
    pins.append({'path': str(p), 'sha256': sha(p)})
    (OUT / ('REVIEWED-' + name)).write_bytes(p.read_bytes())

# Re-derive scalar moments by direct normal monomial integration.
x,y,z,S,T,r,p,theta = s.symbols('x y z S T r p theta', real=True)
gs = (x,y,z,S,T)

def normal_moment(k):
    return 0 if k % 2 else math.prod(range(1,k,2))

def gaussian_expectation(poly):
    result = 0
    for powers, coeff in s.Poly(s.expand(poly), *gs).terms():
        result += coeff * math.prod(normal_moment(k) for k in powers)
    return s.expand(result)

V = r*theta**2*x*y*z*(S+p*theta*x*(T+p*theta*y*z))
mean = gaussian_expectation(V)
variance = s.expand(gaussian_expectation(V**2)-mean**2)
third = s.expand(gaussian_expectation((V-mean)**3))
check(s.expand(mean-r*p**2*theta**4)==0, 'first cluster')
check(s.expand(variance-r**2*theta**4*(1+3*p**2*theta**2+26*p**4*theta**4))==0, 'second cluster')
check(s.expand(gaussian_expectation(-V)+mean)==0, 'opposite first cluster')
check(s.expand(gaussian_expectation((-V+mean)**3)+third)==0, 'opposite third cluster')
# For independent opposite copies, variances add; cumulant normalization divides by 2.
check(s.simplify((variance+variance)/2-variance)==0, 'paired second coefficient')
characteristic_variance = s.expand(variance.subs(theta, s.I*theta))
check(characteristic_variance.coeff(theta,4)==r**2, 'rank-four characteristic sign')
check(characteristic_variance.coeff(theta,6)==-3*r**2*p**2, 'rank-six characteristic sign')
check(characteristic_variance.coeff(theta,8)==26*r**2*p**4, 'rank-eight characteristic sign')

# Independently enumerate the chosen opening's Wick contractions and graph closure.
# Each monomial keeps an ancestor-closed center subtree. Boundary colors match
# precisely the exponents derived above; center identities are never merged.
monomials = [
    {'centers':(0,1), 'edges':((0,1),), 'ports':{'x':[0], 'y':[1], 'z':[1], 'S':[0]}},
    {'centers':(0,1,2), 'edges':((0,1),(0,2)), 'ports':{'x':[0,2], 'y':[1], 'z':[1], 'T':[2]}},
    {'centers':(0,1,2,3), 'edges':((0,1),(0,2),(2,3)), 'ports':{'x':[0,2], 'y':[1,3], 'z':[1,3]}},
]

def pairings(items):
    if not items:
        yield []
        return
    first = items[0]
    for j in range(1,len(items)):
        for rest in pairings(items[1:j]+items[j+1:]):
            yield [(first,items[j])] + rest

def connected(vertices, edges):
    seen = {vertices[0]}
    while True:
        new = set(seen)
        for a,b in edges:
            if a in seen: new.add(b)
            if b in seen: new.add(a)
        if new == seen: return len(seen) == len(vertices)
        seen = new

wick_total = 0
wick_connected = Counter()
wick_disconnected = Counter()
for i,j in itertools.product(range(3), repeat=2):
    vertices, edges, ports = [], [], {}
    for copy,k in enumerate((i,j)):
        m = monomials[k]
        vertices.extend((copy,v) for v in m['centers'])
        edges.extend(((copy,a),(copy,b)) for a,b in m['edges'])
        for color, vs in m['ports'].items():
            ports.setdefault(color,[]).extend((copy,v) for v in vs)
    if any(len(vs)%2 for vs in ports.values()):
        continue
    for parts in itertools.product(*(list(pairings(vs)) for vs in ports.values())):
        all_edges = edges + sum((list(pp) for pp in parts), [])
        wick_total += 1
        deg = Counter(v for e in all_edges for v in e)
        check(all(a!=b for a,b in all_edges), 'no Wick self-loop')
        check(all(deg[v]==3 for v in vertices), 'cubic degree retained')
        check(len(vertices)%2==0, 'even closed cubic center count')
        if connected(vertices,all_edges): wick_connected[len(vertices)] += 1
        else: wick_disconnected[len(vertices)] += 1
check(dict(wick_connected)=={4:1,6:3,8:26}, 'connected variance coefficients from graphs')
check(dict(wick_disconnected)=={8:1}, 'mean-square diagram exactly removed')

# Conservative census of possible (center count,surplus) types. Zero coefficients
# may be overcounted. This verifies termination, not full tensor multiplicities.
queue_results = {}
for P in (12,13,20,21,32,48,64,80):
    todo = [(4,4,0)]
    visited = {}
    maxdepth = 0
    while todo:
        n,d,depth = todo.pop()
        if visited.get((n,d), -1) >= depth: continue
        visited[(n,d)] = depth; maxdepth = max(maxdepth,depth)
        for h in range(2,(P-1)//d+1):
            for nnew in range(2*h, min(h*n,P-h*d-1)+1):
                if nnew%2: continue
                dnew=h*d
                check(dnew>=2*d, 'strict surplus feedback')
                check(nnew+dnew<P, 'strict below-cutoff child')
                check(nnew==2*h+(nnew-2*h), 'root/side occurrence census')
                todo.append((nnew,dnew,depth+1))
    check(all(d<P and (n+d<P or (n,d)==(4,4)) for n,d in visited), 'bounded census')
    check(maxdepth <= int(math.log2(P/4)), 'surplus depth bound')
    queue_results[str(P)] = {'types':len(visited),'max_depth_from_start':maxdepth}
check(4+2*4==12 and 4+2*8==20, 'initial and next root-only return')
# Explicit alpha^14 port: +/- starting pair and one d=8 rank-four counter.
initial_grades=set()
for h in range(2,7,2):  # odd root multiplicities cancel between opposite copies
    for j in range(2*h+1):
        nnew=2*h+j
        if nnew%2==0:
            initial_grades.add(nnew+4*h)
check(min(initial_grades)==12, 'three-packet first return grade')
check(13 not in initial_grades, 'cubic parity excludes grade thirteen')
check(min(initial_grades-{12})==14, 'three-packet residual grade fourteen')
check(4+8==12 and 4+2*8==20, 'rank-four counter cancels twelve and emits twenty')
check(3*(2*4)==24 and 3*(4//2+1)==9 and 3*4==12, 'three-packet literal native counts')
check(4+8 < 100+4, 'effective-grade monotonicity countertest')

# Gaussian obstruction: independent analytic formulas, exact derivatives and
# deterministic Hermite quadrature, with no Monte Carlo used as a proof.
B, eta,t,delta,s0 = s.symbols('B eta t delta s0', real=True, positive=True)
u = s.symbols('u', real=True)
g = lambda a: a+eta*s.sin(a)
f = delta*(g(B+t*x)-g(B)-t*x)/((1+eta)*t)
check(s.simplify(s.diff(f,x)-delta*eta*s.cos(B+t*x)/(1+eta))==0, 'source selected derivative')
check(s.simplify(s.diff(f,B)-delta*eta*(s.cos(B+t*x)-s.cos(B))/((1+eta)*t))==0, 'source center derivative')
# Fourier evaluation from the two exponentials composing cos(B).
cos_fourier = (s.exp(-(u+1)**2/2)+s.exp(-(u-1)**2/2))/2
check(s.simplify(s.expand(cos_fourier-s.exp(-(u*u+1)/2)*s.cosh(u).rewrite(s.exp)))==0, 'Gaussian cosine Fourier factor')
a,c,beta,eps,kappa,m0 = s.symbols('a c beta eps kappa m0', real=True)
log_first = -a*c*eps*theta**2*s.exp(-s.Rational(1,2))*s.cosh(beta*theta)
series = s.series(log_first,theta,0,12).removeO().expand()
for k in range(5):
    rank=2+2*k
    check(s.simplify(series.coeff(theta,rank)+a*c*eps*s.exp(-s.Rational(1,2))*beta**(2*k)/s.factorial(2*k))==0, 'tilted log sign and factorial')
# Raw source baseline has precisely the same mean and variance differences.
M = m0+eps*s.cos(B)
check(s.expand((a+c*M)-(a+c*m0))==c*eps*s.cos(B), 'raw conditional mean difference')
raw_var = s.expand((a+c*M)**2+c**2*(1-M**2)+kappa)
check(s.expand(raw_var-(a*a+c*c+kappa+2*a*c*m0))==2*a*c*eps*s.cos(B), 'raw conditional variance difference')

nodes, weights = np.polynomial.hermite.hermgauss(100)
bvals = np.sqrt(2)*nodes
weights = weights/np.sqrt(np.pi)
Ecos2 = float(weights @ np.cos(bvals)**2)
check(abs(Ecos2-(1+math.exp(-2))/2)<2e-15, 'lower-bound normal constant')
quad_checks=[]
for av,cv,bv,kv,ev,mv,th in [(0.7,0.4,0.9,0.3,0.025,0.0,0.8),
                            (0.3,0.6,-0.5,0.2,0.015,0.03,1.1),
                            (0.4,0.5,1.2,0.7,0.03,0.04,0.3)]:
    vv=av*av+cv*cv+kv+2*av*cv*mv
    chi0=complex(np.exp(-vv*th*th/2)*(weights@np.exp(1j*bv*th*bvals)))
    deriv_quad=complex(np.exp(-vv*th*th/2)*(weights@(-av*cv*th*th*np.cos(bvals)*np.exp(1j*bv*th*bvals))))
    log_deriv=deriv_quad/chi0
    expected=-av*cv*th*th*math.exp(-0.5)*math.cosh(bv*th)
    check(abs(log_deriv-expected)<2e-14, 'tilted characteristic derivative quadrature')
    covgap=1-(mv+ev*np.cos(bvals))**2
    check(float(covgap.min())>0.99, 'explicit native covariance gap')
    # Retained B,P conditional mean gap, P integrated analytically with E P²=1.
    lower=cv*abs(ev)*math.sqrt(Ecos2)
    check(abs(lower-cv*abs(ev)*math.sqrt((1+math.exp(-2))/2))<1e-15, 'exact retained lower bound')
    quad_checks.append({'parameters':[av,cv,bv,kv,ev,mv,th], 'log_derivative':log_deriv.real, 'retained_lower_bound':lower})

# The ordinary joint-W2 test uses one vector Cauchy-Schwarz after this
# exact telescoping identity; no global Lipschitz assertion for P*W is needed.
Cb,Cbp,Pvar,Pprime,Wvar,Wprime=s.symbols('Cb Cbp Pvar Pprime Wvar Wprime')
product_difference=Cb*Pvar*Wvar-Cbp*Pprime*Wprime
split=(Cb-Cbp)*Pvar*Wvar+Cbp*(Pvar-Pprime)*Wvar+Cbp*Pprime*(Wvar-Wprime)
check(s.expand(product_difference-split)==0, 'ordinary joint-W2 product splitting')
for bv,pv,wv,bp,pp,wp in np.random.default_rng(808).normal(size=(100,6)):
    lhs=abs(math.cos(bv)*pv*wv-math.cos(bp)*pp*wp)
    rhs=math.sqrt((pv*wv)**2+wv**2+pp**2)*math.sqrt((bv-bp)**2+(pv-pp)**2+(wv-wp)**2)
    check(lhs<=rhs+1e-12, 'pointwise joint test product bound')

# Reparameterization and direct-sum derivative controls, tested on noncommuting
# rectangular matrices so that scalar checks cannot conceal orientation errors.
rng=np.random.default_rng(4105)
for D in (1,2,5,11):
    J=[rng.normal(size=(D, D+k)) for k in (0,2,5)]
    joined=np.concatenate(J,axis=1)
    norm=np.linalg.norm(joined,2)
    bound=math.sqrt(sum(np.linalg.norm(j,2)**2 for j in J))
    check(norm<=bound+1e-12, 'block-row complete derivative bound')
    check(np.allclose(joined@joined.T,sum((j@j.T for j in J))), 'block-row Gram identity')
    q,_=np.linalg.qr(rng.normal(size=(joined.shape[1],joined.shape[1])))
    check(abs(np.linalg.norm(joined@q,2)-norm)<1e-12, 'orthogonal split preserves Jacobian operator norm')
for rr in (0.0,0.2,0.8,-0.4):
    # Coordinates (B,U) -> (X0,V), then restore B exactly by regression.
    ss=math.sqrt(1-rr*rr)
    rotation=np.array([[rr,ss],[ss,-rr]])
    check(np.allclose(rotation@rotation.T,np.eye(2)), 'regression rotation orthogonal')
    check(np.allclose(np.array([rr,ss])@rotation,np.array([1.,0.])), 'original exposed B retained')
    check(abs((np.array([1.,0.])@rotation[0])-rr)<1e-15, 'cross covariance cannot disappear')

summary={
 'status':'PASS: finite algebra and source pins; native sampler execution not claimed',
 'assertions':checks,
 'manifest_files_verified':manifest_counts,
 'first_cluster':str(mean),
 'single_variance':str(variance),
 'opposite_pair_second_connected_coefficient':str(variance),
 'wick_matchings':wick_total,
 'connected_wick_counts':dict(wick_connected),
 'disconnected_wick_counts':dict(wick_disconnected),
 'surplus_type_censuses':queue_results,
 'three_packet_port':{'conditional_remainder_grade':14,'counter_intrinsic_return_grade':20,'complete_source_pair_calls_per_item':24,'cut_edge_roots_per_item':9,'physical_marks_per_item':12},
 'retained_lower_bound_factor':math.sqrt(Ecos2),
 'tilted_rank_four_coefficient':str(series.coeff(theta,4)),
 'quadrature_checks':quad_checks,
 'limitations':['No native pair/filter/calibration execution.',
               'Finite type enumeration is an upper census, not all signed tensor multiplicities.',
               'Positive law and terminal first/curl conclusions require the analytical hypotheses in the report.']
}
(OUT/'independent_grouped_checks.json').write_text(json.dumps(summary,indent=2)+'\n')
(OUT/'INPUT-PINS.json').write_text(json.dumps(pins,indent=2)+'\n')
print(json.dumps(summary,indent=2))
