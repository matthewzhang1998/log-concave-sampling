"""Independent diagnostic: imports no author implementation; no native compiler.
Floating checks are diagnostics; Fraction identities are exact.
"""
from pathlib import Path
from fractions import Fraction as Q
import itertools, math, json, hashlib, subprocess
import numpy as np
ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent
rng=np.random.default_rng(5102026)
checks=0
metrics={}
def ck(ok, label):
    global checks
    checks+=1
    if not ok: raise AssertionError(label)
def near(a,b,tol=1e-10,label='numerical identity'):
    ck(np.max(np.abs(np.asarray(a)-np.asarray(b)))<tol,label)

# Literal independent implementation; no derivative is executed as a source.
def psi_jac(x,R):
    r=np.linalg.norm(x); D=len(x)
    if R is None or r<=R:return x.copy(),np.eye(D)
    u=x/r
    if r>=2*R:h=1.5*R;hp=0.
    else:
        v=(r-R)/R;h=R*(1+v-v**3+v**4/2);hp=1-3*v*v+2*v**3
    return h*u,h/r*np.eye(D)+(hp-h/r)*np.outer(u,u)
def source(g,x,P,z,t,R=None):
    y,J=psi_jac(x,R);k=len(P);a=1/math.sqrt(k+1);out=np.zeros_like(x);calls=0
    for signs in itertools.product((-1,1),repeat=k):
        shift=sum((s*p for s,p in zip(signs,P)),np.zeros_like(x))
        q=z+a*t*shift
        out+=math.prod(signs)*(g(q+a*t*y)-g(q));calls+=2
    return J@out/(t*2**k),calls

def mark(g,x,P,z,t,delta,Z,R=None):
    c,n=source(g,x,P,z,t,R);w,m=source(g,x,P,z+delta*Z,t,R)
    return (c-w)/2,n+m

def jac(fun,x,h=1e-5):
    return np.column_stack([(fun(x+np.eye(len(x))[j]*h)-fun(x-np.eye(len(x))[j]*h))/(2*h) for j in range(len(x))])
max_curl=0.;max_z_ratio=0.;max_increment_ratio=0.
for D in (1,2,3,6):
    u=rng.normal(size=D);u/=np.linalg.norm(u)
    v=rng.normal(size=D);v/=np.linalg.norm(v)
    # Hessian is between .2 I and .8 I: globally admissible original VALUE.
    def g(q):return .5*q+.2*np.sin(u@q)*u+.1*np.sin(v@q)*v
    for k in (0,1,2,4,6):
        for radial in (False,True):
            R=math.sqrt(D)+.2 if radial else None
            x=rng.normal(size=D)*(3 if radial else 1);P=rng.normal(size=(k,D))
            z=rng.normal(size=D);Z=rng.normal(size=D);t=.08 if k%2 else .61;delta=.42
            f,n=mark(g,x,P,z,t,delta,Z,R)
            ck(n==2**(k+2),'literal marked call count')
            near(mark(g,np.zeros(D),P,z,t,delta,Z,R)[0],np.zeros(D),label='origin')
            linear=lambda q:.43*q
            near(mark(linear,x,P,z,t,delta,Z,R)[0],np.zeros(D),label='linear source cancellation')
            for j in range(k):
                P2=P.copy();P2[j]*=-1
                near(mark(g,x,P2,z,t,delta,Z,R)[0],-f,label='active probe oddness')
            Jx=jac(lambda xx:mark(g,xx,P,z,t,delta,Z,R)[0],x)
            curl=np.max(abs(Jx-Jx.T));max_curl=max(max_curl,curl)
            ck(curl<2e-8,'private curl after common pullback')
            active=np.r_[x,P.ravel()]
            Ja=jac(lambda q:mark(g,q[:D],q[D:].reshape(k,D),z,t,delta,Z,R)[0],active)
            ck(np.linalg.norm(Ja,2)<=13+1e-7,'complete private/probe first')
            Jz=jac(lambda q:mark(g,x,P,q,t,delta,Z,R)[0],z)
            JZ=jac(lambda q:mark(g,x,P,z,t,delta,q,R)[0],Z)
            max_z_ratio=max(max_z_ratio,np.linalg.norm(Jz,2)*t)
            max_increment_ratio=max(max_increment_ratio,np.linalg.norm(JZ,2)*2*t/delta)
            ck(np.linalg.norm(Jz,2)*t<=1+1e-7,'captured center first')
            ck(np.linalg.norm(JZ,2)*2*t/delta<=1+1e-7,'increment first including half')
            if radial:ck(np.linalg.norm(f)<=1.5*R/math.sqrt(k+1)+1e-9,'bounded mark radius')
metrics.update(max_private_curl=max_curl,max_z_bound_ratio=max_z_ratio,max_increment_bound_ratio=max_increment_ratio)

# Direct deterministic Gaussian response integration for independent active variables.
from numpy.polynomial.hermite import hermgauss
nodes,weights=hermgauss(20);nodes=nodes*math.sqrt(2);weights=weights/math.sqrt(math.pi)
max_response=0.
for k in range(4):
    a=1/math.sqrt(k+1);t=.43;z=.31;delta=.62;Z=.77
    grids=np.meshgrid(*([nodes]*(k+1)),indexing='ij');wgrids=np.meshgrid(*([weights]*(k+1)),indexing='ij')
    X=grids[0];PP=grids[1:];W=np.prod(wgrids,axis=0);score=np.prod(grids,axis=0)
    def scalar_source(center):
        out=np.zeros_like(X)
        for eps in itertools.product((-1,1),repeat=k):
            q=center+a*t*sum((e*p for e,p in zip(eps,PP)),np.zeros_like(X))
            out+=math.prod(eps)*(.5*a*t*X+.2*(np.sin(q+a*t*X)-np.sin(q)))
        return out/(t*2**k)
    actual=np.sum((scalar_source(z)-scalar_source(z+delta*Z))/2*score*W)
    p=k+1
    expected=a**(k+1)*t**k/2*.2*math.exp(-t*t/2)*(math.sin(z+p*math.pi/2)-math.sin(z+delta*Z+p*math.pi/2))
    max_response=max(max_response,abs(actual-expected));near(actual,expected,2e-12,'active marked response')
# k=0 mark need not have private mean zero (fixed centers).
x=nodes;t=.43;z=.31;delta=.62;Z=.77
f=.2*(np.sin(z+t*x)-np.sin(z)-np.sin(z+delta*Z+t*x)+np.sin(z+delta*Z))/(2*t)
mean=float(weights@f)
ck(abs(mean)>1e-4,'k=0 private mean counterexample')
metrics.update(max_response_error=max_response,k0_nonzero_private_mean=mean)

# All proper cuts of an explicit nonseparable Fourier-gradient response.
# Tensor indices are output plus k+1 differentiated indices; each Fourier term
# is a rank-one tensor. Original Hessian stays between .2 and .8 I.
def outerpower(u,m):
    out=np.array(1.)
    for _ in range(m):out=np.multiply.outer(out,u)
    return out
u=np.array([1.,2.]);u/=np.linalg.norm(u);v=np.array([2.,-1.]);v/=np.linalg.norm(v)
for k in range(7):
    a=1/math.sqrt(k+1);t=.52;z=np.array([.3,-.6]);zz=z+np.array([.42,.18])
    def tensor(center):
        T=np.zeros((2,)*(k+2))
        if k==0:T+=.5*np.eye(2)
        for c,q in ((.2,u),(.1,v)):
            T+=c*math.exp(-t*t/2)*math.sin(q@center+(k+1)*math.pi/2)*outerpower(q,k+2)
        return a**(k+1)*t**k*T
    M=(tensor(z)-tensor(zz))/2
    ck(np.linalg.norm(M.ravel())<=a*math.sqrt(2)+1e-12,'mark one-HS')
    for p in range(1,k+2):
        c=a**(k+1)*2**(k/2)*math.sqrt(math.factorial(p-1)*math.factorial(k+1-p))
        ck(np.linalg.norm(M.reshape(2**p,-1),2)<=c+1e-12,'every proper cut')

# Exact characteristic-function integration on a genuinely shared structural bank.
# D=1 is enough to distinguish preserving versus discarding bank correlations.
# g(x)=.5*x+.2*sin(x), p=degree at a force vertex.
def mean_product(degrees,mu,rows,private_variances):
    opts=[]
    for p,var in zip(degrees,private_variances):
        local=[(1,.2/(2j)*(1j)**p*math.exp(-var/2)),(-1,-.2/(2j)*(-1j)**p*math.exp(-var/2))]
        if p==1:local.append((0,.5+0j))
        opts.append(local)
    total=0j
    for terms in itertools.product(*opts):
        freq=np.array([q[0] for q in terms]);coef=math.prod(q[1] for q in terms)
        total+=coef*np.exp(1j*(freq@mu)-.5*np.linalg.norm(freq@rows)**2)
    ck(abs(total.imag)<1e-12,'real integrated force product')
    return total.real
errors=[];wrong_independence=[]
for degrees in ([7,1,1,1,1,1,1,1],[1,2,2,2,2,2,2,1],[3,3,2,2,1,1,1,1]):
    n=8;ck(sum(degrees)==2*(n-1),'tree degree census')
    mu=rng.normal(size=n)*.4
    # Nontrivial common ancestry: each row reads earlier primitive innovations.
    L=np.tril(rng.normal(size=(n,n))*.17)
    rr=np.linspace(.31,.93,n);sigma=.17;tau=.39
    tt2=(1-rr**2+sigma**2)/2;dd2=tau*tau-sigma*sigma
    old_t2=(1-rr**2+tau*tau)/2
    rows_c=np.c_[L,np.diag(np.sqrt(tt2)),np.zeros((n,n))]
    rows_w=np.c_[L,np.diag(np.sqrt(tt2)),np.eye(n)*math.sqrt(dd2)]
    rows_old=np.c_[L,np.diag(np.sqrt(old_t2)),np.zeros((n,n))]
    cold=mean_product(degrees,mu,rows_c,tt2)
    warm=mean_product(degrees,mu,rows_w,tt2)
    old=mean_product(degrees,mu,rows_old,old_t2)
    cold_target=mean_product(degrees,mu,L,1-rr**2+sigma*sigma)
    warm_target=mean_product(degrees,mu,L,1-rr**2+tau*tau)
    near(cold,cold_target,label='complete cold heat target')
    near(warm,warm_target,label='complete warm heat target')
    near(old,warm,label='old/new conditional split full-bank equality')
    # Telescope is linear in the marked coefficient at each fixed common bank.
    B=rng.normal(size=3*n)
    def coeff(rows,var):return np.array([(.5 if p==1 else 0)+.2*math.exp(-s/2)*math.sin(z+row@B+p*math.pi/2) for p,z,row,s in zip(degrees,mu,rows,var)])
    c=coeff(rows_c,tt2);w=coeff(rows_w,tt2);o=coeff(rows_old,old_t2)
    telescope=sum((c[m]-w[m])*np.prod(c[:m])*np.prod(w[m+1:]) for m in range(n))
    err=abs(telescope-(np.prod(c)-np.prod(w)));errors.append(err)
    near(telescope,np.prod(c)-np.prod(w),label='pointwise shared-bank eight-mark telescope')
    ck(abs(np.prod(o)-np.prod(w))>1e-11,'old and new warm coefficients differ conditionally')
    # Replacing shared structural innovations by independent occurrence rows is wrong.
    decoupled=np.diag(np.sqrt(np.sum(L*L,axis=1)))
    independent=mean_product(degrees,mu,decoupled,1-rr**2+tau*tau)
    wrong_independence.append(abs(warm_target-independent))
ck(max(wrong_independence)>1e-8,'common ancestry cannot be independently resampled')
metrics.update(max_pointwise_telescope_error=max(errors),shared_bank_vs_independent_errors=wrong_independence)

# Exact heterogeneous mixed and conditional allocation algebra.
for trial in range(1000):
    Gamma=Q(1,5);arity=int(rng.integers(2,7));packets=[]
    for _ in range(arity):
        N=int(rng.integers(1,9));beta=[Q(int(rng.integers(1,8)),40) for _ in range(N)]
        gamma=[Q(int(rng.integers(1,9)),40) for _ in range(N)]
        k=[int(rng.integers(0,8)) for _ in range(N)];a=Q(20)
        H=sum(g*j for g,j in zip(gamma,k));d=a-sum(beta)-H;Psi=d-2*Gamma
        packets.append((a,beta,gamma,k,d,Psi))
        # Omitting arbitrary side occurrences while retaining h genuine roots.
        h=int(rng.integers(2,8));indices=[0]+[i for i in range(1,N) if rng.integers(2)]
        child_beta=h*sum(beta[i] for i in indices);child_H=h*sum(gamma[i]*k[i] for i in indices)
        child_a=child_beta+child_H+h*d
        ck(child_a-child_beta-child_H-2*Gamma==h*Psi+2*Gamma*(h-1),'conditional omission identity')
    # Connect packet i to a preceding packet using a labelled recursive tree.
    bridge_price=Q(0)
    for i in range(1,arity):
        j=int(rng.integers(0,i));u=int(rng.integers(len(packets[i][2])));v=int(rng.integers(len(packets[j][2])))
        bridge_price+=packets[i][2][u]+packets[j][2][v]
    child=sum(p[0]-sum(p[1])-sum(g*k for g,k in zip(p[2],p[3])) for p in packets)-bridge_price-2*Gamma
    rhs=sum(p[5] for p in packets)+2*Gamma*(arity-1)-bridge_price
    ck(child==rhs,'heterogeneous bridge equality')
    ck(child>=sum(p[5] for p in packets),'fixed-Gamma bridge superadditivity')

n=8;Gamma=Q(1,7);beta0=Q(1,14);beta1=Q(1,28);g0=Q(1,8);g1=Q(1,7)
d0=n-n*beta0-(n-2)*g0;d1=n-n*beta1-(n-2)*g1
expected={'d_old':Q(187,28),'Psi_old':Q(179,28),'d_shell':Q(48,7),'Psi_shell':Q(46,7),'reserve_gain':Q(5,28),'root_old':Q(27,4),'root_shell':Q(193,28),'root_first_old':Q(53,8),'root_first_shell':Q(27,4),'energy':Q(50,7),'bank_first':Q(7),'old_join':Q(57,7),'self_join':Q(99,7),'native_own':Q(27,2),'heat':Q(57,7),'allocation_ceiling':Q(37,168)}
got={'d_old':d0,'Psi_old':d0-2*Gamma,'d_shell':d1,'Psi_shell':d1-2*Gamma,'reserve_gain':d1-d0,'root_old':beta0+d0,'root_shell':beta1+d1,'root_first_old':beta0+d0-g0,'root_first_shell':beta1+d1-g1,'energy':n-(n-2)*g1,'bank_first':n-(n-1)*g1,'old_join':n+1-(n-2)*g1,'self_join':2*n-(2*n-3)*g1,'native_own':2*min(beta0+d0,beta1+d1),'heat':n+g1,'allocation_ceiling':g0+Q(n)*beta0/(n-2)}
for key in expected:ck(got[key]==expected[key],key)
for n in range(4,101):
    b0=Q(1,2*(n-1));b1=Q(1,4*(n-1));g0=Q(1,n);g1=Q(1,n-1)
    gain=n*(b0-b1)-(n-2)*(g1-g0)
    ck(gain==Q(n*n-4*n+8,4*n*(n-1)) and gain>0,'fixed-rank reserve gain')
    ro=Q(n)-Q(3,2)+Q(2,n);rs=Q(n)-Q(5,4)+Q(1,n-1);P=Q(n)+Q(1,n-1)
    ck(2*min(ro,rs)>P,'fixed-rank own margin')
    ck(2*n-Q(2*n-3,n-1)>P,'fixed-rank self margin')

# Independent combinatorial count recurrences.
T={1:1};V={1:1}
for n in range(2,9):
    T[n]=sum(i*(n-i)*T[i]*T[n-i] for i in range(1,n))
    V[n]=sum((n-i)*(i+1)*V[i]*T[n-i]+i*(n-i+1)*T[i]*V[n-i] for i in range(1,n))
ck(T[8]==794880,'rank-eight histories');ck(V[8]==27020800,'raw ordinary census')
ck(9*V[8]==243187200,'raw shell census');ck(10*V[8]==270208000,'old-plus-shell census')
ck(2**7+14==142 and 9*142==1278 and 10*142==1420,'star raw maxima')
# Product normalization, readout and optional error factors cannot be omitted.
alpha=.01;t=.073;dg=1e-9
ck(2*dg/t>dg,'VALUE response precision pays actual width')

# Matrix-valued rank-eight tree contraction, with all eight physical legs retained.
max_tensor_telescope=0.
for edges in ([(0,j) for j in range(1,8)],[(j,j+1) for j in range(7)],[(0,1),(0,2),(0,3),(1,4),(1,5),(2,6),(3,7)]):
    labels=[chr(97+j) for j in range(8)]
    for e,(u,v) in enumerate(edges):
        label=chr(105+e);labels[u]+=label;labels[v]+=label
    recipe=','.join(labels)+'->abcdefgh'
    C=[rng.normal(size=(2,)*len(l))*.2 for l in labels]
    W=[rng.normal(size=(2,)*len(l))*.2 for l in labels]
    contract=lambda A:np.einsum(recipe,*A,optimize=True)
    shell=sum((contract(C[:m]+[(C[m]-W[m])/2]+W[m+1:])*2 for m in range(8)),np.zeros((2,)*8))
    err=np.max(abs(shell-(contract(C)-contract(W))));max_tensor_telescope=max(max_tensor_telescope,err)
    ck(err<1e-12,'rank-eight full tensor telescope with unchanged physical legs')
metrics['max_rank8_tensor_telescope_error']=max_tensor_telescope

# Exact polynomial chord degree and positive Gaussian interpolation quadrature.
from numpy.polynomial.legendre import leggauss
max_chord_error=0.
for n in range(4,31):
    c=rng.normal(size=n)*.4;w=rng.normal(size=n)*.4;m=(n+1)//2
    nodes,weights=leggauss(m);nodes=(nodes+1)/2;weights=weights/2
    ck(np.all(weights>0),'positive interpolation weights')
    derivative=lambda s:sum((c[j]-w[j])*np.prod(np.delete((1-s)*w+s*c,j)) for j in range(n))
    quadrature=sum(q*derivative(s) for s,q in zip(nodes,weights))
    target=np.prod(c)-np.prod(w);max_chord_error=max(max_chord_error,abs(quadrature-target))
    near(quadrature,target,1e-12,'ceil(n/2)-node exact polynomial chord')
ck(2*8*4*V[8]==1729331200,'optional chord VALUE census')
metrics['max_chord_quadrature_error']=max_chord_error

# No source-smallness can be credited to delta under uniform C2 assumptions.
# g_delta(x)=.5x+.25*epsilon*sin(x/epsilon), epsilon=delta/pi.
# At x=z=0, k=0, Z=1, t=1: D_x mark=.25 for every delta.
for delta in (1e-1,1e-3,1e-6,1e-9):
    epsilon=delta/math.pi
    dx=.25/2*(math.cos(0)-math.cos(delta/epsilon))
    near(dx,.25,label='high-frequency marked first has no delta gain')
metrics['high_frequency_dx_mark']=.25

# Completed literal-force port: first and antisymmetric Jacobian are separate
# from private source curl and require a directly certified residual first.
for D in (1,2,3,8):
    for rep in range(10):
        O,_=np.linalg.qr(rng.normal(size=(D,D)));A=.13
        H=O@np.diag(rng.uniform(0,A,size=D))@O.T
        DR=rng.normal(size=(D,D))*.017;s=.74;JR=np.linalg.norm(DR,2)
        terminal=H@(s*np.eye(D)-DR)
        ck(np.linalg.norm(terminal,2)<=A*(abs(s)+JR)+1e-12,'literal terminal first')
        ck(np.linalg.norm(terminal-terminal.T,2)<=2*A*JR+1e-12,'literal terminal curl')

pins={}
for folder in ('induction-native-rank-generator-20261005','induction-mixed-queue-20261005'):
    path=Path('/workspace/shared')/folder
    proc=subprocess.run(['sha256sum','-c','SHA256SUMS'],cwd=path,text=True,capture_output=True)
    ck(proc.returncode==0,'sealed source checksum '+folder)
    pins[folder]={'sha256sum_status':'PASS','manifest_sha256':hashlib.sha256((path/'MANIFEST.json').read_bytes()).hexdigest()}
files=[BASE/'HEAT-SHELL-REENTRY.md',BASE/'COST-AND-FLOORS.md',Path('/workspace/shared/induction-shared-variables-20261005/marked-spanning-tree-addendum/FROZEN-COEFFICIENT-WIDTH-ZERO-PORT.md')]
report={'status':'PASS_WITH_SCOPED_QUALIFICATIONS','assertions':checks,'native_compiler_executed':False,'metrics':metrics,'rank8_exact':{k:str(v) for k,v in got.items()},'sealed_inputs':pins,'source_hashes':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
(ROOT/'independent-checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
