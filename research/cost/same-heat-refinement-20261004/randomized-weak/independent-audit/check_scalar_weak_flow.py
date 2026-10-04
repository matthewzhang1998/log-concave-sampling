"""Independent deterministic fixtures, not a proof by numerical experiment."""
from pathlib import Path
from fractions import Fraction as F
import hashlib, json, math
import numpy as np
from scipy.integrate import solve_ivp, simpson
from scipy.special import ndtr, roots_hermitenorm, roots_legendre, gammaln

ROOT = Path(__file__).resolve().parent
REPORT = ROOT / 'SCALAR-WEAK-FLOW-REVIEWED-14ff69c491ab.md'
T = math.pi/2
CHECKS = 0

def check(test, message):
    global CHECKS
    CHECKS += 1
    assert bool(test), message

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def U(z,a): return a*(z*z/4 + .2*(1-np.cos(z)))
def b(z,a): return a*(z/2 + .2*np.sin(z))
def bp(z,a): return a*(.5+.2*np.cos(z))

check(digest(REPORT)=='14ff69c491ab5cbcc3dfcfa196156cefced1c685f4905e5df663a6f8e48f8f27','reviewed report snapshot hash')

# Exact harmonic perturbation, keep cancellation, density change of variables.
quadratic=[]
for a in (.2,.1,.03,.01,.003):
    for lam in (.25,.5,1.):
        w=math.sqrt(1+a*lam); c=math.cos(T*w); s=math.sin(T*w)/w
        target=1/w**2
        check(abs(c*c*target+s*s-target)<1e-14,'invariant quadratic variance')
        for sigma in (.5,.25,.1):
            delta=sigma*sigma; k=math.sqrt(1-delta)
            vmu=c*c*target+k*k*s*s
            vy=vmu+delta
            check(abs(vy-target-delta*(1-s*s))<1e-14,'exact keep excess variance')
            ratio=(math.sqrt(vy)-math.sqrt(target))/(a*delta)
            check(0<ratio<1,'quadratic weak keep coefficient')
            z=np.linspace(-4,4,19)
            den=k*k+delta*c*c
            rn=den**(-.5)*np.exp(-delta*(w*math.sin(T*w)*z)**2/(2*den))
            direct=math.sqrt(target/vmu)*np.exp(z*z*(1/target-1/vmu)/2)
            check(np.max(np.abs(rn-direct))<1e-12,'inverse phase density D1')
            # Same-record endpoint difference is (s-1) sigma K, not sigma K.
            z0,Z,K=.7,-.8,1.3
            y=c*z0+s*k*Z+sigma*K
            exact=c*z0+s*(k*Z+sigma*K)
            check(abs((y-exact)-(1-s)*sigma*K)<1e-14,'same-record keep identity')
            quadratic.append(dict(a=a,lambda_=lam,sigma=sigma,weak_w2_over_a_delta=ratio))

# Nonquadratic D1 and D2 density fixtures. Potential Hessian is in [.3a,.7a].
z=np.linspace(-9,9,361)
gamma=np.exp(-z*z/2)/math.sqrt(2*math.pi)
g,gw=roots_hermitenorm(28); gw=gw/math.sqrt(2*math.pi)
density=[]
for a in (.02,.08):
    zz,pp=np.meshgrid(z,g,indexing='ij'); n=zz.size
    initial=np.r_[zz.ravel(),pp.ravel()]
    def vector_field(t,y):
        q=y[:n]; p=y[n:]
        return np.r_[p,-q-b(q,a)]
    sol=solve_ivp(vector_field,(0,-T),initial,method='DOP853',rtol=2e-11,atol=2e-12)
    check(sol.success,'inverse nonlinear ODE success')
    p0=sol.y[n:,-1].reshape(zz.shape)
    check(np.max(np.abs(p0-zz)/(1+np.abs(zz)+np.abs(pp)))<2*a,'inverse first perturbation')
    h=np.exp(-U(z,a)); Z1=simpson(gamma*h,x=z); pi=gamma*h/Z1
    for delta in (.01,.0625,.25):
        k=math.sqrt(1-delta)
        r=np.sum(gw[None,:]/k*np.exp(-delta*p0*p0/(2*k*k)),axis=1)
        gd=np.exp(-delta*z*z/(2*k*k))/k
        M=simpson(pi*gd,x=z)
        mu=pi*r; nu=pi*gd/M
        check(abs(simpson(mu,x=z)-1)<2e-9,'D1 exact normalized density')
        ch1=simpson((mu-nu)**2/nu,x=z)
        check(ch1/(a*a*delta*delta)<20,'D1 scaled chi-square fixture')
        A=np.sum(gw[None,:]*np.exp(-U(k*k*z[:,None]+k*math.sqrt(delta)*g[None,:],a)),axis=1)
        Zk=float(np.sum(gw*np.exp(-U(k*g,a))))
        conv=gamma*A/Zk
        check(abs(simpson(conv,x=z)-1)<2e-10,'D2 convolved density normalization')
        ch2=simpson((conv-pi)**2/pi,x=z)
        check(ch2/(a*a*delta*delta)<20,'D2 scaled chi-square fixture')
        density.append(dict(a=a,delta=delta,mu_mass=float(simpson(mu,x=z)),chi2_mu_nu_over_a2_delta2=float(ch1/(a*a*delta*delta)),chi2_convolution_pi_over_a2_delta2=float(ch2/(a*a*delta*delta))))

# Deterministic integration of conditional stratified variances.
x,w=roots_legendre(48); w=w/2
strata=[]
a=.07; sigma=.3; k=math.sqrt(1-sigma*sigma); z0=.8; Z=-1.1
for N in (1,2,4,8,16,32,64):
    h=T/N; nodes=(np.arange(N)[:,None]+(x[None,:]+1)/2)*h
    q=np.cos(nodes)*z0+k*np.sin(nodes)*Z
    vals=np.cos(nodes)*b(q,a)
    means=vals@w
    variances=(vals*vals)@w-means*means
    vd=h*h*float(variances.sum())
    # Conservative uniform bound for |F'|.
    R=abs(z0)+abs(Z); L=2*a*R
    check(vd<=N*L*L*h**4/12*(1+1e-9),'stratum variance Lipschitz bound')
    lip=math.sqrt(N)*h*h*L/math.sqrt(2*math.pi)
    check(math.sqrt(vd)<=lip*(1+1e-8),'Gaussian Poincare conditional bound')
    strata.append(dict(N=N,var_D=vd,var_D_times_N3_over_a2=vd*N**3/a**2,buffered_upper=lip*math.sqrt(vd)/(2*sigma)))

# Exact polynomial current J and its rank-three Taylor remainder, one clock.
S=(x+1)*T/2
q=np.cos(S)*z0+k*np.sin(S)*Z
Fv=np.cos(S)*b(q,a); mean=float(w@Fv); D=-T*(Fv-mean)
B=k*Z-T*mean
psi=lambda v:v**4+6*sigma*sigma*v*v+3*sigma**4
lhs=float(w@(psi(B+D)-psi(B)))
quad=float((w@(D*D))/2*(12*B*B+12*sigma*sigma))
t=(x+1)/2
j=0.; r3=0.
for ti,wi in zip(t,w):
    j+=wi*(1-ti)*float(w@(D*D*(12*(B+ti*D)**2+12*sigma*sigma)))
    r3+=wi*(1-ti)**2/2*float(w@(D**3*24*(B+ti*D)))
check(abs(w@D)<1e-14,'conditional rank-one cancellation')
check(abs(lhs-j)<1e-13,'exact conditional current J')
check(abs(lhs-quad-r3)<1e-13,'rank-three remainder coefficient')
current=dict(lhs=lhs,rank_two=quad,rank_three=r3,integral_current=j)

# Explicit C1 cap and nonnegative C1 cubic Hermite filter.
def cap(z,B):
    az=abs(z)
    if az<=B: return z,1.
    if az>=2*B: return math.copysign(1.5*B,z),0.
    t=(az-B)/B
    return math.copysign(B*(1+t-t*t/2),z),1-t

def kernel(r,e):
    if r<=-e: return 0.,0.
    if r>=e: return math.sin(r),math.cos(r)
    t=(r+e)/(2*e)
    v=(-2*t**3+3*t*t)*math.sin(e)+(t**3-t*t)*2*e*math.cos(e)
    d=((-6*t*t+6*t)*math.sin(e)+(3*t*t-2*t)*2*e*math.cos(e))/(2*e)
    return v,d
filters=[]
for e in (.1,.03,.01,.003):
    rr=np.linspace(-e,e,1001); hd=np.array([kernel(float(r),e) for r in rr])
    truth=np.where(rr>=0,np.sin(rr),0)
    check(hd[:,0].min()>-1e-14,'Hermite nonnegativity')
    check(np.max(np.abs(hd[:,1]))<2,'Hermite first uniform bound')
    area=simpson(np.abs(hd[:,0]-truth),x=rr)
    check(area<e*e,'kernel integrated epsilon squared')
    filters.append(dict(epsilon=e,absolute_error_over_epsilon2=float(area/e**2)))
for zc in np.linspace(-3,3,121):
    cv,cd=cap(float(zc),.8)
    check(abs(cv)<=1.2+1e-14 and 0<=cd<=1,'cap bound and first')

# Literal stored-layer graph with one Hessian action per actual force query.
# Tangents include node motion, not just differentiation of force values.
def graph(inputs, counts, direction=None):
    M=len(counts); d=np.zeros_like(inputs) if direction is None else direction
    zz,ZZ,KK=inputs[:3]; dz,dZ,dK=d[:3]
    C0,C0p=cap(zz,.9); CZ,CZp=cap(ZZ,.9)
    aa=.06; sig=.25; kap=math.sqrt(1-sig*sig); eps=.017
    queries=0; hvps=0; stored=[]; off=3
    def eval_q(layer,t,dt):
        val=math.cos(t)*C0+kap*math.sin(t)*CZ
        der=math.cos(t)*C0p*dz+kap*math.sin(t)*CZp*dZ+(-math.sin(t)*C0+kap*math.cos(t)*CZ)*dt
        if layer>0:
            s,ds,fs,dfs=stored[layer-1]; hh=T/len(s)
            for sj,dsj,fj,dfj in zip(s,ds,fs,dfs):
                H,Hp=kernel(t-sj,eps)
                val-=hh*H*fj
                der-=hh*(Hp*(dt-dsj)*fj+H*dfj)
        return val,der
    for layer,N in enumerate(counts):
        G=inputs[off:off+N]; dG=d[off:off+N]; off+=N
        ss=(np.arange(N)+ndtr(G))*T/N
        ds=np.exp(-G*G/2)/math.sqrt(2*math.pi)*dG*T/N
        fs=[]; dfs=[]
        for sj,dsj in zip(ss,ds):
            qj,dq=eval_q(layer,float(sj),float(dsj))
            fs.append(float(b(qj,aa))); queries+=1
            dfs.append(float(bp(qj,aa))*dq); hvps+=int(direction is not None)
        stored.append((ss,ds,np.array(fs),np.array(dfs)))
    ss,ds,fs,dfs=stored[-1]; hh=T/len(ss)
    val=kap*ZZ+sig*KK-hh*np.sum(np.cos(ss)*fs)
    der=kap*dZ+sig*dK-hh*np.sum(-np.sin(ss)*ds*fs+np.cos(ss)*dfs)
    return float(val),float(der),queries,hvps
rng=np.random.default_rng(74109)
dag=[]
for counts in ((3,),(3,5),(3,5,7),(2,3,4,5)):
    for repetition in range(6):
        inputs=rng.normal(size=3+sum(counts)); direction=rng.normal(size=len(inputs))
        val,der,nq,nh=graph(inputs,counts,direction)
        h=2e-6
        fp=graph(inputs+h*direction,counts)[0]; fm=graph(inputs-h*direction,counts)[0]
        fd=(fp-fm)/(2*h)
        check(nq==sum(counts),'additive stored gradient query count')
        check(nh==sum(counts),'one HVP per stored gradient per first sweep')
        check(abs(der-fd)<2e-7*(1+abs(der)),'full clock/source directional derivative')
        dk=np.zeros_like(inputs); dk[2]=1
        check(abs(graph(inputs,counts,dk)[1]-.25)<1e-14,'exact unread K derivative')
        dag.append(dict(counts=counts,queries=nq,hvps=nh,first_error=abs(der-fd)))

# Exact rational all-layer exponents, including absent prior bank for M=1.
exponents=[]
for twice_R in range(5,31):
    R=F(twice_R,2); M=max(1,math.ceil(R-F(3,2)))
    alpha=(R-F(3,2))/2
    betas=[max(F(0),F(2,3)*(R-(M-j+F(3,2)))) for j in range(1,M)]
    beta_final=max(F(0),(R+alpha-F(5,2))/3)
    expected=max(F(0),(2*R-5)/3,R/2-F(13,12))
    actual=max([F(0),beta_final]+betas)
    check(actual==expected,'optimized all-layer exponent B')
    check(F(M)+F(3,2)>=R,'Picard tail target')
    for j,beta in enumerate(betas,1):
        check(M-j+F(3,2)+F(3,2)*beta>=R,'previous-layer exponent target')
    check(F(5,2)-alpha+3*beta_final>=R,'final buffered exponent target')
    check(F(3,2)+2*alpha==R,'keep exponent target')
    exponents.append(dict(R=str(R),M=M,sigma_exponent=str(alpha),earlier_counts=[str(v) for v in betas],last_count=str(beta_final),max_count=str(actual)))

# Near-kink mean defect: C2 potentials, no uniform b' modulus.
kinks=[]
for ratio in (1.,.1,.01,.001,.0001):
    aa=.07; eta=.2; dd=.1; eps=ratio*dd
    defect=aa*eta*(math.sqrt(dd*dd+eps*eps)-eps)
    scaled=defect/(aa*dd)
    check(0<scaled<=eta,'near-kink centered mean defect')
    if ratio<=.01: check(scaled>.98*eta,'near-kink defect retains first order')
    kinks.append(dict(epsilon_over_D=ratio,defect_over_a_absD=scaled))

# Isotropic random-scale fixture, with N=1 and exactly integrated clocks.
mixtures=[]
S=(x+1)*T/2
A=T*np.cos(S)**2; Bn=T*np.sin(S)*np.cos(S)
check(abs(float(w@Bn)-.5)<1e-13,'quadratic B expectation')
check(float(w@(Bn*Bn))-(float(w@Bn))**2>0,'quadratic B variance positive')
for aa in (.1,.01,.001):
    lam=.7; sig=.25; kap=math.sqrt(1-sig*sig)
    V=(aa*lam*A)**2/(1+aa*lam)+kap*kap*(1-aa*lam*Bn)**2+sig*sig
    ss=np.sqrt(V); mean=float(w@ss); var=float(w@((ss-mean)**2))
    # Huge d chosen to separate radial fluctuation from fixed Gaussian radius.
    dim=round(1e6/aa**4)
    # Rigorous coarser moment inequalities avoid huge-d gamma cancellation:
    # Gaussian Poincare gives Var(chi_d)<=1, hence E chi_d>=sqrt(d-1).
    lower=max(0,math.sqrt(dim-1)*math.sqrt(var)-math.sqrt(1/(1+aa*lam)))
    claimed_scale=math.sqrt(dim)*aa*aa/sig
    check(lower>0,'positive dimension-mixture lower bound')
    mixtures.append(dict(a=aa,d=dim,sd_sqrtV=math.sqrt(var),sd_over_a=math.sqrt(var)/aa,radial_lower=lower,lower_over_sqrtd_a2_over_sigma=lower/claimed_scale))
check(mixtures[-1]['lower_over_sqrtd_a2_over_sigma']>50*mixtures[0]['lower_over_sqrtd_a2_over_sigma'],'dimension-uniform weak factor cannot hold')

out=dict(status='PASS: deterministic algebra and numerical fixtures; analytic and source boundaries are in the audit',check_count=CHECKS,report_sha256=digest(REPORT),script_sha256=digest(Path(__file__)),quadratic=quadratic,nonquadratic_density=density,stratified_variance=strata,current=current,filters=filters,stored_dag=dag,exponents=exponents,kink=kinks,dimension_mixture=mixtures)
(ROOT/'scalar_weak_flow_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(status=out['status'],check_count=CHECKS,report_sha256=out['report_sha256'],nonquadratic_density=density,dimension_mixture=mixtures),indent=2))
