#!/usr/bin/env python3
"""Positive rational sector-Gauss builder and exact certificate diagnostics.

Only high-precision *proposals* use mpmath. Root signs, disjointness,
weight enclosures, positivity, mass and first moments use Fraction.
This certifies scalar quadrature arithmetic, not native source execution.
"""
from fractions import Fraction as F
import math, json
import mpmath as mp


def legendre_coeffs(n):
    p0,p1=[F(1)],[F(0),F(1)]
    if n==0:return p0
    for k in range(1,n):
        p=[F(0)]*(k+2)
        for i,v in enumerate(p1):p[i+1]+=(2*k+1)*v/(k+1)
        for i,v in enumerate(p0):p[i]-=k*v/(k+1)
        p0,p1=p1,p
    return p1


def peval(p,x):
    z=F(0)
    for c in p[::-1]:z=z*x+c
    return z


def mul(a,b):
    q=[x*y for x in a for y in b]
    return min(q),max(q)


def ieval(p,iv):
    z=(F(0),F(0))
    for c in p[::-1]:
        z=mul(z,iv);z=(z[0]+c,z[1]+c)
    return z


def sq(iv):
    a,b=iv
    return (F(0) if a<=0<=b else min(a*a,b*b),max(a*a,b*b))


def certified_legendre(n,bits=100):
    if n<1:raise ValueError('degree must be positive')
    p=legendre_coeffs(n)
    dp=[(i+1)*p[i+1] for i in range(n)]
    mp.mp.dps=max(40,math.ceil((bits+120)*math.log10(2)))
    roots,_=mp.gauss_quadrature(n,'legendre')
    D=1<<bits
    pos=[]
    for xx in list(roots)[(n+1)//2:]:
        j=int(mp.floor(xx*D))
        a,b=F(j-1,D),F(j+2,D)
        assert -1<a<b<1 and peval(p,a)*peval(p,b)<0
        pos.append((a,b))
    assert len(pos)==n//2
    brackets=[(-b,-a) for a,b in pos[::-1]]
    if n%2:brackets.append((F(0),F(0)))
    brackets+=pos
    assert len(brackets)==n
    assert all(brackets[i][1]<brackets[i+1][0] for i in range(n-1))
    for a,b in brackets:
        assert (a==b==0 and peval(p,F(0))==0) or peval(p,a)*peval(p,b)<0
    xs=[];wraw=[];wivs=[]
    for a,b in brackets:
        x=(a+b)/2
        xx=sq((a,b))
        d2=sq(ieval(dp,(a,b)))
        den=mul((1-xx[1],1-xx[0]),d2)
        assert den[0]>0
        wi=(2/den[1],2/den[0])
        xs.append(x);wivs.append(wi);wraw.append(sum(wi)/2)
    # Mirror rounding exactly. Mathematical Legendre symmetry justifies it.
    for i in range(n//2):
        j=n-1-i
        assert xs[i]==-xs[j]
        lo=max(wivs[i][0],wivs[j][0]);hi=min(wivs[i][1],wivs[j][1])
        assert lo<=hi
        wivs[i]=wivs[j]=(lo,hi)
        wraw[i]=wraw[j]=(lo+hi)/2
    total=sum(wraw)
    ws=[2*w/total for w in wraw]
    assert all(w>0 for w in ws)
    assert sum(ws)==2 and sum(x*w for x,w in zip(xs,ws))==0
    node_err=max((b-a)/2 for a,b in brackets)
    weight_err=sum(max(abs(w-lo),abs(w-hi)) for w,(lo,hi) in zip(ws,wivs))
    return {'x':xs,'w':ws,'brackets':brackets,'weight_intervals':wivs,
            'node_error':node_err,'weight_l1_error':weight_err}


def ceil_log2_positive(q):
    q=F(q)
    if q<=1:return 0
    j=max(0,q.numerator.bit_length()-q.denominator.bit_length())
    while F(1<<j)<q:j+=1
    while j and F(1<<(j-1))>=q:j-=1
    return j


def exact_B(N,K,eta,W=F(1)):
    eta=F(eta);W=F(W)
    assert 0<eta<=1 and W>=0
    fact=math.factorial(K);root=math.isqrt(fact)
    if root*root<fact:root+=1
    return W*(1<<(N+(3*K+1)//2))*root*eta**(-((K+1)//2))


def exact_parameters(N,K,P,eta,epsilon,W=F(1),caller_weight=None):
    """W is full-family mass; default N*W caller mass assumes rows <=1.
    Supply the actual caller_weight for larger or nonuniform caller rows.
    """
    eta=F(eta);epsilon=F(epsilon);P=max(F(1),F(P));W=F(W)
    assert 0<eta<=F(1,2) and 0<epsilon<1
    W1=F(caller_weight) if caller_weight is not None else N*W
    assert W>=0 and W1>=0
    B=1+exact_B(N,K,eta,W)+exact_B(N,K+1,eta,W1)
    J=max(1,ceil_log2_positive(32*B/epsilon))
    m=max(1,(J+1)//2)
    h=eta/(32*P);M=math.ceil(1/h)+J
    assert F(1,1<<J)<=epsilon/(32*B)
    assert F(1,1<<(2*m))<=epsilon/(32*B)
    return {'N':N,'K':K,'P':P,'eta':eta,'epsilon':epsilon,
            'B_star_rational':B,'J':J,'m':m,'h':h,
            'interval_bound_per_axis':M,'node_bound_per_tree':2*m*m*M*M,
            'node_tolerance':epsilon*eta/(512*B*P),
            'weight_tolerance':epsilon/(512*B)}


def certified_legendre_to_tolerance(n,node_tolerance,weight_tolerance):
    nd,wt=F(node_tolerance),F(weight_tolerance)
    assert nd>0 and wt>0
    bits=max(40,ceil_log2_positive(1/min(nd,wt))+12*n+32)
    while True:
        try:
            rule=certified_legendre(n,bits)
            if rule['node_error']<=nd and rule['weight_l1_error']<=wt:
                rule['certification_bits']=bits
                return rule
        except (AssertionError,ValueError):
            pass
        bits*=2


def prepare_certified_rule(N,K,P,eta,epsilon,W=F(1),caller_weight=None):
    """Total scalar preparation; returns a lazy potentially enormous node rule."""
    params=exact_parameters(N,K,P,eta,epsilon,W,caller_weight)
    base=certified_legendre_to_tolerance(params['m'],params['node_tolerance'],params['weight_tolerance'])
    panels=gap_panels(params['J'],params['h'])
    return params,base,panels

def gap_panels(J,h):
    """Coordinate intervals with 1-x in dyadic panels down to 2^-J."""
    ans=[]
    for j in range(J):
        lo=F(1)-F(1,1<<j);hi=F(1)-F(1,1<<(j+1))
        length=hi-lo;num=math.ceil(length/h)
        ans += [(lo+i*length/num,lo+(i+1)*length/num) for i in range(num)]
    assert sum(b-a for a,b in ans)==1-F(1,1<<J)
    assert len(ans)<=math.ceil(1/h)+J
    return ans


def map_rule(base,a,b):
    mid=(a+b)/2;half=(b-a)/2
    rr=[(mid+half*x,half*w) for x,w in zip(base['x'],base['w'])]
    assert sum(w for x,w in rr)==b-a
    assert sum(x*w for x,w in rr)==(b*b-a*a)/2
    assert all(a<x<b and w>0 for x,w in rr)
    return rr


def sector_nodes(base,panels):
    one=[item for a,b in panels for item in map_rule(base,a,b)]
    for r,wr in one:
        for s,ws in one:
            w=wr*ws*r
            assert w>0 and w<=(1-r)*(1-r*s)
            yield r,r*s,w


def B_log(N,K,eta,W=1):
    return math.log(W)+(N+1.5*K)*math.log(2)+.5*math.lgamma(K+1)-.5*K*math.log(eta)


def theorem_census(N,K,P,eta,epsilon,W=1,caller_weight=None):
    # Default caller_weight=N*W assumes actual caller rows <=1.
    # Logarithms avoid enormous constant overflow. Float ceiling here is a
    # displayed estimate; production uses outward intervals at integer walls.
    W1=caller_weight if caller_weight is not None else N*W
    logs=[0,B_log(N,K,eta,W),B_log(N,K+1,eta,W1)]
    lm=max(logs);lb=lm+math.log(sum(math.exp(x-lm) for x in logs))
    z=math.log(32)+lb-math.log(epsilon)
    J=max(1,math.ceil(z/math.log(2)));m=max(1,math.ceil(z/math.log(4)))
    M=math.ceil(32*max(1,P)/eta)+J
    return {'N':N,'K':K,'P':P,'eta':eta,'epsilon':epsilon,
            'log_B_star':lb,'J':J,'m':m,'interval_bound_per_axis':M,
            'node_bound_per_tree':2*m*m*M*M,
            'note':'sufficient conservative bound; all nodes genuinely cost execution'}


def run_checks():
    results=[]
    def ck(name,test,details=None):
        assert test,name
        results.append({'check':name,'pass':True,**({'details':details} if details else {})})
    for k in range(101):
        ck(f'local_factorial_{k}',(k+1)**k<=4**k*math.factorial(k))
    for K in (7,8,10,11):
        for l in range(31):
            ck(f'price_factorial_K{K}_l{l}',math.factorial(K+2*l)<=2**(K+4*l)*math.factorial(K)*math.factorial(l)**2)
    for n in (1,2,3,4,8,12):
        rule=certified_legendre(n,100)
        ck(f'positive_rational_gauss_{n}',sum(rule['w'])==2 and all(w>0 for w in rule['w']),
           {'node_error':float(rule['node_error']),'weight_l1_error':float(rule['weight_l1_error'])})
        for k in range(2*n):
            actual=sum(w*x**k for x,w in zip(rule['x'],rule['w']))
            exact=F(0) if k%2 else F(2,k+1)
            bound=rule['weight_l1_error']+2*k*rule['node_error']
            ck(f'gauss_moment_{n}_{k}',abs(actual-exact)<=bound)
    pars=exact_parameters(9,7,27,F(1,2),F(1,100))
    ck('exact_parameter_chooser',F(1,1<<pars['J'])<=pars['epsilon']/(32*pars['B_star_rational']))
    adap=certified_legendre_to_tolerance(8,F(1,10**35),F(1,10**35))
    ck('adaptive_precision_budget',adap['node_error']<=F(1,10**35) and adap['weight_l1_error']<=F(1,10**35),{'bits':adap['certification_bits']})
    rule=certified_legendre(4,90)
    # Small illustrative atlas; intentionally not the worst-case theorem h.
    panels=gap_panels(4,F(1,8))
    nodes=list(sector_nodes(rule,panels))
    keep=F(15,16)
    ck('positive_sector_mass',2*sum(w for t1,t2,w in nodes)==keep**3)
    ck('positive_sector_count',len(nodes)==(len(panels)*4)**2)
    for cell1 in panels:
        for cell2 in panels:
            xlo,xhi=1-cell1[1],1-cell1[0]
            ylo,yhi=1-cell2[1],1-cell2[0]
            f=lambda x,y:x+y-x*y
            ck('cell_gap_comparability',f(xhi,yhi)<=2*f(xlo,ylo))
    # Diagnostic analytic fixture on full triangle, not a tensor-oracle test.
    # F(tlarge,tsmall)=exp(.7*tlarge+.2*tsmall).
    mp.mp.dps=60
    truth=mp.quad(lambda r:mp.quad(lambda s:r*mp.exp(mp.mpf('.7')*r+mp.mpf('.2')*r*s),[0,1]),[0,1])
    num=sum(float(w)*math.exp(.7*float(t1)+.2*float(t2)) for t1,t2,w in nodes)
    tail_bound=2*float(1-keep)*math.exp(.9)
    ck('scalar_analytic_fixture',abs(float(truth)-num)<=tail_bound,
       {'exact_one_sector':float(truth),'quadrature_truncated':num,'absolute_error':abs(float(truth)-num),'declared_tail_bound':tail_bound})
    census=[theorem_census(9,7,27,2**-q,10**-p) for q in (1,3,5) for p in (2,8,20)]
    return {'status':'PASS','checks':results,'census':census,
            'scope':'Exact rational scalar checks. No executed native sampler, uniform Sobolev oracle, or numerical evidence substituted for the proof.'}

if __name__=='__main__':
    import pathlib
    out=run_checks()
    path=pathlib.Path(__file__).with_name('checks.json')
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'checks':len(out['checks']),'output':str(path)},indent=2))
