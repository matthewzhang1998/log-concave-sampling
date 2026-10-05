#!/usr/bin/env python3
"""Finite positive scalar atlas. Never evaluates the target tensor or OU means.

The copied gauss_rational_core is the pinned prior scalar root isolator.
Every node/weight output is rational; the iterator is lazy.
"""
from fractions import Fraction as F
from itertools import product
import math, json, pathlib, hashlib
import gauss_rational_core as core

H=128

def parameters(N,K,D,eta,epsilon,W=F(1),caller_weight=None,caller_weights=None):
    eta,epsilon,W=F(eta),F(epsilon),F(W)
    assert isinstance(D,int) and D>=1
    assert 0<eta<=F(1,2) and 0<epsilon<1 and W>=0
    if caller_weights is None:
        W1=N*W if caller_weight is None else F(caller_weight)
        caller_weights={0:W,1:W1}
    else:
        caller_weights={int(q):F(w) for q,w in caller_weights.items()}
        assert 0 in caller_weights and caller_weights[0]==W
    assert all(q>=0 and w>=0 for q,w in caller_weights.items())
    dsqrt=math.isqrt(D)
    if dsqrt*dsqrt<D:dsqrt+=1
    envelope=sum(core.exact_B(N,K+q,eta,w) for q,w in caller_weights.items())
    M=1+D*dsqrt*envelope
    J=max(1,core.ceil_log2_positive(32*M/epsilon))
    m=max(1,(J+1)//2)
    assert 32*M*F(1,1<<J)<=epsilon
    assert 32*M*F(1,1<<(2*m))<=epsilon
    return {'N':N,'K':K,'D':D,'eta':eta,'epsilon':epsilon,
            'M_rational':M,'J':J,'m':m,'subdivisions_per_panel':H,
            'interval_count_per_axis':H*J,
            'node_count_per_tree':2*(H*J*m)**2,
            'node_tolerance':epsilon/(64*M),
            'weight_tolerance':epsilon/(64*M),
            'caller_weights':caller_weights}

def gap_panels(J,H=H):
    assert J>=1 and H>=128
    ans=[]
    for j in range(J):
        g=F(1,1<<(j+1));lo=1-2*g
        for k in range(H):
            a=lo+k*g/H;b=lo+(k+1)*g/H
            assert b-a<=g/128 and 1-b>=g
            ans.append((a,b))
    assert len(ans)==J*H
    assert sum(b-a for a,b in ans)==1-F(1,1<<J)
    return ans

def prepare_certified_rule(N,K,D,eta,epsilon,W=F(1),caller_weight=None,caller_weights=None):
    par=parameters(N,K,D,eta,epsilon,W,caller_weight,caller_weights)
    base=core.certified_legendre_to_tolerance(par['m'],par['node_tolerance'],par['weight_tolerance'])
    panels=gap_panels(par['J'])
    return par,base,panels

def sector_nodes(base,panels):
    """One sector, with its Jacobian already in weight. Swap edge roles for other."""
    yield from core.sector_nodes(base,panels)

def serializable(obj):
    if isinstance(obj,F):return str(obj)
    if isinstance(obj,dict):return {str(k):serializable(v) for k,v in obj.items()}
    if isinstance(obj,list):return [serializable(v) for v in obj]
    return obj

def raw_call_census():
    trees=[[(0,1),(0,2)],[(0,1),(1,2)],[(0,2),(1,2)]]
    out={}
    for name,groups,centers,old in (
        ('triple_cubic',[range(0,3),range(3,6),range(6,9)],[1,4,7],[()]),
        ('three_three_six',[range(0,3),range(3,6),range(6,12)],[1,4,7,10],
         list(product(range(6,9),range(9,12))))):
        records=[]
        for oldhits in old:
            for edges in trees:
                slots=[groups[i] for edge in edges for i in edge]
                for hits in product(*slots):
                    kk=[0]*sum(map(len,groups))
                    for c in centers:kk[c]+=1
                    for c in oldhits+hits:kk[c]+=1
                    records.append({'raw':sum(2**(k+1) for k in kk),
                                    'degree':sum(kk),'max_order':max(kk)})
        out[name]={'histories_including_old_hit_pairs':len(records),
                   'raw_main_value_max':max(r['raw'] for r in records),
                   'max_order':max(r['max_order'] for r in records),
                   'degree_values':sorted(set(r['degree'] for r in records))}
    return out

def run_checks():
    checks=[]
    def ck(name,test,detail=None):
        assert test,name
        checks.append({'check':name,'pass':True,**({'details':detail} if detail is not None else {})})
    # Exact arithmetic identities for the atlas and all representative cells.
    panels=gap_panels(3)
    ck('atlas_exact_size_and_mass',len(panels)==384 and sum(b-a for a,b in panels)==F(7,8))
    for a,b in panels:
        ck('local_radius_ratio',64*(b-a)/(1-b)<=F(1,2))
    base=core.certified_legendre(4,110)
    one=[p for a,b in panels for p in core.map_rule(base,a,b)]
    keep=F(7,8)
    # Tensor-product total mass computed algebraically; no huge materialization.
    ck('positive_two_sector_mass',2*sum(x*w for x,w in one)*sum(w for x,w in one)==keep**3)
    ck('first_moment_exact',sum(x*w for x,w in one)==keep**2/2)
    # Sample all parent panel combinations, and first/last representative subcells.
    for ix in (0,127,128,255,256,383):
        for iy in (0,127,128,255,256,383):
            a,b=panels[ix];c,d=panels[iy]
            xlo,xhi=1-b,1-a;ylo,yhi=1-d,1-c
            gg=lambda x,y:x+y-x*y
            ck('heat_comparability',xhi<=2*xlo and yhi<=2*ylo and gg(xhi,yhi)<=2*gg(xlo,ylo))
            for r,wr in core.map_rule(base,a,b):
                for s,ws in core.map_rule(base,c,d):
                    ck('literal_node_weight',0<wr*ws*r<=(1-r)*(1-r*s))
    # Numerical checks below are diagnostics, not the analytic proof.
    import cmath
    boundary=[]
    for a in (0,.001,1/64,.03,1/16,.1,.25,.5,.9,.999999):
        best=1e10
        for j in range(721):
            z=a+(1-a)/64*cmath.exp(2j*math.pi*j/720)
            margin=-math.log(abs(z))-abs(cmath.phase(z))
            best=min(best,margin)
        ck('ou_disk_diagnostic',best>0,{'center':a,'minimum_margin':best})
    # Scalar Hermite identity E H_a(X1) H_b(X2) H_c(X3), known pairings.
    def pair_moment(a,b,c,r,s):
        x=(a+b-c)//2;y=(a+c-b)//2;z=(b+c-a)//2
        if (a+b+c)%2 or min(x,y,z)<0:return 0
        return math.factorial(a)*math.factorial(b)*math.factorial(c)/(math.factorial(x)*math.factorial(y)*math.factorial(z))*r**x*(r*s)**(y+z)
    for a,b,c in product(range(7),repeat=3):
        r=.31;s=.73
        base_m=pair_moment(a,b,c,1,1)
        rhs=(math.sqrt(r)**(a+b)*(s*math.sqrt(r))**c)*base_m
        ck('sector_hermite_factorization',abs(pair_moment(a,b,c,r,s)-rhs)<=1e-10*(1+abs(rhs)))
    census=raw_call_census()
    ck('cubic_raw_census',census['triple_cubic']=={'histories_including_old_hit_pairs':243,'raw_main_value_max':44,'max_order':3,'degree_values':[7]})
    ck('mixed_raw_census',census['three_three_six']=={'histories_including_old_hit_pairs':5832,'raw_main_value_max':72,'max_order':4,'degree_values':[10]})
    pars=parameters(9,7,100,F(1,1<<20),F(1,10**10))
    ck('exact_parameter_error',4*pars['M_rational']*(F(1,1<<pars['J'])+F(1,1<<(2*pars['m'])))<=pars['epsilon']/4)
    adaptive=core.certified_legendre_to_tolerance(8,F(1,10**25),F(1,10**25))
    ck('rational_adaptive_budget',adaptive['node_error']<=F(1,10**25) and adaptive['weight_l1_error']<=F(1,10**25))
    # Eta and D affect J logarithmically. Show exact counts, not asymptotic slogans.
    report=[]
    for N,K in ((9,7),(12,10)):
        for q in (1,10,100):
            for D in (1,100,10**6):
                p=parameters(N,K,D,F(1,1<<q),F(1,10**8))
                report.append({k:serializable(p[k]) for k in ('N','K','D','eta','J','m','node_count_per_tree')})
    return {'status':'PASS','checks':checks,'raw_census':census,'census':report,
            'scope':'Scalar exact-arithmetic and census tests plus labelled diagnostics; native source execution is specified, not run.'}

if __name__=='__main__':
    out=run_checks();path=pathlib.Path(__file__).with_name('checks.json')
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'check_count':len(out['checks']),'raw_census':out['raw_census'],'output':str(path)},indent=2))
