#!/usr/bin/env python3
"""Exact/finite diagnostics. Does not execute imported native source programs."""
from pathlib import Path
import importlib.util, itertools, json, math
from collections import Counter
import numpy as np
import sympy as sp
HERE=Path(__file__).resolve().parent
SPEC=importlib.util.spec_from_file_location('rank_generator','/workspace/shared/rank-indexed-positive-returns-20261005/generate_rank_terms.py')
gen=importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(gen)

def hit(ast, vertex, edge):
    if ast[0]=='g':
        return ('g',ast[1],ast[2]+(edge,)) if ast[1]==vertex else None
    if ast[0]=='r':
        c=hit(ast[2],vertex,edge)
        return None if c is None else ('r',ast[1]+1,c)
    a,b=hit(ast[1],vertex,edge),hit(ast[2],vertex,edge)
    assert (a is None)!=(b is None)
    return ('mul',a,ast[2]) if a is not None else ('mul',ast[1],b)

def double_star():
    a=('r',2,('mul',('r',2,('g',1,(1,))),('r',2,('g',2,(1,)))))
    coefficient=2
    for n in range(3,9):
        chosen=1 if n<=5 else 2
        da=hit(a,chosen,n-1)
        assert da in gen.derivative(a,n-1)
        a=('r',n,('mul',da,('r',2,('g',n,(n-1,)))))
        coefficient*=n
    return coefficient,a

def geometry(ast, values):
    rows={}; sig={}; d={}; columns=[]; structural=[]; primitive=[]; orders={}
    def walk(a, caller, cy, address=()):
        if a[0]=='mul':
            walk(a[1],caller,cy,address+(0,));walk(a[2],caller,cy,address+(1,));return
        assert a[0]=='r'
        r=float(values(a,address)); child=a[2]
        row={j:r*v for j,v in caller.items()}
        if child[0]=='g':
            i=child[1];name=('W',i);columns.append(name);primitive.append(name)
            sigma=math.sqrt((1-r*r)/2);row[name]=sigma
            rows[i]=row;sig[i]=sigma;d[i]=r*cy;orders[i]=len(child[2]);return
        name=('G',address);columns.append(name);structural.append(name)
        row[name]=math.sqrt(1-r*r)
        walk(child,row,r*cy,address+(2,))
    walk(ast,{},1.)
    ids=sorted(rows);A=np.array([[rows[i].get(c,0.) for c in columns] for i in ids])
    return dict(A=A,sigma=np.array([sig[i] for i in ids]),d=np.array([d[i] for i in ids]),
                structural=[columns.index(c) for c in structural],columns=columns,orders=orders)

def floor_certificate(A,sigma,structural):
    n=A.shape[0]; best=float(min(sigma*sigma));choice={'method':'primitive_diagonal'}
    for s in range(1,min(n,len(structural))+1):
        for S in itertools.combinations(range(n),s):
            R=tuple(i for i in range(n) if i not in S)
            for J in itertools.combinations(structural,s):
                T=A[np.ix_(S,J)]
                smallest=float(np.linalg.svd(T,compute_uv=False)[-1])
                if smallest<1e-10:continue
                eta=1/smallest
                if not R:ell=1/(eta*eta)
                else:
                    u=float(np.linalg.norm(A[np.ix_(R,J)],2))
                    ell=min(float(min((sigma[list(R)])**2))/(1+2*eta*eta*u*u),1/(2*eta*eta))
                if ell>best:best=ell;choice={'method':'structural_minor','S':S,'J':J}
    return best,choice

def geometry_checks():
    rng=np.random.default_rng(20261005);ratios=[];errors=[];improvements=[]
    for rank in range(2,7):
        histories=gen.histories(rank)
        for _ in range(12):
            _,ast=histories[int(rng.integers(len(histories)))]
            def clocks(a,address):
                return float(np.sqrt(1-10**(-rng.uniform(.02,5))))
            g=geometry(ast,clocks);A=g['A'];s=g['sigma'];d=g['d'];G=A@A.T
            assert np.max(np.abs(np.diag(G)+s*s+d*d-1))<1e-12
            ell,choice=floor_certificate(A,s,g['structural'])
            eig=float(np.linalg.eigvalsh(G)[0]);ratios.append(eig/ell)
            improvements.append(ell/min(s*s))
            assert eig>=ell*(1-1e-7)
            c=.37;h=math.sqrt(1-c);delta=ell/2
            vals,vecs=np.linalg.eigh(G-delta*np.eye(rank));assert min(vals)>0
            B=(vecs*np.sqrt(vals))@vecs.T;t2=s*s+(1-c)*delta
            within=c*G+(1-c)*B@B.T+np.diag(t2)
            error=np.max(np.abs(within-(G+np.diag(s*s))));errors.append(float(error))
            assert error<2e-12
            assert np.max(c*np.diag(G)+(1-c)*np.sum(B*B,axis=1))<=1+1e-12
    return {'cases':len(ratios),'minimum_eigenvalue_over_certificate':min(ratios),
            'maximum_floor_improvement':max(improvements),'maximum_reconstruction_error':max(errors)}

def obstruction_checks():
    coeff,ast=double_star();orders=sorted(len(es) for _,es in gen.vertices(ast))
    assert coeff==40320 and orders==[1]*6+[4,4]
    rows=[];vinf=.25+math.exp(-4)/8-math.exp(-2)/4
    a,b=.5,.25
    for K in [4,8,16,32,64,128,256]:
        x=K**-2
        g=geometry(ast,lambda node,address:math.sqrt(1-2*x) if node[2][0]=='g' and node[2][1] in (1,2) else .5)
        A=g['A'];G=A@A.T;v=G[0,1]/(1-2*x);V=K*K*(1-2*x)*v
        expected=np.array([[v*(1-2*x)+x,v*(1-2*x)],[v*(1-2*x),v*(1-2*x)+x]])
        assert np.max(np.abs(G[:2,:2]-expected))<1e-13
        e=np.zeros(8);e[0]=1;e[1]=-1
        assert abs(e@G@e-2*x)<1e-13
        assert np.linalg.matrix_rank(A[:2,g['structural']],tol=1e-10)==1
        variance=.25+math.exp(-4)*(1+math.exp(-8*V))/8-math.exp(-2)*(1+math.exp(-4*V))/4
        mean=math.exp(-1)*(1-math.exp(-2*V))/2
        prefactor=b*b*a**6*K*K*math.exp(-1)
        varL=prefactor*prefactor*variance
        lower_diff=(1-math.exp(-2))**2/8
        assert variance>=lower_diff
        rows.append({'K':K,'primitive_half_heat':x,'parent_variance':v,'Var_S':variance,'mean_S':mean,
                     'leading_Var_L':varL,'Var_L_over_K4':varL/K**4,
                     'normalized_center_weight':x*x/(x**1.5*x**1.5)})
    assert abs(rows[-1]['Var_S']-vinf)<1e-12
    # Low-frequency difference mode is visible on a fixed interior bridge interval.
    nodes,weights=np.polynomial.legendre.leggauss(128)
    bridge=sum(float(w)*(.5-.25)/2 *
        (2*float(r)*.25*math.exp(-2*(1-float(r)**2))*(1-math.exp(-4*float(r)**2)))
        for r,w in zip(.25+(.5-.25)*(nodes+1)/2,weights))
    assert bridge>0
    # Equal t=O(sqrt(x)) is unavoidable. Test a proposed O(1) extraction fails.
    x=1/4096;B=np.ones((2,2))*.9+x*np.eye(2);h2=.5
    eig=float(np.linalg.eigvalsh(h2*B-.01*np.eye(2))[0]);assert eig<0
    # Exact coefficient depends polynomially on parent old Gaussian variables, all observers independent.
    return {'history_coefficient':coeff,'force_degrees':orders,'structural_clock_count':7,
            'center_primitive_orders':[5,5],'ast':ast,'variance_limit':vinf,
            'difference_mode_variance':(1-math.exp(-2))**2/8,
            'difference_mode_bridge_mass_rho_025_to_05':bridge,
            'illegal_constant_extraction_min_eigenvalue':eig,'rows':rows,
            'scope':'One-node leading scalar coefficient; actual six leaf factors have uniform exponentially small correction.'}

def allocation_checks():
    tests=[]
    for n in [2,3,4,8,12]:
        m=2*n;beta=(m-2)/(m-1)
        for K in [1e-4,1.,1e8,1e30]:
            C=100.;gap=.01
            alpha=.1*min(1.,math.sqrt(gap/(C*K)),(gap/C)**(1/beta),gap/(C*K))
            root=K*alpha**2;other=alpha**beta
            assert C*root<=gap and C*other<=gap and C*K*alpha<=gap
            logproduct=math.log(root)+(m-1)*math.log(other)
            assert abs(logproduct-(m*math.log(alpha)+math.log(K)))<1e-10
            tests.append({'n':n,'K':K,'alpha':alpha,'beta':beta,'root_radius':root,'nonroot_radius':other})
    return tests

q,p,t=sp.symbols('q p t')
def gauss(poly,var):
    poly=sp.Poly(sp.expand(poly),var);out=0
    for (d,),c in poly.terms():
        if d%2==0:out+=c*(sp.factorial2(d-1) if d else 1)
    return sp.expand(out)
def R(poly):
    expr=sp.Poly(sp.expand(poly),q);remaining=expr.as_expr();out=0
    for n in range(expr.degree(),0,-1):
        c=sp.Poly(remaining,q).coeff_monomial(q**n)
        h=sp.hermite_prob(n,q);remaining=sp.expand(remaining-c*h)
        out+=c*sp.hermite_prob(n-1,q)
    return sp.expand(out)
def noise_moment(m,mu,var):
    return sp.expand(sum(sp.binomial(m,2*j)*mu**(m-2*j)*var**j*(sp.factorial2(2*j-1) if j else 1) for j in range(m//2+1)))
def second_riesz_checks():
    checks=[]
    for k in range(4):
        F=(q+q*q-1)*sp.hermite_prob(k,p);Fq=sp.diff(F,q)
        C=gauss(F*F,q);tau=R(F)*Fq;assert sp.expand(gauss(tau,q)-C)==0
        J=R(tau-C)*Fq;base=p+p*p/sp.Integer(10);mu=base+t*F;var=1+(1-t*t)*C
        for m in [3,4,5]:
            lhs=sp.diff(gauss(gauss(noise_moment(m,mu,var),q),p),t)
            rhs=t*t*sp.factorial(m)/sp.factorial(m-3)*gauss(gauss(J*noise_moment(m-3,mu,var),q),p)
            assert sp.expand(lhs-rhs)==0
            checks.append({'physical_degree':k,'test_power':m,'exact_identity':True})
    return {'checks':checks,'scope':'Polynomial identity only; this unbounded-Q polynomial is not a witness for the uniform derivative hypothesis.'}

def main():
    report={'geometry':geometry_checks(),'rank8_obstruction':obstruction_checks(),
            'finite_allocation':allocation_checks(),'second_riesz':second_riesz_checks()}
    (HERE/'checks.json').write_text(json.dumps(report,indent=2)+'\n')
    (HERE/'rank8_ast.json').write_text(json.dumps(report['rank8_obstruction']['ast'],indent=2)+'\n')
    print(json.dumps({'geometry':report['geometry'],'rank8_variance_constant':report['rank8_obstruction']['variance_limit'],
                      'bridge_interior_mass':report['rank8_obstruction']['difference_mode_bridge_mass_rho_025_to_05'],
                      'allocation_cases':len(report['finite_allocation']),
                      'second_riesz_identity_cases':len(report['second_riesz']['checks'])},indent=2))
if __name__=='__main__':main()
