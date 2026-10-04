#!/usr/bin/env python3
"""Independent tests of the finite ordered-refresh graph, not a CW7 implementation."""
import json, math, os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from fractions import Fraction as F
from pathlib import Path
import numpy as np
OUT=Path(__file__).parent
MARKS=[F('2.4'),F('2.7'),F('2.8'),F('2.9'),F(3),F('3.1'),F('3.2'),F('3.3'),F('3.4'),F('3.5'),F('3.55')]
checks=[]
def ck(name,ok,**data):
    if not ok:raise AssertionError((name,data))
    checks.append(dict(name=name,passed=True,**data))
def grad(x):
    x=np.asarray(x);r=np.sqrt(np.abs(x));v=r*r-2*r+2*np.log1p(r)
    sm=r<1e-3
    v=np.where(sm,sum(2*(-1)**(m+1)*r**m/m for m in range(3,12)),v)
    return .5*x+.2*np.sign(x)*v

def hvp(x):
    r=np.sqrt(np.abs(x));return .5+.2*r/(1+r)

# Expressions are the normalized physical state Cw + sum_j a_j p_j.
def expr(c,e=None):return (np.array(c,dtype=float),{} if e is None else dict(e))
def add(x,y,s=1.):
    e=x[1].copy()
    for k,v in y[1].items():e[k]=e.get(k,0)+s*v
    return x[0]+s*y[0],e

def scale(x,a):return x[0]*a,{k:v*a for k,v in x[1].items()}
class Graph:
    def __init__(self,a,N=3,M=2):
        self.a=a;self.N=N;self.M=M;self.cs=[];self.es=[];self.ds=[];self.stages=[];self.n=len(MARKS)+2
        e=np.eye(self.n);x=expr(e[0]);x=self.flow(x,e[1],math.pi/2,-1)
        for i,J in enumerate(MARKS):x=self.flow(x,e[i+2],a**float((J-1)/3),i)
        self.terminal=len(self.cs);self.cs.append(x[0]);self.es.append(x[1]);self.ds.append(1.);self.stages.append(len(MARKS))
        m=len(self.cs);d=np.array(self.ds);self.C=np.array(self.cs)*d[:,None];self.A=np.zeros((m,m));self.d=d
        for i,row in enumerate(self.es):
            for j,v in row.items():self.A[i,j]=d[i]*v/d[j]
        self.P=x[0].copy()
    def flow(self,x,momentum,T,stage):
        t=np.arange(self.N+1)*T/self.N;co=np.cos(t);si=np.sin(t)
        if stage==-1:co[-1]=0;si[-1]=1
        base=[add(scale(x,co[i]),expr(si[i]*momentum)) for i in range(self.N+1)]
        states=base
        for _ in range(self.M):
            ids=[]
            for j in range(self.N):
                ids.append(len(self.cs));self.cs.append(states[j][0]);self.es.append(states[j][1]);self.ds.append(self.N**-.5);self.stages.append(stage)
            states=[]
            for i in range(self.N+1):
                row=expr(base[i][0],base[i][1])
                for j in range(i):
                    w=math.cos(t[i]-t[j+1])-math.cos(t[i]-t[j])
                    row[1][ids[j]]=row[1].get(ids[j],0)-self.a*w
                states.append(row)
        return states[-1]
    def val(self,w,C=None):
        C=self.C if C is None else C;p=np.zeros(len(self.d))
        for i in range(len(p)):
            z=C[i]@w+self.A[i]@p
            p[i]=self.d[i]/math.sqrt(self.a)*grad(math.sqrt(self.a)/self.d[i]*z)
        return p
    def protected(self,i,all_later=True):
        C=self.C.copy()
        for row,stage in enumerate(self.stages):
            if row!=self.terminal and (stage>=i if all_later else stage==i):C[row,i+2]=0
        return C

rng=np.random.default_rng(517901)
projection=[];principal=[]
for a in [.05,.02,.008]:
    g=Graph(a);m=len(g.d);w=rng.normal(size=g.n);w0=w.copy()
    ck('physical_final_coisometry',abs(g.P@g.P-1)<1e-13,a=a)
    for i,J in enumerate(MARKS):
        C=g.protected(i);bad=g.protected(i,False);s=g.P[i+2];gam=s*s;h=a**float((J-1)/3);eps=a**float(J-2)
        beta=(1/gam+1/eps**2)**-.5;B=beta*np.r_[np.eye(g.n)[i+2]/s,-1/eps]
        Ct=np.block([[C,np.zeros((m,1))],[C,eps*np.eye(m)[:,[-1]]]])
        S=np.block([[g.A+g.A.T,g.A.T],[g.A,np.zeros_like(g.A)]])
        e=np.zeros(2*m);e[g.terminal]=1
        ck('annihilator_all_later_force_rows',np.max(abs(Ct@B-beta*e))<1e-12,a=a,mark=str(J))
        ck('signed_coisometry',abs(B@B-1)<1e-12,a=a,mark=str(J))
        if i<len(MARKS)-1:
            ck('stage_only_projection_has_later_leak',np.linalg.norm(bad[:-1,i+2])>1e-8,a=a,mark=str(J))
        full=g.val(w)[-1];prot=g.val(w,C)[-1]
        price=abs(full-prot) # physical r=1
        bound=a*h*sum(a**float(2*(JJ-1)/3) for JJ in MARKS[i:])*abs(w[i+2])
        ck('finite_complete_projection_bound',price<=4*bound+2e-12,a=a,mark=str(J),ratio=price/max(bound,1e-300))
        wi=w.copy();wi[i+2]=0
        ck('projection_exact_zero_on_protected_coordinate',abs(g.val(wi)[-1]-g.val(wi,C)[-1])<2e-13,a=a,mark=str(J))
        projection.append(dict(a=a,mark=str(J),ratio=price/max(bound,1e-300),beta=beta,nominal_beta=a**float(max((J-1)/3,J-2))))
        # Literal simultaneous Picard VALUE and first recurrence. No derivative of a Hessian.
        q=np.linalg.norm(abs(S),2);cc=np.linalg.norm(Ct,2)
        ck('signed_absolute_contraction',q<.375,a=a,mark=str(J),q=q)
        z=np.r_[w,.7];p=np.zeros(2*m);D=np.zeros((2*m,g.n+1));dd=np.r_[g.d,g.d];sg=np.r_[np.ones(m),-np.ones(m)]
        P=np.r_[g.P,0.]
        for it in range(1,7):
            at=Ct@z+S@p;phys=math.sqrt(a)*at/dd
            pnew=sg*dd/math.sqrt(a)*grad(phys);Dnew=(sg*hvp(phys))[:,None]*(Ct+S@D)
            DG=Ct.T@Dnew/beta;physder=B@DG
            pb=float(B@DG@B)
            ck('finite_principal_block',abs(pb)<=beta/(1-q)+2e-12,a=a,mark=str(J),iteration=it,principal=pb,bound=beta/(1-q))
            expected=hvp(phys[g.terminal])*(P+g.A[-1]@D[:m])
            ck('finite_recorded_physical_row',np.max(abs(physder-expected))<2e-11,a=a,mark=str(J),iteration=it)
            curl=np.outer(P,physder)-np.outer(physder,P)
            cb=2*np.linalg.norm(g.A[-1])*cc/(1-q)
            ck('recorded_curl_from_terminal_first',np.linalg.norm(curl,2)<=cb+2e-11,a=a,mark=str(J),iteration=it)
            if it==6:principal.append(dict(a=a,mark=str(J),q=q,principal_ratio=abs(pb)*(1-q)/beta))
            p,D=pnew,Dnew

for i,J in enumerate(MARKS):
    x=(J-1)/3;chi=max(x,J-2)
    ck('projection_mark_exact',1+3*x==J,mark=str(J))
    ck('twin_mark_exact',2+(J-2)==J,mark=str(J))
    ck('later_refreshes_no_projection_order_loss',all(1+x+2*(K-1)/3>=J for K in MARKS[i:]),mark=str(J))
ck('pure_fourth_weighted_body',F(3,4)-max((MARKS[0]-1)/3,MARKS[0]-2)-F(1,4)==F(1,30))
chain=MARKS[1:]
for older,newer in zip(chain[:-1],chain[1:]):
    u=(F(327,50)-older)/4;e=newer-F(7,2)+u+F(1,100)
    for J in [older,newer]:
        chi=max((J-1)/3,J-2);gates=[F(3,2)-u+chi-e,F(3,2)-u-chi+e,F(3,2)-u]
        ck('exact32_transition_strict_gates',min(gates)>0,older=str(older),newer=str(newer),child=str(J),gates=list(map(str,gates)))
chi=F(31,20);u=F(3,4);e=F(81,100)
final_gates=[F(3,2)-u+chi-e,F(3,2)-u-chi+e,F(3,2)-u]
ck('exact32_final_strict_gates',min(final_gates)==F(1,100),gates=list(map(str,final_gates)))
for p in [2,3,4,8,16]:
    c=F(13);dmax=F(16);b=F(7,3);delta=F(1,200);E=F(20)
    B=math.ceil((E+b+c/2+(c+dmax)/(4*p)+1)*2*p/delta)
    exponent=delta*B/(2*p)-b-c/2-(c+dmax)/(4*p)
    ck('ambient_high_moment_order_choice',exponent>=E+1,p=p,B=B,exponent=str(exponent))
ck('wrong_order_can_lose_projection_grade',1+(F('3.55')-1)/3+2*(F('2.4')-1)/3<F('3.55'))
summary=dict(status='PASS',checks=len(checks),projection_cases=len(projection),principal_cases=len(principal),
             max_projection_ratio=max(x['ratio'] for x in projection),max_principal_ratio=max(x['principal_ratio'] for x in principal))
(OUT/'outer_refresh_port_checks.json').write_text(json.dumps(dict(summary=summary,projection=projection,principal=principal,checks=checks),indent=2)+'\n')
print(json.dumps(summary,indent=2))
