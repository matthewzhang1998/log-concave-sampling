"""Diagnostics for actual VALUE graph, first/curl, coherent current, and pins.
Analytic Hessians occur only in diagnostic checks, never in shifted_atom.py.
"""
import hashlib,json,math
from pathlib import Path
import numpy as np
from shifted_atom import shifted_atom,coherent_pair,marked_current
P=Path(__file__).resolve().parent
rng=np.random.default_rng(516881)
checks={}
def ok(name,condition):
    if not bool(condition): raise AssertionError(name)
    checks[name]=checks.get(name,0)+1
class Source:
    def __init__(self,A,D):
        self.A=A; self.D=D; self.calls=[]
        r=np.vstack([np.eye(D),rng.normal(size=(D+1,D))]);r/=np.linalg.norm(r,axis=1)[:,None]
        self.rows=r
    def __call__(self,x):
        x=np.asarray(x);self.calls.append(x.copy())
        return self.A*(.2*x+.4*np.tanh(x@self.rows.T)@self.rows/len(self.rows))
    def jac(self,x):
        z=np.asarray(x)@self.rows.T
        return self.A*(.2*np.eye(self.D)+.4*self.rows.T@np.diag(1-np.tanh(z)**2)@self.rows/len(self.rows))
def differentiate(g,sites,C,a,b,c):
    J=[g.jac(x) for x in sites]
    d0=J[0]@C[0]
    d1=J[1]@(C[1]-a*d0);d1b=J[2]@C[1]
    du=J[3]@(C[2]-b*d1);dv=J[4]@(C[2]-b*d1b)
    de=J[5]@(C[3]-c*du)-J[6]@(C[3]-c*dv)
    return du,dv,du-dv,de
commutator_max=0.
for D in [2,3,5]:
 for A in [.03125,.125,.5]:
  g=Source(A,D);n=4*D
  K,_=np.linalg.qr(rng.normal(size=(n,n)))
  sigma=np.sqrt(.75)
  C=[sigma*K[j*D:(j+1)*D] for j in range(4)]
  P3=C[3]/sigma
  for trial in range(8):
   a,b,c=rng.uniform(.1,1,size=3)
   Y=rng.normal(size=D);omega=rng.normal(size=n)
   x=[row@omega+.5*Y for row in C]
   g.calls=[]
   e,rec,sites=shifted_atom(g,*x,a,b,c)
   ok('exact_seven_queries',len(g.calls)==7)
   ok('original_site_record',all(np.array_equal(t,s) for t,s in zip(g.calls,sites)))
   ok('mark_energy',np.linalg.norm(rec.delta)<=a*b*A**3*np.linalg.norm(x[0])+1e-14)
   ok('response_energy',np.linalg.norm(e)<=a*b*c*A**4*np.linalg.norm(x[0])+1e-14)
   du,dv,dr,de=differentiate(g,sites,C,a,b,c)
   ok('u_first',np.linalg.norm(du,2)<=A*(1+b*A*(1+a*A))+1e-13)
   ok('v_first',np.linalg.norm(dv,2)<=A*(1+b*A)+1e-13)
   ok('mark_first',np.linalg.norm(dr,2)<=9*A/4+1e-13)
   ok('response_first',np.linalg.norm(de,2)<=21*A/8+1e-13)
   lift=P3.T@de
   ok('full_bank_curl',np.linalg.norm(lift-lift.T,2)<=13*A*A/2+1e-13)
   zeta=rng.normal(size=n);eps=2e-5
   eplus=shifted_atom(g,*[xx+eps*cc@zeta for xx,cc in zip(x,C)],a,b,c)[0]
   eminus=shifted_atom(g,*[xx-eps*cc@zeta for xx,cc in zip(x,C)],a,b,c)[0]
   ok('jacobian_full_replay',np.linalg.norm((eplus-eminus)/(2*eps)-de@zeta)<1e-7)
   # Exact current identity on the same realized mark, no law approximation.
   keep=rng.normal(size=D);theta=.38;q=.81;v=.14;t=rng.normal(size=D)
   mark,z,current=marked_current(rec,keep,v,q,theta)
   zplus=marked_current(rec,keep,v,q,theta+eps)[1]
   zminus=marked_current(rec,keep,v,q,theta-eps)[1]
   derivative=(np.exp(1j*t@zplus)-np.exp(1j*t@zminus))/(2*eps)
   expected=1j*(t@current)*np.exp(1j*t@z)
   ok('positive_current_identity',abs(derivative-expected)<1e-8)
   # Full affine row pure term is exactly symmetric.
   Jdiff=g.jac(sites[5])-g.jac(sites[6])
   pure=P3.T@Jdiff@C[3]
   ok('full_row_symmetry',np.linalg.norm(pure-pure.T)<1e-13)
   commutator_max=max(commutator_max,np.linalg.norm(g.jac(sites[1])@g.jac(sites[3])-g.jac(sites[3])@g.jac(sites[1])))
   # Mean-value chord matrices computed only as a diagnostic.
   nodes,weights=np.polynomial.legendre.leggauss(32);nodes=(nodes+1)/2;weights=weights/2
   def secant(base,shift):return sum(w*g.jac(base-s*shift) for s,w in zip(nodes,weights))
   J1=secant(x[1],a*rec.h0)
   J2=secant(x[2]-b*rec.h1b,b*(rec.h1-rec.h1b))
   J3=secant(x[3]-c*rec.v,c*rec.delta)
   chain=-a*b*c*J3@J2@J1@rec.h0
   ok('shifted_ordered_chain',np.linalg.norm(e-chain)<2e-12)
ok('noncommuting_original_hessians',commutator_max>1e-9)
# Linear exact witness, multiple dimensions and nontrivial edges.
for D in [1,2,7]:
 for A in [.03125,.25,.5]:
  x=rng.normal(size=(4,D));a,b,c=.7,.4,.9
  e,rec,_=shifted_atom(lambda z:A*z,*x,a,b,c)
  ok('linear_exact_response',np.allclose(e,-a*b*c*A**4*x[0],atol=1e-15))
  ok('linear_exact_mark',np.allclose(rec.delta,a*b*A**3*x[0],atol=1e-15))
# Empirical block is a CHECK of the true covariance algebra, not an executing tensor.
D=3;A=.4;g=Source(A,D);Y=rng.normal(size=D);N=14000
xs=[.5*Y+np.sqrt(.75)*rng.normal(size=(N,D)) for _ in range(3)]
g.calls=[];rec=coherent_pair(g,*xs,.9,.8)
ok('vectorized_pair_five_original_sites',len(g.calls)==5)
W=np.concatenate([rec.v,rec.delta],axis=1);W-=W.mean(0)
B=W.T@W/N
ok('whole_block_psd',np.linalg.eigvalsh(B).min()>-1e-13)
for theta in [0,.2,.5,1]:
 L=np.concatenate([np.eye(D),theta*np.eye(D)],axis=1)
 u=(1-theta)*rec.v+theta*rec.u;u-=u.mean(0)
 ok('coherent_block_readout',np.linalg.norm(u.T@u/N-L@B@L.T)<1e-13)
# Rational native term powers after physical sqrt(u).
from fractions import Fraction as F
r=(F(1),F(-1,2));delta=(F(3),F(0));ac=(F(1),F(0));mu=r
add=lambda *xs:tuple(sum(x[i] for x in xs) for i in range(2))
mul=lambda k,x:tuple(k*t for t in x)
physical=(F(0),F(1,2))
terms=[add(physical,mul(3,r),delta),add(physical,mul(2,r),delta,mu),add(physical,mul(2,r),ac,delta),add(physical,mul(4,r),delta,mul(F(-1,2),mu)),add(physical,mul(4,r),delta)]
expected=[(F(6),F(-1)),(F(6),F(-1)),(F(6),F(-1,2)),(F(13,2),F(-5,4)),(F(7),F(-3,2))]
ok('native_five_exact_powers',terms==expected)
for beta in [F(0),F(1),F(3,2),F(19,10),F(2)]:
 leading=F(6)-beta
 ok('native_row_domination',all(a+b*beta>=leading for a,b in terms))
pins={
 '/workspace/shared/next-history-m4-port-20261005/MANIFEST.json':'ac9845847d49f85da1ce7cb9728e18c2f0a6ed8f7da4de1befe685878ca7250f',
 '/workspace/shared/quartic-corrected-mean-join-20261005/MANIFEST.json':'dde7c7e3220b18f65d13039587af6cd9690db0264f10f87c00fcbcb6870ccad9',
 '/workspace/shared/v9-curation-work/frozen/prerequisites/research-source/High Acc Ideas/ai-bucket/30_low_acc.tex':'7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8'}
for path,digest in pins.items():ok('sealed_input_pin',hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest)
result={'status':'PASS','assertions':sum(checks.values()),'checks':checks,'max_noncommuting_hessian_commutator':commutator_max,'scope':'Actual seven-VALUE local source/current diagnostics and exact native scaling; not a true F2/F3 genealogy or native compiler execution.'}
(P/'checks.json').write_text(json.dumps(result,indent=2)+'\n')
(P/'INPUT-PINS.json').write_text(json.dumps(pins,indent=2)+'\n')
print(json.dumps(result,indent=2))
