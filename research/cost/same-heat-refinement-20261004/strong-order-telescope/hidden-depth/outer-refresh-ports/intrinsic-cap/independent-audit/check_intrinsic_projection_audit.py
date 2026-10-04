#!/usr/bin/env python3
"""Independent deterministic linear-algebra audit; not a simulation of the full cap grammar."""
import hashlib, json
from fractions import Fraction as F
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
BASE=HERE.parent
RNG=np.random.default_rng(202610041828)
checks=0; worst=0.; reverses=[]; intermediate=[]; cases=[]
def require(x, label):
 global checks
 checks += 1
 if not x: raise AssertionError(label)
def close(x,y,label,tol=3e-10):
 global worst
 e=float(np.linalg.norm(np.asarray(x)-np.asarray(y)))
 worst=max(worst,e)
 require(e<=tol, f'{label}: {e}')
def op(x): return float(np.linalg.norm(x,ord=2))
def root(x):
 e,v=np.linalg.eigh((x+x.T)/2)
 require(float(e.min())>-1e-12,'positive root gap')
 return (v*np.sqrt(np.maximum(e,0)))@v.T
for n,d in [(5,1),(5,2),(9,1),(9,3),(17,2),(23,4)]:
 for trial in range(6):
  U=np.linalg.qr(RNG.normal(size=(n,d)))[0]
  P=U.T; Q=U@U.T; R=np.eye(n)-Q
  B=RNG.normal(size=(d,n)); B*=.19/op(B)
  M=U@B
  close(P@P.T,np.eye(d),'coisometry')
  close(Q@M,M,'fixed left support')
  J=RNG.normal(size=(d,n)); J=U@(J/op(J)); C=J-J.T
  require(np.linalg.matrix_rank(J,tol=1e-10)<=d,'first rank')
  require(np.linalg.matrix_rank(C,tol=1e-10)<=2*d,'curl rank')
  require(np.linalg.norm(J)<=np.sqrt(d)*op(J)+1e-10,'first HS')
  require(np.linalg.norm(C)<=np.sqrt(2*d)*op(C)+1e-10,'curl HS')
  x=RNG.normal(size=n); cov=np.eye(n)-M@M.T
  close(Q@(M@x),M@x,'conditional pair mean')
  close(Q@cov@Q+R,cov,'conditional pair covariance')
  close(M.T@Q,M.T,'retained incoming cross covariance')
  close(R@root(cov),R,'exact complementary root block')
  close(root(cov)@R,R,'root cross block')
  close(Q@root(cov)-Q,root(cov)-np.eye(n),'supported forward root perturbation')
  for j in range(1,7):
   L=np.linalg.matrix_power(M@M.T,j)
   close(R@L,np.zeros((n,n)),'forward whole root-square support')
   require(np.linalg.norm(L)<=np.sqrt(d)*op(L)+1e-11,'whole root-square HS')
  rev=float(np.linalg.norm(R@M.T)); reverses.append(rev)
  require(rev>1e-4,'reverse source transpose must leak in hostile case')
  revroot=np.linalg.matrix_power(M.T@M,2)
  require(np.linalg.norm(R@revroot)>1e-8,'reverse root-square may leak')
  sym=(M@M+(M@M).T)/2
  s2=.4; Sself=s2*np.eye(n)-sym
  projected_self=Q@Sself@Q+s2*R
  leak=float(np.linalg.norm(projected_self-Sself)); intermediate.append(leak)
  require(leak>1e-5,'individual self target not fixed')
  correction=sym-M@M.T
  close(-sym+correction,-M@M.T,'self orientation cancellation')
  close(Q@(-sym+correction)@Q,-M@M.T,'complete supported covariance defect')
  # A whole conditional polynomial tensor, rather than separate cyclic contractions.
  T=U@RNG.normal(size=(d,n*n))
  require(np.linalg.norm(T)<=np.sqrt(d)*op(T)+1e-10,'singleton fixed-S tensor cut')
  T3=T.reshape(n,n,n)
  close(np.linalg.norm(T3),np.linalg.norm(T3.transpose(1,0,2)),'slot permutation HS invariance')
  require(np.linalg.norm(R@T3.transpose(1,0,2).reshape(n,-1))>1e-4,'slot permutation need not retain range')
  # Literal first, adjoint, and all-private-zero algebra of completed wrapper.
  H=RNG.normal(size=(n,n+3)); y0=RNG.normal(size=n); w=RNG.normal(size=n+3); z=RNG.normal(size=n)
  value=lambda a,b:Q@(y0+H@a)+R@b
  full=np.concatenate([Q@H,R],axis=1)
  direction=RNG.normal(size=2*n+3); h=1e-5
  fd=(value(w+h*direction[:n+3],z+h*direction[n+3:])-value(w,z))/h
  close(fd,full@direction,'literal wrapper first',tol=3e-9)
  output=RNG.normal(size=n)
  close(full.T@output,np.concatenate([H.T@Q@output,R@output]),'literal wrapper adjoint')
  close(value(np.zeros(n+3),np.zeros(n)),Q@y0,'literal wrapped zero')
  close(P@(Q@y0),P@y0,'physical zero preservation')
  require(op(full)<=max(op(H),1)+1e-10,'wrapper full first radius')
  cases.append({'n':n,'d':d,'trial':trial,'reverse_leak':rev,'individual_self_change':leak})
# Exact working cap versus ambient loss; this only checks arithmetic, not 269 source rows.
B=F(52,5); u=F(3,4)
require(B*u==F(39,5),'working cap exponent')
require(F(83,5)-2<F(16),'ambient c=2 obstruction')
require(F(83,5)>16,'physical cap row restored')
require(F(1607,100)>16 and F(141,10)>14 and F(22603,1250)>18,'reported three minima strict')
# Public adaptation: verify the bundled baseline, retain unbundled pins as
# explicitly unverified external inputs. No finite algebra check is removed.
research_root = next(p for p in Path(__file__).resolve().parents if p.name == 'research')
source_specs=[('exact-slack/exact32-proof-v2.tex','c7487ca24ffa44ea73b52e42e4b6f715ad964171fa5bd0a61cf9323eb646cbff')]
external_pins={'LOW30':'7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8','LOW31':'3de62d349338304cce43af4e49c8862e8c260b1f5816234d63d6abfbc1d004ba'}
pins={}
for path,expected in source_specs:
 digest=hashlib.sha256((research_root/path).read_bytes()).hexdigest(); pins[path]=digest
 require(digest==expected,'pinned native source')
result={'status':'PASS','checks':checks,'cases':len(cases),'max_equality_residual':worst,'min_reverse_transpose_leak':min(reverses),'min_individual_self_target_change':min(intermediate),'source_pins':pins,'external_source_pins_not_reverified':external_pins,'candidate_sha256':hashlib.sha256((BASE/'INTRINSIC-MARKED-CAP.md').read_bytes()).hexdigest(),'limitations':['No full 10849-level family is simulated.','The complete local-type ledger and proper-cut contracts are proof inputs; these tests do not establish them.','The Exact32 three minima are arithmetic checks of the stated numbers, not a fresh 269-row census.'],'case_details':cases}
(HERE/'intrinsic_projection_audit_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('case_details','source_pins')},indent=2))
