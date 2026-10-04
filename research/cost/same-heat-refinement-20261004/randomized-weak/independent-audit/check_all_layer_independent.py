"""Independent fixtures for paired currents and a smooth capped path-gap witness."""
from pathlib import Path
import json, math, hashlib
from fractions import Fraction
import numpy as np
from scipy.special import roots_legendre,roots_hermitenorm,ndtr
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parent
checks=0
def check(test,msg):
 global checks
 checks+=1
 assert bool(test),msg

for name,pin in [('ALL-LAYER-REVIEWED-dbf2b9b04c72.md','dbf2b9b04c72cd83970273cb487b2ac800ff9358dec945e1ef489f054b829d6b'),('RESERVE-REVIEWED-02db073eb1dc.md','02db073eb1dca83c37ac8a8924cc34798550bb02f572e427a2053d053b2037c1')]:
 check(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==pin,'reviewed note snapshot pin')

g,w=roots_legendre(72); t=(g+1)/2; tw=w/2
u,uw=roots_legendre(20); uw=uw/2
u1,u2=np.meshgrid(u,u); probs=(uw[:,None]*uw[None,:]).ravel()
C=np.array([[.16,.07],[-.09,.13],[.05,-.11]])
D=np.stack((u1.ravel(),u2.ravel()),axis=1)@C.T
q=np.array([.2,-.17,.4])
check(np.linalg.norm(probs@D)<1e-15,'centered smooth uniform-source vector')
def F(x):
 return np.column_stack((.3*x[:,0]+x[:,1]**2+.2*x[:,0]*x[:,2],np.tanh(x[:,0]-.4*x[:,2])+.2*x[:,1]))
def DF(x):
 J=np.zeros((len(x),2,3)); J[:,0,0]=.3+.2*x[:,2]; J[:,0,1]=2*x[:,1]; J[:,0,2]=.2*x[:,0]
 e=1-np.tanh(x[:,0]-.4*x[:,2])**2
 J[:,1,0]=e;J[:,1,1]=.2;J[:,1,2]=-.4*e
 return J
def phi(y): return np.cos(y[:,0]+.7*y[:,1])+.2*y[:,0]*y[:,1]+.4*y[:,0]**2
v=np.array([1.,.7]); H0=np.array([[.8,.2],[.2,0.]])
def grad(y): return -np.sin(y@v)[:,None]*v+y@H0
def Hess(y): return -np.cos(y@v)[:,None,None]*np.outer(v,v)+H0
qq=np.repeat(q[None,:],len(D),axis=0); f0=F(qq)
baseD=np.einsum('nij,nj->ni',DF(qq),D)
direct=float(probs@(phi(F(qq+D))-phi(f0))); rank1=rank2=0.
for ti,wi in zip(t,tw):
 x=qq+ti*D; y=F(x); jd=np.einsum('nij,nj->ni',DF(x),D)
 rank1+=wi*float(probs@np.sum((jd-baseD)*grad(y),axis=1))
 rank2+=wi*(1-ti)*float(probs@np.einsum('ni,nij,nj->n',baseD,Hess(y),jd))
check(abs(direct-rank1-rank2)<2e-13,'vector paired-current identity and orientation')
current=dict(direct=direct,rank_one=rank1,rank_two=rank2,error=direct-rank1-rank2)

def cap(v,B=1.2):
 v=np.asarray(v); az=np.abs(v);tt=np.clip((az-B)/B,0,1)
 values=np.where(az<=B,v,np.sign(v)*B*(1+tt-tt*tt/2))
 first=np.where(az<=B,1.,np.where(az>=2*B,0.,1-tt))
 return values,first

# S=s_-+(s_+-s_-)Phi(G) is uniform within a fixed interval away from every kink.
times=np.array([.35,.8,1.2]); s0,s1=.5,.7;h=.4;lam=.5
S=s0+(s1-s0)*t
GH,GHw=roots_hermitenorm(100);GHw/=math.sqrt(2*math.pi)
path=[]
for a in (.1,.03,.01,.003,.001):
 vX=1/(1+a*lam)
 kap=np.where(times[None,:]>S[:,None],np.sin(times[None,:]-S[:,None]),0.)
 A=np.cos(times)[None,:]-a*lam*h*kap*np.cos(S[:,None])
 B=np.sin(times)[None,:]-a*lam*h*kap*np.sin(S[:,None])
 abar=tw@A;bbar=tw@B
 normal=np.cross(abar,bbar);normal/=np.linalg.norm(normal)
 varX=float(GHw@(cap(math.sqrt(vX)*GH)[0]**2));varZ=float(GHw@(cap(GH)[0]**2))
 gram=varX*np.outer(abar,abar)+varZ*np.outer(bbar,bbar)
 dA=A-abar;dB=B-bbar
 cov=np.einsum('n,ni,nj->ij',tw,varX*dA,dA)+np.einsum('n,ni,nj->ij',tw,varZ*dB,dB)
 transverse=float(normal@cov@normal)
 check(abs(float(normal@gram@normal))<1e-13,'capped smooth-clock conditional mean plane')
 check(transverse>0,'smooth-clock transverse variance positive')
 check(float(normal@(gram-cov)@normal)<0,'unavailable covariance subtraction reserve')
 check(.002<math.sqrt(transverse)/a<.01,'strong path gap coefficient')
 # Actual directional first includes cap roots and smooth clock motion.
 inp=np.array([.8,1.6,-.3]);direction=np.array([.7,-.2,.4])
 def evaluate(inp,tangent=None):
  X,Xp=cap(math.sqrt(vX)*inp[0]);Z,Zp=cap(inp[1]);ss=s0+(s1-s0)*ndtr(inp[2])
  kk=np.where(times>ss,np.sin(times-ss),0.);kp=np.where(times>ss,-np.cos(times-ss),0.)
  source=math.cos(ss)*X+math.sin(ss)*Z
  val=np.cos(times)*X+np.sin(times)*Z-a*lam*h*kk*source
  if tangent is None:return val
  dX=Xp*math.sqrt(vX)*tangent[0];dZ=Zp*tangent[1]
  ds=(s1-s0)*math.exp(-inp[2]**2/2)/math.sqrt(2*math.pi)*tangent[2]
  dsource=math.cos(ss)*dX+math.sin(ss)*dZ+(-math.sin(ss)*X+math.cos(ss)*Z)*ds
  der=np.cos(times)*dX+np.sin(times)*dZ-a*lam*h*(kp*ds*source+kk*dsource)
  return val,der
 val,der=evaluate(inp,direction);eps=1e-6
 fd=(evaluate(inp+eps*direction)-evaluate(inp-eps*direction))/(2*eps)
 check(np.max(np.abs(fd-der))<1e-9,'smooth capped clock complete directional first')
 path.append(dict(a=a,transverse_variance=transverse,joint_W2_lower=math.sqrt(transverse),lower_over_a=math.sqrt(transverse)/a,first_error=float(np.max(np.abs(fd-der)))))

# N-stratum extension of the transverse fixture at fixed mean coefficient plane.
stratified_path=[]
a=.03; varX=float(GHw@(cap(math.sqrt(1/(1+a*lam))*GH)[0]**2));varZ=float(GHw@(cap(GH)[0]**2))
for N in (1,2,4,8,16,32,64):
 nodes=s0+(s1-s0)*(np.arange(N)[:,None]+t[None,:])/N
 kern=np.where(times[None,None,:]>nodes[:,:,None],np.sin(times[None,None,:]-nodes[:,:,None]),0.)
 CX=kern*np.cos(nodes[:,:,None]);CZ=kern*np.sin(nodes[:,:,None])
 meansX=np.einsum('p,npk->nk',tw,CX);meansZ=np.einsum('p,npk->nk',tw,CZ)
 abar=np.cos(times)-a*lam*h*np.mean(meansX,axis=0)
 bbar=np.sin(times)-a*lam*h*np.mean(meansZ,axis=0)
 normal=np.cross(abar,bbar);normal/=np.linalg.norm(normal)
 rx=CX-meansX[:,None,:];rz=CZ-meansZ[:,None,:]
 var=(a*lam*h/N)**2*(varX*np.einsum('p,np->',tw,(rx@normal)**2)+varZ*np.einsum('p,np->',tw,(rz@normal)**2))
 scaled=float(var*N**3/a**2)
 check(1e-5<scaled<3e-5,'N-stratum transverse variance has a2 N^-3 scale')
 stratified_path.append(dict(N=N,variance=float(var),variance_times_N3_over_a2=scaled))

# Exact power arithmetic for a hypothetical improved reserve law price tau^p.
reserve_exponents=[]
for R in (Fraction(5,2),Fraction(7,2),Fraction(11,2),Fraction(15,2)):
 for depth in (1,2,3):
  for power in (1,2,3):
   beta=max(Fraction(0),Fraction(2,3*power)*(R-(depth+Fraction(1,2))-power))
   # tau scales like the residual a*N^-3/2 in this covariance-gap allocation.
   tau_power=1+Fraction(3,2)*beta
   check(Fraction(1,2)+depth+power*tau_power>=R,'reserve formal target exponent')
   reserve_exponents.append(dict(R=str(R),remaining_force_depth=depth,power=power,count_exponent=str(beta)))

# A conditional buffer cost certificate cannot beat its original strong scale.
reserve=[]
for ratio in (.01,.1,.3,1.,3.,10.,100.):
 cost=ratio+1/ratio
 check(cost>=2,'reserve AM-GM cost envelope')
 reserve.append(dict(tau_over_e=ratio,total_cost_over_e=cost))
# Near-kink Gaussian reserve has a conditional first-order mean price.
GH,GHw=roots_hermitenorm(240);GHw/=math.sqrt(2*math.pi)
reserve_bias=[]
for ratio in (.1,.01,.001):
 eta=.2
 scaled=eta*math.sqrt(2/math.pi)*quad(lambda x:(math.sqrt(x*x+ratio*ratio)-ratio)*math.exp(-x*x/2),0,math.inf,epsabs=1e-12,epsrel=1e-12)[0]
 check(.13<scaled<.18,'conditional reserve first-order mean defect')
 reserve_bias.append(dict(smoothing_over_tau=ratio,defect_over_a_tau=scaled))
out=dict(status='PASS: vector paired-current identity, smooth capped Gaussian-clock path-gap fixture, and conditional reserve cost envelope',check_count=checks,vector_current=current,smooth_capped_path_gap=path,stratified_path_gap=stratified_path,reserve_exponents=reserve_exponents,reserve_envelope=reserve,reserve_bias=reserve_bias,scope='No posterior-law or general-algorithm lower bound. The reserve envelope audits the stated conditional C2 certificate, not all possible reserve constructions.',script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(ROOT/'all_layer_independent_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
