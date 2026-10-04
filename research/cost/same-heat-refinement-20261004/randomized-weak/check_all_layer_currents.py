import json, math
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad

g,w=leggauss(150)
du=g*.13
dw=w/2
t=(g+1)/2
tw=w/2
a,c,eta,q=.05,.2,.07,.27
def b(x): return a*(x/2+c*(np.sqrt(x*x+eta*eta)-eta))
def bp(x): return a*(.5+c*x/np.sqrt(x*x+eta*eta))
def F(x): return -b(x)
def Fp(x): return -bp(x)
def phi(x): return np.cos(3*x)+x+x*x
def phip(x): return -3*np.sin(3*x)+1+2*x
def phipp(x): return -9*np.cos(3*x)+2
direct=float(dw@(phi(F(q+du))-phi(F(q))))
vcur=0.; qcur=0.
for ti,wi in zip(t,tw):
    yt=F(q+ti*du)
    vcur+=wi*float(dw@((Fp(q+ti*du)-Fp(q))*du*phip(yt)))
    qcur+=wi*(1-ti)*float(dw@(Fp(q)*Fp(q+ti*du)*du*du*phipp(yt)))
assert abs(direct-vcur-qcur)<1e-12

kink=[]
for ratio in [1e-1,1e-2,1e-3,1e-4]:
    eps=.01; et=ratio*eps
    mean_abs_smooth=.5*(math.sqrt(eps*eps+et*et)+et*et/eps*math.asinh(eps/et))
    norm=c*(mean_abs_smooth-et)/eps
    kink.append({'eta_over_epsilon':ratio,'mean_defect_over_a_epsilon':norm})
assert abs(kink[-1]['mean_defect_over_a_epsilon']-c/2)<3e-5

times=np.array([.35,.8,1.2]); nodes=.6+.1*g; nodeweights=w/2; h=.4; lam=.5
paths=[]
for aa in [.1,.03,.01,.003,.001]:
    A=[]; B=[]
    for s in nodes:
        k=np.where(times>s,np.sin(times-s),0.)
        A.append(np.cos(times)-aa*lam*h*k*math.cos(s))
        B.append(np.sin(times)-aa*lam*h*k*math.sin(s))
    A=np.array(A); B=np.array(B); abar=nodeweights@A; bbar=nodeweights@B
    normal=np.cross(abar,bbar); normal/=np.linalg.norm(normal)
    vx=1/(1+aa*lam)
    gram=vx*np.outer(abar,abar)+np.outer(bbar,bbar)
    sigma=sum(ww*(vx*np.outer(x-abar,x-abar)+np.outer(z-bbar,z-bbar)) for ww,x,z in zip(nodeweights,A,B))
    base=float(normal@gram@normal); trans=float(normal@sigma@normal)
    assert abs(base)<1e-14 and trans>0
    assert float(normal@(gram-sigma)@normal)<0
    paths.append({'a':aa,'base_normal_variance':base,
                  'clock_transverse_variance':trans,
                  'joint_w2_lower_bound':math.sqrt(trans),
                  'lower_bound_over_a':math.sqrt(trans)/aa})

reserve=[]
for ratio in [.1,.03,.01,.003,.001]:
    # Adaptive quadrature resolves the small smooth kink around zero.
    mean=math.sqrt(2/math.pi)*quad(lambda x:(math.sqrt(x*x+ratio*ratio)-ratio)*math.exp(-x*x/2),0,math.inf,epsabs=1e-12,epsrel=1e-12)[0]
    reserve.append({'eta_over_tau':ratio,'mean_bias_over_a_c_tau':mean})
assert abs(reserve[-1]['mean_bias_over_a_c_tau']-math.sqrt(2/math.pi))<.005

out={'status':'PASS: finite current identity, fixed-caller feedback, and zero-Gram path fixture',
     'current_identity':{'direct':direct,'rank_one_feedback':vcur,'rank_two':qcur,
                         'residual':direct-vcur-qcur},
     'fixed_caller_feedback':kink,'retained_path_gap':paths,'intermediate_reserve_mean':reserve,
     'scope':'No posterior-law or arbitrary-algorithm lower bound is asserted.'}
Path(__file__).with_name('all_layer_current_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
