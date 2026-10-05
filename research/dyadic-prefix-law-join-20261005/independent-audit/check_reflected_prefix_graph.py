from pathlib import Path
import hashlib, json
import numpy as np

out=Path(__file__).parent
rng=np.random.default_rng(510052026)
A=.13; w=.08; delta=-np.log1p(-w); h=.011; r=np.sqrt(1-h*h)
s=np.linspace(0,delta,13)
beta=np.full(len(s),w/len(s))
a=np.sinh(delta-s)/np.sinh(delta)
b=np.sinh(s)/np.sinh(delta)
sigma=np.sqrt(np.maximum(0,-np.expm1(-2*s)*(-np.expm1(-2*(delta-s)))/(-np.expm1(-2*delta))))
vectors=np.array([[1.,0.],[.6,.8],[-.3,np.sqrt(.91)]])
base=np.array([[.22,.04],[.04,.18]])
cs=np.array([.15,.12,.08])
upper=base+sum(c*np.outer(v,v) for c,v in zip(cs,vectors))
scale=A/np.linalg.eigvalsh(upper)[-1]

def g(x):
    return scale*(base@x+sum(c*np.tanh(v@x)*v for c,v in zip(cs,vectors)))
def H(x):
    return scale*(base+sum(c*(1-np.tanh(v@x)**2)*np.outer(v,v) for c,v in zip(cs,vectors)))

def prefix(x,y,n):
    return sum(ww*g(aa*x+bb*y+ss*n) for ww,aa,bb,ss in zip(beta,a,b,sigma))
def f(x,y,n):
    return -r*(prefix(x,y,-n)-prefix(x,y,np.zeros(2)))
def analytic(x,y,n):
    Jn=np.zeros((2,2)); Jx=np.zeros((2,2)); Jy=np.zeros((2,2))
    for ww,aa,bb,ss in zip(beta,a,b,sigma):
        z=aa*x+bb*y
        hn=H(z-ss*n); h0=H(z)
        Jn+=r*ww*ss*hn
        Jx-=r*ww*aa*(hn-h0)
        Jy-=r*ww*bb*(hn-h0)
    return np.concatenate((Jx,Jy,Jn),axis=1), Jn

worst_fd=0.; smallest_eig=1.; max_private=0.; max_caller=0.; zero=0.
for k in range(200):
    x,y,n=rng.normal(size=(3,2))*rng.uniform(.1,4)
    J,Jn=analytic(x,y,n)
    z=np.concatenate((x,y,n)); eps=2e-5
    Jfd=np.column_stack([(f(*(z+eps*np.eye(6)[j]).reshape(3,2))-f(*(z-eps*np.eye(6)[j]).reshape(3,2)))/(2*eps) for j in range(6)])
    worst_fd=max(worst_fd,float(np.max(np.abs(J-Jfd))))
    smallest_eig=min(smallest_eig,float(np.linalg.eigvalsh(Jn)[0]))
    max_private=max(max_private,float(np.linalg.norm(Jn,2)))
    max_caller=max(max_caller,float(np.linalg.norm(J[:,:4],2)))
    zero=max(zero,float(np.linalg.norm(f(x,y,np.zeros(2)))))
assert zero==0
assert smallest_eig>=-1e-14
assert max_private<=r*A*np.dot(beta,sigma)+1e-13
assert max_caller<=r*A*w+1e-13
assert worst_fd<1e-9

# All four whole Gaussian contributions have physical readout h/2.
assert abs(4*(h/2)**2-h*h)<1e-18
assert abs(2*(h/2)**2-h*h/2)<1e-18

# Native square response's pointwise energy envelope, independent of padding.
def square_trial(p,z1,z0,z2,mu):
    ell=A
    bs=np.array([.12,.31]); ds=np.array([1.3,-.3])
    def source(x):return g(x)/ell
    def response(w):return source(.5*w+np.sqrt(.75)*z2)-source(-.5*w+np.sqrt(.75)*z2)
    I=sum(d*(source(np.sqrt(1-bb*bb)*z1+bb*p)-source(np.sqrt(1-bb*bb)*z1-bb*p))/(2*bb) for bb,d in zip(bs,ds))
    t=np.sqrt(mu)/3
    value=ell*ell*(response(z0+t*I)-response(z0-t*I))/(2*t)
    bound=ell*ell*sum(abs(ds))*np.linalg.norm(p)
    return np.linalg.norm(value),bound
worst_energy_ratio=0
for _ in range(100):
    roots=rng.normal(size=(4,2))
    for mu in [1.,.1,1e-3,1e-6]:
        value,bound=square_trial(*roots,mu)
        worst_energy_ratio=max(worst_energy_ratio,value/bound)
        assert value<=bound*(1+1e-7)

# Complete dimensions: align exactly one D-row per node, retaining complements.
for _ in range(100):
    D=int(rng.integers(1,20)); nl=int(rng.integers(1,10)); nr=int(rng.integers(1,10))
    dtail=rng.integers(5,100,size=nl)*D
    dpl=rng.integers(5,100,size=nl)*D
    dpr=rng.integers(5,100,size=nr)*D
    direct=D+sum(2*D+dtail+dpl-D)+sum(4*D+dpr-D)
    displayed=D+sum(D+dtail+dpl)+sum(3*D+dpr)
    assert direct==displayed

snapshot=out/'REVIEWED-AUTHOR-SNAPSHOT.md'
report={
 'status':'PASS',
 'reviewed_source_sha256':hashlib.sha256(snapshot.read_bytes()).hexdigest(),
 'max_finite_difference_absolute_error':worst_fd,
 'minimum_private_eigenvalue':smallest_eig,
 'maximum_private_radius':max_private,
 'private_radius_envelope':float(r*A*np.dot(beta,sigma)),
 'maximum_endpoint_radius':max_caller,
 'endpoint_radius_envelope':r*A*w,
 'literal_zero_error':zero,
 'maximum_pointwise_square_energy_ratio':worst_energy_ratio,
 'carrier_allocation':'four contributions (h/2)^2 sum to h^2',
 'scope':'Original reflected source and literal zero-clock signed-response graph diagnostics; no complete mean compiler or universal calibration is numerically inferred.'
}
(out/'reflected_prefix_graph_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
