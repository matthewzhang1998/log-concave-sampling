import math, json, hashlib
from pathlib import Path
import numpy as np
import mpmath as mp
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite import hermgauss

OUT=Path(__file__).resolve().parent
checks=0
metrics={}
def check(x,msg):
    global checks
    if not x: raise AssertionError(msg)
    checks+=1

def Kmat(delta,Delta,z):
    L=delta-Delta; q=np.exp(-delta)
    S=np.sqrt(-np.expm1(-2*delta))
    return np.array([[np.exp(-z),2*q*np.sinh(z)/S],
                     [0,np.exp(z-L)*np.sqrt(-np.expm1(-2*Delta))/S]],dtype=complex)

max_det=0.; max_norm=0.
for delta in [1e-6,1e-3,.02,.2,math.log(2)]:
  for k in [1,2,4,9,20]:
    Delta=delta/(2**k); L=delta-Delta
    for theta in [.001,.02,.1,.3,.5,.7,.99,.999]:
      x=L*theta
      for fac in [0,.25,.75,1]:
        y=fac*math.sqrt(x*(L-x))
        M=Kmat(delta,Delta,x+1j*y)
        exact=((1-np.exp(-2*x))*(1-np.exp(-2*(L-x)))
               -4*np.exp(-2*delta)*np.sin(y)**2)/(1-np.exp(-2*delta))
        det=np.linalg.det(np.eye(2)-M@M.conj().T)
        max_det=max(max_det,abs(det-exact))
        excess=np.linalg.norm(M,2)-1
        max_norm=max(max_norm,excess)
        check(abs(det-exact)<3e-10,'complex determinant')
        check(excess<3e-10,'complex contraction')
metrics['complex_det_max_error']=max_det
metrics['complex_contraction_max_excess']=max_norm

# Taylor jet for geometric moments; the implementation does not enumerate atoms.
def geometric_moments(delta,Delta,a,b,order):
    n=b-a+1
    if n==1: return [mp.mpf(1)]+[mp.mpf(0)]*order,mp.exp(-2*Delta*a)
    q=mp.exp(-2*Delta)
    den=[1-q]+[-q*(mp.mpf(1)/(n-1))**j/mp.factorial(j) for j in range(1,order+1)]
    num=[1-q**n]+[-q**n*(mp.mpf(n)/(n-1))**j/mp.factorial(j) for j in range(1,order+1)]
    coeff=[]
    for j in range(order+1):
        coeff.append((num[j]-sum(den[i]*coeff[j-i] for i in range(1,j+1)))/den[0])
    mass=mp.exp(-2*Delta*a)*coeff[0]
    return [coeff[j]*mp.factorial(j)/coeff[0] for j in range(order+1)],mass

def gauss_panel(delta,Delta,a,b,m):
    n=b-a+1
    if n<=m:
      return [mp.mpf(j)*Delta for j in range(a,b+1)],[mp.exp(-2*j*Delta) for j in range(a,b+1)]
    moment,mass=geometric_moments(delta,Delta,a,b,2*m)
    H=mp.matrix(m);H1=mp.matrix(m)
    for i in range(m):
      for j in range(m): H[i,j]=moment[i+j]; H1[i,j]=moment[i+j+1]
    L=mp.cholesky(H); Li=L**-1
    J=Li*H1*Li.T; J=(J+J.T)/2
    nodes,V=mp.eigsy(J)
    weights=[mass*V[0,j]**2 for j in range(m)]
    check(all(0<x<1 for x in nodes),'Gauss nodes inside panel')
    check(all(w>0 for w in weights),'positive Gauss weights')
    normw=[w/mass for w in weights]
    for r in range(2*m):
      err=abs(sum(normw[i]*nodes[i]**r for i in range(m))-moment[r])
      check(err<mp.mpf('1e-45'),'Gauss exact moments')
    return [(a+(b-a)*nodes[i])*Delta for i in range(m)],weights

mp.mp.dps=200
max_moment=mp.mpf(0)
for delta0,N,a,b,m in [(.3,128,8,15,4),(.7,512,32,127,6),(.001,1024,128,255,8)]:
    delta=mp.mpf(str(delta0));Delta=delta/N
    moments,mass=geometric_moments(delta,Delta,a,b,2*m)
    ws=[mp.exp(-2*j*Delta) for j in range(a,b+1)]
    for r in range(2*m+1):
      direct=sum(ws[j-a]*(mp.mpf(j-a)/(b-a))**r for j in range(a,b+1))/sum(ws)
      max_moment=max(max_moment,abs(direct-moments[r]))
      check(abs(direct-moments[r])<mp.mpf('1e-100'),'closed-form normalized moments')
    gauss_panel(delta,Delta,a,b,m)
# This large panel is never enumerated.
delta=mp.mpf('.2'); N=2**22; Delta=delta/N
large_nodes,large_weights=gauss_panel(delta,Delta,2**19,2**20-1,8)
check(len(large_nodes)==8,'large lattice compressed to 8 nodes')
metrics['geometric_moment_max_error']=str(max_moment)
metrics['large_scalar_panel_atoms']=2**19
metrics['large_scalar_panel_nodes']=8

def position_rule(delta,k,m):
    N=2**k;Delta=delta/N
    if N==1:return [mp.mpf(0)],[mp.mpf(1)]
    nodes=[mp.mpf(0),(N-1)*Delta];weights=[mp.mpf(1),mp.exp(-2*(N-1)*Delta)]
    for j in range(k-1):
      a=2**j;b=2**(j+1)-1
      for aa,bb in [(a,b),(N-1-b,N-1-a)]:
        x,w=gauss_panel(delta,Delta,aa,bb,m);nodes+=x;weights+=w
    check(abs(sum(weights)-(1-mp.exp(-2*delta))/(1-mp.exp(-2*Delta)))<mp.mpf('1e-100'),'position exact mass')
    return nodes,weights

def sympower(K,n):
    # Rows: normalized input Hermite monomials. Columns: output monomials.
    a,b,c=float(K[0,0].real),float(K[0,1].real),float(K[1,1].real)
    T=np.zeros((n+1,n+1))
    for j in range(n+1):
      for l in range(n-j+1):
        k=j+l
        T[j,k]=(math.sqrt(math.comb(n,j)/math.comb(n,k))*math.comb(n-j,l)
                  *a**(n-j-l)*b**l*c**j)
    return T

max_operator=0.;max_scaled=0.
for delta0,k,m in [(.02,4,2),(.3,5,3),(.69,6,4),(.001,6,3)]:
    delta=mp.mpf(str(delta0));N=2**k;Delta=delta/N
    nodes,weights=position_rule(delta,k,m)
    for n in range(33):
      exact=sum(math.exp(-2*j*float(Delta))*sympower(Kmat(float(delta),float(Delta),j*float(Delta)),n) for j in range(N))
      approx=sum(float(w)*sympower(Kmat(float(delta),float(Delta),float(x)),n) for x,w in zip(nodes,weights))
      err=np.linalg.norm(exact-approx,2)
      max_operator=max(max_operator,err);max_scaled=max(max_scaled,err/(N*4**(-m)))
      check(err<=10*N*4**(-m),'finite-chaos operator quadrature bound')
metrics['position_operator_max_error']=max_operator
metrics['position_operator_max_scaled_error']=max_scaled

# Full matrix-linear scalar coefficient calibration of all dyadic levels.
def klin(delta):
    return -mp.expm1(-2*delta)/4-delta**2/mp.expm1(2*delta)
def midpointcoef(delta):
    d=delta/2
    return mp.tanh(d)*(d*mp.exp(-d))**2
max_linear=mp.mpf(0)
for delta in [mp.mpf('1e-8'),mp.mpf('.001'),mp.mpf('.2'),mp.log(2)]:
  for K in [1,2,5,10,20]:
    total=mp.mpf(0)
    for k in range(K):
      Delta=delta/(2**k)
      mass=-mp.expm1(-2*delta)/(-mp.expm1(-2*Delta))
      total+=mass*midpointcoef(Delta)
    Delta=delta/(2**K)
    total+=(-mp.expm1(-2*delta)/(-mp.expm1(-2*Delta)))*klin(Delta)
    rel=abs(total-klin(delta))/klin(delta)
    max_linear=max(max_linear,rel)
    check(rel<mp.mpf('1e-100'),'exact linear dyadic covariance calibration')
metrics['linear_covariance_max_relative_error']=str(max_linear)

# Nonlinear scalar midpoint covariance versus the positive OU square identity.
lx,lw=leggauss(100); gh,gw=hermgauss(100); gh=math.sqrt(2)*gh;gw=gw/math.sqrt(math.pi)
rx,rw=leggauss(64); rs=(rx+1)/2; rws=rw*(rs) # 2r dr
max_nonlinear=0.
for Delta,a,b,frequency in [(.02,.3,-.2,1.),(.1,.4,-.5,5.),(.3,1.1,.3,8.),(.6,-.1,.9,12.)]:
    d=Delta/2;alpha=1/(2*math.cosh(d));sigma=math.sqrt(math.tanh(d))
    us=np.r_[(lx+1)*d/2,d+(lx+1)*d/2]
    uw=np.r_[lw*d/2,lw*d/2]
    psi=np.where(us<=d,np.sinh(us)/np.sinh(d),np.sinh(Delta-us)/np.sinh(d))
    v=np.where(us<=d,us,us-d)
    kap2=(-np.expm1(-2*v))*(-np.expm1(-2*(d-v)))/(-np.expm1(-2*d))
    ell=np.where(us<=d,a*np.sinh(d-us)/np.sinh(d),b*np.sinh(us-d)/np.sinh(d))+alpha*(a+b)*psi
    mean=ell[:,None]+sigma*psi[:,None]*gh[None,:]
    gm=.5*(mean+np.exp(-.5*frequency**2*kap2[:,None])*np.sin(frequency*mean)/frequency)
    mv=np.sum((uw*np.exp(-us))[:,None]*gm,axis=0)
    var=float(np.dot(gw,mv*mv)-np.dot(gw,mv)**2)
    square=0.
    for r,wr in zip(rs,rws):
      tau2=kap2+sigma**2*psi**2*(1-r*r)
      z=ell[:,None]+r*sigma*psi[:,None]*gh[None,:]
      H=.5*(1+np.exp(-.5*frequency**2*tau2[:,None])*np.cos(frequency*z))
      jr=sigma*np.sum((uw*np.exp(-us)*psi)[:,None]*H,axis=0)
      square+=wr*np.dot(gw,jr*jr)
    err=abs(square-var)/max(var,1e-30)
    max_nonlinear=max(max_nonlinear,err)
    check(err<2e-7,'nonlinear Gaussian midpoint square identity')
metrics['nonlinear_square_max_relative_error']=max_nonlinear

# Direct source envelopes against increasingly endpoint-resolved positive rules.
max_envelope=0.
for Delta in [1e-5,.001,.05,.5]:
    d=Delta/2
    for n in [8,32,128]:
      x,w=leggauss(n)
      u=(x+1)*d/2; ww=w*d/2
      kap=np.sqrt((-np.expm1(-2*u))*(-np.expm1(-2*(d-u)))/(-np.expm1(-2*d)))
      gamma1=ww*np.exp(-u)*np.sinh(u)/np.sinh(d)
      gamma2=ww*np.exp(-(d+u))*np.sinh(d-u)/np.sinh(d)
      env=float(np.sum((gamma1+gamma2)/kap)/math.sqrt(Delta))
      max_envelope=max(max_envelope,env)
      check(env<5,'inner endpoint first envelope')
metrics['inner_caller_envelope_max']=max_envelope

result={'status':'PASS','assertions':checks,'metrics':metrics,
        'scope':'Scalar geometry, finite-chaos position contractions/quadrature, moment-jet positive Gauss setup, exact linear dyadic decomposition, nonlinear scalar midpoint-square identity, and source-envelope checks. Imported native square action/mean LAW not executed.'}
(OUT/'dyadic_bridge_square_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
