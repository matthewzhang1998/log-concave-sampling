#!/usr/bin/env python3
"""Independent finite diagnostics. Does not execute LOW30 mean compiler."""
from pathlib import Path
from fractions import Fraction
import json
import numpy as np

rng=np.random.default_rng(5102026)
checks=0
maxima={k:0. for k in ['E_mark_ratio','Delta_mark_ratio','E_first_ratio','Delta_first_ratio','full_curl_ratio','current_first_ratio','joint_first_ratio','jacobian_fd_error','current_identity_fd_error','noncommutation']}
def check(x,msg):
    global checks
    checks+=1
    if not bool(x): raise AssertionError(msg)
def norm(x): return np.linalg.norm(x,2)
def vertex(A,a,b,c,C,D,d,O,Y,k=1.):
    calls=0
    v=np.array([1.,0.]); w=np.array([.6,.8]); eye=np.eye(2)
    def query(x):
        nonlocal calls
        calls+=1
        g=A*(.45*x+.2*np.sin(k*(v@x))/k*v+.12*np.sin(k*(w@x))/k*w)
        J=A*(.45*eye+.2*np.cos(k*(v@x))*np.outer(v,v)+.12*np.cos(k*(w@x))*np.outer(w,w))
        return g,J
    x=[C[i]@O+D[i]@Y+d[i] for i in range(4)]
    T=[np.hstack([C[i],D[i]]) for i in range(4)]
    h0,J0=query(x[0]); H0=J0@T[0]
    h1,J1=query(x[1]-a*h0); H1=J1@(T[1]-a*H0)
    h1b,J1b=query(x[1]); H1b=J1b@T[1]
    U,J2=query(x[2]-b*h1); DU=J2@(T[2]-b*H1)
    V,J2b=query(x[2]-b*h1b); DV=J2b@(T[2]-b*H1b)
    rootU,J3=query(x[3]-c*U); rootV,J3b=query(x[3]-c*V)
    E=rootU-rootV; DE=J3@(T[3]-c*DU)-J3b@(T[3]-c*DV)
    return dict(E=E,Delta=U-V,U=U,V=V,DE=DE,DDelta=DU-DV,DU=DU,DV=DV,
                x=x,J3=J3,J3b=J3b,calls=calls,noncomm=norm(J3@J3b-J3b@J3))

Dphys=2; blocks=3; n=blocks*Dphys
for trial in range(400):
    A=float(rng.uniform(.01,.5)); a,b,c=rng.uniform(0,1,3)
    if trial%19==0: a=0.
    if trial%23==0: b=0.
    if trial%29==0: c=0.
    rows=rng.normal(size=(4,blocks))
    rows=rows/np.linalg.norm(rows,axis=1)[:,None]*rng.uniform(.1,1,(4,1))
    C=[np.kron(rows[i:i+1],np.eye(2)) for i in range(4)]
    sigma=np.linalg.norm(rows[3]); P=C[3]/sigma
    D=[rng.uniform(-1,1)*np.eye(2) for _ in range(4)]
    d=rng.normal(size=(4,2)); O=rng.normal(size=n); Y=3*rng.normal(size=2)
    k=float(10**rng.uniform(-1,2))
    out=vertex(A,a,b,c,C,D,d,O,Y,k)
    E,Delta,DE,DD=out['E'],out['Delta'],out['DE'],out['DDelta']
    LU=A*(1+b*A*(1+a*A)); LV=A*(1+b*A)
    check(out['calls']==7,'exact 7 VALUE queries')
    check(norm(P@P.T-np.eye(2))<1e-13,'whole-bank coisometry')
    boundD=a*b*A**3*norm(out['x'][0]); boundE=a*b*c*A**4*norm(out['x'][0])
    check(norm(Delta)<=boundD+2e-13,'pointwise Delta mark')
    check(norm(E)<=boundE+2e-13,'pointwise E mark')
    for section in [slice(0,n),slice(n,None)]:
        check(norm(DE[:,section])<=21*A/8+2e-13,'complete/caller E first')
        check(norm(DD[:,section])<=9*A/4+2e-13,'complete/caller Delta first')
    df=P.T@DE[:,:n]
    pure=sigma*P.T@(out['J3']-out['J3b'])@P
    rest=-c*P.T@out['J3']@out['DU'][:,:n]+c*P.T@out['J3b']@out['DV'][:,:n]
    check(norm(df-pure-rest)<2e-13,'full derivative decomposition')
    check(norm(pure-pure.T)<2e-13,'full pure term symmetric')
    curl=norm(df-df.T)
    check(curl<=2*c*A*(LU+LV)+2e-13,'full square curl')
    origin=vertex(A,a,b,c,C,D,d,np.zeros(n),Y,k)
    check(norm(origin['E'])<=a*b*c*A**4*norm(D[0]@Y+d[0])+2e-13,'captured-origin envelope')
    # Current includes both source and fresh-root derivatives.
    theta=float(rng.uniform()); q=float(rng.uniform()); variance=float(rng.uniform(.02,2.))
    G=rng.normal(size=2); Ut=(1-theta)*out['V']+theta*out['U']
    DUt=(1-theta)*out['DV'][:,:n]+theta*out['DU'][:,:n]
    DZ=np.hstack([-q*DUt,np.sqrt(variance)*np.eye(2)])
    Ltheta=(1-theta)*LV+theta*LU
    LZ=np.sqrt(variance+q*q*Ltheta*Ltheta)
    check(norm(DZ)<=LZ+2e-13,'full direct-current first')
    DJoint=np.vstack([np.hstack([DD[:,:n],np.zeros((2,2))]),DZ])
    LJoint=np.sqrt((9*A/4)**2+LZ**2)
    check(norm(DJoint)<=LJoint+2e-13,'joint mark/current first')
    t=rng.normal(size=2); Z=np.sqrt(variance)*G-q*Ut; eps=1e-5
    mark=-q*Delta
    fd=(np.sin(t@(Z+eps*mark))-np.sin(t@(Z-eps*mark)))/(2*eps)
    rhs=np.cos(t@Z)*(t@mark)
    err=abs(fd-rhs)
    check(err<1e-9,'pointwise current test identity')
    direction=rng.normal(size=n); direction/=norm(direction); h=1e-6
    plus=vertex(A,a,b,c,C,D,d,O+h*direction,Y,k)['E']
    minus=vertex(A,a,b,c,C,D,d,O-h*direction,Y,k)['E']
    fd_err=norm((plus-minus)/(2*h)-DE[:,:n]@direction)
    check(fd_err<1e-7,'independent finite-difference derivative check')
    if boundD>1e-15: maxima['Delta_mark_ratio']=max(maxima['Delta_mark_ratio'],norm(Delta)/boundD)
    if boundE>1e-15: maxima['E_mark_ratio']=max(maxima['E_mark_ratio'],norm(E)/boundE)
    maxima['E_first_ratio']=max(maxima['E_first_ratio'],norm(DE[:,:n])/(21*A/8))
    maxima['Delta_first_ratio']=max(maxima['Delta_first_ratio'],norm(DD[:,:n])/(9*A/4))
    if c>0: maxima['full_curl_ratio']=max(maxima['full_curl_ratio'],curl/(2*c*A*(LU+LV)))
    maxima['current_first_ratio']=max(maxima['current_first_ratio'],norm(DZ)/LZ)
    maxima['joint_first_ratio']=max(maxima['joint_first_ratio'],norm(DJoint)/LJoint)
    maxima['jacobian_fd_error']=max(maxima['jacobian_fd_error'],fd_err)
    maxima['current_identity_fd_error']=max(maxima['current_identity_fd_error'],err)
    maxima['noncommutation']=max(maxima['noncommutation'],out['noncomm'])
check(maxima['noncommutation']>1e-8,'noncommuting Hessian fixture exercised')

# Exact rational linear-chain sign and coefficient, independent of floating-point tests.
for A in [Fraction(1,17),Fraction(1,3),Fraction(1,2)]:
    a,b,c=Fraction(2,3),Fraction(3,5),Fraction(5,7)
    x0,x1,x2,x3=Fraction(13,11),Fraction(-7,3),Fraction(2,9),Fraction(19,5)
    h0=A*x0; h1=A*(x1-a*h0); h1b=A*x1
    U=A*(x2-b*h1); V=A*(x2-b*h1b)
    E=A*(x3-c*U)-A*(x3-c*V)
    check(U-V==a*b*A**3*x0,'exact rational Delta')
    check(E==-a*b*c*A**4*x0,'exact rational four-force sign')

# Exact two-variable exponent bookkeeping: monomial=(A exponent,u exponent).
half=Fraction(1,2)
r=(Fraction(1),-half); delta=(Fraction(3),Fraction(0)); mu=r; acurl=(Fraction(1),Fraction(0)); readout=(Fraction(0),half)
def add(*xs): return tuple(sum(x[i] for x in xs) for i in range(2))
def mul(x,z): return tuple(z*t for t in x)
terms=[add(readout,mul(r,3),delta),add(readout,mul(r,2),delta,mu),
       add(readout,mul(r,2),delta,acurl),add(readout,mul(r,4),delta,mul(mu,-half)),
       add(readout,mul(r,4),delta)]
expected=[(6,-1),(6,-1),(6,-half),(Fraction(13,2),Fraction(-5,4)),(7,Fraction(-3,2))]
check(terms==expected,'all five native physical error powers')
check(add(readout,r)==(1,0),'native linear first')
check(add(readout,mul(r,2),mul(mu,-half))==(Fraction(3,2),Fraction(-1,4)),'native self-reserve first')
# Exact ratios against A^6/u use sqrt(u), (A^2/u)^(1/4), A/sqrt(u).
check(add(terms[2],(-6,1))==(0,half),'third reduction ratio')
check(add(terms[3],(-6,1))==(half,Fraction(-1,4)),'fourth reduction ratio')
check(add(terms[4],(-6,1))==(1,-half),'fifth reduction ratio')

result={'status':'PASS','assertions':checks,'random_cases':400,'seed':5102026,
        'scope':'Independent local finite source, full derivatives/curl, current and exact exponent diagnostics. Does not execute LOW30 mean compiler or establish true-history ancestry.',
        'maxima':maxima,'native_error_exponents':[[str(v) for v in row] for row in terms]}
p=Path(__file__).with_name('independent_shifted_atom_checks.json'); p.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
