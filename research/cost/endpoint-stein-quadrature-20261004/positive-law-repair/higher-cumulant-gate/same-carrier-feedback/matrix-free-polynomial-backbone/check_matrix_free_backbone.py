#!/usr/bin/env python3
"""Author diagnostics for the literal source; no mean compiler implementation."""
import json, math, hashlib
from pathlib import Path
import numpy as np
import sympy as sp
HERE=Path(__file__).resolve().parent
rng=np.random.default_rng(51026)
counts={}; diagnostics={}
def check(ok,group,msg):
    counts[group]=counts.get(group,0)+1
    if not bool(ok): raise AssertionError(f'{group}: {msg}')
def close(x,y,tol,group,msg):check(np.max(np.abs(np.asarray(x)-np.asarray(y)))<=tol,group,msg)
def norm(x):return float(np.linalg.norm(x,2))
def quadrature(n=5):
    t,w=np.polynomial.legendre.leggauss(n); t=(t+1)/2;w=w/2
    return t,w,np.sqrt(1-t*t)
def value_graph(z,G,N,M,g,t,weights,c):
    B=np.zeros_like(z);F=np.zeros_like(z)
    for ti,wi,ci in zip(t,weights,c):
        x=ti*z+ci*G;v=(x+N)/2;s=x/4+N/2+M/4
        b=g(s);u=g(v);q=g(b);T=x-u+q
        F+=wi*g(T);B+=wi*g(x)
    return F,B,F-B

def first_graph(z,G,N,M,g,jac,t,weights,c):
    d=len(z);EG=np.zeros((d,d));EN=EG.copy();EM=EG.copy();EZ=EG.copy()
    BG=EG.copy();BZ=EG.copy()
    for ti,wi,ci in zip(t,weights,c):
        x=ti*z+ci*G;v=(x+N)/2;s=x/4+N/2+M/4;b=g(s);T=x-g(v)+g(b)
        H,V,J,W,H0=jac(T),jac(v),jac(b),jac(s),jac(x)
        P=np.eye(d)-V/2+J@W/4;Q=-V/2+J@W/2;R=J@W/4
        dx=H@P-H0
        EG+=wi*ci*dx;EN+=wi*H@Q;EM+=wi*H@R;EZ+=wi*ti*dx
        BG+=wi*ci*H0;BZ+=wi*ti*H0
    return np.concatenate((EG,EN,EM),axis=1),EZ,BG,BZ

def symbolic():
    k=sp.symbols('k',real=True)
    C=sp.Matrix([[sp.Rational(1,4),sp.Rational(1,4)],[sp.Rational(1,4),sp.Rational(5,16)]])
    R=sp.Matrix([[sp.Rational(1,2),0],[sp.Rational(1,2),sp.Rational(1,4)]])
    check(R*R.T==C,'symbolic','conditional joint factor')
    H=sp.Matrix([[sp.Rational(1,2),sp.Rational(3,8)],[sp.Rational(3,8),sp.Rational(3,8)]])
    mean=sp.Matrix([[sp.Rational(1,2)],[sp.Rational(1,4)]])
    check(C+mean*mean.T==H,'symbolic','unconditional covariance')
    row=sp.Matrix([[k,-k*k]])
    check(sp.expand((row*C*row.T)[0])==k*k/4-k**3/2+5*k**4/16,'symbolic','conditional variance polynomial')
    check(sp.expand((-k/2+k*k/2)**2+(k*k/4)**2)==k*k/4-k**3/2+5*k**4/16,'symbolic','row covariance expansion')
    m=sp.Matrix([[sp.Rational(1,2),sp.Rational(1,2),0],[sp.Rational(1,4),sp.Rational(1,2),sp.Rational(1,4)]])
    check(m*m.T==H,'symbolic','unconditional v,w covariance')
    A=sp.symbols('A',positive=True)
    check(sp.expand(1+A*(2+A)-(1+A)**2)==0,'symbolic','nested VALUE floor')

def quadratic():
    rows=[]
    for d in [1,2,5,16]:
      for A in [0.,.001,.05,.125,.5]:
        Q,_=np.linalg.qr(rng.normal(size=(d,d))); ev=rng.uniform(0,A,size=d)
        K=(Q*ev)@Q.T;L=np.eye(d)-K/2+K@K/4
        for n in [1,3,7]:
          t,w,c=quadrature(n);beta=w@c
          close(w.sum(),1,1e-14,'quadratic','positive mass one')
          close(w@t,.5,1e-14,'quadratic','exact first moment')
          check(np.all(w>0) and np.all(t>0) and np.all(t<1),'quadratic','interior positive rule')
          z,G,N,M=rng.normal(size=(4,d));calls=[0]
          def g(x):calls[0]+=1;return K@x
          F,B,E=value_graph(z,G,N,M,g,t,w,c)
          expected=.5*K@L@z+beta*K@L@G+(-K@K/2+K@K@K/2)@N+(K@K@K/4)@M
          close(F,expected,1e-13,'quadratic','literal polynomial action identity')
          check(calls[0]==5*n,'ledger','full E occurrence source leaves')
          close(E,F-B,1e-15,'quadratic','same occurrence baseline subtraction')
          Cov=(beta*K@L)@(beta*K@L).T+(-K@K/2+K@K@K/2)@(-K@K/2+K@K@K/2).T+(K@K@K/4)@(K@K@K/4).T
          formula=beta**2*K@K@L@L+np.linalg.matrix_power(K,4)@(np.eye(d)-K)@(np.eye(d)-K)/4+np.linalg.matrix_power(K,6)/16
          close(Cov,formula,1e-13,'quadratic','full oriented raw covariance')
          zmean=.5*K@L@z
          exact=.5*(K@z)-.25*(K@K@z)+.125*(K@K@K@z)
          close(zmean,exact,1e-14,'quadratic','direct exact unit-buffer mean')
          close(value_graph(np.zeros(d),np.zeros(d),np.zeros(d),np.zeros(d),g,t,w,c)[0],0,0,'zero','literal total root zero')
          if A>0 and n>1:
            wrong=sum(w*w*c*c)*K@K@L@L+np.linalg.matrix_power(K,4)@(np.eye(d)-K)@(np.eye(d)-K)/4+np.linalg.matrix_power(K,6)/16
            check(norm(Cov-wrong)>0,'quadratic','independent-node covariance differs')
        rows.append({'D':d,'A':A,'spectrum_max':float(ev.max())})
    diagnostics['quadratics']=rows

def nonlinear():
    max_fd=0.;max_ratio=0.;min_commute=0.;rows=[]
    for d in [2,5,9]:
      for A in [.001,.03,.125,.5]:
        Q,_=np.linalg.qr(rng.normal(size=(d,d)))
        K=(Q*rng.uniform(.1*A,.3*A,size=d))@Q.T
        directions=rng.normal(size=(4,d));directions/=np.linalg.norm(directions,axis=1)[:,None]
        b=.5*A/4
        def g(x):return K@x+b*(np.tanh(directions@x)@directions)
        def jac(x):
          u=np.tanh(directions@x)
          return K+b*directions.T@np.diag(1-u*u)@directions
        t,w,c=quadrature(5);beta=w@c;eta=A/2+A*A/4;q=A/2+A*A/2;r=A*A/4
        bound=A*math.sqrt(beta*beta*(1+eta)**2+q*q+r*r)
        curlbound=2*beta*A*eta+A*math.sqrt(q*q+r*r)
        for trial in range(9):
          z,G,N,M=rng.normal(size=(4,d))*(1+trial)
          F,B,E=value_graph(z,G,N,M,g,t,w,c)
          J,Jz,Jbg,Jbz=first_graph(z,G,N,M,g,jac,t,w,c)
          check(norm(J)<=bound+1e-13,'nonlinear_ports','full private first')
          check(norm(Jz)<=A*(1+eta)/2+1e-13,'nonlinear_ports','raw caller first')
          check(norm(Jbg)<=A*beta+1e-13,'nonlinear_ports','baseline private first')
          close(Jbg,Jbg.T,1e-14,'nonlinear_ports','baseline symmetric Jacobian')
          square=np.zeros((3*d,3*d));square[:d]=J;curl=norm(square-square.T)
          check(curl<=curlbound+1e-13,'nonlinear_ports','noncommuting curl')
          check(math.sqrt(2)*bound<=2*A+1e-13,'guards','normalized declared ell')
          check(math.sqrt(2)*curlbound<=3*A*A+1e-13,'guards','normalized declared curl')
          dz,dG,dN,dM=rng.normal(size=(4,d));eps=1e-5
          ep=value_graph(z+eps*dz,G+eps*dG,N+eps*dN,M+eps*dM,g,t,w,c)[2]
          em=value_graph(z-eps*dz,G-eps*dG,N-eps*dN,M-eps*dM,g,t,w,c)[2]
          fd=(ep-em)/(2*eps);chain=Jz@dz+J@np.r_[dG,dN,dM]
          err=norm((fd-chain)[:,None]);max_fd=max(max_fd,err)
          check(err<2e-9,'chain_rule','complete primal ancestry directional derivative')
          zeros=np.zeros(d);E0=value_graph(z,zeros,zeros,zeros,g,t,w,c)[2]
          check(norm(E0[:,None])<=A*A*(.25+A/8)*norm(z[:,None])+1e-13,'origins','caller-origin amplitude')
          J0z=first_graph(z,zeros,zeros,zeros,g,jac,t,w,c)[1]
          check(norm(Jz-J0z)<=A*(1+eta)+1e-13,'origins','anchored caller first')
          n1,n2=jac(z),jac(G)
          min_commute=max(min_commute,norm(n1@n2-n2@n1))
          max_ratio=max(max_ratio,curl/curlbound)
        rows.append({'D':d,'A':A,'private_first_bound':bound,'curl_bound':curlbound})
    check(min_commute>1e-6,'nonlinear_ports','fixture genuinely noncommuting')
    diagnostics['nonlinear']={'max_directional_error':max_fd,'max_curl_to_bound':max_ratio,'largest_Hessian_commutator':min_commute,'rows':rows}

def certificates():
    rows=[]
    for D in [1,2,20,1000,1000000]:
      for A in [1e-6,.001,.05,.125,.5]:
        R=math.sqrt(D)+math.sqrt(8*math.log(1/A));L=A/3
        e=L*math.sqrt(2)*math.exp(-2*math.log(1/A))
        bound=A*(3+2*A)*e;exact=math.sqrt(2)/3*(3+2*A)*A**4
        close(e,math.sqrt(2)*A**3/3,1e-15,'envelope','all Gaussian contraction envelope')
        close(bound,exact,1e-15,'envelope','order-four anisotropic certificate')
        check(bound<=math.sqrt(2)/3*(3+2*A)*A**4*math.sqrt(D)*(1+1e-14),'envelope','all-dimension grade')
        vals=np.linspace(0,A,100)
        check(np.max(1-vals+vals*vals/2)<=1,'envelope','history Gaussian contraction')
        check(math.sqrt(3/8)*A<=1,'envelope','nested action Gaussian contraction')
        rows.append({'D':D,'A':A,'R':R,'e_each':e,'bias_bound':bound})
    diagnostics['certificates']=rows

def main():
    symbolic();quadratic();nonlinear();certificates()
    out={'status':'PASS','assertions':sum(counts.values()),'groups':counts,'scope':'Literal raw source, analytic port formulas and explicit certificate diagnostics; not an implementation of the imported completed mean compiler.','diagnostics':diagnostics}
    (HERE/'matrix_free_backbone_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'assertions':out['assertions'],'groups':counts},indent=2))
if __name__=='__main__':main()
