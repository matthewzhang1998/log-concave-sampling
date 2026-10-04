#!/usr/bin/env python3
"""Finite original-VALUE identities for the terminal-flat affine-feedback escape.

The Gaussian pair below is its analytical conditional reference, not an
implementation of the admitted pair compiler. No oracle derivative is a VALUE.
"""
import json, math
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
T=np.array([[.5,.05],[.05,.5]])
N=np.diag([1.,0.]); e1=np.array([1.,0.])
c,s=.6,.8; delta=.025
I=np.eye(2); count=0; max_fd=0.; max_cov=0.
rng=np.random.default_rng(20261004)

def check(x,y,tol=2e-11):
    global count
    err=float(np.max(np.abs(np.asarray(x)-np.asarray(y))))
    assert err<=tol,(err,x,y)
    count+=1

def sym(x): return (x+x.T)/2
def ps(z):
    u=10*(z+.5)
    return 0. if abs(u)>=1 else (z+.5)*math.exp(1-1/(1-u*u))
def dps(z):
    u=10*(z+.5)
    if abs(u)>=1: return 0.
    b=1/(1-u*u)
    return math.exp(1-b)*(1-2*u*u*b*b)
def sqrtps(A):
    val,vec=np.linalg.eigh(sym(A)); assert min(val)>-1e-11
    return (vec*np.sqrt(np.maximum(val,0)))@vec.T

cases=[]
for a in [1/16,1/32,1/128,1/512]:
  eps=a**.9; q=a*eps; kap=a*a
  def g(y): return T@y+delta*q*ps(y[0]/q)*e1
  def H(y): return T+delta*dps(y[0]/q)*N
  for angle in [0.,.02,math.pi/4,math.pi/2]:
    C0=np.concatenate([c*I,s*math.cos(angle)*I,np.zeros((2,2)),s*math.sin(angle)*I],axis=1)
    C1=C0.copy(); C1[:,4:6]=-q*T
    C0p=np.concatenate([C0,np.zeros((2,2))],axis=1)
    C1p=np.concatenate([C1,np.zeros((2,2))],axis=1)
    B=np.concatenate([T@C0,sqrtps(I-T@T)],axis=1)
    P=np.concatenate([I,np.zeros((2,8))],axis=1)
    K=P.T@B
    check(C0@C0.T,I); check(C0@C1.T,I); check(B@B.T,I)
    def G(w): return kap*(C0p.T@g(C0p@w)-C1p.T@g(C1p@w))
    def JG(w): return kap*(C0p.T@H(C0p@w)@C0p-C1p.T@H(C1p@w)@C1p)
    def E(w): return kap*T@(g(C0p@w)-g(C1p@w))
    def JE(w): return kap*T@(H(C0p@w)@C0p-H(C1p@w)@C1p)
    w0=np.zeros(10); w0[4:6]=e1
    word=-kap*delta*T@N
    DZ=kap*q*T@(T+delta*N)@T
    check(E(w0),kap*q*T@T@T@e1)
    check(JE(w0)[:,:2],c*word)
    check(JE(w0)[:,2:4],s*math.cos(angle)*word)
    check(JE(w0)[:,6:8],s*math.sin(angle)*word)
    check(JE(w0)[:,4:6],DZ)
    O_expected=word@word.T+DZ@DZ.T-c*c*sym(word@word)
    records=[w0]
    for _ in range(30):
        w=rng.normal(size=10)
        # Put some exact argument records into the narrow nonlinear slab.
        x=C0p@w; w[0]+=(q*(-.5+.11*rng.normal())-x[0])/c
        records.append(w)
    for w in records:
      GG=G(w); JJ=JG(w); FF=K@JJ
      check(E(w),B@GG)
      check(JE(w),B@JJ)
      check(JJ,JJ.T)
      check(GG,kap*(C0p.T@(g(C0p@w)-g(C1p@w))+(C0p-C1p).T@g(C1p@w)))
      curl=FF-FF.T
      O=-sym(FF@curl)
      check(P@O@P.T,JE(w)@JE(w).T-sym((JE(w)@P.T)@(JE(w)@P.T)))
      if np.array_equal(w,w0): check(P@O@P.T,O_expected)
      U=np.concatenate([K,-np.eye(10)],axis=1)/math.sqrt(2)
      V=np.concatenate([np.eye(10),K.T],axis=0)/math.sqrt(2)
      block=np.zeros((20,20)); block[:10,:10]=JJ; block[10:,10:]=JJ
      check(U@block@V,curl/2)
      assert np.linalg.norm(U,2)<=1+1e-12 and np.linalg.norm(V,2)<=1+1e-12
      count+=2
      # Full positive reference: conditional means M p and Q p.
      s0=.3; rho=.09; L0=2.; tau=kap*L0
      M=s0*FF
      Q=s0*rho/(2*tau)*curl
      weight=.07; b=.15; d=-weight*tau*2/(2*s0*s0*rho*b)
      vout=1.; keep=vout-b*b-d*d
      assert keep>0 and min(np.linalg.eigvalsh(np.eye(10)-M@M.T))>0
      assert min(np.linalg.eigvalsh(np.eye(10)-Q@Q.T))>0
      cov=b*b*np.eye(10)+d*d*np.eye(10)+keep*np.eye(10)+b*d*(M@Q.T+Q@M.T)
      expected=np.eye(10)-weight*O
      check(cov,expected)
      max_cov=max(max_cov,float(np.max(np.abs(cov-expected))))
      # Flipping the fork sign produces the opposite orientation (nonzero test).
      if np.array_equal(w,w0):
        wrong=np.eye(10)-b*d*(M@Q.T+Q@M.T)
        assert np.linalg.norm(wrong-expected)>1e-16
        count+=1
      dh=q*2e-5
      num=np.column_stack([(G(w+dh*np.eye(10)[j])-G(w-dh*np.eye(10)[j]))/(2*dh) for j in range(10)])
      ferr=float(np.max(np.abs(num-JJ))/kap)
      assert ferr<2e-5,ferr
      count+=1; max_fd=max(max_fd,ferr)
    # Original midpoint baseline has DV N=kap T and DW N=0 at V=0.
    JN=np.zeros((2,10)); JN[:,6:8]=kap*T
    mixed=JE(w0)@JN.T
    check(mixed,s*math.sin(angle)*word@(kap*T))
    if angle:
      assert np.linalg.norm(mixed)>0; count+=1
    cases.append(dict(a=a,angle=angle,word_skew_over_kappa=float(np.linalg.norm(word-word.T)/kap),
                      orientation_22_over_kappa_squared=float(O_expected[1,1]/kap**2),
                      mixed_baseline_norm=float(np.linalg.norm(mixed))))

# A COMMON input rotation preserves the complete field, including canceling
# newly nonzero V columns. It is not the changed-source construction above.
for _ in range(100):
  X=rng.normal(size=(2,8)); Y=rng.normal(size=(2,8))
  X[:,6:8]=0; Y[:,:6]=0
  ang=rng.uniform(-math.pi,math.pi)
  R=np.eye(8)
  R[2:4,2:4]=math.cos(ang)*I; R[2:4,6:8]=math.sin(ang)*I
  R[6:8,2:4]=-math.sin(ang)*I; R[6:8,6:8]=math.cos(ang)*I
  check(R@R.T,np.eye(8)); check((X@R)@(Y@R).T,X@Y.T)

# Positive-width finite-clock identity: independent fine banks per product.
# Exact identities hold for any finite convex average of the actual firsts.
a=1/32; eps=a**.9; q=a*eps; kap=a*a; angle=math.pi/4
for v in [q/64,q/8,q/2,.03]:
  # Closures currently use the final loop's a; reset by evaluating directly.
  C0=np.concatenate([c*I,s/math.sqrt(2)*I,np.zeros((2,2)),s/math.sqrt(2)*I,np.zeros((2,2))],axis=1)
  C1=C0.copy(); C1[:,4:6]=-q*T
  B=np.concatenate([T@C0[:,:8],sqrtps(I-T@T)],axis=1)
  P=np.concatenate([I,np.zeros((2,8))],axis=1); K=P.T@B
  w=np.zeros(10); w[4:6]=e1/math.sqrt(1-v*v)
  fine=rng.normal(size=(128,10)); fine=np.r_[fine,-fine]
  J=sum(kap*(C0.T@H(C0@(math.sqrt(1-v*v)*w+v*z))@C0-C1.T@H(C1@(math.sqrt(1-v*v)*w+v*z))@C1) for z in fine)/len(fine)
  Jf=K@J; C=Jf-Jf.T
  check(C,K@J-J@K.T)
  check(-P@sym(Jf@C)@P.T,(B@J)@(B@J).T-sym((B@J@P.T)@(B@J@P.T)))

# The smaller genuine-gradient V partner, including its known Gaussian body.
smaller_partner_cases=[]
for a in [1/16,1/32,1/128]:
  eps=a**.9; q=a*eps; kap=a*a; angle=math.pi/4
  Cplus=np.concatenate([c*I,s*math.cos(angle)*I,np.zeros((2,2)),s*math.sin(angle)*I,np.zeros((2,2))],axis=1)
  Cminus=Cplus.copy(); Cminus[:,6:8]*=-1
  Cp1=Cplus.copy(); Cp1[:,4:6]=-q*T
  Cm1=Cminus.copy(); Cm1[:,4:6]=-q*T
  Q=2*kap*q*s*math.sin(angle)*T@T
  LG=np.zeros((10,10)); LG[6:8,4:6]=Q; LG[4:6,6:8]=Q
  Bcommon=np.zeros((2,10)); Bcommon[:,:2]=T/c; Bcommon[:,8:10]=sqrtps(I-T@T/c**2)
  check(LG,LG.T); check(Bcommon@Bcommon.T,I); check(Bcommon@LG,np.zeros((2,10)))
  def gradlift(w,Ca,Cb): return kap*(Ca.T@g(Ca@w)-Cb.T@g(Cb@w))
  def jaclift(w,Ca,Cb): return kap*(Ca.T@H(Ca@w)@Ca-Cb.T@H(Cb@w)@Cb)
  for _ in range(60):
    w=rng.normal(size=10)
    if _%2==0:
      w[0]+=(q*(-.5+.09*rng.normal())-(Cplus@w)[0])/c
    nl=gradlift(w,Cplus,Cp1)-gradlift(w,Cminus,Cm1)-LG@w
    jnl=jaclift(w,Cplus,Cp1)-jaclift(w,Cminus,Cm1)-LG
    literal=kap*delta*q*(Cplus.T@(ps((Cplus@w)[0]/q)*e1)-Cp1.T@(ps((Cp1@w)[0]/q)*e1)
                        -Cminus.T@(ps((Cminus@w)[0]/q)*e1)+Cm1.T@(ps((Cm1@w)[0]/q)*e1))
    D=kap*T@(g(Cplus@w)-g(Cp1@w)-g(Cminus@w)+g(Cm1@w))
    check(nl,literal); check(jnl,jnl.T); check(Bcommon@nl,D)
    jD=kap*T@(H(Cplus@w)@Cplus-H(Cp1@w)@Cp1-H(Cminus@w)@Cminus+H(Cm1@w)@Cm1)
    check(Bcommon@jnl,jD)
  # Conditional averages: reflecting fine V makes plane cancellation exact.
  for v in [q/64,q/4,.1]:
    for root_v in [np.zeros(2),np.array([.2*q,-.1*q])]:
      root=np.zeros(10); root[4:6]=e1/math.sqrt(1-v*v); root[6:8]=root_v
      fine0=rng.normal(size=(100,10)); reflected=fine0.copy(); reflected[:,6:8]*=-1
      fine=np.r_[fine0,reflected]
      mats=[]
      for Ca,Cb in [(Cplus,Cp1),(Cminus,Cm1)]:
        A0=sum(H(Ca@(math.sqrt(1-v*v)*root+v*z)) for z in fine)/len(fine)
        A1=sum(H(Cb@(math.sqrt(1-v*v)*root+v*z)) for z in fine)/len(fine)
        BB=kap*T@(A0-A1); FF=kap*q*T@A1@T
        JJ=kap*T@(A0@Ca-A1@Cb)
        mats.append((BB,FF,JJ))
      Bp,Fp,Jp=mats[0]; Bm,Fm,Jm=mats[1]
      db=Bp-Bm; sb=Bp+Bm; df=Fp-Fm; tt=s*s*math.sin(angle)**2
      actual=sym(Jp@(Jp-Jm).T)
      expanded=sym((1-tt)*Bp@db.T+Fp@df.T+tt*Bp@sb.T)
      check(actual,expanded)
      if not root_v.any():
        check(Bp,Bm); check(Fp,Fm); check(actual,2*tt*Bp@Bp.T)
        check((Jp-Jm)[:,:6],np.zeros((2,6)))
      smaller_partner_cases.append(dict(a=a,v=v,root_v=root_v.tolist(),mixed_norm=float(np.linalg.norm(actual)),
                                        away_plane_extra_norm=float(np.linalg.norm(actual-2*tt*Bp@Bp.T))))

out=dict(checks=count,max_normalized_finite_difference_error=max_fd,max_pair_covariance_error=max_cov,
         cases=cases,smaller_partner_cases=smaller_partner_cases,
         note="Finite VALUE, first, selected-word, positive-reference, smaller full-gradient partner, whole mixed field, and rotation algebra only. No finite compiler theorem or nonlinear-baseline joint endpoint is inferred.")
(ROOT/'projected_gradient_escape_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
