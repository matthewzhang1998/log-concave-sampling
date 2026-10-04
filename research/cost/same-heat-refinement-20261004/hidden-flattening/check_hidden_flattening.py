import json,math
from pathlib import Path
import numpy as np

mixtures=[]; decoders=[]
for a in [.2,.1,.05,.025]:
    lam=.5; t=1/(1+a*lam); tau=a**1.5
    truevar=a*t
    for n in [1,4,16]:
        mixvar=truevar+t*t*tau*tau/n
        diff=math.sqrt(mixvar)-math.sqrt(truevar)
        mixtures.append({'a':a,'replicas':n,'target_variance':truevar,'mixture_variance':mixvar,
                         'w2':diff,'n_w2_over_a_2p5':n*diff/a**2.5})
        assert mixvar>truevar
    modecoef=sum((-a*lam)**j for j in range(9))
    b=a/2; tb=1/(1+b*lam); vc=a*(modecoef*modecoef+1)/4
    vtrue=b*tb+tb*tb*vc
    vdrop=b*tb+tb*tb*vc/(1+lam*tb*vc)
    assert vtrue>vdrop
    decoders.append({'a':a,'finite_mode_steps':8,'exact_normalized_path_variance':vtrue,
                     'dropped_partition_variance':vdrop,
                     'w2_over_a_1p5':(math.sqrt(vtrue)-math.sqrt(vdrop))/a**1.5})

whitening=[]
for k in [1,2,4,8,16,32]:
    # Fixed x0=0, normalized V=0 Gaussian decoder; vectors ordered R then X.
    BR=np.zeros((k,k)); BX=np.zeros((k,k))
    for i in range(k):
        for j in range(i+1):
            BR[i,j]=.5**(i-j+1)
            BX[i,j]=.5**(i-j)/math.sqrt(2)
    L=np.block([[np.eye(k),np.zeros((k,k))],[BR,BX]])
    LX=L[k:,:]
    LC=np.zeros((k,2*k))
    for i in range(k):
        LC[i]=.5*L[i]
        if i: LC[i]+=.5*L[k+i-1]
    a=.1; lam=.5; b=a/2
    exact_precision=np.eye(2*k)+a*lam*(LX.T@LX)-a*lam/(1+b*lam)*(LC.T@LC)
    eig=np.linalg.eigvalsh(exact_precision)
    singular=np.linalg.svd(L,compute_uv=False)
    assert singular[0]/singular[-1]<5
    assert eig[0]>.5
    whitening.append({'decoder_steps':k,'latent_scalar_dimension':2*k,
                       'gaussian_map_condition':float(singular[0]/singular[-1]),
                       'normalized_joint_precision_min':float(eig[0]),
                       'normalized_joint_precision_max':float(eig[-1])})

# Matrix-free quadratic mean/noise polynomials have no inverse-heat power
# in the degree for a fixed requested accuracy order.
A=np.array([[.6,.2],[.2,.4]]); c=np.array([.3,-.4]); grad=np.array([.1,-.2])
polys=[]
for b in [.1,.03,.01]:
    exact=np.linalg.inv(np.eye(2)+b*A)
    for m in [2,4,8]:
        term=np.eye(2); P=np.eye(2)
        C=np.eye(2); root=np.eye(2); coeff=1.
        for j in range(1,m+1):
            term=term@(-b*A); P+=term
            C=C@(b*A); coeff*=(-.5-(j-1))/j; root+=coeff*C
        vals,vec=np.linalg.eigh(exact); rex=(vec*np.sqrt(vals))@vec.T
        assert np.linalg.norm(P-exact,2)<=b**(m+1)/(1-b)+1e-14
        polys.append({'b':b,'degree':m,'inverse_error':float(np.linalg.norm(P-exact,2)),
                      'square_root_error':float(np.linalg.norm(root-rex,2))})

out={'status':'PASS: hidden-mean mixture mismatch, decoder normalization tilt, conditioned Gaussian whitening, and quadratic matrix-free actions',
     'scope':'No universal impossibility theorem; affine/quadratic simplification requires a separately certified source graph.',
     'latent_mean_mixture':mixtures,'decoder_partition_drop':decoders,
     'favorable_joint_whitening':whitening,'quadratic_polynomials':polys}
Path(__file__).with_name('hidden_flattening_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
