#!/usr/bin/env python3
"""Algebraic diagnostics for the new independent-node program, not m3 closure."""
import json
from pathlib import Path
import numpy as np
rng=np.random.default_rng(81267)
checks=0
cases=[]
for n_pair in range(1,9):
  for rep in range(12):
    a=rng.uniform(.03,.47,n_pair)
    tau=np.ravel(np.column_stack((a,1-a)))
    wp=rng.uniform(.2,1.0,n_pair); wp/=wp.sum()
    w=np.repeat(wp/2,2)
    c=np.sqrt(1-tau*tau); n=len(tau)
    q=float(rng.uniform(.15,.95))
    x=np.r_[q,np.zeros(n)]
    Y=np.zeros((n,n+1));Y[:,0]=tau*q
    Y[np.arange(n),1+np.arange(n)]=c
    assert abs(w.sum()-1)<1e-14; checks+=1
    assert abs(w@tau-.5)<1e-14; checks+=1
    ix=np.arange(n)[::2]
    score=np.r_[1/q,np.zeros(n)]
    score[1+ix]=-tau[ix]/c[ix]
    cy=Y@score
    target=tau.copy();target[ix]=0
    assert np.max(np.abs(cy-target))<1e-12; checks+=1
    sig=1/q**2+np.sum((tau[ix]/c[ix])**2)
    assert abs(score@score-sig)<1e-10; checks+=1
    precision=1/q**2+np.sum((tau/c)**2)
    covY=Y@Y.T;covxy=x@Y.T
    var=x@x-covxy@np.linalg.solve(covY,covxy)
    assert abs(var-1/precision)<1e-10; checks+=1
    assert np.linalg.matrix_rank(np.vstack((x,Y)))==n+1;checks+=1
    terminal=x-w@Y
    assert np.linalg.matrix_rank(np.vstack((Y,terminal)))==n+1;checks+=1
    private=w*c
    beta2=private@private;beta1=private.sum()**2
    assert beta2<=beta1+1e-14;checks+=1
    assert abs(np.linalg.norm(private)-np.sqrt(beta2))<1e-14;checks+=1
    assert w[np.setdiff1d(np.arange(n),ix)]@tau[np.setdiff1d(np.arange(n),ix)]<=.5+1e-14;checks+=1
    cases.append({'N_inner':n,'conditional_x_variance':float(var),'all_ancestor_precision':float(precision),'beta2_squared':float(beta2),'common_beta_squared':float(beta1)})
res={'status':'PASS algebraic diagnostics only; no mean-current admission','assertions':checks,'cases':len(cases),'scope':['Changed independent-node covariance and raw Gaussian first','Marked-score independence only at the stated subset cut','All-inner-site posterior precision','Full inner-plus-terminal record has no original Gaussian innovation for invertible linear g'],'sample_cases':cases[:5]+cases[-5:]}
Path(__file__).with_name('independent_node_cut_checks.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps({k:v for k,v in res.items() if k!='sample_cases'},indent=2))
