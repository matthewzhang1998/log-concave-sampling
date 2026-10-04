from pathlib import Path
import numpy as np,json
rng=np.random.default_rng(99); count=0; worst=0
for n in [1,2,4,8,12]:
 for k in range(n+1):
  Q,_=np.linalg.qr(rng.normal(size=(n,n))); D=Q[:,:k]@Q[:,:k].T; I=np.eye(n)
  U=np.concatenate([I,-D],axis=1)/np.sqrt(2); V=np.concatenate([D,I],axis=0)/np.sqrt(2)
  R=np.block([[I-D/2,-D/2],[-D/2,(I-D)/np.sqrt(2)+D/2]])
  for a,b in [(R@R.T,np.eye(2*n)-V@V.T),((I-D)/np.sqrt(2)@((I-D)/np.sqrt(2)).T,I-U@U.T)]:
   er=np.linalg.norm(a-b); worst=max(worst,float(er)); assert er<1e-12; count+=1
  H=rng.normal(size=(n,n)); H=(H+H.T)/2
  HB=np.block([[H,np.zeros_like(H)],[np.zeros_like(H),H]])
  er=np.linalg.norm(U@HB@V-(H@D-D@H)/2); worst=max(worst,float(er)); assert er<1e-12; count+=1
out={'status':'PASS','checks':count,'max_error':worst,'scope':'New exact projector-fill covariance and one-pair skew identity checks; finite pair contracts remain imported.'}
Path(__file__).with_name('projector_fill_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
