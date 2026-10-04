# Publication adaptation: all 60,813 mathematical checks retained; 11 public source hashes replace 13 historical input hashes.
from fractions import Fraction as F
from pathlib import Path
import hashlib, json, math
checks=0

def check(x):
 global checks
 assert x
 checks+=1

# Current recurrence endpoint difference and rational trajectory.
c=F(0)
current=[]
for p2 in range(14,401):
 P=F(p2,2); R=P+F(1,2)
 if P==7:
  c=max(c+P-F(49,10),F(2,3)*R*c)
 else:
  left=c+P-F(49,10); right=R/P*(c+P-F(27,5))
  check(right-left==(c-F(27,5))/(2*P))
  c=max(left,right)
 if R>=F(17,2):
  n=int(2*R-17)
  H=sum((F(1,k) for k in range(17,n+17)),F(0))
  check(c/R==F(78,85)+n-F(54,5)*H)
 if R in [F(15,2),F(8),F(17,2),F(10),F(20),F(50),F(100)]:
  current.append({'P':str(R),'c':str(c),'ratio':float(c/R)})

# Granted ideal shared-core ledger (not an admitted source construction).
P=F(15,2); c=F(21,10)
for _ in range(400):
 R=P+F(1,2)
 c=max(R/P*c,R-F(27,5))
 check(c==R-F(27,5))
 P=R

# General shared-core endpoint reduction when packet exponent b<=5.4.
for p2 in range(15,90):
 P=F(p2,2);R=P+F(1,2);S=R/P
 for bi in range(55):
  b=F(bi,10)
  left=b+R-F(27,5)
  right=S*b+max(F(0),R-F(27,5)*S)
  check(right-left==(b-F(27,5))/(2*P))
  check(right<=left)
  for ti in range(11):
   s=1+(S-1)*F(ti,10)
   f=s*b+max(F(0),R-F(27,5)*s)
   check(f<=left)

# No-copy seed rebuild: arbitrary terminal heat dilation times c=0 stays zero.
for S in [F(1),F(16,15),F(5),F(100),F(12345,7)]:
 check(S*F(0)==0)
# A positive incumbent's normalized inherited ratio does not fall.
P=F(15,2); c=F(21,10); ratio=c/P
for _ in range(400):
 R=P+F(1,2); c=R/P*c
 check(c/R==ratio);P=R

# Explicit target-specific small-kappa budget, including special first stage.
for N in range(1,101):
 eps=F(1,N+1)
 k0=F(3,4)*eps
 c=5*k0;P=F(15,2)
 ratio=c/P
 check(ratio==eps/2)
 for i in range(1,N+1):
  k=eps*P/F(2**(i+1))
  R=P+F(1,2);c=R/P*(c+k);P=R
  ratio+=eps/F(2**(i+1))
  check(c/P==ratio)
 check(c/P<eps)

# Additive same-heat corrections: the running max is the exponent.
c=F(13,2);bs=[F(0),F(1),F(9),F(3),F(19,2),F(0)]
for i,b in enumerate(bs):
 c=max(c,b)
 check(c==max([F(13,2)]+bs[:i+1]))

# Fixed nine weights cancel cubic but not quartic cumulants.
a=[F(1,6)]*8+[F(-1,3)]
check(sum(a)==1)
check(sum(x*x for x in a)==F(1,3))
check(sum(x**3 for x in a)==0)
check(sum(x**4 for x in a)==F(1,54))
# Holder lower bound equality for equal weights.
for m in range(1,101):
 check(sum(F(1,m)**4 for _ in range(m))==F(1,m**3))
# Smooth scalar fixture E=A sin(G): explicit nonzero fourth cumulant.
ms2=(1-math.exp(-2))/2
ms4=(3-4*math.exp(-2)+math.exp(-8))/8
k4=ms4-3*ms2*ms2
check(k4<-.25)

# Exact variance decomposition for U+mean(H_i) vs independent U_i+H_i.
for N in range(1,101):
 check(F(1)+F(1,N) >= F(2,N))
 check((F(1)+F(1,N))-F(2,N)==1-F(1,N))

# Analytic Gaussian W2 witness for aliased replicas: no 1/N decay.
alpha=.05;v=.5
alias=math.sqrt(v+alpha*alpha)-math.sqrt(v)
for N in [1,2,10,100,1000]:
 independent=math.sqrt(v+alpha*alpha/N)-math.sqrt(v)
 check(alias>=independent)
 if N>1:check(alias>independent)

# Fixed-seed Gaussian resolvent truncation retains its deterministic bias.
for L in range(1,15):
 for x in [F(1,10),F(1,3),F(1,2)]:
  t=sum(((-x)**j for j in range(L+1)),F(0))
  check(t-1/(1+x)==-(-x)**(L+1)/(1+x))

root=Path(__file__).parent.parent
manifest = {'files': [{'file': 'COST-OPTIMIZATION-PROOFS-RECONSTRUCTED.md', 'sha256': 'ece35741d1cae46c6bd85d82851fe56f18376d5478da4aae74d315483260ca9e'}, {'file': 'GROWTH-AND-SUBCRITICAL-TARGETS-RECONSTRUCTED.md', 'sha256': 'edf398befbc5c6efb131b60824fa64e06c914f67bfebb2670464cee2c55e9e2a'}, {'file': 'NO-COPY-ACCEPTANCE-CRITERION-RECONSTRUCTED.md', 'sha256': 'a69a234b182beccb4f1ae091b57e983dc3a9638c8b9f3b9c37a7f2a2cd7e7d2f'}, {'file': 'check_cost_schedules.py', 'sha256': 'e67a0abe6b38c7fdcd27f2f7c1d10ac6c16b0a5814b058a3240e9b49a0cc7a9a'}, {'file': 'check_extended_cost_schedules.py', 'sha256': '03e4c36781cd9f0d093631362c34f3297217291e4e614d0600d1390371cf6954'}, {'file': 'check_no_copy_acceptance_arithmetic.py', 'sha256': '9c8a992e44dbded925e23348c86bf4f65f9edf2a57cf1eaa58598386bad2efe0'}, {'file': 'check_recovery_additional.py', 'sha256': '80192e85a3e08d512ff286527cfebc94db88aeb7ab4c4d338486d914170659ed'}, {'file': 'cost_schedule_checks.json', 'sha256': 'd386ca0fbd567c62e5ddef6c5bf1dd698093ec927c5d44a6e8b7dee91a21dcc7'}, {'file': 'extended_cost_schedule_checks.json', 'sha256': '50b6c44b67bde1bdb83599d22459d2bf9e55c4e5b9eddba31aeebe0de2003b76'}, {'file': 'no_copy_acceptance_arithmetic.json', 'sha256': 'd373d46e5442942358df428b41f4eec25ae924db21a3802244a109b00e91a78b'}, {'file': 'recovery_additional_checks.json', 'sha256': '11893de074ef3be7ae5579d7c17687c69d9352ba51f064d1ac01e887ed233d7e'}]}
verified=[]
for entry in manifest['files']:
 path=root/entry['file'];actual=hashlib.sha256(path.read_bytes()).hexdigest()
 check(actual==entry['sha256'])
 verified.append({'file':entry['file'],'sha256':actual})

out={
 'status':'PASS','checks':checks,
 'scope':'Exact rational recurrence, budget, variance, and cumulant algebra; numerical evaluation only for stated smooth scalar / Gaussian diagnostics. No new all-order source theorem.',
 'current_trajectory':current,
 'ideal_shared_core_limit_ratio':1,
 'nine_weight_fourth_power_sum':'1/54',
 'sin_gaussian_fourth_cumulant':k4,
 'public_cost_source_hashes_verified':verified,
 'whole_family_no_copy_status':'OPEN',
}
Path(__file__).with_name('recurrence_candidate_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':'PASS','checks':checks,'sin_gaussian_fourth_cumulant':k4,'source_hashes_verified':len(verified)}))
