from fractions import Fraction as F
from pathlib import Path
import json

checks=0
rows=[]
c=F(0)
for j in range(1,27):
    P=7+F(j-1,2);R=P+F(1,2);S=R/F(22,5)
    c=max(c+P-F(49,10),S*c)
    rows.append(dict(P=str(R),late_omission=str(c),ratio=str(c/(2*R))))
    checks+=1
assert rows[0]['late_omission']=='21/10'
assert rows[1]['late_omission']=='47/10'
assert rows[2]['late_omission']=='799/88'

for P2 in range(15,121):
    P=F(P2,2);R=P+F(1,2);S=R/P;T=2*R/3
    for J1000 in [2001,2005,2010,2020,2030,2040,2050,2100,2500,2900,3000,3500,3899]:
        J=F(J1000,1000);h=(2*J+1)/3
        assert S<h<T
        # Child-protection width exponent and clipped twin exponent.
        sch=(h-S)/2
        er=J-2
        ec=max(F(0),er-(S-1)/2)
        effective=ec+(S-1)/2
        chi=max(J-2,(J-1)/3)
        assert sch>0
        assert max(effective,(J-1)/3)==chi
        # Correct generic projection price, in root exponent units.
        assert (S-1)/2+S+3*sch==J
        # Unclipped twin price r_child A*^2 eps_child <= r a^J.
        assert (S-1)/2+2*S+(er-(S-1)/2)>=J
        checks+=5
    for oldc in [F(0),F(21,10),F(47,10),F(27,5),F(39,5),F(100)]:
        Sd=R/F(22,5)
        bound=max(oldc+P-F(49,10),Sd*oldc)
        for n in range(101):
            s=1+(Sd-1)*F(n,100)
            assert s*oldc+max(F(0),R-F(27,5)*s)<=bound
            checks+=1

# Exact scalar whole-block gauge identity, with a nonlinear translated
# original source. This checks algebraic equivalence, not a host theorem.
for ki in range(1,10):
    k=F(ki,10)
    for zi in range(-7,8):
        z=F(zi,11); U=F(zi+2,9); Cw=F(2*zi-1,13)
        M=F(1,17);L=F(3,7);theta=F(zi+4,19)
        f=lambda x: theta+x+F(1,5)*x*x*x
        fscaled=lambda y:k*f(y/k)
        old=Cw+L*(U/k)+M*f(z)
        new=k*Cw+L*U+M*fscaled(k*z)
        assert new==k*old
        checks+=1

out=dict(status='PASS',checks=checks,
         scope='Exact rational cost, clipped-width, generic projection, and whole-block gauge algebra only. Source proofs are in the accompanying audit notes.',
         late_omission_rows=rows)
Path(__file__).with_name('extended_cost_schedule_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print('PASS:',checks,'exact rational checks')
for row in rows[:6]:print(row['P'],row['late_omission'],float(F(row['late_omission'])))
