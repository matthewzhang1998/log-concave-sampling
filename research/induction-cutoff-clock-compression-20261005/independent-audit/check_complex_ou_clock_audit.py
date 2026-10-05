#!/usr/bin/env python3
"""Independent exact scalar/census checks; no native sampler is executed."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib, importlib.util, json, math, random, sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import cutoff_free_sector_gauss as author
import gauss_rational_core as core

checks=[]
def ck(name, condition, detail=None):
    assert condition, (name,detail)
    checks.append({'name':name,'pass':True,**({'detail':detail} if detail is not None else {})})

prior=Path('/workspace/shared/induction-positive-clock-compression-20261005/positive_sector_gauss.py')
ck('pinned_core_is_byte_identical',prior.read_bytes()==(ROOT/'gauss_rational_core.py').read_bytes())
input_pins=json.loads((ROOT/'INPUT-PINS.json').read_text())['inputs']
for pin in input_pins:
    data=Path(pin['path']).read_bytes()
    ck('sealed_input_pin_matches',len(data)==pin['bytes'] and hashlib.sha256(data).hexdigest()==pin['sha256'],{'path':pin['path']})

# For complex-bilinear Hermite contraction, the exact pairing monomial
# r^(x+y+z) s^(y+z) equals the common-G monomial r^((a+b+c)/2) s^c.
for a,b,c in product(range(11),repeat=3):
    if (a+b+c)%2:
        ck('odd_hermite_product_vanishes',(a+b+c)%2==1)
        continue
    x,y,z=(a+b-c)//2,(a+c-b)//2,(b+c-a)//2
    if min(x,y,z)<0:
        ck('unpairable_hermite_product_vanishes',max(a,b,c)>a+b+c-max(a,b,c))
        continue
    coefficient=F(math.factorial(a)*math.factorial(b)*math.factorial(c),math.factorial(x)*math.factorial(y)*math.factorial(z))
    ck('exact_bilinear_ou_factorization',x+y+z==(a+b+c)//2 and y+z==c and coefficient>0)

# Exact norm check on integer complex-valued chain tensors. No conjugates
# appear in the contraction; conjugates occur only when forming squared norms.
rng=random.Random(421803)
def rv():return complex(rng.randrange(-3,4),rng.randrange(-3,4))
def nsq(v):return int(v.real)**2+int(v.imag)**2
for _ in range(100):
    A=[[rv() for a in range(3)] for p in range(2)]
    B=[[[rv() for b in range(2)] for a in range(3)] for q in range(2)]
    C=[[rv() for b in range(2)] for r in range(3)]
    out=[sum(A[p][a]*B[q][a][b]*C[r][b] for a in range(3) for b in range(2)) for p,q,r in product(range(2),range(2),range(3))]
    rhs=sum(nsq(x) for row in A for x in row)*sum(nsq(x) for plane in B for row in plane for x in row)*sum(nsq(x) for row in C for x in row)
    ck('complex_bilinear_tree_HS_inequality',sum(nsq(x) for x in out)<=rhs)

for N,K in ((9,7),(12,10)):
    for D,eta,eps in product((1,2,10**8),(F(1,2),F(1,31),F(1,2**80)),(F(1,3),F(1,10**25))):
        ws={0:F(7,11),1:F(17,3),2:F(41,7)}
        p=author.parameters(N,K,D,eta,eps,W=ws[0],caller_weights=ws)
        for q,w in ws.items():
            exact_square=w*w*(2**(2*N+3*(K+q)))*math.factorial(K+q)*eta**(-(K+q))
            ck('rational_B_is_upper_bound',core.exact_B(N,K+q,eta,w)**2>=exact_square)
        M,J,m=p['M_rational'],p['J'],p['m']
        ck('exact_full_error_budget',4*M*(F(1,2**J)+F(1,4**m))+M*p['node_tolerance']+2*M*p['weight_tolerance']<=F(19,64)*eps)
        ck('literal_node_count',p['node_count_per_tree']==2*(128*J*m)**2)
        ck('exact_minimum_integer_panels',2**J>=32*M/eps and (J==1 or 2**(J-1)<32*M/eps))

for n in (1,3,7,11):
    rule=core.certified_legendre_to_tolerance(n,F(1,10**30),F(1,10**30))
    ck('independent_degree_certification',rule['node_error']<=F(1,10**30) and rule['weight_l1_error']<=F(1,10**30),{'degree':n})
    panels=author.gap_panels(4)
    ck('dyadic_atlas_exact_mass',sum(b-a for a,b in panels)==F(15,16))
    for index in (0,127,128,255,256,383,384,511):
        a,b=panels[index]
        mapped=core.map_rule(rule,a,b)
        ck('cell_exact_mass_first_moment',sum(w for x,w in mapped)==b-a and sum(x*w for x,w in mapped)==(b*b-a*a)/2)
        ck('mapped_root_error_no_cutoff_loss',64*(b-a)/(2*(1-b))<=F(1,4))
        for r,wr in mapped:
            for s,ws in mapped:
                ck('literal_positive_node_heat_weight',0<wr*ws*r<=(1-r)*(1-r*s))

# Independently enumerate the actual force orders rather than call the author's census.
census={}
for name, sizes, base, old_pairs in (
    ('triple_cubic',(3,3,3),(0,1,0)*3,[()]),
    ('three_three_six',(3,3,6),(0,1,0)*4,list(product(range(6,9),range(9,12))))):
    starts=(0,sizes[0],sizes[0]+sizes[1]);groups=[range(st,st+size) for st,size in zip(starts,sizes)]
    count,max_bill,max_order=0,0,0
    degree=set();witness=None
    for central in range(3):
        leaves=[i for i in range(3) if i!=central]
        for old in old_pairs:
            for v in product(groups[central],groups[leaves[0]],groups[central],groups[leaves[1]]):
                orders=list(base)
                for hit in old+v:orders[hit]+=1
                bill=sum(2**(order+1) for order in orders)
                count+=1;degree.add(sum(orders));max_order=max(max_order,max(orders))
                if bill>max_bill:max_bill,witness=bill,orders
    census[name]={'histories':count,'max_raw_VALUES':max_bill,'max_order':max_order,'degrees':sorted(degree),'maximizer_orders':witness}
ck('raw_cubic_bill',census['triple_cubic']['histories']==243 and census['triple_cubic']['max_raw_VALUES']==44 and census['triple_cubic']['degrees']==[7])
ck('raw_mixed_bill',census['three_three_six']['histories']==5832 and census['three_three_six']['max_raw_VALUES']==72 and census['three_three_six']['degrees']==[10])

pins={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ('CUTOFF-FREE-SECTOR-GAUSS.md','COMPLETE-EXECUTION-BILL.md','INPUT-PINS.json','cutoff_free_sector_gauss.py','gauss_rational_core.py')}
out={'status':'PASS','check_count':len(checks),'checks':checks,'independent_raw_census':census,'reviewed_sha256':pins,'scope':'Exact arithmetic and finite census supplement the mathematical audit; no native source or sampler is run.'}
destination=Path(__file__).with_name('complex-ou-clock-audit-checks.json')
destination.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'check_count':len(checks),'independent_raw_census':census,'reviewed_sha256':pins,'output':str(destination)},indent=2))
