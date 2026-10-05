#!/usr/bin/env python3
"""Independent finite exponent audit. Pure combinatorics, not a source oracle."""
from itertools import product
from fractions import Fraction as F
from pathlib import Path
import json, hashlib, math

OUT=Path(__file__).resolve().parent
ROOT=OUT.parent

def local_single(p):
    """(gap power, extra eta power) using one mass at each cubic clock."""
    high=[j for j,v in enumerate(p) if v>2]
    assert len(high)<=1
    if not high:
        return F(0),F(0)
    j=high[0]
    a=F(p[j],2)-1
    # A term in min(x_other,z_j)^(-a) charges exactly one clock.
    losses=[max(F(0),a-1)]
    for v in range(3):
        if v!=j:
            losses.append(max(F(0),F(p[v],2)+a-1))
    return a,max(losses)

three=[]
for d in range(4):
    count=0
    for hits in product(range(3),repeat=d):
        p=[0,1,0]
        for h in hits:p[h]+=1
        a,e=local_single(p)
        assert e==0
        assert a<=max(F(0),F(d-1,2))
        count+=1
    three.append({'d':d,'ordered_patterns':count})

mixed=[]
for d in range(4):
    count=0; max_pc=0;max_pl=0;max_diag=F(0); max_full_local_eta=F(0)
    branches={'center_heat':0,'unheated_squared_clocks':0}
    for old_a,old_b in product(range(3),repeat=2):
        for hits in product(range(6),repeat=d):
            p=[[0,1,0],[0,1,0]]
            p[0][old_a]+=1;p[1][old_b]+=1
            da=[0,0]
            for h in hits:
                i,j=divmod(h,3);p[i][j]+=1;da[i]+=1
            powers=[];losses=[]
            for i in range(2):
                a,e=local_single(p[i]);powers.append(a);losses.append(e)
                assert a<=F(da[i],2)
                assert e<=F(max(0,da[i]-2),2)
            assert sum(powers)<=F(d,2)
            assert sum(losses)<=F(max(0,d-2),2)
            max_full_local_eta=max(max_full_local_eta,sum(losses))
            P=[p[0][j]+p[1][j] for j in range(3)]
            assert sum(P)==4+d
            max_pc=max(max_pc,P[1]);max_pl=max(max_pl,P[0],P[2])
            if P[1]>4:
                branches['center_heat']+=1
                a=F(P[1],2)-2
                # The canonical old floor omits center clock 1.
                assert a<=F(d,2)
                assert a<=2
                assert all(F(P[j],2)+a<=2 for j in (0,2))
                # Integrating the old bridge is the only power loss.
                diag_loss=max(F(0),a-1)
            else:
                branches['unheated_squared_clocks']+=1
                diag_loss=sum(max(F(0),F(P[j],2)-2) for j in range(3))
            assert diag_loss<=F(max(0,d-2),2)
            max_diag=max(max_diag,diag_loss)
            count+=1
    assert count==9*6**d
    mixed.append({'d':d,'ordered_patterns':count,'maximum_center_power':max_pc,
                  'maximum_leaf_power':max_pl,'maximum_diagonal_eta_loss':str(max_diag),
                  'maximum_full_local_eta_loss':str(max_full_local_eta),'proof_branches':branches})

# Weighted Pruefer identity with literal force-hit multiplicities.
pruefer=[]
for weights in [(3,3,3),(3,3,6),(1,2),(1,2,3,4),(2,3,4,5,6),(1,2,3,4,5,6)]:
    n=len(weights);summed=0
    for seq in product(range(n),repeat=n-2):
        deg=[1]*n
        for i in seq:deg[i]+=1
        summed+=math.prod(weights[i]**deg[i] for i in range(n))
    expected=math.prod(weights)*sum(weights)**(n-2)
    assert summed==expected
    pruefer.append({'force_counts':weights,'exact_histories':summed})

# Arbitrary-d local lemma: enumerate aggregate powers rather than ordered hits.
arbitrary_local=[]
for d in range(33):
    count=0;largest_loss=F(0)
    for p1 in range(d+2):
        for pc in range(1,d+2-p1):
            p=[p1,pc,1+d-p1-pc]
            assert min(p)>=0 and sum(p)==1+d
            if max(p)>2:
                j=max(range(3),key=lambda j:p[j])
                a=F(p[j],2)-1
                others=[v for v in range(3) if v!=j]
                exponent_terms=[[F(p[v],2) for v in others]+[a]]
                for charged in others:
                    exponent_terms.append([F(p[v],2)+(a if v==charged else 0) for v in others])
                for exps in exponent_terms:
                    assert sum(exps)==F(d-1,2)
                    loss=sum(max(F(0),z-1) for z in exps)
                    assert loss<=max(F(0),F(d-3,2))
                    largest_loss=max(largest_loss,loss)
                assert a<=max(F(0),F(d-1,2))
            count+=1
    arbitrary_local.append({'d':d,'aggregate_patterns':count,'maximum_eta_loss':str(largest_loss)})

def compositions(total,length):
    if length==1:
        yield (total,)
    else:
        for x in range(total+1):
            for tail in compositions(total-x,length-1):
                yield (x,)+tail

rank_explicit=[]
for n in range(2,9):
    for q in range(6):
        total=n+q-2; half_E=max(F(0),F(n+q-4,2))
        local_count=0; bridge_count=0
        max_local=F(0);max_bridge=F(0)
        # Relax to every nonnegative r_i=d_i-1 with the required total.
        for r in compositions(total,n):
            k=sum(max(F(0),F(x-2,2)) for x in r)
            assert k<=half_E
            max_local=max(max_local,k);local_count+=1
        # Relax to every half-integer assignment to the n-1 bridge edges.
        for twice_b in compositions(total,n-1):
            prefix=0; k=F(0)
            for j,b in enumerate(twice_b,1):
                prefix+=b
                k=max(k,F(prefix,2)-j)
            assert k<=half_E
            max_bridge=max(max_bridge,k);bridge_count+=1
        assert max_local+max_bridge<=max(0,n+q-4)
        rank_explicit.append({'n':n,'q':q,'relaxed_hit_allocations':local_count,
                              'relaxed_bridge_allocations':bridge_count,
                              'maximum_local_loss':str(max_local),
                              'maximum_bridge_loss':str(max_bridge),
                              'claimed_e':max(0,n+q-4)})

files=['THREE-CUBIC-EXACT-HEAT.md','MIXED-THREE-THREE-SIX.md','FINITE-CUBATURE-AND-N-COPY-RULE.md','N-CUBIC-RANK-EXPLICIT-BOUNDS.md']
result={'status':'PASS','three_cubic_patterns':three,'mixed_patterns':mixed,
        'mixed_total_patterns':sum(x['ordered_patterns'] for x in mixed),
        'weighted_pruefer':pruefer,
        'arbitrary_local_aggregate_patterns':arbitrary_local,
        'arbitrary_local_total_patterns':sum(x['aggregate_patterns'] for x in arbitrary_local),
        'rank_explicit_checks':rank_explicit,
        'rank_explicit_total_allocations':sum(x['relaxed_hit_allocations']+x['relaxed_bridge_allocations'] for x in rank_explicit),
        'same_source_contract':'The executed allocation uses m_star=max_j m_j; every analytical hit estimate bounds those same widths and residual-root rows.',
        'reviewed_sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
        'scope':'Exact finite integer/rational exponent checks only. The audit text supplies the mathematical arguments. Native sources were not executed.'}
(OUT/'exact_exponent_ledger_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
