#!/usr/bin/env python3
"""Deterministic diagnostics for the independent modulated-jet audit.
These checks supplement, and do not replace, the mathematical proof.
"""
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np


# Publication-only source locator and original/public pin mapping.
import sys
_ARCHIVE=next(p for p in Path(__file__).resolve().parents if (p/'publication_sources.py').is_file())
sys.path.insert(0,str(_ARCHIVE))
from publication_sources import source_path, source_label, verify_source_pin, pin_report

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
rng = np.random.default_rng(20261005)
checks = []

def check(name, actual, bound, relation='le', tolerance=1e-11):
    actual, bound = float(actual), float(bound)
    ok = actual <= bound + tolerance if relation == 'le' else actual >= bound - tolerance
    checks.append({'name':name,'actual':actual,'bound':bound,'relation':relation,'pass':bool(ok)})
    assert ok, checks[-1]

def max_cut(tensor):
    rank = tensor.ndim
    maximum = 0.0
    for mask in range(1, 2**rank-1):
        left = [i for i in range(rank) if (mask >> i) & 1]
        right = [i for i in range(rank) if not ((mask >> i) & 1)]
        if len(left) > len(right):
            continue  # Transposition preserves the singular values.
        dleft = math.prod(tensor.shape[i] for i in left)
        matrix = tensor.transpose(left+right).reshape(dleft, -1)
        maximum = max(maximum, float(np.linalg.norm(matrix, 2)))
    return maximum

def perturbation(x, sigma, eta=1/8):
    n = len(x)-1
    ans = np.zeros((n+1,n+1))
    co, si = np.cos(x[1:]/sigma), np.sin(x[1:]/sigma)
    ans[0,0] = -eta*sigma*sigma*np.cos(x[0])*co.sum()
    ans[np.arange(1,n+1),np.arange(1,n+1)] = -eta*np.cos(x[0])*co
    ans[0,1:] = ans[1:,0] = eta*sigma*np.sin(x[0])*si
    return ans

def tensor_entries(x, sigma, eta=1/8):
    n = len(x)-1
    d = math.exp(-(1+sigma*sigma)/2)
    out = np.zeros((n+1,)*4)
    out[(0,)*4] = eta*d*sigma**4*np.cos(x[0])*np.cos(x[1:]/sigma).sum()
    for i in range(1,n+1):
        for zeros in range(4):
            val = eta*d*sigma**zeros
            if zeros % 2:
                val *= -np.sin(x[0])*np.sin(x[i]/sigma)
            else:
                val *= np.cos(x[0])*np.cos(x[i]/sigma)
            for positions in itertools.combinations(range(4), zeros):
                key = [i]*4
                for p in positions: key[p] = 0
                out[tuple(key)] = val
    return out

def tensor_fourier(x, sigma, eta=1/8):
    n = len(x)-1
    d = math.exp(-(1+sigma*sigma)/2)
    out = np.zeros((n+1,)*4)
    for i in range(1,n+1):
        for sign in [-1,1]:
            k = np.zeros(n+1); k[0]=1; k[i]=sign/sigma
            out += eta*d*sigma**4/2*np.cos(k@x)*np.einsum('a,b,c,d->abcd',k,k,k,k)
    return out

for n in [3,4,5,7]:
    sigma = n**-0.5
    for rep in range(3):
        x = rng.normal(size=n+1)
        K = perturbation(x,sigma)
        check(f'hessian_op_n{n}_{rep}',np.linalg.norm(K,2),1/4)
        check(f'hessian_min_n{n}_{rep}',np.linalg.eigvalsh(0.5*np.eye(n+1)+K)[0],1/4,'ge')
        T = tensor_entries(x,sigma)
        check(f'fourier_derivative_n{n}_{rep}',np.max(np.abs(T-tensor_fourier(x,sigma))),1e-11)
        check(f'center_all_cuts_n{n}_{rep}',max_cut(T),2)

for n in [3,4]:
    sigma=n**-0.5; r2=math.sqrt(1-2/n); rleaf=1/math.sqrt(2)
    bank=rng.normal(size=(5,n+1)); other=rng.normal(size=(5,n+1))
    t=1-1/n**3
    bankp=t*bank+math.sqrt(1-t*t)*other
    def sites(q):
        Y=math.sqrt(3)/2*q[0]
        X=Y/2+math.sqrt(3)/2*q[1]
        return rleaf*Y+q[2]/2, r2*X+sigma*q[3], rleaf*X+q[4]/2
    a,x,c=sites(bank); ap,xp,cp=sites(bankp)
    T,Tp=tensor_entries(x,sigma),tensor_entries(xp,sigma)
    leafheat=math.exp(-(n+1)/8)
    leaves=[0.5*np.eye(n+1)+leafheat*perturbation(q,sigma) for q in [a,c,ap,cp]]
    gamma=0.75*(1-2/n)
    J=gamma*np.einsum('au,cv,buvh,dw,fz,ewzh->abcdef',leaves[0],leaves[1],T,leaves[2],leaves[3],Tp,optimize=True)
    J0=gamma/16*np.einsum('bach,edfh->abcdef',T,Tp,optimize=True)
    check(f'six_tree_cuts_n{n}',max_cut(J),4)
    check(f'leaf_difference_cuts_n{n}',max_cut(J-J0),4*leafheat)
    check(f'leaf_difference_HS_n{n}',np.linalg.norm((J-J0).ravel()),4*leafheat*math.sqrt(n+1))
    c_n=gamma/1024*math.exp(-(1+1/n))
    for i in range(1,n+1):
        explicit=c_n*(math.cos(x[0])*math.cos(xp[0])*math.cos(x[i]/sigma)*math.cos(xp[i]/sigma)+sigma*sigma*math.sin(x[0])*math.sin(xp[0])*math.sin(x[i]/sigma)*math.sin(xp[i]/sigma))
        check(f'diagonal_contraction_n{n}_i{i}',abs(J0[(i,)*6]-explicit),1e-11)

# Closed-form scalar diagnostics: no Monte Carlo is used for the variance check.
for n in [16,32,64,128,129,256,1000,10000,1000000]:
    v=15/16-7/(8*n); gamma=0.75*(1-2/n)
    cn=gamma/1024*math.exp(-(1+1/n))
    check(f'gamma_lower_n{n}',gamma,0.5,'ge')
    check(f'cn_lower_n{n}',cn,2**-14,'ge',0)
    check(f'cos_square_sd_n{n}',(1-math.exp(-4*v))/math.sqrt(8),1/3,'ge')
    check(f'sd_M_bound_n{n}',1/3-2*n**-1.5,1/4,'ge')
    check(f'conditional_sd_bound_n{n}',1/12-1/(2*n),1/32,'ge')
    for j in range(11):
        t=1-2*n**-3*j/10
        low=math.exp(-v*(1-t)*n)
        high=math.exp(-v*(1+t)*n)
        mu=(low+high)/2; nu=(low-high)/2
        check(f'mu_lower_n{n}_t{j}',mu,1/3,'ge')
        check(f'nu_upper_n{n}_t{j}',abs(nu),0.5)
        # U-V and U+V are independent; the displayed scalar is a linear
        # combination of their cosines. This computes its exact variance.
        varminus=(-math.expm1(-2*v*(1-t)))**2/2
        varplus=(-math.expm1(-2*v*(1+t)))**2/2
        conditional_var=((mu+nu/n)/2)**2*varminus+((mu-nu/n)/2)**2*varplus
        check(f'exact_conditional_sd_n{n}_t{j}',math.sqrt(conditional_var),1/32,'ge')
        check(f'exact_J0_sd_n{n}_t{j}',cn*math.sqrt(n*conditional_var),2**-19*math.sqrt(n),'ge',0)
    if n>=128:
        leaferror=6*math.exp(-(n+1)/8)
        check(f'leaf_error_n{n}',leaferror,2**-20, tolerance=0)
        check(f'final_sd_n{n}',2**-19-leaferror,2**-20,'ge',0)

source_paths=[
    ROOT/'MODULATED-C2-SIX-TREE-COUNTEREXAMPLE.md',
    ROOT.parent/'COVARIANCE-TEST-PORT-REDUCTION.md',
    ROOT.parent/'TIGHT-BRIDGE-LEDGER-AND-REROOT-GATE.md',
    ROOT/'../../../EXACT-FIVE-CLOCK-NATIVE-SOURCE-ADAPTER.md',
    ROOT/'../../../grouped-kernel-response/PAIRED-ROTATION-COVARIANCE-ALL-CUT-LEMMA.md',
    HERE/'INDEPENDENT-MODULATED-JET-AUDIT.md',
    Path(__file__),
]
sources={source_label(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
report={'status':'PASS','description':'Deterministic diagnostics, not a substitute for the proof.','check_count':len(checks),'sources_sha256':sources,'checks':checks}
out=HERE/'modulated_jet_audit_checks.json'
out.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'status':report['status'],'check_count':len(checks),'report':str(out)}))
