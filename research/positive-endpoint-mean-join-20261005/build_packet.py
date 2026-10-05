from pathlib import Path
import hashlib,json
P=Path(__file__).resolve().parent
B=Path('/workspace/shared/recovery-20261004/cost/endpoint-stein-quadrature-20261004')
G=B/'positive-law-repair/higher-cumulant-gate'
paths=[
Path('/workspace/shared/quartic-corrected-mean-join-20261005/MANIFEST.json'),
Path('/workspace/shared/quartic-corrected-mean-join-20261005/ACTUAL-QUARTIC-CORRECTED-MEAN-JOIN.md'),
Path('/workspace/shared/quartic-corrected-mean-join-20261005/independent-audit/INDEPENDENT-JOIN-AUDIT.md'),
G/'THIRD-ORDER-FULL-REVERSE-OU-LAW-AND-ENDPOINT.md',
G/'same-carrier-feedback/SAME-CARRIER-P3-MEAN-AND-FULL-COVARIANCE-GATE.md',
G/'same-carrier-feedback/FULL-COVARIANCE-CLOSURE-AND-RESUMMED-MEAN-FRONTIER.md',
G/'same-carrier-feedback/RESOLVENT-COVARIANCE-STABILITY-AND-FINITE-J-KERNEL.md',
G/'same-carrier-feedback/GAUSSIAN-COVARIANCE-MIXTURE-ONE-ENERGY-LAW.md',
G/'order-reentry/RECTANGULAR-FIRST-COEFFICIENT-MEAN-AND-GRAM-RETURN.md',
G/'same-carrier-feedback/independent-resolvent-audit/INDEPENDENT-RECTANGULAR-FILTER-AND-GRAM-RETURN-AUDIT.md',
G/'same-carrier-feedback/mixed-k-current/MIXED-K-CURRENT-FINITE-CLOSURE.md',
G/'same-carrier-feedback/mixed-k-current/INDEPENDENT-FINITE-C0-CUBIC-RESERVE-AUDIT.md',
G/'rank3-continuation/POSITIVE-QUADRATIC-REFERENCE-FOR-SKEW-GAUSSIANIZATION.md',
G/'rank3-continuation/EXACT-FIVE-CLOCK-NATIVE-SOURCE-ADAPTER.md',
Path('/workspace/shared/shrinking-buffer-skew-join-20261005/MANIFEST.json'),
Path('/workspace/shared/shrinking-buffer-skew-join-20261005/POSITIVE-SKEW-SHRINKING-BUFFER-JOIN.md'),
Path('/workspace/shared/shrinking-buffer-skew-join-20261005/independent-audit/INDEPENDENT-BUFFERED-SKEW-JOIN-AUDIT.md'),
Path('/workspace/shared/fourth-cumulant-return-20261005/comparison/FOURTH-CONDITIONAL-CUMULANT-COMPARISON.md'),
Path('/workspace/shared/fourth-cumulant-return-20261005/comparison/POSITIVE-FOURTH-CUMULANT-BUFFERED-CONSUMER.md'),
B/'ENDPOINT-STEIN-PACKET-AND-POLYLOG-QUADRATURE.md',
B/'positive-law-repair/POSITIVE-SECOND-ORDER-LAW-AND-QUADRATIC-AMPLIFIER.md',
B/'positive-law-repair/independent-audit/INDEPENDENT-POSITIVE-LAW-AUDIT.md',
]
# Input hashes identify source-qualified theorem dependencies, not executable native code.
pins=[]
for f in paths:
    if not f.is_file():
        print('MISSING',f)
        continue
    pins.append({'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
assert pins[0]['sha256']=='dde7c7e3220b18f65d13039587af6cd9690db0264f10f87c00fcbcb6870ccad9'
(P/'INPUT-PINS.json').write_text(json.dumps({'scope':'unchanged sealed sources; literal native compiler guards remain mandatory','inputs':pins},indent=2)+'\n')
graph={
 'stage_h':'1/2','retained':['original exterior caller and versions','fresh standard Z'],
 'banks':[
  {'id':'M','target':'N(m3(Z),I)','readout':'-1/2','baseline_after_readout':'1/4','input_type':'completed own-mean LAW','owned_roots':'entire sealed five-bank mean source plus complete final compiler and all perpendicular/response/filter roots'},
  {'id':'H','target':'N(0,1/4 I+Cov(H_j|Z))','readout':'1/2','baseline_after_readout':'1/16','input_type':'true-orientation Gram LAW','owned_roots':'all coarse, p, filter, J_Q, response, means and keeps'},
  {'id':'K','target':'N(0,1/4 I+2 Sym R2 K(Z))','readout':'1/2','baseline_after_readout':'1/16','input_type':'mixed covariance LAW','owned_roots':'all four-clock nodes, three chronological adapters, coarse/public/filter/fill roots'},
  {'id':'V3','target':'negative cubic reference for -h^3 kappa3(H|Z)/6','readout':'1','baseline_after_readout':'1/8','input_type':'positive completed cubic VALUE packet','owned_roots':'all five-clock coarse roots and C0,C1,C0 banks, originals, shields, publics, complements and fills'},
  {'id':'keep','target':'N(0,1/4 I)','readout':'1','baseline_after_readout':'1/4','input_type':'untouched standard Gaussian','owned_roots':'fresh D-root'}],
 'independence':'M,H,K,V3,keep banks mutually independent after same Z/exterior captures; repeated internal banks also independent unless their native graph specifies a shared root',
 'covariance_restorations':['Cov(H_j|Z)+2 Sym R2 K to Cov(F2|Z): L2 standard-Z HS O(A^4 sqrt(D))','Cov(F2|Z) to Cov(F3|Z): L2 standard-Z HS O(A^4 sqrt(D))'],
 'reference_regression':'consume all service tapes; use analytic reference Gaussian S covariance 1/2 I+h^2 Cov(F3|Z); keep external 1/4 unread; pay A^6 regression and output-slot symmetrization',
 'errors':{'M':'C A^(39/10)+Lambda A^(59/15)+Lambda A^4','covariance':'Lambda A^4','third_target':'C A^4','cubic_prior':'Lambda A^5','cubic_feedback':'Lambda A^6','positive_skew_consumer':'Lambda A^4','P3_to_true_target':'C A^4'},
 'outer_schedule':{'rho':'3/4','r_j_squared':'1-rho^j','Delta':'(1-r_previous^2)/4','stop':'rho^J<=A^(19/25)','terminal':'original positive conditional Stein with physical error A^2 s_J^5','exact_contraction':'r/t','local_grade':'(Delta/t) A^k s^(2k-1)'},
 'forbidden':['terminal LAW as RAW','posterior force mean substituted for canonical m3','reference coupling root treated as actual private root','integrated old private root retained as observer','original HVP or cumulant tensor as producer'],
 'root_dimension':'2D+d_M+d_H+d_K+d_3 per buffered stage, plus full terminal tape',
 'cost':'sum complete native original-VALUE serial counts over J, plus finite modes/anchors/restoration; no unit-cost analytical services',
 'guards':'all actual literal imported fixed-order source radii, complete dimensions, public factors, source versions, covariance gaps, root and clipping amplitudes, precision and caller paths at every alpha_j'
}
(P/'ENDPOINT-GRAPH.json').write_text(json.dumps(graph,indent=2)+'\n')
print('pins',len(pins))
