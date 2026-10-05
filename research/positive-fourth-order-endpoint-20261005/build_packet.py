from pathlib import Path
import hashlib,json
P=Path(__file__).resolve().parent
OLD=Path('/workspace/shared/positive-endpoint-mean-join-20261005')
MEAN=Path('/workspace/shared/dyadic-prefix-law-join-20261005')
assert hashlib.sha256((OLD/'MANIFEST.json').read_bytes()).hexdigest()=='af11e79e0b121f38ae9ce3559e9cd6d813d61e1e13a44514622e96e7ea1638f5'
assert hashlib.sha256((MEAN/'MANIFEST.json').read_bytes()).hexdigest()=='130131317183c6f00ea3566c187296912db7850081904322215a3ec9d5eb37f2'
s=(OLD/'MEAN-COVARIANCE-SKEW-ENDPOINT.md').read_text()
s=s.replace('# A positive grade-39/10 endpoint from the canonical mean','# A positive log-grade-four endpoint from the completed canonical mean')
s=s.replace('2026-10-05. A bounded general-C2 construction, conditional on the literal fixed-order native guards of the pinned sources. The theorem below composes the outer reverse-OU schedule explicitly. No terminal LAW is promoted to a RAW force, and no derivative-valued producer is added.','2026-10-05. A new bounded general-C2 construction using the completed dyadic-prefix canonical m3 mean. The finite mean bank is actually replaced; the complete full-covariance and negative-skew banks and their positive reference proofs are retained. The theorem below composes the reverse-OU schedule through its separately priced positive terminal. All literal fixed-order native/public-log/precision guards are required. Imported native compilers are not numerically executed by this packet.')
a=s.index('Let g=grad U, g(0)=0, and 0<=Dg<=A I.')
b=s.index('The exact analytical carrier',a)
s=s[:a]+'''Let g=grad U, g(0)=0, and 0<=Dg<=A I. Write rho=3/4. In the declared nonempty sufficiently-small-A/public-log window there is a finite positive original-g-VALUE endpoint program E with

    W2(Law(E),mu_U) <= Lambda_ep A^4 sqrt(D)+e_absolute,
    mu_U(dx) proportional to exp(-|x|^2/2-U(x)) dx.       (1)

Lambda_ep>=1 is ONE fixed polynomial of a frozen public-log budget at the declared fixed native orders and complete graph dimensions. It includes every stage and every actual caller, precision, source-version, gap, count and readout parameter required by the imported interfaces. It stays at leading order. There is no numerical-constant C A^4 assertion and no strict-power log-absorption trick.

At this bounded fixed order, original VALUE counts and complete root dimension/D are fixed public-log polynomials under the literal imported model. In the inherited physical posterior normalization, the external sqrt(A) readout gives Lambda_ep A^(9/2) sqrt(D), with all physical mode/numerical floors restored at their actual weights. This is a scaling corollary under that normalization, not a redefinition of the present A or a general extra half-order for arbitrary targets.

The precise new gate is the optional completed own-mean LAW of ACTUAL-DYADIC-PREFIX-LAW-JOIN.md. The previous endpoint certificate stopped at grade 39/10 because its mean bank stopped below four. Its fixed-buffer covariance, skew, full-current, exact bridge and terminal interfaces are reused below; their grade-four substantive debts remain paid. No inference from the mean exponent alone is used.

'''+s[b:]
a=s.index('M. The optional unit-variance own-mean completion')
b=s.index('\nR. The existing full',a)
s=s[:a]+'''M. Execute the optional unit-variance own-mean completion of the NEW sealed dyadic-prefix canonical source, including its entire final raw-split compiler. Its analytical target is N(m3(Z),I), with integrated conditional error

    Lambda_M A^4 sqrt(D)+e_M.                         (M1)

This follows by adding the source's L2 mean error and its separate own-mean completion error. The intermediate target N(E[T_source|Z],I) is translated to N(m3(Z),I) using exactly that L2 standard-Z estimate; neither equality of the two means nor a pointwise-Z error is asserted.

At each normalized source radius a, this literal bank uses

    w=a, h_inner=a^2, eta=a^2 K_log^2,
    epsilon_m=a^2, epsilon_Sigma=a, epsilon_outer=a^3,
    prefix mean order 4, prefix square padding a^(1/2),
    tail mean/coefficient-mean order 6 with FULL native-radius padding,
    corrected negative-fourth Gram paths order 6,
    cubic/mixed-K/true-quartic selected pairs order 5,
    final own-mean order 4 with padding a.              (M2)

K_log=C L^b is a fixed sufficiently large public-log power chosen by the mean packet's finite graph admission proof; b is not the number of reverse-OU stages. For this section a=A. Later a=alpha_j at EVERY invocation. The graph has the independent reflected/anchored prefix mean and dyadic covariance banks with exact total h_inner^2 carrier; covariance restoration is consumed at the original Gaussian endpoint pair before the actual coherent caller shift. Every owned prefix root and its retained-caller source port survive in the executing mean graph. It has the original corrected-Gram/quartic five-bank tail, common-carrier rotation with all perpendicular roots, and complete final own-mean completion. Replacing only a formal target or scalar exponent would not define this M.

The mean proof retains the undominated order-six tail prior: after terminal Lipschitz readout its finite outer sum is a^7 sum omega v^(-5/2) <= C[a^(9/2)+a^4/K_log^3]. It does not use the invalid pointwise order-six prior domination at v near a^2 K_log^2. The RAW endpoint contribution a^4 K_log and all other public factors remain in Lambda_M. Section 6.1 verifies that these admissions and factors are compatible with the complete endpoint schedule.
'''+s[b:]
s=s.replace('C A^p sqrt(D)+Lambda A^q sqrt(D)\n                                +Lambda A^4 sqrt(D)+e_T. (8)','Lambda_T A^4 sqrt(D)+e_T.                         (8)')
a=s.index('    W2(K_new(r,t;z),K_exact(r,t;z))')
b=s.index('\n## 6.',a)
s=s[:a]+'''    W2(K_new(r,t;z),K_exact(r,t;z))
      <= Lambda_j (Delta/t) A^4 s^7 sqrt(D)
                         +b_rt |R_M|+e_(r,t)(z).       (10)

The factor 2 from sqrt(Delta)=2Delta/s is included in Lambda_j: sqrt(v0) alpha^4=2(Delta/t)A^4 s^7. The new mean theorem is applied on the fresh STANDARD internal carrier Z; the exterior z may have an approximate non-Gaussian entering law. Anchoring puts every finite x_M in the same Hessian/first class. All actual caller and finite-precision profiles remain those of f and x_M; no uniform unbounded-caller floor is asserted. Lambda_j is evaluated from the ENTIRE new mean graph plus the unchanged endpoint covariance/skew banks at alpha=A s^2, then bounded by the common Lambda_ep.
'''+s[b:]
a=s.index('Choose r_j^2=1-rho^j')
b=s.index('\n## 7.',a)
s=s[:a]+'''Choose r_j^2=1-rho^j, rho=3/4, and s_j^2=rho^j. Let J be the FIRST integer with

    s_J^5<=A^2, equivalently rho^J<=A^(4/5).            (11)

Initialize z_0=0. At r_0=0 both exact and finite first kernels do not read the incoming state: a_rt=0, and the finite anchored mode is zero. Thus initialization introduces no law debt. Every buffered stage has Delta_j=s_(j-1)^2/4 and uses (9). At r_J->1 execute the ORIGINAL separately audited positive conditional Stein terminal, with its finite mode, anchored original VALUES, positive finite clock rule, actual root tape and floors. Its conditional nonlinear law error is

    C A^2 s_J^5 sqrt(D)<=C A^4 sqrt(D).                (12)

Set the terminal's positive quadrature delta_terminal<=alpha_T=A s_J^2. Its normalized law error C(alpha_T^2+alpha_T delta_terminal) is O(alpha_T^2); multiplying by s_J gives (12). Any additional absolute numerical tolerance is chosen against this same readout using its original polylog clock theorem. Its full cost is M_terminal+1+n_terminal original VALUES (plus a separately requested canonical force only if actually requested), and its two D-dimensional Gaussian roots remain in the endpoint tape. No buffered mean/covariance/cubic native bank is called at h=1 or zero reserve.

The EXACT reference kernel has Wasserstein contraction r/t, including terminal contraction r_J. Conditional coupling and integration of (10) give E_j<=(r_(j-1)/r_j)E_(j-1)+epsilon_j. A buffered local error at stage j therefore receives exact final weight r_j. The intrinsic sum is

    sum_(j=1)^J r_j (Delta_j/r_j) Lambda_j A^4 s_(j-1)^7
      <= Lambda_ep A^4 sum_j Delta_j s_(j-1)^7
      <= Lambda_ep A^4/[4(1-rho^(9/2))].              (13)

A fixed numerical enlargement of Lambda_ep pays this denominator and (12). In particular no extra J loss is needed for the intrinsic errors. Equations (10)-(13), rather than mean improvement by itself, prove the normalized endpoint grade in (1).

Use the existing fixed conditional-mode iteration count M=3 at buffered stages. Its residual obeys |R_M|<=A^4 s^8 r |z|. After b_rt and final reference contraction, its price is at most A^4 Delta s^6 r ||z||_2. The exact OU marginals have second moment at most D. Choose propagated constant numerical floors at most c sqrt(D), and the sum of propagated coefficients of linear entering-caller floors at most a sufficiently small fixed c. The imported finite-precision model permits these budgets. Then ||z_j||_2<=sqrt(D)+E_j, the error recursion and discrete Gronwall (sum Delta_j<=1) give ||z_j||_2<=C sqrt(D), and the total finite-mode debt is O(A^4 sqrt(D)). A finite M=3 terminal mode has its own A^4 s_J^8 r_J||z_J||_2 bill, no larger; retain any separate inherited physical-mode tilt at its original weight. If budgets are left unrestricted, use C(sqrt(D)+e_budget) and retain its resulting mode/caller-floor terms instead of claiming the simplified moment bound.

The first stop obeys

    rho A^(4/5)<s_J^2<=A^(4/5),
    A^(9/5)<alpha_j=A s_(j-1)^2<=A  (1<=j<=J).        (14)

The terminal alpha_T=A s_J^2 lies between rho A^(9/5) and A^(9/5). Thus J=O(log(1/A)) and log(1/alpha_j)<=9/5 log(1/A), with only a fixed additive constant for the terminal. There is no hidden exponentially small clock. The minimal-stop lower bound is stated for the ACTUAL buffered alphas, not for a hypothetical call below the terminal.

### 6.1 One frozen public-log budget and actual admission at every call

Choose L>=1 to bound every imported public logarithm across all buffered stages and the terminal, including D, A, absolute precision, caller-profile parameters, complete source dimensions, finite graph sizes, actual gaps, coefficients, all source/version keys and service counts. Freeze all topologies/modes/orders/clocks/allocations, choose their counts with this budget, and close the fixed-order size inequalities by the mean packet's finite public-log bootstrap. Since there are O(log(1/A)) stages and each fixed-order dimension/count is polynomial in its local logarithms, enlarging L by a fixed multiple absorbs log J, log L and the new exact stage list. This is a pre-enumeration of one finite graph, not growing native order or an adaptive source search.

Use ONE K_log=C L^b for every mean call, replacing a in (M2) by its local alpha_j. The cutoff condition eta_j<=w_j/2 is alpha_j K_log^2<=1/2. It follows from A K_log^2<=1/2, but remains recorded for each literal graph. The smallest literal tail share u>=c alpha_j^2 K_log^2 has main normalized native radius <=C Lambda_0/K_log. The imported tree has finitely many positive powers of that entering radius times fixed public-log factors. Choose b exceeding the maximum fixed log-power/positive-K-decay ratio and C large enough for all these finitely many guard families. The normalized self-reserve first, source/readout/selector/gap/shield/clipping factors, and q-rescaled packet amplitudes are admitted by the same enumeration. Mixed-K has strict alpha_j^(1/3) K_log^(-2/3) margin; the prefix normalized radius has strict alpha_j^(1/2) margin; its actual endpoint first must still be <=1. The final own-mean order-four completion is also separately checked at its FULL actual radius Lambda alpha_j.

The outer endpoint covariance and cubic banks have fixed positive allocations and radii equal to fixed public factors times alpha_j or their documented positive powers. They require their own literal normalized-radius, root-amplitude, covariance-gap and clipping guards; the inner mean admission does not certify them automatically. Their log factors, as well as terminal quadrature and finite-mode floors, are included in L and Lambda_ep. Guard records are required for EVERY source occurrence and native reentry, using its complete active dimension and actual absolute floor. A generic statement that alpha_j<=A does not replace these checks.

This gives a nonempty sufficiently-small-A/public-log window at fixed external dimension/precision/caller parameters. One can choose A small enough for A K_log^2<=1/2, all strict-alpha margins and the original endpoint guards after the fixed b,C are chosen. Arbitrarily shrinking user-requested precision or increasing dimension can shrink that window. No dimension/precision-independent numerical threshold, actual numerical guard certificate for unspecified A,D, or production arbitrary-precision scalar Gauss solver is claimed. All literal imported guards remain hypotheses, while the finite log argument shows they are jointly compatible rather than formally contradictory.

Crucially, native admission and error accounting are separate. The integrated prior alpha_j^(9/2)+alpha_j^4/K_log^3 and the RAW alpha_j^4 K_log contribution remain inside Lambda_j alpha_j^4; no pointwise prior-domination inequality is reintroduced. This is why the leading logarithm cannot be erased by moving to the endpoint.
'''+s[b:]
s=s.replace('all five mean-join service tapes, smoothing/bridge roots, known-row perpendicular roots and the final own-mean compiler tapes','all five corrected tail service tapes, the entire new dyadic prefix mean and covariance banks, smoothing/bridge roots, known-row perpendicular roots and every final own-mean compiler reentry tape')
s=s.replace('The inner mean source\'s own final completion uses its separately admitted actual first/curl/energy/origin certificate.','The inner mean source\'s own final completion uses its new dyadic-prefix actual first/curl/energy/origin certificate, including the direct reflected-source retained first and the restored origin. Its radius is the actual full public-factor radius, not merely A.')
s=s.replace('For Q_M use the entire sealed five-bank mean graph and the final native raw-split own-mean recurrence; every reentered occurrence gets a fresh complete tape.','For Q_M use the NEW complete dyadic-prefix/five-tail mean graph and the final native raw-split own-mean recurrence. In particular its underlying source has d_source=D+sum_LAW[D+d_tail+d_pref]+sum_RAW[3D+d_pref]; d_pref includes BOTH complete prefix banks, and d_tail is recomputed at tail mean order six. Its prefix mean count is bounded by 2 J_m N_P,4+J_m+Q_capture/restore, while its dyadic covariance count is the full sum of 2 N_u Q_sq over every (level,position,r) label plus bottom/capture terms, with Q_sq=2K_f+4 and all native descendants. Every reentered final-compiler occurrence pays the ENTIRE source or source-plus-baseline recurrence and gets a fresh complete d_source tape. The source-level d_source is not substituted for the completed bank d_M. Counts include every positive clock, scalar filter, alias/version record, actual source capture, restore and changed-argument ancestor; no old order-25 count is retained.')
a=s.index('## 9. Boundary and extension')
s=s[:a]+'''## 9. Bounded conclusion and remaining debts

The mean target remains canonical m3 until the analytical P3-to-mu_U comparison is used. A posterior-force mean cannot be substituted into (3). The complete covariance LAW is consumed with its own positive variance; no known-carrier subtraction manufactures a RAW force. The new dyadic prefix covariance theorem is used only at its original stationary endpoint law inside M, before the uniform actual caller shift, exactly as in its sealed proof.

This packet closes a finite positive general-C2 endpoint at log-grade four, and inherited physical log-grade 9/2. It uses no new regularity of g beyond the original Hessian sandwich, no derivative-valued producer, no tensor oracle, no signed probability and no all-order recursion. The actual endpoint graph contains an outer negative cubic packet, not an unexecuted outer fourth packet; the quartic packet INSIDE M remains part of that mean construction. Its full fourth current is paid by the cubic consumer.

The canonical P3 analytical discrepancy, full covariance target/service restoration, H-to-F3 skew restoration, cubic-consumer fourth current, finite M=3 mode restoration and completed mean all retain grade-four allowances. These are certified debts of this specific construction, not lower bounds against other methods. Replacing the mean once more, adding more clocks, or merely adding a positive fourth packet does not erase the remaining debts. No >4 endpoint claim, no arbitrary-order recurrence and no order-uniform cost theorem is made.

All sealed sources remain unchanged. INPUT-PINS.json records this new join's exact source lineage; the separate independent audit reviews the actual replacement, positive reference, outer terminal and guards. The checks exercise algebraic identities and finite graph/accounting diagnostics, not the imported native compilers themselves. No external upload or communication is part of this packet.
'''
(P/'FOURTH-ORDER-POSITIVE-ENDPOINT.md').write_text(s)
# Preserve every original input pin and add exact new source snapshots and dependency pins.
pins={}
for root in [OLD,MEAN]:
 for pin in json.loads((root/'INPUT-PINS.json').read_text())['inputs']:
  assert hashlib.sha256(Path(pin['path']).read_bytes()).hexdigest()==pin['sha256'],pin
  pins[pin['path']]=pin
 for name in ['README.md','MANIFEST.json','INPUT-PINS.json']:
  path=root/name;pins[str(path)]={'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
for path in [OLD/'MEAN-COVARIANCE-SKEW-ENDPOINT.md',OLD/'ENDPOINT-GRAPH.json',OLD/'independent-audit/INDEPENDENT-ENDPOINT-JOIN-AUDIT.md',MEAN/'ACTUAL-DYADIC-PREFIX-LAW-JOIN.md',MEAN/'independent-audit/INDEPENDENT-PREFIX-LAW-AUDIT.md']:
 pins[str(path)]={'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
(P/'INPUT-PINS.json').write_text(json.dumps({'scope':'Exact read-only source snapshots; native compilers imported, not numerically executed; every literal guard retained','inputs':list(pins.values())},indent=2)+'\n')
g=json.loads((OLD/'ENDPOINT-GRAPH.json').read_text())
g['banks'][0]['owned_roots']='complete new dyadic prefix mean/covariance banks, all five corrected tail banks and rotations, plus every final own-mean compiler reentry/root'
g['banks'][0]['source_packet']=str(MEAN/'ACTUAL-DYADIC-PREFIX-LAW-JOIN.md')
g['banks'][0]['native_orders']={'prefix_mean':4,'tail_mean_and_coefficient_means':6,'corrected_negative_fourth_paths':6,'cubic_mixed_true_quartic':5,'final_own_mean':4}
g['banks'][0]['local_parameters']={'w':'alpha','h_inner':'alpha^2','eta':'alpha^2 K_log^2','K_log':'C L^b fixed pre-enumerated public-log power','prefix_square_padding':'alpha^(1/2)','final_padding':'alpha','tail_mean_padding':'FULL native radius'}
g['errors']['M']='Lambda_M alpha^4 plus absolute floors; leading logs retained'
g['outer_schedule']['stop']='rho^J<=A^(4/5), first such J'
g['outer_schedule']['local_grade']='Lambda_j (Delta/t) A^4 s^7'
g['outer_schedule']['buffered_alpha_range']='A^(9/5)<alpha_j<=A'
g['outer_schedule']['weighted_sum']='A^4/[4(1-rho^(9/2))], before harmless local factor 2'
g['root_dimension']='2D+d_M_completed+d_H+d_K+d_3 per buffered stage, plus the full terminal 2D tape; do not use uncompleted d_source for d_M'
g['result']={'normalized':'Lambda_ep A^4 sqrt(D)+absolute floors','physical':'Lambda_ep A^(9/2) sqrt(D)+restored physical floors under inherited sqrt(A) readout','leading_logs_retained':True,'native_compilers_executed':False}
(P/'ENDPOINT-GRAPH.json').write_text(json.dumps(g,indent=2)+'\n')
print('Built main proof, graph and',len(pins),'exact input pins')
print('main_sha256',hashlib.sha256((P/'FOURTH-ORDER-POSITIVE-ENDPOINT.md').read_bytes()).hexdigest())
