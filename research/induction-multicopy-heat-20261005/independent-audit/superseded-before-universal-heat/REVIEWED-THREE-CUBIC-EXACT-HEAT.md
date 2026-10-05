# Exact three-copy old-cubic heat and caller budget

2026-10-05. Constructive coefficient/source theorem. This note introduces no new smoothing of the target. It supplies an exact Gaussian redistribution, all original force-hit histories, scalar root-square and complete captured-caller budgets. Finite quadrature and native compilation have their own charged costs; no unproved public-log node count is asserted.

## 1. Target and retained-variable contract

Let B be the COMPLETE original standard coefficient bank. The original physical source-zero readout and the retained external caller Y have zero B row. Original finite five-clock cubic coefficients L_a(B;Y) have their literal original clocks, source versions, private shields, physical marks and affine rows. The target is

    kappa_B(L^(1), L^(2), L^(3)),  L^(j)=sum_a L_a^(j).

Expand all ordered history triples, including different histories which share B. Repeated clock labels remain repeated. The theorem does not replace a diagonal sum by independent clocks or pretend shared banks are independent. The usual 1/3! is inserted once if this is a log-generating coefficient.

Each cubic history has three original force vertices, initially C0-C1-C0, and three physical marks. The BKAR tree contributes four force-hit slots in total. Every resulting connected graph therefore has

    amplitude grade a=9, force count N=9, local width degree K=7,
    physical rank R=9, internal edges M=8.

There are three labelled replica trees and 3^4=81 ordered force-hit assignments per tree, before the original history and physical-permutation census. Repeated hits are included. A C0 leaf hit becomes a marked C1 vertex. The union of the three old force trees and the two bank edges is a marked spanning tree.

Temporary conditional comparison can retain B and new residual roots. The final bank-mean identity integrates them. It is not pointwise cancellation in the old B. No private tape or copy of B may be appended as an unattenuated observer. Every old/new mixed current remains in the common endpoint ledger.

## 2. Exact n-copy redistribution, of which n=3 is the first application

For a labelled replica tree T and edge parameters t_e in [0,1], let

    K_ii=1, K_ij=min_{e on path i--j} t_e,
    c=min_e t_e, delta_i=1-max_{e incident to i} t_e.

For n=1 this definition is not used. The usual partition identity gives

    K = integral_0^1 C(Pi_u) du,

where Pi_u consists of components after edges t_e<u are deleted. Subtract the full connected block on u in [0,c], and at every larger u subtract each singleton component. What remains is a sum of positive block matrices. Consequently

    H := K-c 11^T-diag(delta_i) >=0.                         (2.1)

This is stronger than merely K>=0. Let R be a known scalar square root of H. With the ACTUAL old B, fresh standard complete banks V_1,...,V_n and E_1,...,E_n, put

    B_i = sqrt(c) B + (R V)_i + sqrt(delta_i) E_i.           (2.2)

All new roots are independent of B and the old physical publics. The covariance of (B_i) is exactly K tensor I, and every marginal is the original standard bank. Original derivative contractions still use the original query injection rows, not sqrt(c) times those rows.

In copy i suppose the original force centers are mu_i(Y)+A_i B_i, with Gamma_i=A_i A_i^T. Choose ANY known diagonal D_i>=0 with Gamma_i-D_i>=0. Replace the final independent E_i contribution by a residual query bank of covariance delta_i(Gamma_i-D_i), and enlarge the independent force shields by

    new_t_iv^2 = old_sigma_iv^2 + delta_i (D_i)_vv.         (2.3)

The common center remains mu_i(Y)+sqrt(c) A_i B+A_i(RV)_i. Equation (2.3) plus the residual covariance exactly reconstructs delta_i Gamma_i plus the ORIGINAL private covariance. All means, within-copy covariances and cross-copy covariances are unchanged. This is redistribution of already present Gaussian noise, not reheating. Different copies may use different D_i, chosen deterministically from their hit history and frozen scalar clocks.

The Gaussian-bank part of the affine center row has squared norm at most the original Gamma_i,vv. The retained Y row is unchanged and must be added separately. For the literal cubic, d_v^2+Gamma_i,vv=(1+r_v^2)/2<=1, so its concatenated Y-plus-bank row, and each of its subrows, have norm at most one. For a general input use its actual recorded concatenated row bound. Clocks are frozen labels; no derivative of a chosen scalar square root or of a quadrature label is being hidden.

## 3. A hit-adapted spectral floor for the literal old cubic

Write the original genealogy, conditional on Y=y,

    Y0=q y+sqrt(1-q^2)G0,
    X=sY0+sqrt(1-s^2)G1,
    q1=r1 Y0+sigma1 W1,
    q2=r2 X +sigma2 W2,
    q3=r3 X +sigma3 W3,
    sigma_v^2=(1-r_v^2)/2.

These q_v are coefficient centers before each force's separate original private shield sigma_v is added. For a chosen force j define

    z_1=1-q^2, z_2=z_3=1-s^2,
    m_j=min({sigma_v^2:v!=j}, z_j).

Then the original center covariance Gamma obeys

    Gamma >= (m_j/16) I_3, for EACH j=1,2,3.               (3.1)

Proof. For j=1, discard the independent G1 part and regress all rows on Y0. The target row coefficient is r1 and the other two are r2 s,r3 s, all at most one. For j=2 or 3, regress Y0 onto X: the regression coefficient s(1-q^2)/(1-s^2 q^2) is in [0,1], and Var(X)>=1-s^2. Discard the independent regression residual. In either case, if r_j>=1/2, writing the common row as z=r_j x_j+sum_{v!=j} a_v x_v yields

    |x|^2 <=13 (sum_{v!=j} x_v^2+z^2).

The retained independent W_v and common variance then give Gamma>=m_j I/13. If r_j<1/2, sigma_j^2>=3/8 and the original independent W bank alone gives Gamma>=3m_j I/4. This proves (3.1).

Use D_i=(m_j/32)I in (2.3), with j chosen as the force with more than two inverse-width powers, if there is one. In the derivative orders needed here there cannot be two such forces. If all powers are at most two, any j is allowed; the estimate can simply ignore the extra heat.

## 4. Literal positive clock assumptions and local lemma

The original positive finite clock rule has endpoint panels of variance scale x. Its total panel mass is at most C x, each node mass is at most C x, and the width squared is comparable to x. The structural q/s rules have the corresponding endpoint mass for z=1-q^2 or 1-s^2. Fixed numbers/comparison factors are retained in C. These are the original five-clock mass assumptions already used by the sealed two-copy proof.

Let eta>0 be the minimum ORIGINAL private variance in the finite old program, decreased to at most 1/2. It is not newly added heat. Let L>=1 dominate one plus every old endpoint-panel count, and 1+log(1/eta). Other pre-existing cutoff errors remain the responsibility of the original target's ledger.

For one cubic copy with d bank/caller hits, let p_v be the exact inverse-width powers: p_center starts at one, each hit adds one at its actual force, and sum p_v=1+d. For 0<=d<=3 the complete weighted old-clock sum obeys

    J_d(delta) <= C L^5 (delta+eta)^(-a_d),
    a_d=max(0,(d-1)/2).                                  (4.1)

There is also a nodewise bound with the same power and no panel-count factor, including its literal original node weight.

Proof. If every p_v<=2, bound t_v>=sigma_v and sum x/(x^(p_v/2)); this is bounded or logarithmic. Otherwise choose j with p_j>2, put a=p_j/2-1, and use the hit-adapted m_j. Dyadic summation gives

    sum_{clock j} w_j/(sigma_j^2+delta m_j/32)^(p_j/2)
       <= C m_j^(-a) (delta+eta)^(-a).                    (4.2)

The inequality follows from t_j^2>=eta and t_j^2>=delta m_j/32, using m_j<=1, and the usual split at delta m_j. At these orders a<=1 and every other p_v/2+a<=1, since sum p_v<=4. Bound

    m_j^(-a)<=sum_{v!=j} sigma_v^(-2a)+z_j^(-a).

Every remaining old endpoint mass then has exponent at most its available one; equality contributes only a panel-count logarithm. The same pointwise argument uses x^(1-p_v/2-a)<=1 and proves the nodewise bound. This proof keeps weights at their actual powers. A repeated-label diagonal may be MAJORIZED by the corresponding nonnegative Cartesian product sum, but this is only an inequality, never a change to the target.

## 5. Three-copy main and complete one-hit budgets

Write the two edge gaps as s1=1-t1 and s2=1-t2. In every three-vertex tree, the central replica has delta_c=min(s1,s2), and its leaves have delta=s1,s2. Their tree degrees are 2,1,1. Equation (4.1) therefore bounds the main integrand by

    C L^15 (min(s1,s2)+eta)^(-1/2).

Its full square integral is bounded by 8/3 before the harmless eta regularization. A complete extra captured-caller/bank hit has two cases:

    central hit: (min(s1,s2)+eta)^(-1),
    leaf hit:    (min(s1,s2)+eta)^(-1/2)(s_leaf+eta)^(-1/2).

The central integral is <=2 log((1+eta)/eta). Splitting the leaf integral into s_leaf<=s_other and the complementary sector bounds it by log((1+eta)/eta)+2. Thus, after the finite tree/hit/history/readout constants,

    sum |main coefficients| <= C L^15,
    sum |coefficients| sum_v 1/t_v <= C L^15[1+log(1/eta)]. (5.1)

Both all proper cuts and HS/sqrt(D) have these bounds by the marked-spanning-tree theorem. For a bank derivative, retain its extra bank slot and use its actual row of norm at most one. The bound is uniform in fixed Y and coefficient centers. Physical Gaussian/Hermite energy and matrix-frame constants remain their fixed-rank factors.

## 6. Positive node rule and scalar root-square bound

Use any positive dyadic bridge rule with each endpoint panel mass <=C s and each node weight u_e<=C s_e, with its own certified quadrature error. Section 5 then holds discretely up to fixed constants, irrespective of the number of nodes inside a panel. The original endpoint s<eta can be one final panel; the exact integrand remains regular because the original private heat is positive.

Let b_h be the complete unsigned scalar node factor INCLUDING the two bridge weights u1 u2, all original clock weights, all exact width powers, and bounded injection scalars. Keep adapter/readout normalizations separately explicit. The local node bounds give

    sup_h b_h <= C,
    sup_h b_h sum_v 1/t_v <= C.                           (6.1)

Indeed u1u2<=C s1s2 pays min(s1,s2)^(-1/2) for the main, min(s1,s2)^(-1) for a central caller, and the stated product of square roots for a leaf caller. The regularization only improves these estimates. Consequently

    sum_h b_h^2 <= C L^15,
    sum_h b_h^2 (sum_v 1/t_v)^2
          <= C L^15[1+log(1/eta)].                       (6.2)

These are genuine squared literal weights, not a substituted independent-clock variance theorem. Known finite rank factors and actual inverse physical readout shares must be multiplied back in. They can depend on the chosen finite census and are not automatically public-log.

## 7. Original-VALUE construction and a precise return port

At each of the nine force vertices execute the active finite C_k original-gradient VALUE source at the exact center/width (2.3), including its same-center anchors, radial preparation if required, selected pair/filter, and frozen numerical floors. Here k<=3 for the main (a twice-hit original center is C3); a caller proof can use analytical k=4 but does not query it. The source is normalized by the literal original A; service amplitude alpha is separate.

Choose a genuine marked physical leaf as root and compile the marked spanning tree. For example put alpha^(1/16) at eight nonroots and b_h times all known scalar normalizations times alpha^(17/2) at the root. The product is exactly the required alpha^9 coefficient. Signs belong in the root. Every source covariance/readout share is positive. The complete captured-caller path uses the actual strict ancestors and (6.1); it is not a derivative of a LAW statement.

Conditional on coefficient callers, the imported width-zero native port has intrinsic error C_h rho_root,h^2 sqrt(D). Equation (6.2) pays the scalar part at alpha^17. The actual C_h, all readout shares, source/frame/gap constants and finite side errors must pass the literal native guard. With a zero bank-readout old endpoint having actual old first O(alpha), the proved same-endpoint Riesz join has old mixed scale alpha^10 times the main-energy constant, and self-bank scale alpha^18 times the one-hit logarithm. This is a fixed-family return under those exact original endpoint hypotheses. No original-target heat mismatch is created here.

If the cubature/node/readout constants bring additional inverse-alpha powers, they must be subtracted from these grades. This note does not declare arbitrary desired numerical tolerance achievable at public-log cost. A conservative executable cubature and its possibly expensive node bill are in FINITE-CUBATURE-AND-N-COPY-RULE.md.

The per-node force bill is sum_v 2^(k_v+1) M_native(k_v,b_v,epsilon_v), plus every original capture, old complete version, known scalar Gaussian roots, physical rows, keep, and affected ancestor/anchor replay. For the main nine-force graph sum k_v=7 and max k_v=3. One scalar (n=3) BKAR root acts on the COMPLETE old bank; each cubic residual uses a known 3-by-3 query root. Those roots and all retained internal tapes are charged. A cache is free only for the identical complete version.
