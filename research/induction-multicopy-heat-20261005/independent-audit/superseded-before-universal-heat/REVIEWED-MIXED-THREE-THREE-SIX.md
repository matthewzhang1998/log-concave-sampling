# Exact mixed three/three/six family, with literal diagonal and cross-history losses

2026-10-05. Positive source construction and heat/clock bounds for kappa_B(L3,L3,L6). This is a genuine old/new mixed bank family. The result keeps the same target and exposes its complete caller loss; it does not assert that the loss is unavoidable or that the whole induction closes.

## 1. Which L6 is meant

L3 is a complete finite old cubic coefficient sum on the full bank. L6 is the finite exact-symmetric-OU covariance counterpacket coefficient, including the old covariance bridge and its unchanged private heat. It has two old cubic copies and one covariance bridge hit on each. Its force count is six, physical rank six, width degree four, and amplitude grade six.

A diagonal six-packet has a single repeated old history/node, hence literal w^2. A FULL covariance of a sum sharing a bank also has cross-history products w_a w_b. These are different clock ledgers. The stronger diagonal estimate below must never be imported for the full sum without checking its exact labels.

The new three-copy cumulant has three labelled trees and these source types:

    a=12, N=12, K=10, R=12, M=11.

If L6 is the central replica, the two new bridges have 6*3 choices each, giving 324 ordered new hit histories. If one of the two L3 replicas is central, the edge to the other L3 has 3*3 choices and the edge to L6 has 3*6 choices, giving 162 histories each. The total is 648 new ordered hit histories for one exact old L6 history. Include its old nine covariance hit pairs, old physical permutations and cross-history census separately, without multiplying an already combined count twice.

The original spanning trees plus all three bank bridges form a twelve-force marked spanning tree. All physical marks survive. Maximum main adapter order is C4 when the twice-hit new-tree center is one of the original C2 centers inside L6.

## 2. Exact Gaussian realization, including the six-force residual

Apply the n-copy rule of THREE-CUBIC-EXACT-HEAT.md to the COMPLETE original joint coefficient bank, including the residual roots already owned by L6. The two cubic replicas use their hit-adapted heat. The diagonal six-packet can retain its actual old widths without extracting new heat. For the full cross-history packet, also extract heat from its recorded residual roots, as follows.

Write s=1-c_old for the distinct OLD covariance-bridge gap. The actual old six-packet centers use common old-Q rows plus two independent residual query blocks with covariances

    s (Gamma_a-m_a^old I/32), s (Gamma_b-m_b^old I/32).

Here m^old is whatever spectral floor was actually used by that recorded packet (the sealed canonical choice omits its center). The original private widths are sigma_v^2+s m^old/32. Since Gamma>=m^old I/16,

    Gamma-m^old I/32 >= Gamma/2 >= m_j I/32

for EVERY hit-adapted m_j from the cubic spectral lemma. In the new replica's independent delta_6 part, each of these two residual blocks is therefore able to provide extra independent query heat

    delta_6 s m_j/64

at all three forces of that constituent cubic. Subtract exactly this diagonal covariance from that constituent's new residual covariance. Both remainders are positive. This is redistribution only inside the NEW independent replica innovation; the recorded common B row and every original old six-packet query covariance remain unchanged.

For a given finite hit history choose j in each constituent cubic according to its final inverse-width powers. Its new private widths satisfy

    new_t_v^2 >= sigma_v^2+delta_6 s m_j/64.               (2.1)

The old private heat is preserved as a summand. This does not change L6 as a function of its original bank before the BKAR expectation: it exactly represents L6(B_6) after disintegrating the new B_6 innovation. Changing an OLD redistribution while claiming the same conditional L6(B) would be invalid; that is not the operation here.

All original old covariance derivative contractions, and all new BKAR contractions, keep their ACTUAL original injection rows. Complete center rows have norm at most their original norm, at most one. No quadrature clock or scalar square root is differentiated as a dynamic input.

## 3. A six-packet derivative ledger

Let J6(d,delta) sum the actual old clocks and old bridge weight for d additional bank/caller hits on a six-packet inside a new replica with independent variance delta. Include every ordered hit assignment. This is an analytical bound, not an executed high-derivative query.

Let eta>0 be a lower bound on all ORIGINAL private and structural endpoint variances entering the old floors, and let L include their old panel counts and 1+log(1/eta). For the FULL common-bank cross-history sum, (2.1) gives

    J6(0,delta) <= C L^C,
    J6(1,delta) <= C L^C (delta+eta)^(-1/2),
    J6(2,delta) <= C L^C (delta+eta)^(-1),
    J6(3,delta) <= C L^C eta^(-1)(delta+eta)^(-1).        (3.1)

The constants are fixed enumerable rank constants. They hide no inverse-alpha powers.

Proof. Distribute the d extra hits as d_a+d_b=d. Each constituent cubic already has one OLD covariance hit. For d_a<=2 the cubic dyadic lemma with effective gap delta*s gives local sum

    C L^5 (delta*s+eta)^(-d_a/2).

For d_a=3, its total inverse-width degree is five. Choose the force with largest exponent p if p>2. Its heat costs m_j^(-(p/2-1)); at any remaining clock this exponent plus that clock's own force exponent/2 is at most 3/2, since the total is five. An original single mass pays one. At most eta^(-1/2) is left. If all force exponents are <=2, no deficit occurs. Thus the safe local bound is

    C L^5 eta^(-max(0,d_a-2)/2)(delta*s+eta)^(-d_a/2).

Multiply the two literal clock sums and integrate the old bridge s. The elementary bounds, valid also as delta tends to zero, are

    integral_0^1 (delta*s+eta)^(-1/2) ds <=2/(delta+eta)^(1/2),
    integral_0^1 (delta*s+eta)^(-1) ds <=2L/(delta+eta),
    integral_0^1 (delta*s+eta)^(-3/2) ds <=2 eta^(-1/2)/(delta+eta).

These prove (3.1). Repeated labels may be bounded by a nonnegative Cartesian-product majorant, but their target is never changed.

For a DIAGONAL packet, retaining its actual repeated weights w^2 gives a stronger estimate independent of the new delta:

    J6,diag(0), J6,diag(1), J6,diag(2) <= C L^C,
    J6,diag(3) <= C L^C eta^(-1/2).                      (3.2)

For example, with the old center-center hit, let x=sigma_c^2 and let j additional hits land at the two centers. Their total inverse-width degree is 4+j. The literal center weight is w_c^2, so

    sum_c w_c^2/(x+s m/32)^((4+j)/2)
       <= C L (s+eta)^(-j/2) m^(-j/2).

The remaining d-j leaf hits are paid by their SQUARED masses. At d<=3 their worst combined endpoint exponent is at most d/2<=3/2, below the available two. Structural clocks also retain squared masses. The old bridge integral then gives at most eta^(-1/2). If both old covariance hits land at the same repeated leaf and all three new hits land there, its combined exponent is FIVE. Its squared mass leaves eta^(-1/2), still within (3.2), while the center has only its baseline exponent two and creates no simultaneous deficit. Otherwise a leaf exponent is at most four and is paid logarithmically by its squared mass. Center exponents decrease when an old hit moves to a leaf. A uniform power proof avoids relying on enumeration: let p_c,p_1,p_3 be the combined powers at the three repeated clock labels. Then p_c+p_1+p_3=4+d<=7. If p_c>4, center summation costs (s+eta)^(-a)m^(-a), a=(p_c-4)/2<=3/2. At any leaf, p_v/2+a<=(7-4)/2=3/2<2, so its squared mass pays every m term. Only the old bridge can then leave eta^(-1/2). If p_c<=4, no m cost is needed; at most one leaf can exceed exponent four, and its maximum is five. This again leaves at most eta^(-1/2). For d<=2, the same inequalities replace seven by six and five by four, yielding only logarithms. This covers all nine old hit pairs and every extra-hit distribution.

The stronger diagonal estimate is not imported for the full common-bank covariance. Neither estimate is a lower bound or a no-go theorem.

## 4. The new three-copy clock budget

Write s1,s2 for the NEW edge gaps. If L6 is central, delta_6=min(s1,s2), while each cubic leaf has delta equal to its incident gap. The main is bounded by (min(s1,s2)+eta)^(-1) from J6(2); its integral is logarithmic. A caller on L6 costs eta^(-1) times the same logarithm. A caller on a cubic leaf adds (s_leaf+eta)^(-1/2), bounded conservatively by eta^(-1/2), so is smaller than the stated eta^(-1) allowance.

If a cubic is central, its main factor is (min(s1,s2)+eta)^(-1/2), and J6(1) adds (s_6+eta)^(-1/2). Their joint integral is the logarithmic two-square-root bound from the three-cubic proof. A caller on the central cubic replaces its exponent 1/2 by 1; a caller on the six-packet replaces the exponent on s_6 from 1/2 to 1; a caller on the other cubic adds one more exponent 1/2. Each is bounded by eta^(-1/2) times an already logarithmic main envelope. Thus all are within eta^(-1).

For the diagonal packet, use (3.2), and the new-edge analysis reverts to the simpler three-cubic integrable/logarithmic cases; the only power loss is eta^(-1/2) when the six-packet receives three hits.

Consequently the complete leading coefficient and extended one-hit bounds are

    all proper cuts:        C L^C alpha^12,
    HS:                     C L^C alpha^12 sqrt(D),
    complete caller cuts:   C L^C alpha^12 eta^(-1),         full sum,
                             C L^C alpha^12 eta^(-1/2),     diagonal.
                                                               (4.1)

The marked-spanning-tree theorem supplies all cuts and the single Hilbert factor. Bounds are uniform in fixed retained callers. Physical Hermite and matrix-frame factors at rank twelve remain explicit fixed-rank constants.

## 5. Literal node squares and a positive native return port

Include the old bridge node weight and both new bridge weights in b_h, together with every actual old weight and inverse width. For a main node, the most singular six-packet local degree is d=2; its old factor is at most C u_old/(delta_6*s+eta), hence at most C/(delta_6+eta). The NEW two bridge weights pay this factor when L6 is central. When a cubic is central, its delta^(-1/2) and the six-leaf delta_6^(-1/2) are paid by the two new bridge weights. The diagonal estimate is stronger. Therefore

    sup b_h <=C,           sum b_h <=C L^C,
    sum b_h^2 <=C L^C.                                    (5.1)

For one complete caller,

    sup b_h sum_v 1/t_v <=C eta^(-chi),
    sum b_h sum_v 1/t_v <=C L^C eta^(-chi),
    sum b_h^2(sum_v 1/t_v)^2 <=C L^C eta^(-2chi),          (5.2)

where chi=1 for the full cross-history sum and chi=1/2 for the certified diagonal family. These are safe conservative bounds. All fixed adapter constants and actual inverse physical readout products are still charged separately.

Compile each twelve-force marked tree with eleven nonroots of amplitude alpha^beta and root amplitude b_h times all exact normalization factors times alpha^(12-11 beta). The product is exactly its original alpha^12 coefficient. The old/new fields and all residual roots have zero source-zero bank readout, and a separate Gaussian keep is reserved.

The frozen-callers width-zero port gives intrinsic own error with scalar allowance

    alpha^(24-22 beta) sum_h C_h b_h^2 sqrt(D).

The complete actual caller path includes the worst root contribution

    C L^C alpha^(12-11 beta) eta^(-chi).

For eta=alpha^(2 gamma) and fixed beta>0, a sufficient source-first smallness exponent is

    12-11 beta-2 gamma chi >1.                            (5.3)

At an old endpoint with actual complete old bank first O(alpha), the same-endpoint bank mean join has old mixed scale alpha^13 times the main energy constant and self scale alpha^24 eta^(-chi), each with one sqrt(D). Thus a sufficient scalar test for grade thirteen is

    2 gamma chi <=11,
    24-22 beta >=13,

together with (5.3), all actual readout/node constants and propagated finite floors. For example beta=1/22 leaves root grade 23/2 and requires gamma<21/4 in the full chi=1 case, or gamma<21/2 for diagonal chi=1/2. These are sufficient finite-window conditions, not a theorem about an arbitrarily small inherited cutoff. No new smoothing mismatch is added.

The old coefficient's ORIGINAL truncation/error ledger remains owed if the ultimate target is a continuum clock object. Choosing an inherited eta after seeing these inequalities cannot erase that original debt. A full feedback queue must carry it explicitly.

The source bill is twelve complete original-gradient C_k programs per new history/node, including all old L6 captures, each exact original bridge/root version, old programs at required precision, and all affected replay. The raw main adapter order is at most four; analytical first may reach five. New quadrature cost and readout shares are not replaced by an assumed constant.
