# Root-seed bridge for the original common-carrier execution

5 October 2026. Draft theorem with a complete proof. The executed high-order native program is unchanged. The proof inserts actual order-two root-seed maps, not globally truncated polynomials. All comparisons stop at the original single group keep. This is an environment-stable group theorem, not a carrier-fixed high-order packet prior.

## 1. Exact hypotheses and quantitative conclusion

Freeze the original coefficient/caller/clock/source-version record E. For each packet h choose the marked-leaf root used in LOW30. Such a root has one child and zero external dynamic probes. Its source is therefore a fixed source f_h conditional on E; its selected mean matrix M_h is independent of all new dynamic tapes. Write its radius as r_h, including its actual selected-source normalization. Absorb fixed root readout coefficients and complete source/path constants into a nonnegative envelope s_h >= C_h r_h. All assertions below are conditional on E with uniform deterministic envelopes, or integrated with their explicit Lp versions. Geometry is independent of changing E coordinates.

The original root output is P_{b_h}(f_h; I_h), where I_h is the complete actual child-subtree output. The packet has the literal form

    X_h^a = c_h P_{b_h}(f_h; I_h) + B_h P_h^pub + G_h^buf
          = L_h P_h + u_h^a,       L_h L_h^T = v I_D.

The root compiler's known standard source-zero carrier G_h^root is independent of the incoming I_h. Its primitive carrier bank, the original marked public readouts, and any internal root readout buffer comprise P_h. All nonroot compiler innovations have zero readout when r_h=0 and belong to the private bank. The source-zero identity here is literal. The common execution uses

    H ~ N(0,v I_D), Pi_h = I-L_h^T L_h/v,
    P_h = L_h^T H/v + Pi_h V_h.

The original raw path theorem supplies a global complete Gaussian first C_h r_h for u_h^a, hence a carrier first C_h r_h/sqrt(v). Its root-linear L2 energy is derived, not attributed verbatim to that theorem: output-rank Gaussian Poincare gives ||u_h^a-Eu_h^a||_2<=C_h r_h sqrt(D). The root conditional prior gives |Eu_h^a|<=|c_h|[||M_h||op |E I_h|+a_h], and |E I_h|<=d_h since the ideal child is marginally standard. At the admitted d_h<=C sqrt(D), r_h<=1 and b_h>=2, the displayed bound for a_h below makes this mean O(C_h r_h sqrt(D)). Enlarge s_h to bound both complete first and full L2 energy, giving Lip_H u_h^a<=s_h/sqrt(v) and ||u_h^a||_2<=s_h sqrt(D). The root-seed maps constructed below have the same bounds after enlarging s_h by fixed/public-log constants. Set S=sum_h s_h, Q2=sum_h s_h^2.

Let d_h be an integrated conditional W2 error for replacing the child subtree by its ideal forward tree while retaining its original marked publics identically. The same d_h bounds the unconditioned distance of I_h to standard Gaussian. These errors are supplied by the native tree proof with its strict-ancestor attenuation, now applied to the child tree; the public-retention claim is justified in Section 4. Let a_h be the root's original conditional prior floor evaluated at I_h, and let beta_h be the root finite tail-filter mean floor there. The exact envelope gives

    a_h <= C_h r_h^{b_h} sqrt(D) + C_h r_h d_h,
    beta_h <= epsilon_h^tail sqrt(D) + C_h r_h d_h.

Define

    E_seed = sum_h |c_h| [a_h + beta_h + C_h r_h^2 sqrt(D)].

The first bracket is a genuine ordinary pre-keep packet comparison, proved below. A sharper literal root-seed bound can replace it. Let F_poly be the sum of the strong ideal-forward-to-original-polynomial floors at the chosen finite versions. Then

    W2(group_actual, group_original_polynomial)
      <= (1+S/sqrt(v)) E_seed
         + (2/sqrt(v)) S^2 sqrt(D)
         + (1/(2 sqrt(kappa))) Q2 sqrt(D)
         + sum_h |c_h| [beta_h + ||M_h||op d_h
                         + C ||M_h||op^2 sqrt(D)]
         + F_poly.                                      (1)

Here the source amplitude is inside M_h, so ||M_h||op <= C_h r_h. The Q2 coefficient uses s_h enlarged to bound the centered seed's complete private first. The covariance restoration term includes the actual fixed c_h and root convention. Formula (1) has no polynomial-in-D multiplier and no inverse packet variance shares. For fixed numerical v,kappa it is O(E_native + E_filter + S^2 sqrt(D) + sum r_h d_h + F_poly), with all actual public-log constants exposed in s_h and the floors.

Equation (1) is first a conditional-on-E bound. If its envelopes depend on E, denote the full displayed right-hand side by B(E) and require ||B(E)||_{L2(E)}<infinity; the integrated retained-E coupling has cost at most that L2 norm. Products such as S(E) E_seed(E) are not replaced by products of L2 norms without the needed higher integrability. The output is translated by the unchanged W_old(E), and both laws use the original one independent keep sqrt(kappa) G_keep. W_old reads no new dynamic tape. Formula (1) compares these whole group outputs; it does not remove or duplicate their keep.

## 2. A literal root seed with the same zero carrier

Use LOW30's order-two pair seed (equations P-seed and P-target, lines 5837–5909). With the selected inclusion S0 and contraction C0, let

    D_x(u)=f_h(C0 u+S0 x)-f_h(u),
    m_tail(x)=E D_x-M_h x.

The finite scalar filter supplies an actual source J_h(x,w) whose private and incoming firsts are O(r_h) and whose mean approximates m_tail at its separately paid floor. All root source-private Gaussian inputs (u,w, and shields used by J_h) are fresh and independent of the old child subtree and of P_h. Execute the proof-only map

    Y_h^s = G_h^root + D_{I_h}(u_h)-J_h(I_h,w_h),
    X_h^s = H + u_h^s,
    u_h^s = c_h[D_{I_h}(u_h)-J_h(I_h,w_h)].              (2)

This is exactly the order-two seed's marginal law after collapsing its independent source-zero Gaussian sum into the existing G_h^root. Its untouched eta Gaussian and the other Gaussian summands belong to that sum; no source-private VALUE variable belongs to it. Replacing the sum by G_h^root preserves the seed's joint law with the child input and original marked publics. It is not subtraction from an abstract LAW. The old high-order root carrier can be used because it is a known Gaussian affine row independent of the incoming passive.

The seed residual is globally Lipschitz in its private roots and child input. Its conditional centered energy is at most its private first times sqrt(D), by Gaussian Poincare and rank(Du)<=D. Its conditional mean is c_h M_h I_h plus the tail-filter discrepancy. Since ||I_h||_2<=sqrt(D)+d_h, the usual admitted choices d_h<=sqrt(D) and beta_h<=C_h r_h sqrt(D) give ||u_h^s||_2<=C_h |c_h| r_h sqrt(D). Otherwise enlarge the literal s_h by these displayed extra factors before using (1); do not assert the root-linear energy for free. Propagating through the unchanged nonroot raw programs yields the same root-linear complete packet radius as the original raw-path proof. In particular, the finite root-seed family, the original family, and every mixed family have deterministic global carrier first bounded by S/sqrt(v). No ideal polynomial appears in this telescope.

Condition on I_h and all original public/caller labels, but not on the root carrier. Original P_b is within a_h of K_{M_h}(I_h). The order-two seed is within beta_h+C_h r_h^2 sqrt(D) of that same conditional kernel. The seed mean-noise estimate uses only its uniform private first; the Gaussian covariance correction costs C ||M_h||op^2 sqrt(D). Couple these conditional outputs and keep the marked publics fixed. This proves the pre-keep ordinary packet error used in E_seed. It does not keep G_h^root or H fixed.

At a mixed replacement stage, lift that packet-output coupling by disintegrating each carrier conditional on its packet output. Both sides retain the same original marginal Gaussian carrier. Their displacement is bounded by

    ||H_a-H_s||_2 <= e_h+2 s_h sqrt(D).

Draw the other packet private banks independently and use them identically on both sides. Their mixed remainder is globally H-Lipschitz with constant S_-h/sqrt(v). Hence the stage cost is

    (1+S_-h/sqrt(v)) e_h
       + 2 s_h S_-h sqrt(D)/sqrt(v).

Summing proves the first two terms in (1). This is the environment-stable native replacement sought by the common-carrier interface. The carrier is allowed to move in this proof coupling, and its effect on every readable unchanged packet is explicitly charged.

## 3. Remove only the seed's private centered noise

Let C contain E,H, every original marked public, all nonroot tapes, and all packet null coordinates that the source-zero bank uses. The fresh root-seed source banks Omega_h are independent conditional on C and are not read outside their own seed. Define

    m_h(I_h)=E_{Omega_h} u_h^s,
    F_h=u_h^s-m_h(I_h),  F=sum_h F_h.

The conditional seed means are analytic proof objects and are never evaluated by the executed program. Gaussian Poincare and the complete private first bounds give

    ||F||_2 <= sqrt(Q2 D),
    ||D_Omega F||op <= sqrt(Q2)

pointwise for the second bound. The first uses conditional independence; an l1 S bound would also suffice. The derivative concatenates the independent Omega_h blocks. Centering m_h differentiates no Omega coordinate.

Use the positive path

    Z_t=W_old+H+sum_h m_h(I_h)+tF+sqrt(kappa)G_keep.

Gaussian Riesz in the whole Omega bank, conditional on C, gives

    d/dt E phi(Z_t)
       = t E[(R_Omega F)(D_Omega F)^T : D^2 phi(Z_t)].

The HS energy is at most Q2 sqrt(D); this is one Hilbert factor times one operator factor. The coefficient and F are independent of G_keep. One integration by parts in that same final keep and the continuity-equation estimate give the cost Q2 sqrt(D)/(2 sqrt(kappa)). This is a grouped current comparison after the one final keep, not a conditional-W2 assertion for F before smoothing.

Next pay the actual tail mean discrepancy strongly:

    m_h(I_h)=c_h M_h I_h+b_h(I_h),   ||b_h(I_h)||_2<=|c_h| beta_h.

No private source bank is left readable by the endpoint after its conditional mean has replaced it. The child input and all marked publics remain readable; they have not been averaged or erased.

## 4. Child-tree replacement with the carrier fixed

A marked-leaf root has no external dynamic probes. Thus M_h depends on E alone, not on I_h, marked publics, the root carrier or any nonroot innovation. This property is essential.

The full kernel construction is given in RETAINED-PUBLIC-CHILD-LEMMA.md. The child tree's original marked publics are primitive leaf passives and marked nonleaf public inputs. They are already read by the full tree's affine marked readout. LOW30's tree-replacement coupling keeps these primitives fixed: its spine composition retains the original passive (lines 5912–5917), conditions on the spine's own publics and side inputs, and consumes only internal spine outputs before collapse. Recursive side-subtree replacement uses independent distinguished passives and subtree tapes. It can therefore be constructed as a kernel conditional on the original marked-public tuple, with integrated cost d_h. No uniform-in-public claim is needed.

For completeness, retaining the other marked publics does not make the independent distinguished Gaussian used in the side-input HS estimate disappear. At a side replacement, its coupling is chosen from that side subtree's own private tapes and publics; it is independent of the enclosing distinguished passive and spine innovations. Those are integrated in the L2 estimate after the conditional coupling is built. The physical public marginals are still the original standard Gaussian bank. Internal side outputs are not retained as observers. This is precisely the spine-before-side proof order, with the existing affine readout labels kept.

The root compiler's primitive source-zero bank is independent of the entire child tree conditional on its marked publics. Under the common-carrier construction, conditioning on H changes the conditional distribution of these publics and generally correlates them with the root carrier; their original Gaussian marginal is unchanged. The child's private kernel given its publics is unchanged. Couple the child outputs using the same publics and the same H; integrate over their unchanged Gaussian marginal. The cost remains d_h. Its propagation through c_h M_h is <=|c_h| ||M_h||op d_h. Across packets the child couplings factor conditional on E,H and the packet-specific null public coordinates, so the usual conditional product construction is valid.

After this replacement the grouped mean endpoint is

    W_old+H+sum_h c_h M_h I_h^ideal+sqrt(kappa)G_keep.

Restore, on identical tapes, the ideal root covariance increment

    c_h[(I-M_h M_h^T)^(1/2)-I]G_h^root.

Each root G_h^root is marginally standard and independent of E, the child tree and its packet's marked publics in the preserved packet marginal. It is correlated with H and is generally correlated with those publics conditional on H, in exactly the intended way. Since M_h is E-measurable and uniformly small, the strong restoration costs C |c_h| ||M_h||op^2 sqrt(D). No independence between different packets is asserted or needed for Minkowski.

The identification is pointwise:

    H+c_h M_h I_h^ideal+c_h[(I-M_h M_h^T)^(1/2)-I]G_h^root
      = c_h[M_h I_h^ideal+(I-M_h M_h^T)^(1/2)G_h^root]
        +B_h P_h^pub+G_h^buf.

Choose the original ideal root innovation to be this same known G_h^root; if necessary extend its normalized affine row to an orthogonal coordinate system on the original root carrier tape, leaving unused coordinates private. No fresh independent replacement innovation is inserted at this step. The result is exactly the original ideal forward forest assembled with the same physical source-zero H. LOW30's strong ideal-forest-to-polynomial comparison can now be made on identical tapes packet by packet. Each original packet marginal Gaussian bank is unchanged by the coisometry, so the sum of its strong floors is F_poly. This finishes (1).

## 5. Readable and erased tape ledger

- E, original coefficient banks and old caller labels: retained throughout the conditional proof; actual external admission still follows the original observer ledger.
- H and packet P_h: internally shared by the actual group; may move during the actual-to-seed coupling, with every mixed environment charged. Kept identical during child replacement and final strong restoration. Never exported separately.
- Marked publics: kept in all subtree comparisons, still read by affine packet readouts. Their marginal Gaussian law remains exact; they are not mistaken for independent conditionally on H.
- Original root-native private tapes: consumed only in the ordinary root conditional law comparison. Other packets read them only through H, whose movement has already been charged. The seed comparator uses fresh source-private banks.
- Seed source banks Omega_h: independent of H and child trees, centered conditionally on all those variables, then consumed by the exact Riesz identity. No surviving caller reads them.
- Nonroot compiler innovations: absent from the source-zero affine readout, hence private after conditioning on marked publics and H; consumed in the retained-public subtree prior. Internal spine outputs are never exported.
- The one group keep: independent, unchanged and never read by any source, caller or Riesz coefficient. Used once for the grouped centered-noise current bound. No per-packet keep is appended.
- W_old: unchanged and reads only E/its old tapes, never a new dynamic bank.

## 6. Triple-cubic consequence and costs

For the original rank-nine groups use r_h proportional to b_h alpha^(17/2), and the eight nonroots alpha^(1/16), exactly as before. Then S<=B alpha^(17/2), so every new seed/hybrid/current/covariance substantive term is O(B^2 alpha^17 sqrt(D)). It is beyond the alpha^10 target without a polynomial-D smallness guard.

The child-prior error has its literal native strict-ancestor formula. A safe unattenuated fixed-order choice d_h<=Gamma_h alpha^(b/16) sqrt(D) costs root-times-child alpha^(17/2+b/16); b>=24 reaches alpha^10, and any fixed b>24 leaves a strict margin. The actual per-node radii, complete caller paths, public logarithms and finite census must still be substituted. This argument does not hide a growing-order theorem. Root prior floors and root tail filters have their separately propagated costs; choose them, the ideal strong polynomial floors, scalar cubature and old-bank join budgets as in the original finite theorem.

No additional source evaluations are executed by this proof: the root seed, its conditional means and all couplings are proof-only. This statement is relative to the original compiled versions that meet the displayed child/root floors. If a caller chooses to raise a native order or tighten an actual source version to meet them, its complete compilation and affected replay must be charged in the usual bill. The actual sampler, Gaussian rows, native repetitions, original VALUE bill, source versions, precision floors, private curl and complete terminal-first ledger are unchanged. No derivative bound is inferred from any W2 floor. The seed comparator's ordinary global first is proved directly from its source map.

The independent-share alternative remains valid but is unnecessary for this native-hybrid port if the marked-leaf root and zero-readout privacy hypotheses above hold. This result does not improve quadrature node counts, inverse-width caller costs, original target debts, or arbitrary-order growth.

## 7. Source anchors and scope checks

LOW30 pinned source: /workspace/scratch/e90bd698ae04/research-source/High Acc Ideas/ai-bucket/30_low_acc.tex, SHA256 7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8.

- Native conditional pair contract and its actual erased tapes: 5837–5859.
- Literal root order-two seed and uniform centered-noise estimate: 5860–5896.
- Composition retaining the original passive: 5898–5926.
- Marked-leaf root has only one child, zero external probes: 6008–6015.
- Spine-before-side tree prior and its independent distinguished Gaussian: 6029–6063.
- Strong ideal finite polynomial comparison: 6210–6307.
- Global raw source-zero carrier, its incoming independence and complete paths: 6488–6517, 6544–6618.

The conclusion uses the globally Lipschitz actual source versions supplied by the admitted active radial-pullback VALUE sources. A root with dynamic external probes or geometry that depends on live coefficient callers needs an additional proof and is not included. The theorem also requires that every nonroot innovation has zero source-zero affine readout; if a proposed implementation exports one, it must be moved to the public ledger and the subtree argument rechecked.
