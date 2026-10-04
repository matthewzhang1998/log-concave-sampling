# Extension audit: weighted graph and relative known-center source admission

2026-10-04. This extends, and does not alter, the independently pinned short-refresh audit.

## Decision

**Relative admission passes.** The weighted deterministic graph has the native small-edge, bounded-affine-row, weak-terminal-row, and protected-row geometry independently of its inverse-heat node count. A polynomial-logarithmic node-count hypothesis is not needed for those facts. The known-center law/source induction has original-query exponent

    k_R <= max(k_P,R-3/2), R=P+1/2,

with the stated actual first, common retained origin, retained state grade, and fixed J=2.9 protected width.

This decision is relative to the completed known-center seed and its complete original-gradient/finite-restoration ports. It does not independently reprove the old seed or promote a marginal-law-only sampler into a source. The revised input port 7 now explicitly includes the genuine-gradient reference and finite restoration obligations spelled out in Section 5 below. A finite non-gradient block cannot be inserted into an exact-gradient embedding merely because its VALUE is close to a gradient.

The child-heat normalization in Section 8 passes for fixed or public-logarithmic comparability. The hidden adapter bill in Section 9 is a separate bill, not the known-center recurrence. Section 10 genuinely uses only known-center posterior calls and correctly includes its outer 1/a factor.

## Inspected version and preservation

The reviewed route is

    ../DETERMINISTIC-KNOWN-CENTER-SOURCE-WITH-PROTECTED-REFRESH.md

with inspected SHA-256

    44cfd4b2d5ef98b86df8b824c1b98217281547eb7cf19c3d3ff05f095d816036.

The original short-refresh audit remains unchanged, SHA-256

    43914ab982e96c73e9538851c4ac4c6bb9ed9a980b10176c3a883404223c5917.

The final revision was re-read after input port 7, full-joint zero-translation protection, and the private/caller dimension distinction were added. These additions satisfy the audit's corresponding signature qualifications. Its new quantitative Section 5 also explicitly proves the m_nodes-independent centered Picard bound, weighted original VALUE-error propagation, finite precision count, and actual complete graph cost. A byte-identical reviewed-source snapshot is retained beside this report.

The reviewed LOW30 source remains the supplied file with SHA-256

    7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8.

Relevant native text includes lines 420–456 (actual finite columns), 741–904 (graph, projection, twin, caller), 7240–7445 (full-center/protected and original-source normalization), and 12679–12734 (decoder law and actual derivatives). The completed retained-source and same-version certificates are imported as seed contracts, not newly reconstructed here.

## 1. The diagonal change is an exact program identity

Let a physical normalized node be f_i(q_i), where f_i is a genuine gradient with ||Df_i||<=1. Give a group of N nodes d_i=N^(-1/2). Define

    fhat_i(z)=d_i f_i(z/d_i),
    p_i=d_i f_i(q_i),
    C_i=d_i Cphys_i,
    A_ij=d_i Mphys_ij/d_j.

Then p_i=fhat_i(C_i W+sum_j A_ij p_j) is exactly the same finite physical program after division by d_i. It introduces no approximation and no source deletion.

If f_i=grad U_i, the new potential is d_i² U_i(z/d_i). Its Hessian is Df_i(z/d_i), so it preserves the unit derivative bound, symmetry, and sign type. A negative primitive remains a signed gradient. An original VALUE evaluation still occurs once at its original physical argument; a directional first/adjoint still requires the corresponding original HVP. Neither d_i nor 1/d_i permits omitting that call.

The actual original physical widths must be retained in the oracle and finite-precision census. The scaled map has bounded first radius even if its displayed internal input width is large. This is sufficient under the global C2 bounded-Hessian original oracle model. For a restricted-domain oracle, all those original query arguments would need a separate domain certificate; normalization alone does not supply one.

## 2. Product-integration rows and columns

Let Delta=T/N, t_i=i Delta, and

    w_ij=cos((i-j-1)Delta)-cos((i-j)Delta), j<i.

For 0<T<=pi/2 these are nonnegative. The row sum is 1-cos(t_i). For a fixed j, the column sum also telescopes:

    sum_(i=j+1)^N w_ij = 1-cos((N-j)Delta).

Thus both row and column sums are at most 1 for the quarter step and at most K_h=1-cos(h) for the short step. Deleting endpoint-only rows only decreases these bounds. Schur's test gives

    ||W_quarter||_2<=1, ||W_short||_2<=K_h.

For equal-size Picard groups the d_i factors cancel, so an internal force block has norm at most a or a K_h respectively. The same inequalities hold for the matrices of absolute block norms. There is no unproved conversion from a row-sum bound alone to a Euclidean norm bound.

The individual weights satisfy

    max_j w_quarter,Nj <= (pi/2)/N,
    max_j w_short,Nj <= h²/N.

These bounds prove the unequal-group and terminal estimates below.

## 3. Cross groups, old-state replication, and terminal geometry

The short harmonic baseline contains cos(u_i) times the quarter endpoint. After eliminating the known affine endpoint skeleton, its quarter-force edge is

    -a cos(u_i) w_quarter,Nj d_i/d_j.

For N_i short nodes and N_j quarter nodes this is bounded by

    C a /sqrt(N_i N_j).

Its Frobenius norm, hence its operator norm, is O(a). This is the correct physical force connection; there is no unit nonlinear edge from a completed quarter-state label once its affine skeleton is eliminated.

If the old retained state's nonlinear readout is R_old p_old with ||R_old||=O(a), replication into a new group is

    (d_i cos(t_i))_i R_old.

The first vector has norm at most one. Therefore this cross block is O(a), without dividing the old block by the new d_i or copying its old transcript N times. All outgoing directions accumulate before an old adjoint sweep on its saved tape.

A terminal row using a physical quadrature vector w is a w_j/d_j. With d_j=N^(-1/2),

    ||a w/d||² = a² N sum_j w_j²
               <=a² N (max_j w_j) sum_j w_j.

This is O(a²) for the quarter readout and O(a² h^4) for the short readout. Their fixed number of groups gives terminal norm O(a). The weak terminal path estimate required by the twin remains O(a), and its squared twin price remains O(a² epsilon).

At each fixed target the number of groups is finite. The old admitted block remains in its admitted coordinates. Combining its strict contraction guard with the O(a) new and cross blocks gives a common strict guard for the complete directed graph and its symmetric twin, after the fixed-target small-a choice. If one insists on a specific numerical guard such as 1/8, the old guard must have that much slack; a generic strict native contraction threshold suffices for this extension.

## 4. Stacked rows and centered profiles

For every group, sum_i d_i²=1. A group of bounded pointwise original Gaussian/caller rows therefore has bounded stacked operator and Frobenius norms, with constants depending on the finite number of Gaussian directions and groups rather than the number of deterministic nodes. In the new quarter/short groups the private direct harmonic rows are coisometric before scaling; the full physical-center port has bounded direct coefficient.

The weighted fresh-L column in one short group satisfies

    sum_i d_i² sin²(u_i)<=h².

It is exactly zero in all preceding groups and original caller labels. This supplies the global O(h) restricted affine-row norm missing from the naive unweighted graph. The direct positive-kernel projection proof in the original audit remains valid, independently of this Hilbert-space argument.

The Gaussian tape dimension is unchanged by a deterministic quadrature node. Each induction step adds the physical-dimensional Z and L, and a proxy adds its separately counted twin bank only when constructed. Internal force-vector dimension, Gaussian input dimension, and original-query count are different quantities.

For centered primitive nodes fhat_i(0)=0, a contraction graph satisfies

    ||p(W)|| <= ||CW||/(1-||A||).

The new weighted affine rows have bounded Hilbert--Schmidt norm on the actual Gaussian input directions: summing their covariance traces gives O(Dim) times fixed/logarithmic source constants, not O(N Dim). The analogous actual-row profiles of inherited blocks must come from their completed source contract. This controls the initial Picard error and fixed moments without inserting a spurious sqrt(N) into a radial energy estimate. Nonzero deterministic anchors are retained or removed only by the lawful full-joint centering/gauge mechanism.

For numerical restoration, N can be a fixed power of 1/a. Then log N=O_R(log(1/a)); inverse d_i are also fixed powers at a fixed R. Enumerating their actual VALUE amplification and choosing finite original precisions afterward changes the logarithmic precision/iteration budget, not automatically the original-query heat exponent. This conclusion requires the existing arbitrary-precision finite restoration port. It is not a consequence of a mere bounded first radius.

## 5. Precise source signature and additional qualifications

The following are the required interface facts. Most are already stated in the reviewed route or in its imported completed native contracts; stating them together prevents a hidden strengthening.

### Required input

- A literal known-center finite program at fixed a,y with the claimed posterior error, fresh first, actual caller first, and actual conditional moments.
- A same-tape retained STATE with its own energy grade and common finite origin, not a state estimate obtained by dividing a force error.
- An old retained original-gradient graph, bounded original Gaussian and full-center rows, an admitted absolute operator guard, and an O(a) nonlinear state readout.
- Every retained whole block used as a primitive in the exact signed embedding has a genuine full-joint gradient reference. If its executable version is only approximately gradient, its finite VALUE and zero restorations remain separate and are charged at their newly enumerated consuming widths.
- Access to sufficiently accurate finite realizations of those reference blocks, with their original-query exponent unchanged by tighter finite precision. One immutable low-precision old version cannot simply be asserted to meet an arbitrarily tighter future restoration budget. Rebuild the affected matched source batch where required, before execution/differentiation.
- Actual moving zeros, anchors, captured-label read sets, and all first/adjoint/origin work. Old finite versions must not be quietly replaced by exact references in the emitted sampler or its same-record identity.

The genuine-gradient-reference and finite-restoration clauses are native source obligations. No additional smoothness beyond the global bounded Hessian is introduced. No polylogarithmic deterministic-node-count premise is introduced.

### Preserved output

The quarter harmonic endpoint cancels the direct old-state dependence. Literal finite differentiation therefore yields actual center first I+O(a), full fresh first O(sqrt(a)), and Gaussian-subtracted fresh first O(sqrt(a)a). The short step preserves those scales and creates the coisometry P_new=(0,cos(h)I,sin(h)I).

The new retained state is obtained by replacing the old initial state on the same tape and retaining every new finite force call. Its error acquires the quarter factor a and has grade K+1. If the input actual and retained zeros coincide, deterministic execution on zero Z,L preserves that equality. Exact diagonal normalization also preserves it.

The physical original-gradient force has first O(r), caller O(r/sqrt(a)), and full recorded square curl O(r a). Protection has grade J=2.9; the twin has the same grade and beta comparable to a^.9. Both physical Jacobian sides of the finite proxy follow from its literal finite recurrence and terminal row. Differentiated finite-versus-exact convergence is not assumed.

### Caller and anchor qualification

Adjoin U=(y-y_ref)/sqrt(a) before embedding, and retain all original caller dependencies. An expression constant in the private tape need not be constant in U. In particular, adding B*F_actual(U;0) does not by itself define a full-joint gradient just because it is linear in the private variables. The actual two-node moving-anchor identity or fixed global-origin/flat-gauge construction is the lawful route. The reviewed Section 5 names these routes; they must remain literal in the finite implementation and numerical zero restoration.

Exact coisometry statements also require consistent known coefficient encodings: use the exact known trigonometric row or an explicitly normalized numerical pair, and charge any resulting sampler coefficient approximation. Independently rounded sine/cosine values cannot be called exactly coisometric without that step.

## 6. Child heat and full-host normalization

For the last decoder child at eta=a/2, its own new private L has physical coefficient sqrt(eta)sin(h_eta), h_eta=eta^(19/30). Normalizing the final hidden state by sqrt(a) gives

    gamma_hidden=(eta/a)sin²(h_eta).

Write eta=a ell. Then gamma_hidden is comparable to

    a^(19/15) ell^(34/15).

For ell bounded above/below by admitted public logarithmic factors, this changes only the declared log constants. No new a exponent appears. The coefficient of the last fresh child state is one; it is not multiplied by the inherited 1/2 center contraction.

The child projection's state price is sqrt(eta)eta^J. The host force readout r/sqrt(a) gives

    r sqrt(eta/a)eta^J = r a^J ell^(J+1/2).

This is the claimed O(r a^J) up to logs. A local child force uses its own r/sqrt(eta) and has price r eta^J. Mixing these two normalizations would be an error; the reviewed version keeps them separate.

For eta=A zeta/(1+zeta), an inverse-zeta public-log bound yields the same result. A genuinely smaller power heat eta=A^q with q>1 would change the source query exponent and must be charged explicitly. The new construction does not introduce such a call.

Earlier statistic/observation/decoder records precede the last child's L. Its entering center was exposed before L was drawn. These establish chronology and the last-child protected row only. The rest of a hidden host still needs its supplied statistic, original full graph, caller profiles, and error certificates.

## 7. Cost distinction and the outer application

The deterministic known-center routine performs one complete old call, fixed-depth quarter and short original-gradient groups, and actual finite zeros/anchors/modes. Its proxy, when requested, repeats this complete graph a logarithmic number of times at fixed target. Thus

    Q_R(a)<=Lambda_R[Q_P(a)+a^(-(R-3/2))],
    k_R<=max(k_P,R-3/2).

From k_7=0, the displayed half-order induction has k_R=R-3/2 for R>=7.5. This is a linear query-exponent family, not a sublinear one. Dense matrix arithmetic and memory retain their separate actual bills; no conclusion about their smaller exponent is made.

The existing weak-E empirical arm has normalized error proportional to a^4.9/N. The statistic's physical sqrt(a) insertion yields a^5.4/N. Hence target physical error a^R gives the stated nominal count a^(-(R-5.4))_+. Its multiplication by a complete supplied hidden/force-source cost Q_H remains real. Section 9 correctly does not relabel that hidden recurrence as the standalone known-center recurrence.

LOW30's sufficient-statistic decoder is precisely a chain of fresh known-center Q_(a/2,(z+R_i)/2) calls. Its reference observations share one statistic before Gaussian comparison. The inspected decoder theorem provides C delta_T+C e_R+2^(1-k)sqrt(ad); actual finite caller derivatives, rather than exact-kernel contraction, control common-statistic transport.

For the separate outer proximal Gibbs chain Y=X+sqrt(a)G followed by Q_(a,Y), exact posterior contraction gives rho=(1+a alpha)^(-1). Integrated per-call error delta gives the usual geometric error bill delta/(1-rho). Its required uniform or integrated caller allowance must hold along the actual chain; the reviewed route explicitly states that premise. No inaccessible center is invoked.

With alpha=1/kappa, k_R=R-3/2 and a^(R-1) balanced against epsilon/(kappa sqrt(Dim)), the original-query cost, apart from fixed-order/logarithmic factors, is

    kappa^(2+1/[2(R-1)])
       (Dim/epsilon²)^((R-1/2)/[2(R-1)]).

The dimension/accuracy exponent tends to 1/2. The outer 1/a step count has been included. This does not claim an improvement over the previously supplied denominator-32 algorithm or close the sublinear-order objective.

## 8. Checks and final boundary

`check_weighted_join.py` writes `weighted_join_checks.json`. It verifies:

1. Equal-grid exact row and column bounds and terminal-weight bounds.
2. Four unequal-grid complete graphs, including old-state replication and quarter-to-short connections, with bounded absolute edge norm, affine norm, terminal norm, and protected-L row norm as node count increases.
3. Exact equality of the scaled finite program and its unscaled physical version.
4. The protected/twin coisometry identities and a converged signed reference's gradient symmetry with the normalized nonlinear primitives.
5. Fixed and log-comparable child heats, and twelve steps of the known-center exponent recurrence.

These are diagnostics in addition to the analytic inequalities above. They do not audit every old seed line or instantiate every future caller/finite-precision census.

**Final classification:** the polynomial-node graph obstruction is resolved by the explicit weights; the protected width is independent of R; and the revised relative linear known-center source construction is admitted under its now-explicit complete old genuine-gradient-reference/finite-restoration ports and literal caller/anchor implementation. The hidden-family bill and any stronger outer application remain separate. The original short-refresh audit and its pins are unchanged.
