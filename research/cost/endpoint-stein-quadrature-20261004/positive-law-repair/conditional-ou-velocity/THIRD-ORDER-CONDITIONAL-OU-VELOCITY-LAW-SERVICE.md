# A native conditional OU-velocity correction and its positive mean-law consumer

2026-10-04. New source-qualified construction, pending independent review. This is a conditional buffered-force law service, not yet a same-endpoint posterior sampler or a sublinear complete-cost recurrence.

## Result first

Fix an exposed endpoint caller z and a known OU node `0<=r<1`, with `s=sqrt(1-r^2)`. The exact reverse-OU velocity from the audited positive-law proof is

    m_r(z)=E[g(rz+sG) exp(-U(rz+sG))]/E[exp(-U(rz+sG))].

There is a literal original-gradient VALUE source whose mean approximates m_r through order three, and whose actual split matches the already admitted gradient-mean / near-gradient-mean consumer:

    B(u)=g(x_M+s u),
    E(u,v)=g(x_M+s[u-H_alpha(u,v)])-g(x_M+s u).

Here x_M is a finite conditional mode and H_alpha is the SAME two-root endpoint Stein packet applied to the anchored conditional potential, with Hessian bound `alpha=A s^2`. Every force call is an original gradient VALUE. The shared roots u,v inside E are intact.

The new term E is an actual near-gradient source with

    ||E||_Lp <= C_p A^2 s^3 sqrt(D),
    Lip_(u,v) E <= C A s,
    Curl(P_u^* E) <= C A^2 s^3,
    E(0,0)=0 exactly with recorded-anchor reuse.

Its first need not be O(A^2); its SMALL CURL, rather than an invented small first, is what matches the consumer. B is a genuine gradient in its private u variable and has first As.

Using independent COMPLETE mean banks for B and E, positive output variance shares 1/2+1/2, and the admitted padding mu=A, one obtains an actual positive law service M_r with

    W2(Law(M_r|z), N(m_r(z),I))
        <= Lambda_k sqrt(D)[A^3 s^5 + A^4 s^3] + numerical/mode floors. (R)

Here Lambda_k retains the imported fixed-order public-log factors. The raw source facts hold for A<=1/2; the buffered consumer additionally requires the normalized-radius, covariance-gap and finite-clock guards stated in Section 6. Within those guards this is third order in A for the force/velocity, up to the displayed public-log factors, and has a fixed known unit Gaussian buffer. It is a law return, not a strongly accurate random estimate of m_r. For each fixed requested consumer order and admitted clock nodes, its complete query and private-tape counts have no inverse-A exponent. All original first/adjoint sweeps remain at recorded original-gradient VALUE sites.

The unresolved next step is to compose such noisy conditional services into the exact OU probability-flow endpoint while preserving the endpoint Gaussian variance and all nonlinear current descendants. This note does not call the unit buffer free or subtract an unobserved Gaussian reference from the actual output.

## 1. Inputs, source pins and conditional target

Use the original setting

    g=grad U, g(0)=0, 0<=Dg<=A I, 0<A<=1/2.

The source g is the original mode-anchored gradient VALUE `sqrt(A)[grad V(b_M+sqrt(A)x)-grad V(b_M)]`, with its already priced physical finite-mode residual. Freeze its version and anchor before all subsequent private sampling.

Inputs imported with their actual scope:

1. The new positive endpoint theorem, `../POSITIVE-SECOND-ORDER-LAW-AND-QUADRATIC-AMPLIFIER.md`, SHA a085dfdef3f66f69208612034c17e4e18045fc00a21be3aed64ccfcb5e9b7085, with independent PASS at `../independent-audit/INDEPENDENT-POSITIVE-LAW-AUDIT.md`, SHA df32f63c88ec075180b231bd997f13524792296aaa161b1879c987bfcdd28512.
2. `no-copy-rank-20261004/host-observer-audit/COMPLETE-MEAN-FIRST-ORDERING-AND-REENTRY-RECIPE.md`, SHA a7cde258f1c86defbfdc681f473076293cfe65253479493c6ef5baa8b08a88a0, Sections 2 and 4: the existing complete gradient-mean / near-gradient-mean composition and actual replay rule.
3. Its independent audit SHA d9b4d45c276e72adc4b2590c1c17140bed2fe94dd2a6de5d8f7cc1b3a0a7d21a. This does not supply a generic all-rank decomposition; the exact new B+E decomposition is proved below.

At fixed z,r the conditional law of q is

    nu_(r,z)(dq) proportional to exp(-|q-a|^2/(2s^2)-U(q)) dq,
    a=rz.

It is strongly log-concave. Its exact mean force is m_r(z). The mode x_* solves

    x_*=a-s^2 g(x_*).

No conditional density, potential value, mean or Hessian action is queried.

The target m_r throughout is for the ANCHORED U. If the original global finite mode leaves the already recorded linear tilt ell, the true total OU drift is ell plus the conditional expectation of g under U(q)+ell dot q. Adding the known ell to the returned service restores its constant part; the remaining conditional mean price is at most A s^2 |ell|, because the two conditional q laws have W2 at most s^2 |ell|. This extra original-mode floor is retained whenever the service is used for that true finite-mode posterior.

At r=1, s=0, simply return `g(z)+N` for a fresh unit Gaussian N if a unit-buffer law service is requested. This endpoint requires one original VALUE and has exact target mean. The division by s below is analytical and is not used at s=0.

## 2. Finite conditional mode and literal inside packet

Execute a frozen number M of iterations

    x_0=a,
    x_(k+1)=a-s^2 g(x_k), k=0,...,M-1.

Record g(x_M) once. The contraction is `alpha=A s^2<=1/2`. With

    R_M=x_M-a+s^2 g(x_M)=x_M-x_(M+1),

one has

    |R_M| <= s^2 alpha^M |g(a)|,
    ||D_z x_M|| <= r/(1-alpha),
    |x_M| <= |a|/(1-alpha).

R_M is an analytical residual or an extra explicitly counted last VALUE; no adaptive stopping is required. The exact anchored conditional gradient is

    gbar(xi)=s[g(x_M+s xi)-g(x_M)],
    gbar(0)=0, 0<=D gbar<=alpha I.

Its original query at xi is the actual point `b_M+sqrt(A)[x_M+s xi]`. This is an original-gradient VALUE with the same global and conditional anchors.

Choose the audited positive quadrature Q_alpha with error delta_alpha<=alpha. Draw fresh independent standard u,v in R^D after z,r,x_M are exposed. Execute

    H_alpha(u,v)=sum_i w_i gbar(r_i u+sqrt(1-r_i^2)v),
    Y=u-H_alpha(u,v),
    Q_fin=x_M+sY.

All nodes use the same v. The inner actual bounds are

    ||H_alpha||_Lp<=C_p alpha sqrt(D),
    ||D_(u,v)H_alpha||<=alpha,
    0<=D_u H_alpha<=(alpha/2)I,
    H_alpha(0,0)=0.

The numerical zero is literal only when every repeated site reuses its recorded anchor value. The positive endpoint theorem and linear-tilt restoration give

    W2(Law(Q_fin|z),nu_(r,z))
       <= C A^2 s^5 sqrt(D)+|R_M|+e_num,state.               (1)

Indeed the standardized source error is `C alpha^2 sqrt(D)`, its physical conditional scaling is s, and the true standardized conditional potential has the additional linear term `R_M dot xi/s`. Strong convexity charges exactly |R_M| after scaling. The finite conditional-mode profile must not be folded into the intrinsic sqrt(D) term for an unbounded fixed caller.

## 3. The second-order native force term and its exact mean target

Execute on the same u,v record

    F(u,v)=g(Q_fin),
    B(u)=g(x_M+s u),
    E(u,v)=F(u,v)-B(u).

This is a literal pathwise identity F=B+E. There is one terminal original VALUE at Q_fin and one baseline VALUE at x_M+s u. By (1) and Lip(g)<=A,

    |E_(u,v) F-m_r(z)|
       <= C A^3 s^5 sqrt(D)+A|R_M|+e_num,force.              (2)

This is a mean difference between two laws, proved by their W2 coupling. It does not assert that F, E, or a completed buffered mean is strongly close to the exact deterministic m_r.

The exact Gaussian first-order velocity is `P_r g(z)=E_u g(a+s u)`. Consequently the complete native correction to it is

    [g(x_M+s u)-g(a+s u)] + E(u,v).

The first bracket is a genuine-gradient translation chord. Its mean plus E's mean approximates `m_r(z)-P_r g(z)` with the SAME error in (2). This identifies the actual intended second-order current. It is not merely an arbitrary small residual.

For a fixed potential amplitude expansion, the corresponding analytical second-order coefficient is

    -Cov_gamma(U(a+s u),g(a+s u))
      =-s^2 integral_0^1 E[Dg(a+s u)
                g(a+s(tu+sqrt(1-t^2)v))] dt.                (3)

Equation (3) is only a coefficient identity, obtained from Gaussian covariance integration by parts. No HVP is executed. The finite VALUE graph above instead matches the FULL conditional mean m_r via (2), including its finite mode and shared-root nonlinear shifts. Thus no pointwise Hessian Taylor remainder or dimension-dependent small-displacement premise is used to justify (R).

## 4. Full private first, curl, energy and exact zero

Let B_1=Dg(Q_fin), B_0=Dg(x_M+s u), H_u=D_u H_alpha, H_v=D_v H_alpha. These are analytical actual firsts at the recorded VALUE points. Then

    D_u E=s[(B_1-B_0)-B_1 H_u],
    D_v E=-s B_1 H_v.                                      (4)

Both B_i are symmetric in [0,A I], and H_u is symmetric. Their products are not commuted. Thus

    ||D_(u,v)E||<=As(1+alpha),
    ||D(P_u^*E)-D(P_u^*E)^*||<=C A^2 s^3,                 (5)
    P_u=(I,0), P_u P_u^*=I.

Explicitly the upper-left skew block is `s(H_u B_1-B_1 H_u)`, and the off-diagonal blocks are `-s B_1 H_v` and its negative transpose. Each is bounded with the actual two factors. The leading term B_1-B_0 is symmetric and costs no curl. This is why the source can have first O(A) while having curl O(A^2) with no Hessian-continuity modulus.

The actual VALUE mark is even simpler:

    |E|<=A s |H_alpha|,
    ||E||_Lp<=C_p A^2 s^3 sqrt(D).                          (6)

This has ONE physical Gaussian energy, not a product of two sqrt(D) factors. At all private zeros, H_alpha=0, Q_fin=x_M and the baseline site is also x_M. Recorded same-site reuse therefore gives E(0,0)=0 exactly. Its centered energy and anchored energy are both bounded by (6); there is no hidden origin debt.

B is a genuine gradient on u: for s>0 it is the gradient of `U(x_M+s u)/s`, and `D_u B=s Dg(x_M+s u)` is symmetric. Its centered/anchored version B-B(0) has first As, energy at most C_p As sqrt(D), and exact zero; B(0)=g(x_M) is the recorded deterministic caller term. Embedding it as `(B,0)` on a padded tape preserves the exact full-gradient property. This is not a claim that the entire rectangular F or E on (u,v) is a gradient.

## 5. Caller and query geometry

Treat z as the exposed normalized caller. No private u,v variable is part of z. The conditional mode has `||D_z x_M||<=2r`. Differentiating only once,

    D_z gbar(xi)=s[Dg(x_M+s xi)-Dg(x_M)] D_z x_M,
    ||D_z H_alpha||<=2sAr.

Thus `||D_z B||<=2Ar` and, by subtracting the two terminal firsts as in (4),

    ||D_z E||<=C Ar.                                      (7)

No derivative of this small VALUE mark is declared to be smaller than that.

For the original physical caller y, retain its finite b_M record. The original normalized g has `||D_y g(x)||<=C sqrt(A)` at fixed normalized x. The conditional mode recurrence gives `||D_y x_M||<=C s^2 sqrt(A)` (in addition to whatever explicitly captured labels have been declared). Differentiating gbar once then gives `||D_y H_alpha||<=C s sqrt(A)` and `||D_y B||,||D_y E||<=C sqrt(A)`. These are the ordinary canonical-force caller bounds. They do not need a third derivative or a modulus of the original Hessian. Requested first/adjoint sweeps use original HVPs only at stored VALUE points.

The private sites satisfy a caller-uniform displacement envelope

    |x_M+s(r_i u+s_i v)-x_M|
        <=s sqrt(|u|^2+|v|^2),
    |Q_fin-x_M|<=C s sqrt(|u|^2+|v|^2).

Together with `|x_M|<=2r|z|`, the original physical site is always `b_M+sqrt(A)` times this normalized profile. The finite-mode residual retains `|g(a)|<=A r|z|`, rather than being called uniform in an unbounded caller. No algebraic D-versus-A guard is used.

## 6. Match the admitted positive mean consumers

Freeze z,r,x_M and all versions first. Choose numerical constant variance shares `v_B=v_E=1/2`. Use the existing exact-gradient VALUE mean at fixed order k>=4 for the anchored B source, and the existing physical near-gradient VALUE mean for the anchored E source, with declared scale A, curl parameter kappa=O(A^2), and padding mu=A. Use COMPLETE independent execution banks between these two completed mean branches:

    M_B=B(0)+sqrt(v_B) N_k((B-B(0))/sqrt(v_B)),
    M_E=sqrt(v_E) N^v(E/sqrt(v_E)),
    M_r=M_B+M_E.

These are the positive programs from the imported consumer, not formal Gaussian draws. The join has the following EXPLICIT imported-domain gates, in addition to the raw source bounds:

- The normalized gradient radius is rho_B=sqrt(2)As. Require rho_B<=r_*(k,Delta_B,Lambda_k), the literal b27:compiler:mean threshold for its actual complete active dimension and public-log factor. A<=1/2 alone does not imply this.
- The normalized near-gradient square lift has ell_E<=sqrt(2)As(1+alpha). Require ell_E<=1/4 (or an explicitly supplied no-weaker covariance-gap constant in the chosen eta,zeta allocation) and mu=A in (0,1]. Then the nine-main covariance satisfies Cov(S)<=ell_E^2 I/3, so the path I-t^2 Cov(S) has gap at least 47/48.
- Every inherited finite clock, filter, numerical-radius and source-restoration guard is applied at its actual normalized field and positive output share. Their fixed-order public-log factors are denoted Lambda_k throughout the mean-law and returned-program bounds. No algebraic D-power heat guard is added.

For any fixed consumer order these are ordinary small-A/public-log restrictions, rather than an extra empirical bank count. Their raw sources otherwise meet the input ports exactly:

- B is an actual full-gradient VALUE source at private first As, with its captured origin and bounded caller.
- E is an actual near-gradient VALUE source with full first O(A), full square-lift curl O(A^2), coisometry P_u, mark (6), exact zero, and the caller (7).
- The SAME u,v tape inside each E occurrence contains its full H_alpha ancestors and terminal/base queries. No root is replaced by an independently sampled force inside E.
- Independence is imposed only between COMPLETE mean branches and where the imported mean circuit itself requires independent complete source occurrences.

The gradient branch has error at most

    Lambda_k sqrt(D)(As)^k + e_num,B.

The imported near-gradient law estimate is its literal

    Lambda_k e_E {A(A+mu)+A^3(1+mu^(-1/2))}+e_num,E.

For mu=A, (6) makes this at most `Lambda_k A^4 s^3 sqrt(D)+e_num,E`. With k>=4, `(As)^k<=A^4 s^3`. Conditional convolution of the two independent complete mean laws therefore yields

    W2(Law(M_r|z),N(EF,I))
          <=Lambda_k A^4 s^3 sqrt(D)+e_num,mean.           (8)

Combining (2) and (8) proves (R). The covariance I in the reference is exact because the two target variance shares sum to one. It is not the covariance of the raw F or E.

The caller z, r, conditional finite-mode record and all deterministic parameters are retained in this conditional comparison. The mean circuit's private fine/source/clock banks are integrated out. There is no licence to append one of those private banks, an old H_alpha value from that bank, or its source-zero Gaussian carrier as an additional observer. A downstream program may use the complete M_r output and the captured caller, as provided by the law-only host rule; it may not identify an unobserved coupling Gaussian with one of the actual private Gaussian roots.

The imported actual first/caller program bounds retain their Lambda_k factors: O(Lambda_k A) for the carrier-subtracted correction and O(Lambda_k Ar) in z, O(Lambda_k sqrt(A)) in y, with their actual known Gaussian carrier and origin. The mean law (8) itself is not a proof of those derivatives: their proof is the consumer's finite graph applied to (4)-(7). Any stronger native retained/protected/proxy rebuild still needs that consumer's separate graph-return premises and is not admitted here by a law coupling.

## 7. Complete original work, zeros, and finite versions

The conditional mode and anchor depend only on the captured z,r and original anchored potential. Cache only that actual caller-only subgraph. It costs M+1 original g VALUES, counting the recorded last anchor. A raw B occurrence needs one additional original VALUE at x_M+s u. A raw E occurrence needs

    n_alpha original inner packet node VALUES,
    one terminal original VALUE at Q_fin,
    one baseline original VALUE at x_M+s u,

hence at most `n_alpha+2` private original VALUES. Its captured anchor g(x_M) is reused. Terminal/private coincidences may only be cached under an exact complete semantic key. Every changed raw source argument inside a compiler is a new full execution.

If the admitted fixed-order gradient and near-gradient programs have complete private occurrence counts N_B and N_E, the safe total is

    Q(M_r)<=M+1+N_B+(n_alpha+2)N_E+actual numerical/known-row work. (9)

A first/adjoint direction has the corresponding original HVP site ledger. Captured-mode outgoing directions may be accumulated before one reverse sweep through its stored graph. If its primal record is discarded, its replay is paid. The literal private Gaussian dimension is O(D(N_B+N_E) plus the consumer's clock/fill blocks); it is NOT just 2D, even though each raw E call has two D-dimensional roots.

The dyadic inside quadrature has `n_alpha=O(log^2(1/(A s^2)))`. To call (9) zero inverse-heat exponent, restrict to the supplied finite outer clock whose nonzero widths have `|log s|<=C_J(log(1/A)+L)` and whose node/share counts are fixed-order public-log polynomials. These are the same explicit clock hypotheses as the imported consumers. An arbitrary r with exponentially tiny 1-r has its literal larger logarithmic bill; it is not covered by a free uniform node-count assertion.

All M, quadrature versions, output shares, mean orders and numerical tolerances are frozen before differentiation. Coherent original-gradient numerical errors are propagated by the actual absolute quadrature weights and every consumer readout/inverse width. They are absolute floors, never divided by a possibly zero actual E energy. At all private zeros the raw E cancellation is exact; the gradient branch retains the original captured B(0)=g(x_M). If the entire original caller and all roots are zero, the anchored conditional mode and all g values are zero by same-site anchor reuse. The completed programs retain their imported source-zero carrier convention; no Gaussian Lp estimate is evaluated at zero.

## 8. What has and has not advanced

The new native term is the actual shared-root force chord E, together with the explicit conditional-mode translation of the genuine-gradient branch. Equation (2) directly ties it to the true OU velocity m_r, so it is a meaningful nonlinear second-order correction. Equations (4)-(6) explain its exact match to an already admitted positive VALUE mean consumer without an HVP-valued source, a Hessian continuity rate, or a missing rank-one covariance.

The conditional completed output has a unit Gaussian buffer. Scaling it to a much smaller variance changes every raw radius, curl, padding error, caller and numerical factor; those factors must be recalculated rather than calling (R) a small-noise or strong mean oracle. Even a perfect conditional Gaussian law does not by itself produce the correct final posterior endpoint when inserted into a nonlinear flow. The baseline Gaussian variance, same endpoint, all stage/root descendants, and every time quadrature error still need a positive composition theorem.

Thus this note closes a reusable conditional velocity-law SUBPROBLEM with a fixed order gain. It does not close the general all-order posterior-law amplifier or show `c(P)/P -> 0`.
