# Independent audit: positive fourth correction of Gram covariance noise

2026-10-05. Addendum to the retuned native audit. The analytical higher-current argument below is verified independently. The finite realization imports the path-only selected-pair interface in the bounded fourth-packet source; no generic native-forest closure is used.

## Verdict and resolved variance clarification

**PASS, source-qualified, for the smoothed construction and for the sharper unheated path-only specialization in Section 8 of the reviewed main note.** The same-law positive interpolation is mathematically sound, including its factors `t/4`, `t^2/8`, the negative cubic coefficient `-1/(8s^3)`, its one-energy sixth-current remainder, and its one-energy cubic feedback. The second Riesz tensor requires no derivative of the first tensor. All stated physical scalings are correct.

The initial draft required one exact bookkeeping clarification, now resolved explicitly in its Sections 4 and 5. The imported positive fourth packet has an internal independent Gaussian keep inside its allocated `u_P`. If its cubic reference is written

    sP+h(P)+sqrt(r_P)Z_P,

then use its actual `s^2+r_P=u_P`, with both prescribed fractions positive. In the interpolation set

    k^2=u_keep+r_P,
    u_G+s^2+k^2=u.

Both `s^2` and `k^2` are comparable to u under the fixed allocations. Writing `s^2=u_P` while also keeping an additional internal packet Gaussian would count variance twice. A further regression onto the whole packet carrier, paid using the external keep, is an alternative, but should not be silently assumed.

With the corrected exact allocation and the imported complete path/native guards, the corrected integrated-Y W2 bill is

    Lambda sqrt(D) [A^4/u + A^5/u^2 + A^5/u^(3/2)
                     + A^6/u^(5/2) + A^8/u^(7/2)] + absolute floors.

For `u<=1` and the actual public-log smallness of `alpha=A/sqrt(u)`, this is at most

    Lambda sqrt(D) [A^4/u + A^5/u^2] + absolute floors.

This is a substantive replacement of the former unchanged-mixture law debt. It is not obtained by a formal cumulant match alone.

## 1. The second Riesz tensor is legitimate at first regularity

Freeze the exposed label. Let H be the Hilbert space of real symmetric D-by-D matrices, and let

    E=C-E C, ||D C||_(Gaussian-root -> H)<=L,
    ||E||L2(H)<=e,
    R=D(-L_OU)^(-1)(I-E_G).

The Hilbert-valued spectral calculation gives `||R F||L2(HS)<=||F-EF||L2` for every F in L2. Define

    T1_ab = sum_l (R E_a)_l D_l C_b,
    T2_abc = sum_l (R(T1_ab-E T1_ab))_l D_l C_c.

Pointwise Hilbert-operator ideal inequalities and the spectral bound give

    ||T1||L2(H tensor H)<=L e,
    E T1=Cov_H(C),
    ||T2||L2(H tensor H tensor H)<=L||T1-E T1||L2<=L^2 e.

Only `D C` is needed. R acts on the centered L2 class of T1. No derivative of T1, of `D C`, or of an original Hessian is used. Contraction against a smooth function of C is justified by Gaussian divergence duality, differentiating the callback rather than T1. Standard smooth/Sobolev limits preserve the displayed L2 bounds.

The Frobenius Hilbert space embeds isometrically into two physical tensor indices. Its tensor square and cube therefore become physical rank-four and rank-six tensors without a dimension multiplier. For the actual Gram field, `L=O(A^2)`, `e=O(A^2 sqrt(D))`, hence the rank-six coefficient has L2-HS norm `O(A^6 sqrt(D))`.

## 2. Every required proper cut has the claimed size

Set `R_0=||E||op,infinity`, and let `K0=E[E tensor E]` in physical indices. Its full Hilbert-Schmidt norm satisfies

    ||K0||HS<=L e.

This follows either by averaging T1 or by the Hilbert cross-covariance inequality. The relevant physical cuts are separately justified:

1. Original 2|2: for Frobenius-unit T, Gaussian Poincare gives `Var(Tr(T E))<=L^2`, so the covariance operator norm is at most L².
2. Crossed 2|2: the map, up to a transpose, is `T -> E[E T E]`. Its HS-to-HS norm is at most `R_0^2`.
3. A 1|3 cut evaluated on v is `Cov(Ev,E)` up to permutation. Gaussian cross-covariance duality bounds its HS norm by `L ||Ev||L2<=L R_0 |v|`. In this application E is centered already. The transposed cut has the same bound.

Full symmetrization K=Sym4 K0 is an average of permutations. Thus

    ||K||HS<=L e,
    all 2|2 cuts <= max(L^2,R_0^2),
    all 1|3 and 3|1 cuts <= L R_0.

Since `L,R_0=O(A^2)`, every proper cut is `O(A^4)` and the sole HS energy is `O(A^4 sqrt(D))`. After normalization by u² these become `O(alpha^4)` and `O(alpha^4 sqrt(D))`.

This argument does not infer crossed cuts from only the matrix-space covariance norm, and it does not use an abstract bound for an internal trace.

## 3. Exact factors and cancellation at the same law

Use the actual positive shares discussed above. Let P,Z be standard independent physical Gaussians, independent of the covariance bank. Define

    h_i(P)=Q_iabc H3_abc(P), Q=-K/(8s^3),
    W_t=X_t+sP+t^2 h(P)+kZ,
    X_t|C ~ N(0,u_G I+Cbar+t(C-Cbar)).

The covariance remains positive, since it equals `u_G I+(1-t)Cbar+tC`. Let phi be a smooth compactly supported test function initially. The covariance derivative is

    (1/2) E[E:D^2 phi(W_t)].

The first covariance callback derivative costs `t D^2/2`, and the second costs the same. Applying the two Riesz identities gives exactly

    t/4 Cov_H(C):E D^4 phi(W_t)
      +t^2/8 E[T2:D^6 phi(W_t)].

The first tensor may be symmetrized because D4 phi is fully symmetric. The packet derivative is `2t E[h.D phi(W_t)]`. Three integrations by parts in P identify its leading contribution as

    2t s^3 Q:E D^4 phi(W_t) = -t/4 K:E D^4 phi(W_t).

It cancels at the very same W_t. In particular the cubic shift is not set to zero in either leading term, and X_t is not replaced by an independent matched Gaussian during cancellation.

The coefficient `-1/8` also agrees with cumulant normalization: the mixture's fourth cumulant is `3 Sym4 Cov(C)`, while a small cubic `Q:H3` added to `sP` contributes leading fourth cumulant `24s^3 Sym4 Q`. Hence the requested opposite current is precisely Q above. This is only a consistency check; the interpolation is the actual law proof.

## 4. Sixth-current conversion spends five keep derivatives

T2 can remain correlated with X_t and with its entire owned coarse bank. It is independent of the untouched Z. Integrating five of the six test derivatives against kZ yields a physical vector coefficient involving the fifth Gaussian Hermite tensor. Conditional Gaussian Wick isometry gives

    ||velocity_6||L2 <= C k^(-5) ||T2||L2(HS)
                     <= C A^6 u^(-5/2) sqrt(D).

One physical test derivative remains. There is no extra dimension factor: one contracts the five Hermite indices first, leaving a physical vector. The covariance floor in X_t ensures the Price callbacks are regular, but does not replace the separately retained Z used here.

## 5. Explicit cubic-feedback tensors avoid a hidden trace estimate

The feedback argument can be made completely explicit. Work first at normalized variance one, so s and k are bounded positive numerical constants. Write `a=t^2` and

    J=sI+a Dh,
    Dh_ij=3 Q_ijkl H2_kl(P),
    D2h_ijk=6 Q_ijkl P_l,
    D3h_ijkl=6 Q_ijkl.

Three-fold chain rule after packet integration by parts leaves the following test currents, after the leading `s^3 Q:D4 phi` is removed:

- Rank four: Q contracted with `J tensor J tensor J - s^3 I tensor I tensor I`.
- Rank three: three permutations of Q contracted with `a D2h` and J.
- Rank two: Q contracted with `a D3h`.

Let q be a common proper-cut bound for Q and H=`||Q||HS`. Here `q=O(alpha^4)` and `H=O(alpha^4 sqrt(D))`.

For rank four, tensor multiplication in the three input indices gives the pointwise bound H times a sum containing at least one `||Dh||op`. The fixed Gaussian-chaos matrix bound is

    ||Dh||Lp(op)<=Lambda_p q.

One can establish it without separating an unbounded trace: for a symmetric quadratic Gaussian matrix Q(P,P)-E Q(P,P), Jensen and an independent copy P' reduce its Lp norm to that of the difference. Orthogonal Gaussian rotation identifies the difference as twice the decoupled bilinear matrix Q(U,V). Two applications of the Gaussian matrix-series estimate use the actual proper tensor cuts and introduce only public logarithms. Thus the rank-four coefficient energy is at most

    Lambda H [a q + a^2 q^2 + a^3 q^3].

For rank three, the two-index contraction is explicitly a 2|2 Gram product:

    U_(i,l),(a,m) = sum_jk Q_ijkl Q_ajkm.

Its HS norm is at most `||Q||_(2|2,op)||Q||HS<=q H`. Contracting its m index with P has L2-HS norm exactly `||U||HS`. Fixed-degree Hilbert hypercontractivity plus the operator moment bound for J therefore bounds this row by `Lambda a q H(1+a q)`.

For rank two, the three-index contraction is explicitly a 1|3 Gram product:

    V_ia = sum_jkl Q_ijkl Q_ajkl.

Its HS norm is at most `||Q||_(1|3,op)||Q||HS<=q H`, so this row costs `C a q H`.

These are open Gram products with an output index, not an unsupported internal self-trace contraction. All trace subtractions stay inside their exact Wick polynomials. Transfer respectively three, two and one remaining derivatives to the independent keep. The total continuity-velocity norm is

    Lambda q H (1+Lambda q)^2.

With the actual public-log smallness of alpha, this is `Lambda alpha^8 sqrt(D)` in normalized units, hence

    Lambda A^8 u^(-7/2) sqrt(D)

in physical units. This verifies the draft's feedback claim without requiring an unproved all-tensor closure theorem.

## 6. From the interpolation to W2

The independent keep makes the positive curve of laws smooth. Each current above has an L2 continuity velocity. The dynamic transport bound therefore gives W2 at most the integral of their L2 norms. The t coefficients are bounded and integrable. At t=0 the law is exactly `N(0,uI+Cbar)`; at t=1 it is the covariance mixture plus its independent negative cubic packet and keep. Consequently

    W2 <= C A^6 u^(-5/2) sqrt(D)
           +Lambda A^8 u^(-7/2) sqrt(D).

Finite moments, covariance gaps, and the displayed L2 bounds support cutoff/Sobolev approximation and endpoint limits. The cubic map itself need not be invertible, monotone, or have a globally small derivative. Positivity is ordinary positivity of the pushforward law, while the independent Gaussian supplies the smoothing.

## 7. Original-g finite genealogy and coefficient replacements

The separate four-path target is consistent with this coefficient and orientation:

- Differentiating `B_g(X)^2` gives two choices on each covariance side, hence four trees.
- Their edges are `(1,2),(3,4),(hit_left,hit_right)`. The two hit vertices are the two degree-two centers, and the two unhit vertices are the endpoints. Thus there are two C1 and two C0 vertices and no star.
- The two coarse ancestors remain correlated by the covariance-Riesz clock, with the exact common caller Y. The four primitive shields remain independent conditional on that bank.
- The positive coefficient is `a_k^2 (1-q_k^2) product(w_i r_i) r_hit_left r_hit_right`, including squared coarse weights and both derivative-chain correlation factors.
- Full physical symmetrization is either part of the positive path comparison or separately paid by its same-law Hermite interpolation; it is not inferred from a single cumulant contraction.

Replacing the actual j kernel by pure g changes the fourth coefficient by `C A^5 sqrt(D)` in integrated standard-Y HS norm. This uses one Hilbert covariance inequality with a uniformly Lipschitz covariance field and the L2 difference as its sole energy. It does not differentiate a small difference.

The coefficient-only contraction `r_i -> sqrt(1-tau^2) r_i` creates `B_tau=B(P_s_tau g)` and a private shield at least `tau/sqrt(2)`. Its fourth-coefficient bias is `C A^4 tau sqrt(D)`. With `tau=A/sqrt(u_P)` and comparable u_P,u, cubic Hermite isometry converts these biases respectively to

    C A^5 u^(-3/2) sqrt(D),
    C A^5 u^(-2) sqrt(D).

Each actual path root is at an unhit leaf. Its signed amplitude `-omega alpha_P/(8R4 sigma_hit_left sigma_hit_right)` produces the coefficient `-K_tau/8`, with the original adapter's normalization constants included in R4. Distinct primitive hit shields have the required positive dyadic endpoint sums. In particular their inverse square-root factors are independently integrable; an additional caller inverse at one hit produces at most the recorded logarithmic sum. Complete roots, correlated ancestors, native banks, anchors, captures and packet keeps remain in the executed graph.

The imported path-only fixed-order comparison supplies `Lambda sqrt(u_P) alpha_P^5 sqrt(D)=Lambda A^5 u^(-2) sqrt(D)` plus floors, and actual residual/caller first `Lambda A`. This audit checks compatibility of the new genealogy, scalar normalization and buffer use with that interface. It does not numerically execute the selected-pair graphs or replace their actual guards and counts.

## 8. Full return, scope and nonclaims

The complete native Gram-to-mixture comparison costs `Lambda A^4/u` times its inherited one-energy caller profile. Add the independent complete path packet and preserve the independent keep throughout this replacement. Next replace the entire packet by its positive cubic reference, paying its native error and coefficient substitutions. Finally apply the same-law interpolation above. These are conditional W2 compositions, with all complete owned banks consumed and only genuine exposed labels retained.

The field bounds and analytical interpolation are uniform in the frozen exposed Y. The finite source substitutions retain their integrated standard-Y scope. Their combination is therefore an integrated-Y law theorem, exactly as stated. No observer of the owned coarse G, packet P or native roots survives their LAW comparisons.

The new first-port and root/cost assertions come from the literal sum of the actual original-VALUE graphs. They cannot be inferred from the interpolation. Raising the mean order and adding the packet require corresponding fresh complete root dimensions, expanded native counts, anchor copies, and replay costs. The exact independent variance identity must use the actual internal packet keep, as specified at the start.

The proof establishes this rank-four correction only. It does not establish an arbitrary-rank recurrence, a general network feedback theorem, or a zero-buffer covariance source.


## 9. Independent review of the unheated path-only sharpening

Section 8 of the main note clears the audit under the same literal native/finite-filter/path guards. It does not merely set tau=0 in the old seven-clock heat proof. Its crucial distinction is that the original primitive coefficient clocks of B_g are already finite and strictly interior. Only the new covariance-Riesz operator is approximated, before expansion.

The original primitive shield is `sigma_i=sqrt(1-r_i^2)/sqrt(2)`. For a positive dyadic panel of mass O(Delta), the node shields have `sigma_i^2` comparable to Delta. Summing node weights first gives

    sum w/sigma = O(1),
    sum w/sigma^2 = O(number of dyadic panels),
    sum w^2/sigma^2 = O(1),
    sum w^2/sigma^3 = O(1).

For the squared sums, `sum_panel w_i^2 <= (sum_panel w_i)^2 = O(Delta^2)`, so the panel costs are O(Delta) and O(sqrt(Delta)). Their geometric series converge even as the finite lower cutoff tends to zero. The unsquared second-power sum has an O(1) contribution per panel and is only logarithmic in the cutoff. The terminal midpoint obeys the same estimates.

The two hit vertices are distinct degree-two centers; root the path at either unhit degree-one endpoint. The root amplitude includes one inverse shield at each hit. A direct root caller adds an inverse shield at a third, unhit slot, whose first-power weighted sum is integrable. A hit caller adds one inverse to its own hit slot and is multiplied by at least the root and its first ancestor amplitude: this produces at worst the unsquared `sigma_hit^(-2)` logarithmic sum. Every other slot retains a separately integrable first inverse. Shared coarse/Riesz Gaussian ancestors have bounded affine injection coefficients; their common ownership requires adding actual paths, not changing the primitive shield powers.

For whole-coarse feedback the path amplitude is `p=product rho_i`, and its paid bound is `p^2(1+sum_i sigma_i^(-1)) sqrt(D)`. Each hit slot therefore has at worst a squared-weight third inverse, controlled above. Nonhit slots have only squared-weight first inverses. Squaring all coarse, Riesz and primitive weights is essential and is explicit in the draft. Thus the feedback remains `Lambda alpha_P^8 sqrt(D)` with no inverse minimum-gap power.

The complete order-six prior can be summed directly. A root prior has factors `alpha_P^6 omega^6/(sigma_hit_left^6 sigma_hit_right^6)`; panel masses supply sixth powers before the third powers of their gaps are divided out. These sums are bounded. A descendant prior includes a single root weight ratio and extra strict-ancestor powers of alpha_P, so its remaining inverse hit shields are integrable. Its complete total is `Lambda alpha_P^6 sqrt(D)`. The actual complete residual/caller first remains `Lambda A` by the preceding unsquared path calculation. Every actual amplified radius still needs its native order-six admission check.

For the additional Riesz quadrature, the full actual pure-g covariance field has `Lip_HS(C_g)=O(A^2)` and centered HS energy `O(A^2 sqrt(D))`, independently of the smallest individual original clock gap. Approximating R on the unexpanded Hilbert-valued covariance field with operator error epsilon_R therefore gives covariance tensor error

    C epsilon_R A^4 sqrt(D).

Choose `epsilon_R<=A^2/u`. The cubic-reference map costs `u^(-3/2)` per coefficient-HS error, so this is `C A^6 u^(-5/2) sqrt(D)`. The new positive Riesz clock has no inverse gap in an original private shield. There is no need for an unjustified uniform approximation of D2B, a derivative of an original Hessian, or the old tau-dependent telescoping estimate.

Removing coefficient heat removes its `A^5/u^2` bias. The remaining pure-g replacement costs `A^5/u^(3/2)`, and the complete unsmoothed result is consequently

    Lambda sqrt(D) [A^4/u + A^5/u^(3/2)
                     + A^6/u^(5/2) + A^8/u^(7/2)] + floors
      <= Lambda sqrt(D) [A^4/u + A^6/u^(5/2)] + floors.

The absorptions are exactly `A^5/u^(3/2)=(A/sqrt(u))(A^4/u)` and `A^8/u^(7/2)=(A^2/u)(A^6/u^(5/2))`. The A-Lipschitz caller therefore pays `Lambda sqrt(D)[A^5/v+A^7/v^(5/2)]` at comparable variance shares. All retained-label and complete-bank scope remains as in the base theorem.

The final native-mean order formula is also correct. At `u>=c A^beta`, the gradient term is dominated by `A^4/u` when `beta<=2(b_M-4)/(b_M-3)` for fixed `b_M>3`. The other retuned errors and first ports permit every fixed `beta<2` subject to their literal public-log guards. This does not assert uniform cost in a growing b_M.

## 10. Independent evidence and source qualification

`check_positive_gram_correction.py` independently verifies the exact Price/Riesz/Hermite identities with degree-two, -four, -six and -eight polynomial tests, including the `-1/8` factor and the positive cubic covariance feedback. Its unbounded polynomial covariance fixture is explicitly only an algebra test, not an example satisfying the theorem interface. Separate bounded smooth positive covariance fixtures check all physical proper cuts and the two explicit feedback Gram products. The checker also verifies every four-path derivative placement, exact visible/internal/external variance shares, scalar powers, and dyadic shield sums. The JSON output records all reviewed hashes.

The fourth-packet dependency was independently PASS at main source SHA256 `48c23f0c8966a28eae3b6849a0c63348e62f971340d2f5016902580aeb46fee4`. Its independent audit and manifest were inspected. The present verdict uses only its path portion and verifies compatibility of the new covariance genealogy. Neither checker numerically executes the complete native selected-pair/filter program or independently replaces those imported source-qualified interfaces.
