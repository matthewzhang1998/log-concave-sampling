# Independent audit: exact next-history targets and producer obstruction

2026-10-05. This audit verifies analytical target identities and the scope of the existing ports. It does not construct a new native original-VALUE producer.

## 1. Coherent force and centering conventions

All forces use the same stationary OU history:

    F0,t = 0,
    Fk,t = ∫₀∞ e^(−s) g(X_(t+s)−F_(k−1,t+s)) ds.

At a fixed retained state Y, set

    V = F2 − E[F2|Y],
    R = Δ3 − E[Δ3|Y],
    Δ3 = F3 − F2,
    U = V+R = F3 − E[F3|Y].

Let e(Y)²=E[|R|²|Y]. The coherent contraction gives ||e||_(L²Y)≤||Δ3||₂≤C A³√D. Let L bound the complete conditional Gaussian-bank firsts of F2 and F3. Then L≤CA, whereas only first(R)≤2L is available. Small energy is not small first. Finite Gaussian-bank approximation is understood before passage to true histories.

## 2. Exact tensor increments

All expectations below are conditional on this same Y. Define normalized symmetrization Sym_r(T)=(1/r!)∑_(π∈S_r)πT, retaining physical tensor slots. For repeated vectors this convention gives the usual binomial coefficients.

Mean:

    E[F3|Y]−E[F2|Y] = E[Δ3|Y].

Covariance: with C=E[V⊗V], B=E[V⊗R+R⊗V], E=E[R⊗R],

    Cov(F3|Y)−Cov(F2|Y) = J = B+E.

Equivalently J=Cov(F3,R)+Cov(R,F2), with consistent physical output/input orientation. The R⊗R term is mandatory.

Third centered moment (equal to third cumulant):

    κ3(F3)−κ3(F2)
      =3 Sym₃ E[V⊗V⊗R]
       +3 Sym₃ E[V⊗R⊗R]
       +E[R⊗R⊗R].

Fourth centered moment:

    Q4(F3)−Q4(F2)
      =4 Sym₄ E[V⊗V⊗V⊗R]
       +6 Sym₄ E[V⊗V⊗R⊗R]
       +4 Sym₄ E[V⊗R⊗R⊗R]
       +E[R⊗R⊗R⊗R].

For matrices C,D, define the ordered bilinear pairing

    P(C,D)_(ijkl)=C_ij D_kl+C_ik D_jl+C_il D_jk.

This map is intentionally not separately symmetrized. Since κ4(F)=Q4(F)−P(CovF,CovF), the exact fourth cumulant increment is

    κ4(F3)−κ4(F2)
      = Q4(F3)−Q4(F2)
        −P(C,J)−P(J,C)−P(J,J).

The two cross-pairing terms cannot be replaced by one without a declared symmetrization convention. An equivalent and usually cleaner formula uses mixed cumulants:

    κ4(V+R)−κ4(V)
      =4 Sym₄ κ(V,V,V,R)
       +6 Sym₄ κ(V,V,R,R)
       +4 Sym₄ κ(V,R,R,R)+κ4(R).

These are exact identities, not producer interfaces. For delayed displacement −qF3, the mean increment is −qEΔ3, the covariance increment q²J, the cubic increment −q³δκ3, and the quartic increment +q⁴δκ4.

## 3. Dimension-safe stability and the multiple-error-slot trap

The Gaussian cross-covariance inequality gives, for centered R and an L-Lipschitz Gaussian map W,

    ||E[R⊗W]||HS ≤ L ||R||₂.

Apply it to the two-term covariance telescoping identity to obtain

    ||J||_(L²Y;HS) ≤ 2L ||e||_(L²Y) ≤ C A⁴√D.

Applied separately, it only gives

    ||Cov(R|Y)||HS ≤ first(R)e(Y),
    tr Cov(R|Y) = e(Y)².

Therefore the squared-error covariance is O(A⁴√D) in the available dimension-safe integrated HS ledger, although its integrated trace is O(A⁶D). The latter is not an O(A⁶√D) bound. No extra power of A may be assigned to an additional R slot solely from its small total energy.

The imported mixed-quadratic one-energy lemma gives

    ||E[R⊗W1⊗W2]||HS ≤2L²||R||₂

when the other two centered maps have first ≤L. Telescope using U and V in the other slots, rather than expanding all seven terms. There are three terms, hence

    ||κ3(F3|Y)−κ3(F2|Y)||_(L²Y;HS)
       ≤6L²||e||_(L²Y) ≤C A⁵√D.

The sealed mixed Wick-cubic estimate gives

    ||κ(R,W1,W2,W3)||HS ≤6L³||R||₂.

Four-slot telescoping with other slots only U or V yields

    ||κ4(F3|Y)−κ4(F2|Y)||_(L²Y;HS)
       ≤24L³||e||_(L²Y) ≤C A⁶√D.

These arguments never differentiate R. In expanded binomial formulas, an additional R slot has first ≤2L, so every term with at least one R has the same grade5/grade6 guarantee under this method. The estimates do not justify a grade hierarchy according to the number of R occurrences.

## 4. Exact linear test

Take g(x)=Ax in any physical dimension, independently in each coordinate. Set

    Z_j = ∫₀∞ e^(−s) s^(j−1)/(j−1)! X_s ds.

Then

    Fk = ∑_(j=1)^k (−1)^(j−1) A^j Z_j,
    E[Z_j|X0=x]=2^(−j)x,
    Δ3=A³Z3,
    m4(x)−m3(x)=−A⁴x/16.

The centered Z_j has Brownian kernel

    h_j(r)=√2 e^(−r)∑_(ℓ=0)^(j−1) r^ℓ/(ℓ! 2^(j−ℓ)).

The conditional covariance matrix of (Z1,Z2,Z3), per physical coordinate, is exactly

    [ 1/4    1/4     3/16  ]
    [ 1/4    5/16    9/32  ]
    [ 3/16   9/32   19/64  ].

Consequently,

    Cov(F2|x)=[A²/4−A³/2+5A⁴/16] I,
    Cov(F3|x)=[A²/4−A³/2+11A⁴/16−9A⁵/16+19A⁶/64] I,
    J=[3A⁴/8−9A⁵/16+19A⁶/64] I,
    Cov(F2,Δ3|x)=[3A⁴/16−9A⁵/32] I,
    Cov(Δ3|x)=(19A⁶/64) I.

This proves genuine nonzero order-A⁴ mean-history and covariance obligations. Keeping only the Δ3 covariance loses the leading 3A⁴/8 term in J. Both third and fourth cumulants vanish in this linear test, so it proves no nonzero cubic or quartic obligation.

Here J is positive: its coefficient is A⁴(24−36A+19A²)/64, and the quadratic has negative discriminant. The linear test cannot prove sign-indefiniteness of the general increment.

## 5. What positivity can and cannot mean

A covariance increment need not be PSD for general coherently coupled Gaussian maps satisfying precisely the available first/energy bounds. The one-dimensional example

    V=A G, U=(A−A³)G, R=−A³G

has first(U),first(V)=O(A), ||R||₂=A³, but

    Var(U)−Var(V)=−2A⁴+A⁶<0 for 0<A<1.

This is a counterexample to inferring PSD from the exported numerical ports. It is NOT asserted to arise as F2,F3 of one convex-gradient force hierarchy. No same-g negative example for the specific F3−F2 increment is proved in this audit; the F4−F3 increment below supplies a genuine negative example one history later.

Thus an independent additive bank whose covariance equals J is not justified by these ports. An independent R-only bank supplies Cov(R), discards the cross covariance, and fails even in the exact linear example. A signed-covariance correction to a sufficiently buffered positive law is a different possibility; it requires its own gap and producer proof. Positivity of the executing law does not require positivity of every analytical covariance difference.

A universally valid positive analytical representation uses the whole joint block. Put Z=(V,R) in R^(2D) and let P_s denote the OU semigroup on its same complete conditional Gaussian bank. Then

    Cov(Z)=2∫₀∞ E[D P_s Z (D P_s Z)*] ds ≥0,
    Cov(F3)=[I I] Cov(Z) [I I]*.

The off-diagonal blocks are essential. This is a positive full-block Gram identity, not a positive increment identity. Implementing it requires a finite original-VALUE marked-history supplier that preserves the V/R common tape and all roots and callers. Derivatives in this display remain analytical; the identity alone does not authorize derivative leaves or establish a native graph.

## 6. No relative comparison follows from existing marginal contracts

The sealed canonical-mean source exports its own mean, first/curl/energy bounds, and an optional complete-bank marginal LAW. It does not export coupling to F3, or coupling to a previous source strong enough to identify mixed targets. Common-known-carrier alignment preserves each existing whole joint Gaussian law; it does not create the missing V/R correlation.

Even exact equal marginal Gaussian laws can disagree after known-carrier subtraction. For independent G,H,

    T1=G, T2=cos(A)G+sin(A)H

both have exactly N(0,I) law, but residuals after subtracting G have covariance respectively zero and 2(1−cos A)I. Therefore zero marginal LAW error cannot imply a small relative residual error or its needed cross-covariance.

Similarly, a common baseline plus small remainder and exact own mean does not establish a target-relative coupling. Sources S±=A G±A³G share the baseline A G and zero own mean, obey even stronger remainder bounds than required, yet have opposite leading covariance increments relative to A G. Additional paired/common-tape contracts are necessary.

The existing target-restoration comparisons can lawfully pay J, δκ3 and δκ4 at their stated grades, with the original standard-Y integrated scope. They cannot improve those grades or execute them merely by asserting a relative common tape. No pointwise arbitrary-Y promotion is obtained.

## 7. Audit conclusion

The next-history m4 packet may honestly isolate the exact marked-history/mixed-target producer obligations and the direct analytical stability lemmas. The sealed grade-39/10 approximation to m3 also approximates m4 at that same grade by paying ||m4−m3||₂≤C A⁴√D, but this does not produce a new accuracy grade or eliminate the nonzero order-A⁴ target change. A claimed precision beyond grade four must address that target change explicitly.

The narrow obstruction is missing coherent mixed-target production and closure, not an impossibility theorem for positive laws or arbitrary-order constructions.

## Sources checked

- /workspace/shared/quartic-corrected-mean-join-20261005/README.md and ACTUAL-QUARTIC-CORRECTED-MEAN-JOIN.md.
- /workspace/shared/law-only-reentry-ceiling-20261005/REUSABLE-PORT-AND-DEPTH-RECURRENCE.md, Sections 1–2 and 5.
- /workspace/shared/fourth-cumulant-return-20261005/comparison/FOURTH-CONDITIONAL-CUMULANT-COMPARISON.md, Sections 1–6.

No imported source was modified. The exact linear covariance formulas were checked by symbolic integration of the Brownian kernels.

## 8. A genuine same-g obstruction one history later

For the same g(x)=Ax, the coherent fourth-history increment is Delta4=−A⁴Z4. The exact covariance entries are Cov(Z1,Z4)=1/8, Cov(Z2,Z4)=7/32, Cov(Z3,Z4)=17/64 and Var(Z4)=69/256. Therefore

    Cov(F4|x)−Cov(F3|x)
      =[−A⁵/4+7A⁶/16−17A⁷/32+69A⁸/256]I.

For 0<A≤1/2, divide the coefficient by A⁵ and write it as

    −1/4+7A/16+A²(−17/32+69A/256).

The first two terms are at most −1/32, and the bracket in the last term is negative. Thus this true-history increment is strictly negative definite. A repeatable architecture that keeps the old covariance and only adds independent positive covariance increments cannot be correct at every depth. This does not obstruct replacing the whole covariance bank, using a full joint block law, or implementing a buffered signed covariance update. It also does not change the positive sign of the specific F3−F2 linear increment under current audit.
