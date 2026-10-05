# The exact covariance-test port behind the missing width gain

2026-10-05. General tensor lemma, followed by an explicitly unproved source-specific hypothesis. This does not construct a positive sampler.

## 1. An arbitrary fixed-rank tensor lemma

Let F be a square-integrable real random rank-r tensor on D-dimensional physical slots. Assume every proper matrix flattening of F has operator norm at most K pointwise. Put Fbar=F-EF and h^2=E||Fbar||HS^2. Suppose the scalar test inequality

    Var <H,F> <= v^2 ||H||HS^2                      (V)

holds for every deterministic tensor H of the same rank. Equivalently the covariance operator of vec(F) is bounded by v^2 I.

Then the connected rank-2r coefficient T=E[Fbar tensor Fbar] satisfies

    every proper cut(T) <= max(4K^2, 2Kv, v^2),
    ||T||HS <= v h.                                (C)

For r>=2, h<=K sqrt(D) by isolating one physical index in F, so v=O(K) gives one dimension-sized Hilbert factor and dimension-free cuts. For r=1 the proper-cut hypothesis on F is empty; the displayed covariance conclusions still use the separately declared h, and no h<=K sqrt(D) inference is made.

### Proof

Every proper cut of Fbar is at most 2K. Fix an external cut of the two old tensor blocks. Regard the two reshaped copies Fbar(omega) as matrix families A_omega,B_omega. Operator Cauchy-Schwarz gives

    ||E[A tensor B]||op
       <= ||E[A A*]||op^(1/2) ||E[B* B]||op^(1/2),

or the reverse orientation.

If the chosen cut splits both old blocks internally, both local factors are proper cuts and the bound is 4K^2. If exactly one old block lies wholly on one side, use that block as a vector. Its E[A A*] is the covariance operator controlled by (V); the other factor is a proper cut, giving 2Kv. If each old block lies wholly on opposite sides, the global cut is precisely the covariance operator, giving v^2. Putting both old blocks entirely on the same side would not be a proper global cut. These cases exhaust all cuts.

As an operator on the whole old tensor space, T is positive semidefinite. Therefore

    ||T||HS^2 = Tr(T^2) <= ||T||op Tr(T) <= v^2 h^2.

This proves (C). No separate scalar trace is bounded by D, and no second Hilbert factor is multiplied in.

The hypothesis (V) is also necessary for the balanced old-block-versus-old-block cut in (C), with its corresponding constant. Thus it identifies the exact missing variance control for that cut; it is not an arbitrary stronger regularity requirement.

## 2. What this would repair in the tight current

For the center-center six-tree Price history, write its leading conditional coefficient as

    C(Q) = A^7 beta J(Q),    beta=w^2/sigma_2^4,

with the complete old/Price Gaussian bank Q, all six physical slots, and the original private shields unchanged. The source-qualified native target already has pointwise proper cuts of J bounded by fixed-order/public-log constants and HS at most Lambda sqrt(D).

A genuinely new source-specific theorem of the form

    Var <H,J(Q)> <= Lambda ||H||HS^2,

UNIFORMLY in the inherited original shields and the declared Price clock, would imply

    Cov(C) all proper cuts <= Lambda A^14 beta^2,
    ||Cov(C)||HS <= Lambda A^14 beta^2 sqrt(D).

The missing sigma_2^-2 factor in the generic derivative/Price bound would disappear, and the original squared clock weights would then be logarithmically summable.

This proposed source-specific variance theorem is OPEN. The tensor lemma in Section 1 does not prove it. The native Jacobian matrix-test theorem controls linear tests of one Jacobian uniformly in heat; its product/selector extension uses an additional root-width ledger. It cannot simply be relabeled as a theorem for this six-slot product of higher heat jets.

The elementary coefficient F(q)=sin(q_1) sum_i e_i^(tensor 3) shows why proper cuts alone are insufficient: all its proper cuts are at most one, but for H=D^-1/2 sum_i e_i^(tensor 3), the variance of <H,F> equals D Var(sin G). Such a field is not automatically a source-qualified original-gradient six-tree. A counterexample to the source-specific hypothesis would need to obey that actual genealogy and its bounded original Hessian.

## 3. A distinct producer obligation remains

Even a proof of this variance inequality would be an analytical coefficient certificate. The direct pointwise Riesz-expanded normalization still contains its private inverse shields. The required positive grouped original-VALUE consumer must realize the connected coefficient without reintroducing those losses into its root radius or its intrinsic descendants. Large LAW-only caller factors may price absolute numerical floors, not this intrinsic current.

Thus the continuation has two separate tests:

1. A source-qualified uniform scalar-test covariance estimate for the complete grouped coefficient (or another actual one-Hilbert/common-heat estimate).
2. A positive sign-programmable grouped VALUE/current producer preserving that estimate, the actual observer boundary, and every history-specific descendant/cost.

Neither is supplied by ordinary Gaussian-terminal transposition, external-probe rerooting, or the bounded auxiliary-cycle construction. No general recurrence or eventual-sublinear conclusion is asserted here.
