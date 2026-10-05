# Independent normalized-ancestry audit

2026-10-05. Read-only audit of `NORMALIZED-ANCESTRY-RESUMMATION.md` and the sealed stage-two port derivation. No sealed inputs edited. Native completion has not been executed or independently certified.

## Verdict

The claimed radial bound, sine-source lower bound, ordered raw derivatives, private first/curl bounds, conditional energy bound, and original-VALUE error propagation are valid with the explicit Gaussian and normalized-rule assumptions below. The source is an inexpensive radial-obstruction escape; it does not achieve a uniform higher-order approximation. No native own-mean completion theorem has been proved by this audit.

## Necessary corrections and scope qualifications

1. In the radial lemma, “Y has covariance q² I” must read “Y is a centered Gaussian with covariance q² I.” Covariance alone does not imply the radial concentration estimate. All actual applications in the proof meet the Gaussian assumption.
2. The outer-rule hypotheses must include omega_i≥0, sum_i omega_i=1 and t_i∈[0,1], in addition to sum_i omega_i t_i=1/2. The moment condition alone does not imply beta≤√3/2 or the stated first/energy bounds. With these assumptions Jensen gives beta=sum omega_i√(1-t_i²)≤√3/2.
3. The imported bracket stated in sealed `SHIFTED-BIAS-AND-PORTS-DERIVATION.md`, §7, is ell_E(a_seed+mu)+ell_E³(1+mu^(-1/2)). At ell_E=2A, a_seed=2A, mu=A it is exactly 6A²+8A³+8A^(5/2). The displayed 16A²+8A³+8A^(5/2) is a valid upper bound, not the exact substitution.
4. The sine obstruction rules out the three-VALUE graph and the specified fourth-probe variant. It does not prove that every finite variance-normalized nesting fails. The one-dimensional sine source is itself radial, so the conclusion should say that this construction loses nonlinear history information, without asserting that the missing information must be nonradial.

## Radial proof checked

The retained Gaussian rows give Var(v)=I/2, Var(w)=3I/8, Cov(v,w)=3I/8 and therefore Var(v-lambda w)=s²I. The history integrals use the same genuine OU coupling. Linearizing at lambda(1) gives the stated errors for F1,a; then linearizing the Gaussian X_a-lambda H1,a at radius q gives (4). The executed graph first replaces h by lambda w and only then applies the Gaussian radial lemma to (v-lambda w)/r. Thus no Gaussian identity is applied to an actual nonlinear descendant.

For 0≤lambda≤A≤1/2, qhat≥√7/4>1/2, q+qhat>1, and q-qhat≤lambda/2. Hence |lambda(q)-lambda(qhat)|≤A lambda≤A². The terminal Lipschitz factor gives

    ||R1 psi_tilde-m3||₂
      ≤ √2 A²[q+s+A(1+k)] + A³ s√D
      ≤ 3.554334 A² + A³√D/√2.

This proves the advertised weaker 4A² bound. Quadratic exactness is algebraic for every symmetric matrix K: the raw terminal is Kx-K²v+K³w. The optional radial probe removes the coefficient mismatch and retains the same dimension-free constant; its external source/anchor dependence is correctly excluded from the frozen-source private-port claim.

## Sine obstruction checked

For f(x)=(x+sin x)/2, Gaussian integration by parts gives exactly the displayed C(u). The linear-in-u contributions cancel between rC(r) and integral_0^1 C(u)du. The remaining d is negative, approximately -0.0001347198783440454. Its j=2 series contribution is -exp(-1)/2880; all later contributions are negative.

The first-chaos remainders are bounded by the two displayed expressions: use Cauchy–Schwarz for the nested-force difference and ||X||₄=3^(1/4) for the quadratic terminal remainder. Their summed coefficient is at most 3.031253 for A≤1/2, so 4A³ is conservative. R1 contributes its exact first-chaos factor 1/2. At A≤10^(-5),

    exp(-1)/5760 - 2A ≥ 0.00004386795 > 1/30000.

The claimed lower bound follows. The sign of the terminal-mean difference is -d A² at leading order; taking its absolute value makes the packet's bound unaffected.

A useful strengthening is immediate by tensorization: g(x)_j=A(x_j+sin x_j)/2 yields the same independent coordinatewise graph and target, hence ||R1 psi_tilde-m3||₂≥A²√D/30000 for every D and the same A range. This is an anchored smooth convex gradient with 0≤Dg≤AI.

## Raw native-port audit

The exact ordered derivative is H3[x_y-H2 v_y+H2 H1 w_y]. No factor commutation or derivative of an HVP is used. Subtracting g(x) gives E_x=(H3-H0)-H3H2/2+H3H2H1/4, E_N=-H3H2/2+H3H2H1/2 and E_M=H3H2H1/4. The PSD interval implies ||H3-H0||≤A, establishing (8).

For the square lift on (G,N,M), the leading G block is symmetric. Triangle inequality on its remaining skew part plus the off-diagonal block norm yields exactly the stated curl estimate. At A=1/2 the raw first coefficient is 1.198550557 and raw curl coefficient is 1.842877072. Under the stated half-variance scaling they become 1.695006453 and 2.606221748, so declarations 2A and 3A² are safe. Curl≤ell_E a_seed is also safe with a_seed=2A.

Pointwise |E|≤A²(|v|+A|w|). Conditional on Z=z, the centered Gaussian row norms of v and w are at most r and k, respectively; their mean terms integrate to |z|/4 and |z|/8. This proves (10) and the nonzero-caller origin bound. The private-root dimension is 3D, with four original VALUES per residual outer node and three per raw terminal. First/adjoint sweeps must include every original source site actually used, including the residual baseline and captures.

The uniform original-VALUE terminal floor is (1+Ar+A²k)nu. The residual baseline and independently executed capture each retain their own allowance. The printed complete positive-law bill remains only a conditional interface: full native guards, independent complete banks, captures, replay, outer-rule accuracy and all absolute floors must be verified separately.
