# Independent audit of the actual-mark split and retained orientation estimate

2026-10-04. Result: PASS within the fixture-specific, finite-clock and stated ownership scope. Audited source: TERMINAL-FLAT-GRADIENT-SPLIT-AND-ACTUAL-MARK.md in this directory, including its equation (15).

## 1. Remainder source and actual energy

The exact remainder identities are correct. With u=g−T id and eta=u(S)−u(S+epsilon Z), the original feedback point is x1=xbar1+a eta, and

    E−E0=raT[g(xbar1)−g(x1)]+r[u(t0)−u(t1)].

The first term is bounded by C r a² q, which is C a times the actual-E energy under the proved lower bound e_actual>=c r a² epsilon sqrt(d). The second is zero unless at least one of its two endpoint first coordinates is in the bump interval of length .2q. Its pointwise magnitude is at most C r a² epsilon |Z|.

Conditioning on S2,U,Z leaves S1 standard. Both terminal first-coordinate maps have derivative at least a fixed numerical positive constant: the t0 derivative is 1+ac H11, and the t1 derivative is

    1+ac H11(x1)+a² e1*H(x1)[H(S)−H(S+epsilon Z)]e1
       >=1−2a²(.6)².

The inverse image of the bump interval therefore has length O(q) and conditional probability O(q), uniformly in the conditioned variables. This justifies keeping the |Z|^p moment before integration and proves the fixed-p gain q^(1/p). It does not require differentiating a small energy bound or estimating third derivatives.

The whole-gradient lift energy follows directly from its orthogonal companion decomposition, not from the selected physical row. The two-query gradient, its coisometry, and the seven-query joint E/E0 count are valid. The fixed-angle rotation leaves the exact affine term kappa q T³Z unchanged, so the corresponding localized four-slab estimate is also valid. The unchanged-tuple baseline coupling costs the stated VALUE remainder, without an unearned kappa factor.

## 2. Proof of the retained-field estimate without a dimension trace

Let P be the fixed physical S coisometry. At one literal Gaussian clock write h=P*H, r=P*R and let J_h,J_r be their conditional heat Jacobians. Set

    C_h=J_h−J_h*, C_r=J_r−J_r*.

The physical orientation is the P readout of −Sym(J C). Since E=H+R,

    O(E)−O(H)
       =−P Sym[J_r C_r+J_h C_r+J_r C_h] P*.          (A)

The source contracts give ||J_h||op<=C kappa, ||C_h||op<=C kappa and ||C_r||op<=C kappa. The last bound uses the original E curl bound and H's small full first, not R's larger O(A) first. Also ||C_r||HS<=2||J_r||HS. Thus every term in (A) satisfies a pointwise Hilbert-Schmidt bound C kappa ||J_r||HS. No raw Hilbert-Schmidt bound for J_h is used.

Conditional Gaussian integration by parts and first-chaos Bessel give

    ||J_r||_(L2 coarse roots;HS)
        <=||R−ER||2/v<=||R||2/v.

This is the full Gaussian source law after integrating the retained coarse roots. It need not be uniform in a fixed root. Summing the literal positive nodes by the triangle inequality and using sum w/v<=Lambda proves

    ||O_J(E)−O_J(H)||_(L2 roots;HS)
        <=C Lambda kappa ||R||2.

It applies with a shared coarse root or the actual supplied root tuple, and never replaces the field by its expectation. It therefore proves equation (15) as stated. For p=2, kappa=A², e_actual=Theta(A^3.9), and q=A^1.9, the resulting field bound has nominal exponent 6.85, versus 5.9 for kappa e_actual.

## 3. Positive reference and retention qualifications

If both actual conditional Gaussian references have an already supplied common positive spectral gap, their same-root W2 distance is bounded by a numerical gap factor times the L2 Hilbert-Schmidt covariance difference. Positivity and the upstream conditional Gaussian source contract remain premises. The inequality alone does not manufacture an endpoint Gaussian law.

Likewise, the pathwise source comparison retains the whole original nonlinear baseline and all original W,V exactly. The later pair comparison may retain only the roots and passive input actually named by its theorem. It cannot be enlarged to keep source-fine/private records which that theorem integrates. The source note states this distinction correctly.

The companion projected-gradient construction in this directory gives an exact all-root positive skew selector for the H orientation, with its finite pair-law and absolute numerical floors charged. Combining these two results establishes a genuine fixture-specific local orientation improvement with fixed-rank extra query exponent zero under the supplied polynomial-logarithmic clock/compiler gates. It does not prove a recursive gain for the returned non-gradient remainder, or a generic all-rank theorem.
