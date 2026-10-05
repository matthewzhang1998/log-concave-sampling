# Independent bounded audit: m3 raw source and exact currents

Date: 2026-10-05. Reviewer: independent m3-current reviewer.

## Verdict

**PASS for the literal three-root raw-source ports/counts and its guarded OWN-mean compilation; PASS for the stated source-qualified exact analytical current identities. The target-bias comparison and native fourth-order m3 current consumer remain OPEN.**

The two separator results also PASS within their stated scopes, with a separately pinned independent review imported below. Neither separator refutes the actual same-g single-I covariance formula.

This is not a full order-four endpoint PASS, an all-order theorem, or a finite current producer. The exact unfilled target is still

    m3(Z) = R1 psi2(Z),
    psi2(x) = E[g(x-F2) | original Gaussian endpoint x].

A separately completed full Cov(F2|Z) service does not fill this mean target.

Two qualifications must accompany the component PASS:

1. The raw note's Section 4 imports the earlier nested-mean theorem. Its A^3 mean estimate requires that theorem's actual Hermite-multiplier clock tolerances, for example both inner and middle tolerances delta<=A^2. Mass one and first moment one half alone do not supply this estimate. Those weaker clock conditions do suffice for the raw ports and OWN-mean target in Section 2.
2. The resummed note's generic Gaussian Riesz statement should be read for Gaussian Sobolev functions with integrable first derivatives, not arbitrary square-integrable C1 functions. All the actual I_H/I_Q sources satisfy the stronger condition, since their complete first derivatives are bounded. Thus this wording qualification does not invalidate any source-qualified identity reviewed here.

The independent script passes 1,993 assertions. Its finite diagnostics supplement the analytical reasoning below; none estimates or proves the missing uniform A^4 trace remainder.

## 1. Frozen targets and imported review

The following SHA256 pins were checked before and after the audit:

- `../P3-MEAN-RAW-SOURCE-AND-SECOND-DECOUPLING-GATE.md`: `87b61a6a4ca03f81f7d910cdaa83b3783cd4792b3a428b89c8f0f629062e4467`.
- `../SINGLE-HISTORY-EXACT-CENTERED-RESUMMED-CURRENT.md`: `276d5873cc7d3f5ca8650955c08abe128cec7eda96680319f5c2a00b45118e9a`.
- `../COHERENT-SHIFT-SEPARATES-UNSHIFTED-M3-FEEDBACK.md`: `cb04bb4e1cd3886ed378bb45593ea855b800e1cd1e532916c5f238cd3141e285`.
- `../GRADIENT-PORTS-DO-NOT-CONTROL-A-CONTRACTED-TRACE.md`: `e6050915d148b7248e97bd26647f5948f14ef50d47dbd365d9c668b18790faa1`.
- Imported separator review, `../independent-mean-separator-audit/INDEPENDENT-ANALYTICAL-AUDIT.md`: `072f1672c7591305cb6d72fb78746ab5ee11e99733c011a40a4e0cbf0de73024`.

The single-history-current author is an independent single-history-current analyst; that author is independent of both separator constructions. I inspected both separator sources and the imported review, and independently checked the main algebra and scope. The present checker was written independently, without importing the author checkers.

Relevant mathematical imports were also inspected: the conditional-innovation localization and the third-order nested-force mean theorem with its independent review. Their finite mean law is an integrated Gaussian-carrier statement, not a strong conditional-force value oracle.

## 2. Literal raw three-root source

Assume the inherited small-radius setting, anchored g=grad U, U in C2, 0<=Dg<=A I. Fix all positive interior nodes, weights, versions, shares, and precision parameters before taking firsts. Denote

    beta_o = sum_i w_i c_i,
    beta_m = sum_j v_j d_j,
    beta_n = sum_k u_k e_k.

Each beta is at most one. The private tape is exactly (G,H,J), dimension 3D. G is shared across every outer occurrence; H across every middle occurrence; J across every innermost occurrence. The exposed Z is outside that tape. The finite raw program is not identified with a Markov path.

At an outer site x let y_j=tau_j x+d_j H and L_j=sum_k u_k g(sigma_k y_j+e_k J). Then

    ||D_y L_j|| <= A/2,   ||D_J L_j|| <= A beta_n,
    ||D_x I2|| <= a_x := A(1+A/2)/2,
    ||D_H I2|| <= A beta_m(1+A/2),
    ||D_J I2|| <= A^2 beta_n.

Every displayed derivative is the ordinary first of the full recorded VALUE graph. In particular D_x I2 need not be symmetric. No Hessian is differentiated.

For the outer chord, write B1=Dg(x-I2), B0=Dg(x). Then

    D_G E_i = c_i[(B1-B0)-B1 D_x I2],
    D_H E_i = -B1 D_H I2,
    D_J E_i = -B1 D_J I2,
    D_Z E_i = r_i[(B1-B0)-B1 D_x I2].

Because B1 and B0 both lie between zero and A I, their difference is symmetric and has norm at most A. Hence valid explicit aggregate upper bounds are

    ||D_G E|| <= beta_o(A+A a_x),
    ||D_H E|| <= A^2 beta_m(1+A/2),
    ||D_J E|| <= A^3 beta_n,
    ||D_Z E|| <= (A+A a_x)/2.

The full first is at most the sum of the first three bounds, and is O(A). For the square lift (E,0,0), a sufficient skew bound is

    ||D(P_G*E)-D(P_G*E)*||
       <= 2 beta_o A a_x
          + sqrt([A^2 beta_m(1+A/2)]^2+[A^3 beta_n]^2)
       = O(A^2).

The first Hessian difference contributes no G-G skew. This is the complete curl calculation, including the off-diagonal H and J blocks. It does not infer a quadratic first from a quadratic VALUE amplitude.

The checker uses genuinely noncommuting positive Hessians and compares all four complete derivatives against finite differences. Its largest directional discrepancy is below 1.9e-11. A separate scalar oscillatory family with bounded positive Hessian has raw chord first bounded below by a fixed multiple of A while its VALUE is O(A^2). Thus an O(A^2) first is not available from these ports.

### Energy and origins

Positive weights, the Lipschitz chord bound, and the actual Gaussian row norms give

    ||E(Z,G,H,J)||_(Lp private)
       <= C_p A^2(|Z|+sqrt(D)).

This uses one physical vector energy. Each z_ijk is a normalized affine Gaussian row in the independent variables Z,G,H,J; conditional means and variances are bounded by the displayed caller profile. No sum over the number of roots or nodes appears in the energy constant.

B is a genuine gradient in G with first at most A beta_o. An analytical potential is sum_i(w_i/c_i)U(r_i Z+c_i G); it is never executed. At fixed Z the two private origins are generally nonzero. Useful bounds are

    |B(Z,0)| <= A|Z|/2,
    |E(Z,0,0,0)| <= A^2(1+A/2)|Z|/4.

They must be executed, captured under the same complete key/version, subtracted, and restored. Subtraction preserves private first/curl and gives the same energy profile up to a constant. The origin's actual captured-Z first is O(A); subtracting it does not improve that order. At Z=G=H=J=0, the exact zero follows recursively from the same-site g(0)=0 anchor. It is not a Gaussian expectation argument.

### Counts, native evaluation, and floors

A raw baseline occurrence uses N_out original g VALUES. A raw chord occurrence uses

    N_out [N_mid(N_in+1)+2]

VALUES: N_in innermost sites and one middle terminal per middle node, followed by one outer terminal and one baseline per outer node. The safe total bill is therefore exactly the source's symbolic bill

    Q_captured + N_B Q_B + N_E Q_E
       + known/numerical/replay work.

The completed B branch has its own D-dimensional raw bank; the completed E branch has its own 3D-dimensional raw bank. The complete compiler dimensions also include their actual filters/clocks/replays. This is not a license to count every raw occurrence as one original call, or to keep a private compiler root after it has been integrated.

Changing any raw argument replays every affected original-g ancestor. First/adjoint sweeps refer only to original HVPs at the recorded VALUE sites. A stored HVP is neither a producer nor a differentiable leaf. The checker explicitly verifies that every first site has an exact recorded primal.

If every original g VALUE has deterministic norm error at most epsilon, a simple safe one-occurrence propagation is

    raw F3 floor <= (1+A+A^2)epsilon,
    raw B floor <= epsilon,
    raw E floor <= (2+A+A^2)epsilon.

Origins and restoration add their own absolute contributions. These are absolute bounds, independent of realized E energy. The subsequent mean compiler must still propagate them through its actual coefficients, callers, modes, and saved versions.

### Guarded own-mean compilation only

The earlier independently admitted gradient/near-gradient mean theorem applies after origin subtraction, with fixed positive shares, actual normalized radii, actual curl, padding mu=A, gradient order at least four, and every inherited active-dimension/filter/clock/precision/caller guard. For example the explicit first bounds above, divided by the branch's square-root variance share, give a sufficient source-radius check. A<=1/2 by itself is not that check.

Under those guards, the source's completed own-mean target and A^4 error are valid. The node counts can be logarithmic at this bounded grade when chosen by the admitted positive compression theorem; arbitrary frozen rules do not acquire that cost property merely from their moments.

Nothing here compares that own mean to m3 at A^4. Section 4 equation (9) is explicitly UNPROVED and receives no target-bias PASS.

## 3. Innovation current: original history and correct transpose

Condition on the original Gaussian endpoint x. The earlier localization defines, on the same future OU history,

    F2 = H_j+S,     E[S|x]=0,
    S = integral L_s dB_s,
    ||L_s|| <= A^2 s exp(-s)/sqrt(2).

Future conditional centering, the conditional Clark-Ocone formula, and stochastic Fubini justify this adapted representation using only first derivatives of g. In particular the covariance constant is

    integral_0^infinity A^4 s^2 exp(-2s)/2 ds = A^4/8.

This is a martingale-integrand bound, not a small full/caller first for S.

For F_theta=H_j+theta S, the ordinary first-order fundamental theorem gives raw equation (5) exactly:

    m3-R1 psi_j
       = -R1 integral_0^1 E[Dg(x-F_theta)S|x] dtheta.

Bounded Dg and finite conditional moments justify the identity under C2. Its direct norm bound is only O(A^3 sqrt(D)).

The complete Malliavin derivatives satisfy

    ||D_s H_j|| <= Lip(j) exp(-s)/sqrt(2),
    ||D_s F2|| <= [A+A^2(s+1/2)]exp(-s)/sqrt(2),
    Lip(j) <= A(1+A/2).

For s<t and s>=t, respectively, the derivative of the inner I_t has the common bound A exp(-|s-t|)/sqrt(2). Integrating the latter bound against exp(-t) gives (s+1/2)exp(-s), establishing the claimed ancestor contribution without differentiating a Hessian.

For a smooth regularization, adapted stochastic integration by parts gives, componentwise,

    E[partial_j g_i(x-F_theta) S_j]
       = -E[sum_(j,l,k) partial_j partial_l g_i(x-F_theta)
                             (L_s)_(j,k)(D_s F_theta)_(l,k)] ds.

The coefficient is therefore L_s(D_s F_theta)*, in exactly the source's order. This proves the positive sign and transpose in raw equation (8). No derivative of L_s is introduced by this step.

A useful explicit deterministic bound is

    ||T_theta||op <= A^3/8+3A^4/16,
    ||T_theta||HS <= sqrt(D)(A^3/8+3A^4/16).

The scalar integrals and the current sign were independently checked. This is a coefficient bound only. It supplies no norm estimate for its contraction with D2g. The identity under C2 is defined by equation (5), not by convergence of individual higher-derivative factors. The original endpoint and every live coefficient ancestor must survive any later score transfer.

## 4. Centered/resummed single-history identities

### Gaussian Riesz orientation and covariance

For F in Gaussian W^(1,2), let B_rho=rho B+sqrt(1-rho^2)B'. The Riesz representation gives

    D(-L)^(-1)(F-EF)(B)
       = integral_0^1 E_[B'] DF(B_rho) d rho.

Pairing its ith row with D phi(F(B)) proves

    tau_F = integral_0^1 E_[B'] DF(B_rho) DF(B)* d rho,
    E[(F_i-EF_i)phi(F)] = E[sum_j(tau_F)_ij partial_j phi(F)],
    E tau_F = Cov(F).

The product order is correct. The generic reversed order fails: for F=(X,X^2-1), tau=[[1,2X],[X,2X^2]]. With phi(F)=F1 F2, the first-row correct expectation is 2 and the transposed one is 1. The checker detects this distinction. The D2g contraction itself sees only the symmetric part, but the underlying Stein identity must still use the correct orientation.

The actual I_H and I_Q functions have bounded full Gaussian firsts, so Gaussian Sobolev closure is available. For I_H the product of Brownian rows is c(t,u)I, yielding J_t^rho J_u in the source's order. For I_Q the corresponding product is d_i d_j I. Thus equations (4)-(5) have the right genealogy, auxiliary rotations, and indices.

Positive kernel mass gives ||tau_H||op<=A^2/4 and ||tau_Q||op<=A^2 beta_Q^2. The conditional covariance kernel has integral mass 1/4; it is not the unconditioned stationary covariance. Physical I and centered xi energies are bounded by A sqrt(D). These estimates remain independent of analytical approximation dimension.

### Exact convolution homotopy and C2 meaning

At each fixed x use two independent complete centered banks, but keep each tau_a correlated with its own xi_a. Differentiating

    E g(x-mu_theta-sqrt(theta)xi_H-sqrt(1-theta)xi_Q)

gives the stated derivative-free equation (7). Applying the preceding Stein identity within one bank, conditional on the other bank, gives the plus tau_H and minus tau_Q terms of equation (1). The covariance expectation cannot be substituted for the random Stein coefficient inside that contraction.

The apparent 1/sqrt(theta) and 1/sqrt(1-theta) endpoint singularities are integrable against the source's L2 energies. Anchored convex mollification preserves the uniform Hessian bound, and the regularized g,Dg converge locally with uniform linear-growth majorants. Therefore the derivative-free form, its integral, and endpoint difference converge in L2. This proves the whole-current C2 interpretation claimed, without supplying pointwise D2g or separately convergent trace factors.

After the L2-contractive R1, the coherent-mean term is at most A||Delta_mu||2, giving exactly delta A^2 sqrt(D) under the stated mean-quadrature guarantee. The old A^3 total estimate in equation (9) is a legitimate import of the already audited decoupling comparison. Neither estimate promotes the current to A^4.

The independent deterministic Gaussian diagnostics use a two-clock Markov bank versus a different one-clock cheap rule, so Delta_mu is genuinely nonzero. They check the complete random Stein coefficients, their covariance expectations, the derivative-free homotopy, and the combined current. The largest Stein-current identity discrepancy is below 9.7e-10.

### Same-clock covariance interpolation

For common clocks, the Gaussian covariance increment Delta_k has zero diagonal because both laws have variance 1-tau_i^2 at node i. Thus all terms differentiating one inner VALUE twice have zero coefficient. The remaining mixed-site derivatives give exactly

    (1/2) sum_(i,j) v_i v_j Delta_k(i,j)
        E[D2g(x-I_lambda):(J_i J_j)].

Individual J_i are symmetric, but different J_i need not commute. The ordered sum and symmetric D2g contraction are the legitimate reasons the expression is well defined. Interpolating E[I I*] yields the covariance equation (11), since the conditional mean is constant along this same-clock interpolation. All the finite-bank correlations remain inside these formulas.

The continuum kernel masses 1/4, pi^2/16, and (pi^2-4)/16 are correct. The scalar pointwise inequality k_M<=k_C is not a Loewner order of clock covariance matrices: even a two-node difference has zero diagonal and opposite-sign eigenvalues. No such Loewner inference is made by the source.

The finite interpolation diagnostics have error below 9.0e-11. Passing from finite analytical approximations to continuous clock integrals requires only the stated strong L2 limits and bounded-first majorants. It does not give a finite executed Markov quadrature with an order-four rate.

### Precise open boundary

The single-history current compares I_H to I_Q. Its identities do not themselves compare the raw two-level I2_Q to true F2, and do not close m3. The separated unshifted and centered leading terms in Section 4 are not silently granted separate C2 convergence at the frozen pin. The uniform one-energy A^4 trace estimate (15), its needed limits, and a finite original-VALUE consumer remain unproved by this artifact.

Any later fixed-D distributional existence supplement would narrow only a limiting-current issue; it would not imply a dimension-uniform A^4 bound. Such a supplement is outside these frozen four pins.

## 5. Separator imports and exact scope

The independent review at the pin in Section 1 validates both separator proofs. My algebra checks agree:

- In the coherent-shift family, the diagonal remainder in the projected small-feedback formula is O(A^3), while the rank-one term survives at A^2. Conditioning leaves the scalar projected OU path, and injectivity of R1 preserves a nonzero limiting function. This refutes replacing the shifted derivative by Dg(x) in the small-feedback m3-m2 formula. It does not refute the single-I covariance formula, whose relevant scaling is different.
- For the gradient-port separator, the uncut Hermite gradient has exact energy A^4(6D+12) and Laplacian mean A^2(2D+4)e1. The cutoffs preserve anchoring and O(A) first, with exponentially small Gaussian-tail error. Since R1 preserves expectation, W=A^2 I produces a current lower bound c A^4 D. This defeats a port-only O(A^4 sqrt(D)) trace bound.
- The second construction chooses E and W independently. If W came from the quadratic same-g source, its actual chord would not be that E. Therefore the result is port insufficiency, explicitly not a same-g m3 counterexample.

Neither a generic two-score/Hilbert bound nor gradient/curl/energy ports alone close the missing trace estimate. A successful theorem needs the actual shared genealogy, a genuine new cancellation estimate, and then a separately audited finite native consumer.

## 6. Reproducibility and final status

Run `python check_m3_gate_independent.py` in this directory. The script validates exact pins, literal raw counts and recorded first sites, complete derivatives and curl, actual nonzero origins, absolute VALUE floors, matrix-quadratic mean calibration, innovation constants/sign, Riesz orientation, unequal-mean centered convolution, same-clock interpolation, and the imported separators' Hermite algebra.

Its 1,993 assertions and deterministic numerical errors are recorded in `m3_gate_independent_checks.json`. Numerical integration is used only for supplementary fixed finite smooth examples. Continuum limits, uniform C2 scope, and the guarded mean-consumer import are justified analytically above.

Final gate status:

- Literal raw F3_Q graph, first/curl/energy/origin/count contracts: PASS.
- Conditional OWN-mean compiler: PASS subject to actual inherited guards and absolute floors.
- Raw exact first current (5): PASS.
- Raw regularized whole current (8): PASS with C2 interpretation through (5).
- Raw leading-covariance candidate (9): UNPROVED / OPEN.
- Centered Stein/convolution and covariance-interpolation identities: PASS for the actual Sobolev sources.
- Separator claims: PASS only in their explicitly bounded scopes.
- Uniform order-four m3 target-bias bound and native mean-current consumer: OPEN.
- Full order-four endpoint / all-order closure: NOT ESTABLISHED by this audit.
