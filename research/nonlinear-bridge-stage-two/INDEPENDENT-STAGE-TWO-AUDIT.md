# Independent audit of the bounded-block stage-two repair

2026-10-05. This is a new audit. No sealed first-repair file was edited.

## Verdict and scope

**PASS for the new mathematical ingredients:** the true conditional ordered-pair construction, its uniform covariance margin, positive finite polylogarithmic pair quadrature for retained-variable L2 error, the shifted order-four cancellation, the executed source ports and counts, and the stated absolute source floors.

This verifies **one finite order-four step** under a fixed orthogonal block decomposition with block dimensions at most b, symmetric 0 <= Dg <= A I, g(0)=0, Lip(Dg_block)<=BA, and continuous D²g_block with Lip(D²g_block)<=CA. Constants may depend on b,B,C. It is not an induction theorem, a general-C2 uniform rate, a large-block dimension-free theorem, or a strong approximation conditional on the private roots.

The final positive own-mean completion remains an import of the same pinned guarded native compilers used by the sealed first repair, at the newly declared dimensions and bounds. This audit verifies the new interfaces and the resulting order, not a new implementation or independent reproof of those compilers. All their guards and absolute floors remain conditions of the result.

Reviewed inputs:

- TWO-TIME-POSITIVE-QUADRATURE.md, including its explicit linear-growth specialization and moment-corrected tensor rule.
- SHIFTED-BIAS-AND-PORTS-DERIVATION.md, including Sections 8–9.
- certify_pair_covariance.py and its exact-arithmetic output.
- check_stage_two_source.py and its explicitly limited diagnostics.
- The sealed first-repair theorem and bridge proof. Its SHA256SUMS rechecked successfully during this audit.

## 1. Genuine joint ancestry and covariance certificate

I independently derived the two Markov-innovation regression rows:

    R1 = 2a/sqrt(exp(2a)-1) * (1,a-1),
    R2 = exp(-a) 2h/sqrt(exp(2h)-1) * (1,2a+h-1).

These follow by subtracting exp(-h) times the first-time OU row from the second-time row before normalization. Consequently K=I-RR^T is the actual conditional innovation covariance given x,N,M. Its scalar source-independent square root reproduces the true joint pair, not two independent one-time bridges.

The rational program uses valid positive Taylor-tail enclosures for exp, monotonic enclosures for (exp(2s)-1)/(2s), and interval Sylvester tests for (9/10)I-R^TR. I reran it: all 2,792 boxes pass, maximum depth 19. The script's final numerical rational inequalities alone do not prove its tail formulas; I checked those formulas separately:

- For s>=8, the first-row squared norm is at most 8s^4 exp(-2s)<1/64.
- F(h)<=2 and h²F(h)<=6 follow from the quadratic and fourth positive exponential-series terms.
- For a>=8, the second-row squared norm is at most exp(-2a)(16a²+14)<=1038(2/5)^16<1/64.
- For h>=8, exp(-2a)[1+(2a+h-1)²] is decreasing in a: its derivative is exactly -2exp(-2a)(2a+h-2)².

The compact and tail estimates give K>=(1/16)I globally. The symmetric 2-by-2 square-root formula has a uniformly positive denominator. All five-input scalar history rows have norm one.

Using one common pair of residual roots across finite nodes does not create a multi-time OU process. The proof never uses such a process. It uses correct pair laws within each node and linearity of expectation across nodes.

## 2. Retained-variable L2 and nonanalytic source values

The main possible failure here was to confuse a joint L2 density pairing with an L2 function of the retained x,N,M. The written proof does not do that. It establishes a pointwise conditional-density bound before taking the retained-variable L2 norm.

I checked its explicit chain at eta=2^-24:

    ||R(z,w)-R0|| <=384 eta,
    ||L0^(-1)L(z,w)-I|| <=24 eta,
    ||K(z,w)-K0|| <=1152 eta,
    ||L0^(-1)C(z,w)L0^(-T)-K0|| <=8192 eta,
    ||G-I|| <=1/128,
    |beta| <=2^-13 |Y|.

The larger-polydisk bounds ||R||<=24 are valid. In particular |f(z)|<=4 and |z f(z)|<=8 follow directly by bounding the exponential denominator with its positive quadratic and fourth terms. The mean normalization correctly retains both the x contribution and the N,M regression contribution; no caller is frozen.

The complex density bound

    2^d (2pi)^(-d) exp(-|r|²/8+8|beta|²)

follows from the stated covariance and mean bounds. The nonlinear function F_Y(u,v), including S(Y), is evaluated only at real arguments. Density domination gives a polynomial times exp(2^-23|Y|²), whose square is integrable against gamma_(3d). This proves the required Banach-valued holomorphy without analyticity of g.

The explicit linear-growth constant also checks:

    C_b=12*8^b*(1-2^-21)^(-(3b+2)/4).

The factor 8^d is the density-mass ratio against N(0,4I_(2d)); its first radial moment gives 2E|r|<=6sqrt(d). The exact tilted Gaussian second moment then supplies the displayed exponent. Dependence on b is allowed and is not hidden dimension independence.

## 3. Positive finite quadrature and both clock distributions

The relative complex disks contain the stated Bernstein ellipses on dyadic panels. Positive weighted Gaussian quadrature has exact panel mass, so head/tail atoms and the panels form an actual finite probability rule. Two marginal rules telescope to the pair bound in L2; the nonlinear conditional integrand need not be an operator on g alone.

The now-explicit constants are sufficient: with internal epsilon=delta/32, an uncorrected marginal costs at most epsilon M; exponential-moment repair has mixing weight alpha<=6epsilon and corrected error at most 13epsilon M; the corrected tensor error is at most 26epsilon M<delta M. For the normalized source class, use delta/C_b as the internal target.

The ancestry correction uses Exp1(a) times Exp1(h). The quadratic correction uses Exp2(a) times Exp1(h), obtained by splitting the independent Exp1 clock square into its two time orderings. Q is symmetric in its last two arguments. The diagonal has zero measure. Both clocks and both weights are necessary; replacing the quadratic distribution by Exp1 times Exp1 would be incorrect.

Marginal node counts are O_b(log²(1/delta)); pair counts are O_b(log⁴(1/delta)). The constants, especially the explicit safe holomorphy radius, are conservative and potentially very large. This is a finite-existence and fixed-grade original-VALUE exponent statement, not a practical small-node claim.

## 4. Shifted cancellation and quantitative remainder

The retained displacement identity is exactly x-F2=S+d+e. The terminal expansion is performed at S. I checked each remaining contribution:

- H[d,e]+H[e,e]/2 gives the B[d0 e0+(A/2)e0²] part.
- Taylor's third-order bound gives (C/6)u_A³.
- The one-future-sample substitution is a weak error. Its first-order term cancels only after the future-clock average, at the center X_a-F1,a. The variance identity gives (B/2)A³ before multiplication by Dg(S), hence (B/2)A⁴ afterward.
- Centered C cancels its quadratic Taylor term, leaving (CA/6)|a|³. This applies at every actual positive single-rule node.
- The /8 polarized stencil is exactly (1/8) times the rectangle integral of D²g. Its leading term is (1/2)D²g(S)[a,b], and its remainder is at most (CA/4)|a||b|(|a|+|b|).

These yield exactly the displayed K_(b,B,C,A). L4 and L6 Gaussian moments, conditional Jensen, and orthogonal block summation give sqrt(D), with the b-dependent factors shown. There is no commutation of source derivatives and no assumption A sqrt(D)<<1.

The sealed one-time operator rule is used only for the linear part Dg(S)d. Its nonlinear centered remainder is bounded nodewise. The pair theorem applies directly to the combined actual C and Q stencil values, so no derivative producer or hidden expectation leaf is executed.

The B-free pair envelopes are valid:

    H_ancestry=sqrt(11/8)A³,
    H_quadratic=A²/2.

Thus sufficient budgets are delta1<=A², delta_ancestry<=A, delta_Q<=A², and delta_out<=A³. Every quadrature term is then O(A⁴ sqrt(D)). The optional B-dependent quadratic envelope permits a looser quadratic pair tolerance but is not needed.

## 5. Actual source interfaces, counts, and completion

The source needs 5D private coordinates after including G. Its residual occurrence cost is

    (5+3J1+5JL+6JQ) N_out

original VALUES before complete-key sharing. The two pair supports generally differ; the reduced shared-support bill cannot be assumed without verifying identical nodes.

I rederived the noncommuting derivative products, including the difference Dg(U)-Dg(U-g(V)) in the U derivative of e. All factors stay in their executed order. The PSD interval [0,A I] yields the stated norms A/2 for the centered S coefficient and A/4 for each polarized coefficient. The leading caller/G-block coefficient is symmetric, so it does not create a spurious order-A curl.

The proposed bounds are safe. At A=1/2:

    normalized first²/A² <=547655/8192<81,
    normalized curl/A² <=(169sqrt(6)+sqrt(126874))/32<25.

Therefore ell_E=9A and a_seed=3A are legitimate conservative declarations, and A<=1/36 suffices for the single ell_E<=1/4 guard. It does not discharge any other guard.

I checked the exact-moment and moment-free caller-energy constants. The captured origin is the literal graph at zero private roots with the caller still live. Its subtraction doubles the raw caller first bound, and its full replay cost remains charged. There is no residual-energy normalization or source-dependent stopped parameter.

Source energy is O(A²)(|z|+sqrt(D)). With mu=A, the imported near-gradient bracket is O(A²), so its integrated allowance is O(A⁴ sqrt(D)), subject to the native guards and restored floors. Complete banks must be independent conditional on the caller and exterior labels. The claim concerns the completed source's own mean and the positive completed law; raw private roots are not added as observers.

The matrix-quadratic reduction and exact own-mean identity check algebraically. The extra ancestry stencil replaces the old K³w term by the weighted K³V term. The linear pair product moment 1/4 and the single moment 1/2 recover the exact canonical mean. This does not assert a raw covariance identity.

## 6. Absolute numerical errors

The stated conservative uniform raw F VALUE floor

    [7/2+14A+(11/2)A²] nu

is valid. Baseline subtraction and anchoring add their separately stated errors. The pair row bounds

    ancestry: A²(epsilon_U+A epsilon_V)sqrt(D),
    quadratic: (A²/4)(epsilon_U+epsilon_V)sqrt(D)

follow from the exact source derivatives at fixed S. Probability-weight floors also follow from the actual combined-stencil sizes. No floor is divided by A, residual energy, a conditional variance, or a small covariance eigenvalue.

Actual numerical rows, probability masses, moments, and positivity still require their declared checks. Exact real-arithmetic quadrature identities must not be silently transferred to rounded weights. Changes of caller, coefficients, numerical version, or discarded primal require the full original-VALUE/HVP replay ledger. No HVP is differentiated.

## 7. Independent executable checks

The independent_stage_two_checks.py script uses a separate covariance implementation and reports 5,107 passing assertions. It checks exact rational scalar port constants and tails, true pair covariance identities, normalized complex covariance/mean estimates on an independent sample, and the shifted cubic-jet identity with arbitrary noncommuting first derivatives and symmetric second-derivative tensors. Maximum cubic-identity discrepancy was about 1.01e-14.

These checks supplement the analytic proof. They do not implement the native mean compilers, certify the entire holomorphy domain by sampling, measure a strong pathwise rate, or turn a small diagnostic quadrature into the proved high-accuracy finite rule.

No substantive defect remains in the reviewed new claims. Any assembled final theorem must preserve the stated qualification, native-guard conditions, tolerances, complete bills, and absolute-floor scope.

## 8. Final integrated theorem and resolved wording

I also reviewed STAGE-TWO-POSITIVE-VALUE-THEOREM.md after assembly. Its main error bound uses the correctly normalized B-free pair errors

    sqrt(11/8) deltaL A³ + (deltaQ/2) A²,

with the internal rule tolerance delta/C_b. Its explicit marginal and pair node counts, guarded completion bracket 36A²+729A³+729A^(5/2), numerical floors, and full source bill agree with the audited component proofs.

Two presentation issues were identified and corrected before approval:

1. A fixed pair node has E[V_j|x]=exp(-(a_j+h_j))x. The zero extra quadratic-source bias uses only the weighted identity sum q_j E[V_j|x]=x/4. The final main statement now makes this distinction explicit.
2. For f(y)=s[g(a+s y)-g(a)] and alpha=A s², Df=s²Dg and D²f=s³D²g, so the structural constants are B_f=Bs and C_f=Cs². These calculations are valid, but they do not automatically certify new live caller/anchor/scale interfaces. The final statement now explicitly requires their recertification, the actual alpha guards, fixed numerical rule versions for a requested sweep, and the full original-g implementation bill.

Both issues are resolved. I reran the author raw-graph diagnostics as well: 8,170 assertions PASS, 55 VALUES per test residual, maximum Jacobian discrepancy 1.48e-10, and maximum quadratic-mean discrepancy 1.45e-16. The exact covariance certificate and the 5,107 independent checks also pass. The integrated theorem is approved within the scope stated at the start of this audit.
