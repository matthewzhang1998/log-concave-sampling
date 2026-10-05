# Independent audit: multilevel positive conditional bridge covariance

2026-10-05. Author files were not modified. The audit concerns the exact pinned snapshot in this directory, not an unpinned later revision.

## Verdict

**PASS, with the stated imported native-action guards, for the positive public-log conditional-covariance approximation and its actual original-VALUE action ports.** This is a constructive theorem, not a completed arbitrary-precision production implementation. The conclusion is stronger than an identity audit: it includes the literal anchored producer, owned pair and coarse roots, caller derivatives, complete action-private integration, one-HS dimension factor, level/node counts, and a polynomial-bit scalar setup argument.

The audited author file is `DYADIC-POSITIVE-BRIDGE-COVARIANCE.md`, SHA-256 `9bc692e910de080b46690bd698e460e386ce6ee1a4de6c6c7d39cddb293ff40e`. Its exact copy is `AUTHOR-SNAPSHOT.md` here.

This passes the missing covariance-clock port. It does **not** pass the remaining complete prefix mean-LAW/caller join, the final own-mean compiler, arbitrary-order recurrence, or an application at non-Gaussian deterministic retained endpoints. It does not assert that a realized action is a PSD matrix. The PSD object is its analytical target covariance; the executing object is a vector-valued original-VALUE graph and its reserve is a positive Gaussian-root pushforward.

## 1. Exact target, conditioning, and dyadic identity

The target is

    Sigma(X,Y) = Cov(integral_0^delta exp(-u) g(X_u) du | X_0=X, X_delta=Y),

for the stationary OU process with covariance exp(-|u-v|) I and noise coefficient sqrt(2). The original retained pair has the joint Gaussian law

    X ~ N(0,I),  Y = exp(-delta) X + sqrt(1-exp(-2delta)) Z,

where Z is independent standard Gaussian. Every covariance-restoration norm is L2 of this exact pair with HS matrix norm. This law must be used before any later caller displacement. In contrast, the action calibration, energy and derivative bounds below are uniform in the retained values.

At interval length Delta, its midpoint conditional on local endpoints (a,b) is

    M = (a+b)/(2 cosh(Delta/2)) + sqrt(tanh(Delta/2)) xi.

After revealing all midpoints at one level, the residual bridges in disjoint child intervals are conditionally independent. The conditional expectation increments from distinct parent intervals have zero cross covariance given the previous skeleton. Consequently total covariance gives the displayed positive level sum and positive residual exactly; no covariance of correlated increments was dropped.

For each bottom bridge, the full isonormal private first is bounded by

    A integral exp(-u) sqrt(Var(X_u | endpoints)) du <= C A Delta^(3/2).

Conditional Gaussian Poincare gives its covariance <= C A² Delta³ I, uniformly in its endpoints. At depth K the sum over 2^K intervals, including exp(-2jDelta) weights, is <= C A² delta³ 4^(-K) I. Thus K=ceil(log_4(C/epsilon)) suffices. This exponential skeleton exists in the proof only.

## 2. Fixed-span complex contraction and position compression

In the author’s standardized local/retained coordinates the 2 by 2 matrix K_z is correct. Independent whitening into symmetric/antisymmetric endpoint coordinates gives an orthogonally equivalent matrix. Its determinant squared is independent of position, and the complex trace margin is proportional to

    sinh(x) sinh(L-x) - exp(-(delta+Delta)) sin²(y).

This agrees with the author’s determinant formula for I-K_z K_z*. The lower-right diagonal entry is positive for Delta<delta and 0<=x<=L. The wedge y²<=x(L-x) therefore implies I-K_z K_z*>=0. Tensor second quantization yields a dimension-free L2 operator contraction, including matrix-valued HS functions. No complex-density determinant estimate is being smuggled into the norm.

The endpoint atoms are retained exactly. Dyadic integer-index panels have width at most their distance to either global boundary. Their parameter-2 Bernstein ellipses remain in the contraction wedge. An m-node Gauss rule for each positive discrete panel measure is exact on degree 2m-1, has positive weights and correct mass. Polynomial approximation of the bounded operator-valued analytic function therefore gives error C 4^(-m) times the panel mass. This argument works for a discrete measure just as for Lebesgue measure; it requires no smoothness of that measure.

At level k its full mass is at most 2^k. Multiplication by ||C_Delta,Q||_(L2;HS)<=C A² Delta³ sqrt(D) gives C epsilon A² delta Delta² sqrt(D). These errors sum geometrically over k; an extra factor K in accuracy is unnecessary. The k=0 case has Delta=delta, L=0 and exactly one retained pair; no complex operator or singular zero-span whitening is used there.

## 3. Scalar setup really is polynomial in public logarithms

The author’s geometric generating function produces all panel moments through degree 2m using finite Taylor-jet multiplication/division. A separate independent recurrence is

    S_0=(1-r^M)/(1-r),
    S_p=[r sum_(l=0)^(p-1) binom(p,l) S_l - M^p r^M]/(1-r),

where S_p=sum_(j=0)^(M-1) j^p r^j and r=exp(-2Delta). An index shift and affine rescaling give the actual panel moments. This is O(m²) scalar arithmetic, plus polynomial-cost exponentials and integer powers. Neither formula visits M lattice points. Log M<=O(k), so representing the index endpoints is inexpensive.

Here is a quantitative precision argument supplying the detail behind the draft’s Vandermonde sentence. Normalize a panel to [0,1], normalize its mass to one, and let M be its atom count. Since delta<=log 2, all normalized atom weights are at least c/M. If M<=m, retain those atoms exactly. Otherwise choose any m distinct support points. Their separation is at least 1/(M-1). For the m by m monomial evaluation matrix V,

    |det V| >= (M-1)^(-m(m-1)/2),  ||V||op<=m.

The moment Hankel matrix H_m dominates (c/M) VV^T. Consequently

    lambda_min(H_m) >= exp[-C(m² log M + m log m)].

The same estimate applies to all smaller principal matrices used in a Stieltjes/Jacobi construction. Their arithmetic and square roots can therefore be certified with polynomially many bits in m,k and the requested output precision. Cancellation in the moment recurrence costs at most polynomially many extra bits, since every division by 1-r loses only O(log(1/delta)+k) bits and there are O(m) such stages. Taylor jets give the same conclusion. Unit-panel normalization prevents powers of delta from becoming unreported arithmetic-operation counts.

For completeness, root separation and positive weights do not require exponentially many scalar operations. The exact Gaussian rule satisfies H_m=V_G diag(w) V_G^T, sum w=1. Thus each weight is >=lambda_min(H_m)/m. If two Gauss nodes were separated by h, two columns of V_G would differ by at most C m^(3/2) h, giving lambda_min(H_m)<=C m³ h². Hence distinct Gauss nodes are separated by at least sqrt(lambda_min(H_m)/(C m³)). Both are exponential-in-polynomial lower bounds, which require polynomial bit precision. Certified polynomial root isolation or symmetric-Jacobi eigensolving is therefore polynomial in these parameters and the requested floor.

Nodes can be encoded within their real panel; endpoint atoms are already exact. Positive rounded weights may be renormalized to the recorded panel mass. The moment/quadrature error, root/weight encodings, pair-root encodings, and any renormalization error receive separate absolute floors. Structural positivity and actual coefficient envelopes, rather than a nonexistent Hessian modulus, preserve the first bounds of the encoded graph.

There are O(km) position nodes at level k and O(K²m) in total. With K,m=O(log(1/epsilon)), this is O(log³(1/epsilon)). Inner bridge and Gaussian covariance clocks each use O(log²(1/epsilon)) nodes. The unexpanded square filter contributes its own fixed-depth logarithmic count. Counting every original source replay gives a conservative log-eighth pattern. Scalar coefficient computation has polynomial bit/operation complexity in these logs, log(1/delta), and the assigned absolute precision. No inverse-Delta or 2^K enumeration is hidden in this claim.

This complexity proof is not a claim that a certified production moment/root solver has been delivered. The author expressly discloses that implementation boundary.

## 4. Local fixed field, one-HS error, and literal producer

For the midpoint field m_Delta, let sigma²=tanh(Delta/2). The derivative is

    J_Delta(xi)=sigma integral exp(-u) psi(u) E[H(ell_u+sigma psi(u)xi+kappa_u N)] du,

where H=Dg is analytical notation. The exact scalar integral is

    c_Delta=(Delta/2) exp(-Delta/2).

Hence J is symmetric PSD and bounded by ell I, ell=A sigma c_Delta<=C A Delta^(3/2). Positive half-bridge quadrature, with strictly interior replacements for short endpoint strips, preserves these properties and gives

    ||J_Q-J||_(L2(a,b,xi);HS)<=C epsilon ell sqrt(D),
    sum gamma_i/kappa_i <= C sqrt(Delta).

The endpoint law for this quadrature is the stationary triple (a,M,b). It is exactly the image of the stationary local pair and independent xi under the midpoint formula. Conditioning further on the unused third endpoint does not alter a within-half bridge given that half’s endpoints. Thus the imported one-time bridge theorem is applied at its correct law, rather than being asserted uniformly in arbitrary endpoints.

The Gaussian covariance identity

    Cov_xi(m_Delta)=2 integral_0^1 r E[(P_r J_Delta)²] dr

is correct. For fixed local endpoints, every chaos component of the symmetric matrix field J_Q is symmetric. Its matrix energies A_n=E[J_n²] are PSD and sum to E[J_Q²]<=ell² I. A uniform scalar multiplier error epsilon therefore yields a Loewner error between -epsilon ell²I and +epsilon ell²I, hence HS error epsilon ell² sqrt(D). This avoids a product of two dimension-sized energies. The inner J approximation is also one-HS: B²-C²=B(B-C)+(B-C)C, with both other factors bounded in operator norm by ell.

The actual source at fixed r,xi,a,b is exactly

    F(N)=sigma sum_i gamma_i/tau_i [g(z_i+tau_i N)-g(z_i)],
    tau_i²=kappa_i²+sigma² psi_i²(1-r²),
    z_i=ell_i(a,b)+r sigma psi_i xi.

It is a genuine D-dimensional gradient in N, has literal zero at N=0, and has no expectation or HVP-valued executing leaf. Its derivative is sigma sum gamma_i H(z_i+tau_iN), so its radius is ell and its pointwise norm <=ell|N|. A common N across inner nodes is legal because only their one-time expectations identify E DF=P_r J_Q.

The caller derivative contains the difference of the shifted and anchor Hessians. Since both lie in [0,A I], that difference has operator norm <=A. Therefore tau_i>=kappa_i and the displayed coefficient sum give

    first_(a,b) F <= C A Delta,
    first_xi F <= C A Delta^(3/2),

uniformly in all captured values and N. These are actual graph derivatives, not derivatives of an integrated clock error or LAW approximation.

## 5. Imported action, exact pair roots, and complete assembly

The original LOW30 `t30:lem:actions`, including equations `t30:eq:first-response` and `t30:eq:square-value`, was inspected directly. Its pinned source hash is listed in INPUT-PINS.json. Its zero-clock action is pointwise odd in p, targets (s0 E DF)²p, and supplies the four displayed calibration/energy/private/captured bounds. The prior independently audited square-clock specialization agrees with that statement. This audit imports that native theorem and its actual finite guards; it does not infer the theorem from numerical tests.

Applying it to the literal F gives local energy/calibration scales A²Delta³[mu]sqrt(D), complete private/p first A²Delta³/sqrt(mu), and local endpoint first A²Delta^(5/2)/sqrt(mu), with declared native logarithmic factors. The owned xi path uses its own captured-first bound and is absorbed into the full private tape.

At each compressed position, the actual local pair is a *joint* bridge draw. One explicit legal implementation is first draw X_t conditional on (X,Y), then draw X_(t+Delta) conditional on (X_t,Y):

    b = [sinh(delta-t-Delta)/sinh(delta-t)] a
        + [sinh(Delta)/sinh(delta-t)] Y + S G2,
    S²=(1-exp(-2Delta))(1-exp(-2(delta-t-Delta)))
        /(1-exp(-2(delta-t))).

The first site uses its usual bridge row and an independent G1. At t=0 reuse X; at t+Delta=delta reuse Y. At level zero both endpoints are reused and no root is manufactured. Expanding the two rows shows full private pair norm <=C sqrt(delta) and original retained-pair coefficient norm <=C. This is not independent marginal sampling of a and b.

The source receives independent COMPLETE pair/coarse/native banks for distinct labels and one COMMON incoming p. All pair and xi roots are integrated before claiming its conditional mean matrix. The sums

    sum_k (delta/Delta_k) Delta_k³ <= C delta³,
    sum_k (delta/Delta_k) Delta_k^(5/2) <= C delta^(5/2),
    sqrt(delta) sum_k (delta/Delta_k) Delta_k^(5/2) <= C delta³

pay, respectively, direct private/p, original caller, and owned-pair-root paths. Weighted triangle inequalities are valid for both disjoint banks and the shared p, and give exactly the claimed four action ports. They introduce neither sqrt(number of nodes) nor sqrt(full private dimension).

No unknown matrix square root is executed. All covariance fields remain analytical certificates. The basic six D-roots per bank (two pair, one xi, three native) and shared p are correct; optional bottom actions need no xi. Every native replay expands each F occurrence into original shifted VALUES and same-caller anchors. The conservative count 2 N_u (2K_f+4) per basic uncached bank is correct. Same-caller capture can reduce it only under the identical caller/source/version key. All changing endpoints, pair roots, xi and action inputs require complete affected replay and fresh anchors. Recorded HVPs may be used for requested first/adjoint sweeps, not as producers.

## 6. Bottom proxy, reserve ledger, and stopping boundary

The closed-form k_Delta in the optional bottom proxy agrees with the exact scalar conditional variance, including k_Delta=Delta³/6+O(Delta⁴). Its source has radius A sqrt(k_Delta) and local endpoint first C A Delta. Thus it has the same weighted ports, and replacing the small positive residual by this small positive proxy changes the error only by C A²delta³4^(-K)sqrt(D).

For any symmetric PSD matrix B and g(x)=Bx, the midpoint field equals sigma²c_Delta²B² and the proxy equals k_Delta B². Every finite inner/position clock preserves the relevant constant and exact mass. The full covariance is therefore restored, not just a scalar trace or leading asymptotic coefficient. Numerical encodings and native floors remain separate.

For eta=zeta=h/2, c=r0²/(2eta), the centered complete-private action comparison costs

    C c² Lip_private(D_Q) energy(D_Q)/zeta
      <= Lambda A⁴delta⁶ sqrt(D)/(h³sqrt(mu)).

Action calibration costs Lambda A²delta³mu sqrt(D)/h. The Gaussian reference has the paid additional covariance r0⁴ Sigma_Q²/(4eta²), costing C A⁴delta⁶ sqrt(D)/h³. Clock restoration costs C epsilon A²delta³sqrt(D)/h. These reproduce equation (24), with all inverse-h factors present.

At delta comparable to A^(4/5), h=A^(9/5), mu=A^(2/5), multiplying by the terminal A gives grades 4,21/5,22/5. Choosing epsilon<=A^(2/5) puts the clock bill at grade >=4. Reserve residual private first is of grade 12/5 and original endpoint first grade 2. These are actual graph bounds. The old full h carrier must be removed if the two h²/2 services are used; no variance may be spent twice.

This does not license extra observers of integrated private roots, transport the covariance certificate to shifted endpoints without payment, or supply the other complete mean-LAW branch. The complete prefix/caller/final compiler assembly remains a separate audit. The covariance action by itself is not certified as a genuine-gradient source for arbitrary re-entry.

## 7. Independent executable diagnostics

`check_multilevel_covariance_independently.py` uses independent whitened coordinates and scalar formulas. It checks 2,000 complex-wedge instances, the exact affine conditional-variance decomposition through nine levels at four horizon scales, the tail’s 4^(-K) normalization, and degree-12 geometric moments for panels up to 2^24 atoms without enumerating them. The largest observed singular value is 0.9999999979007945. The affine decomposition errors are below 7e-73 at 70-digit arithmetic; the large-panel moment relative error is below 3e-59.

`check_literal_value_graph.py` independently executes the anchored original-VALUE primitive and native zero-clock signed response on two-dimensional noncommuting bounded-Hessian fixtures. It checks exact zero, private PSD/radius, finite-difference first, captured scale, pointwise action oddness, original VALUE counts, and full non-diagonal matrix affine calibration. Maximum relative derivative discrepancy is 1.19e-10; exact parity discrepancy is zero; matrix calibration discrepancy is below 2.6e-14. The uncached 48-inner-node, six-filter-node action uses exactly 1,536 original VALUES, matching the census.

These finite checks support the explicit proof and graph inspection. They do not establish universal calibration by Monte Carlo, execute a final mean-LAW compiler, or constitute the certified arbitrary-precision moment solver.
