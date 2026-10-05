# Independent audit: finite-depth shrinking-buffer reentry

2026-10-05. This is a mathematical/source-contract audit, not a new algorithm or a numerical execution of the native compiler.

## Verdict and pins

**PASS under the declared imported native guards, for the finite sequence of actual sources constructed from F_Q.** The corrected note proves a genuine finite-depth own-mean reentry recurrence. It does not improve the certified grade 16/5, return an unbuffered covariance-matched force, or establish an endpoint sampler. All public-log factors, complete source banks, retained callers, and new replay floors remain live at every depth.

Reviewed artifact: `../REUSABLE-PORT-AND-DEPTH-RECURRENCE.md`, SHA256 `1de572cbddecd40ed5c990df98638a5551d714f74cbc0d50dfd5200c19a3d58e`.

Verified pinned inputs:

- `../../law-only-variance-join-20261005/LAW-ONLY-SHRINKING-BUFFER-JOIN.md`: `40891db38221e3b02f2e657a4ef93c3e4fd06c69b29cc03e7db9c5761fe2fc0a`.
- Its `MANIFEST.json`: `673bdfcc6dc54c00ea380c942039add5ecda9f32686e6af0398243923dc6f3a7`.

I also checked the relevant literal hypotheses and physical-first conclusion of LOW30 `b27:raw:self-reserve` and `b27:raw:split`. Their remainder hypothesis is a **centered** Lp energy bound, so the centered-energy port in the reentry note has the right type. Caller origins are still evaluated, subtracted and restored as actual graphs; they do not disappear by centering.

## 1. Genuine targets and higher-depth restoration

Minkowski and the A-Lipschitz source give `||F_k||_2 <= A sqrt(D)+A||F_(k-1)||_2` and `||F_k-F_(k-1)||_2 <= A||F_(k-1)-F_(k-2)||_2`. Starting from `||F1||_2<=A sqrt(D)` proves every displayed contraction estimate, uniformly in finite k for A<=1/2. The same complete Gaussian-bank chain rule gives `L_k<=A(1+L_(k-1))`. A small target bias does not make this first smaller.

For centered Delta at fixed y, write its Gaussian Riesz derivative as `R Delta=D(-L_OU)^(-1)Delta`. Its L2 Hilbert-Schmidt norm is at most `||Delta||_2`. The conditional covariance identity therefore gives

    ||Cov(Delta,F_j|y)||HS
      <= ||R Delta||_(2,HS) ||D F_j||_(infinity,op)
      <= L_j ||Delta||_2.

This is dimension-free. Expanding the covariance difference into two cross terms and integrating y proves the O(A^4 sqrt(D)) restoration. No first bound on Delta is used. The positive Gaussian gap then prices the restoration at O(A^5/sqrt(u) sqrt(D)), safely below the displayed A^5/u^(3/2) service term for u<=1.

The added third-cumulant lemma is also valid. For centered a,b, Gaussian Poincare and their covariance operator bounds imply `Var(a^T M b)<=4 L_a^2 L_b^2 ||M||HS^2`. Duality gives the asserted centered-scalar-L2 to matrix-HS operator bound; applying it to each component of Delta and summing squares introduces no dimension factor. The three-slot centered tensor telescope gives O(A^5 sqrt(D)) for F_k versus F2, and O(A^4 sqrt(D)) for F2 versus F1. These are target restorations only; no skew producer follows from them.

## 2. Actual finite-depth input/output typing

The input is the literal split graph S_k=B_k+E_k, not a completed mean LAW and not a newly anchored convex gradient. The native split theorem accepts the gradient baseline and near-gradient remainder at their actual first, curl, centered-energy and caller constants. A signed baseline remains a genuine gradient, so the bulk input `-q S_k/sqrt(u)` is admissible under the declared guards.

At each bulk node the three **whole** service outputs have independent complete banks conditional on the same retained y. Retaining z adds no private observer. Their references add to the intended mean and covariance, with the true F_k covariance difference charged as above. The original F_Q near branch continues to have numerical O(A) complete-bank first and O(A^3 sqrt(D)) mean discrepancy from every m_k. It therefore avoids propagating Lambda_k into the near-RAW weak bill.

The short-prefix bound is uniform in depth because `||F_(k-1)||_2=O(A sqrt(D))`. The replacement by w g(X0) costs A w^(3/2)+A^2 w before the terminal g; the coherent move into `g(u-wg(u))` costs A^3 w afterward.

The output ports follow from the actual aligned graph. Its preterminal carrier-subtracted displacement has full first Lambda A; source-zero anchors and the inherited caller paths give its one-energy profile. The terminal g difference consequently has energy Lambda A^2. The possibly O(A) Hessian difference on the common carrier multiplies a scalar identity and is symmetric; nonsymmetric/off-carrier paths contain the small displacement first, giving curl Lambda A^2. This is the pinned direct graph proof applied inductively. No derivative is inferred from a LAW comparison. The theorem should continue to be read as this inherited, anchored finite construction, not permission to omit origin/path data from an arbitrary black-box mean estimator.

## 3. Error recurrence and certified ceiling

The recurrence's terms are correctly charged:

- Bulk mean mismatch: at most A e_k.
- Prefix: A^2 w^(3/2)+A^3 w.
- Fresh F_Q endpoint: A^3 sqrt(w), plus its O(A^3) mean mismatch over mass O(w), giving A^4 w.
- Bulk true-force Gaussianization: A^4/w.
- Finite fixed-target outer quadrature: A^4.
- Variable mean/Gram and mixed-K services: Lambda[A^5/w^(3/2)+A^6/w^(7/6)].

The exact positive rule is applied to the fixed true delayed target before its nodewise replacements. No t-dependent surrogate is silently promoted to one OU semigroup source. Additional covariance and mixture restorations are dominated by the displayed terms.

For w=A^alpha, the two unavoidable ledger lines `2+3alpha/2` and `4-alpha` already bound the best minimum by 16/5, uniquely at alpha=4/5. Every other displayed line is larger there. The affine recurrence and its finite geometric expansion are correct. Also `||m_infinity-m3||_2<=A^4/(1-A)sqrt(D)`, so S3 already has the same certified limiting-mean grade. This is an upper-error-ledger ceiling, not an error lower bound.

The prospective 7/2 ledger after replacing the Gaussianization term by A^5/w^(3/2), and the prospective 10/3 prefix-only ledger, are arithmetically correct. Neither is a constructed improved producer.

## 4. Same-source RAW covariance separator

The shared outer carrier is essential. In the linear example it gives coefficient `b_Q=sum omega_j sqrt(1-t_j^2)`. The positive degree-two certificate implies `b_Q>=2/3-epsilon>=7/12`. The centered residual energy bounds the cross covariance by O(Lambda_k A^3), so the integrated actual variance is at least `A^2 b_Q^2-O(Lambda_k A^3)`.

For the correctly time-parametrized force `F1=A integral_0^infinity exp(-s)X_s ds`, its variance is A^2/2 and its conditional mean is Az/2; its conditional variance is A^2/4. Coherent contraction changes this by O(A^3), uniformly in depth. Thus `13/144` is a valid positive leading variance gap, with Lambda_k A sufficiently small. The integrated expectation gap lower-bounds the L2_y gap. This does not rely on taking a quadrature limit.

The completed-law carrier-subtraction example is also correct, with m(0)=0 as now stated. Equal marginal Gaussian laws can have zero versus `2(1-cos A)` carrier-subtracted variance. A marginal LAW therefore cannot supply the missing coupling.

## 5. Costs, roots, radii and precision

The count recurrence expands every S_k and B_k occurrence, including origin/capture/restoration and discarded-primal replays. The inherited Gram and mixed-K services are charged at their full pinned counts. The baseline at the next depth costs one original VALUE per outer node. Fixed depth preserves a public-log polynomial count, but its degree and coefficients can grow with depth.

Before alignment a bulk node has its Y root plus all three complete service dimensions; a near node has four D-roots. Sharing one carrier row across all nodes removes `(N_out-1)D`, leaving exactly

    D + sum_bulk(d_M+d_H+d_K) + 3D N_near.

Each reentered occurrence carries a complete fresh d_k tape. This is a graph dimension identity, not dimension estimated from a LAW target.

At u comparable to w, the physical raw-split first is `Lambda[A+A^(3/2)/sqrt(u)]`. The displayed bound certifies the same O(Lambda A) source class on w>=A; mere smallness of normalized native radii would allow a larger range and would not suffice for this closure. The normalized radius powers are `1-alpha/2` and `3/2-alpha`; the mixed-K power is `1-alpha/3`. The claimed endpoint values 3/5, 7/10, 11/15 and 1/2, 1/2, 2/3 are correct. All share constants, actual Lambda values, dimensions, positive gaps and native orders still require admission at each call.

Inherited target floors acquire A in the mean recurrence, while fresh primitive/replay floors are budgeted anew using complete downstream paths. No intrinsic forcing term is being hidden in numerical precision.

## Diagnostics and corrections

`check_depth_ports.py` independently passes **692 exact rational assertions** covering raw-split scale exponents, physical source closure, complete root-count algebra and the finite affine recurrence. Results and the reviewed hash are in `depth-port-checks.json`. The author's separate 4,721-assertion diagnostic was inspected; its scalar and linear-Gaussian scope is appropriately limited. Neither diagnostic executes the native compiler or proves the nonlinear source theorem.

Two notation issues found during review were corrected before the pinned hash above: the linear F1 formula now uses OU time consistently, and the carrier-subtraction example explicitly takes m(0)=0. No substantive mathematical correction remains outstanding within this stated finite-depth scope.

## Separate audit: independent outer cutoff and weighted closure

**PASS under the stated actual native guards, including the strict aggregate extension beta<3/2.** Reviewed addendum: `../INDEPENDENT-OUTER-CUTOFF-ADDENDUM.md`, SHA256 `ddbe52d0ccf04956f0948d8a4429f5aae135185c141026eaa2d54ebd086c63ee`. The main note's reviewed hash above is unchanged.

This addendum correctly distinguishes the force-tail delay w from the endpoint cutoff eta. In particular, the earlier prospective prefix-only 10/3 bound applies to the tied choice eta=w. It is not a two-parameter obstruction.

### Conditional variance and integrated debts

For `q=1-w`, the identity

    1/v(t)=q^2/(1-q^2)+1/(1-t^2)

is exact. The variance decreases with t. At the largest bulk correlation `t=1-eta`, put `c^2=eta(2-eta)` and `sigma^2=w(2-w)`. Since eta<=w, c^2<=sigma^2 and

    v_min = c^2 sigma^2/(sigma^2+q^2 c^2)
          >= c^2/2 >= (3/4)eta.

Also v<=sigma^2<=2w. The smallest of the actual mean/Gram/K unscaled shares is at least `3 eta/16`; no service is assigned a zero Gaussian gap.

For p>=1, split the two positive terms in 1/v using the stated convexity inequality. The first term contributes O(w^(-p)); the endpoint term is compared with `integral_eta^1 h^(-p) dh`, because `1-t^2` is comparable to `1-t`. This proves the bulk p=1, 3/2 and 7/6 bounds. For p=1/2 use subadditivity. The whole endpoint integral is finite (`integral_0^1 (1-t^2)^(-1/2)dt=pi/2`), proving the bulk O(w^(-1/2)) bound; restriction to a near interval of length eta gives `O(eta/sqrt(w)+sqrt(eta))`.

The positive finite rule has the same estimates. On a dyadic endpoint-gap panel h, its mass is O(h) and node gaps are comparable to h, giving O(h^(1-p)) per panel. Cutting an existing panel at eta only decreases its mass; its node gaps remain comparable to the original panel's scale. The finite early midpoint lies in the near branch when its width is chosen below eta, and contributes O(sqrt(width)) to the square-root singularity. No smoothness or semigroup identity for the surrogate is used. The fixed true target is still quadrature-approximated first.

The generalized recurrence charges the correct singular integral to every original debt. Endpoint mean mismatch has mass O(eta); covariance restoration uses the integrable p=1/2 bound; all higher mixture terms remain dominated by the mean/Gram term. Replacing the branch sets changes literal counts but creates no inverse-eta replication. All native radii and absolute precision paths use the actual smallest u, not w.

### Aggregate ports when eta<A

Section 1's eta>=A assumption is sufficient for O(A) residual first at each node. Section 5 correctly proves that it is unnecessary for the **final positive source**. For the actual per-node residual majorant

    R_j <= Lambda[A+A^(3/2)/sqrt(v_j)],

positive summation and the finite p=1/2 estimate give

    sum_bulk omega_j R_j
      <= Lambda[A+A^(3/2)/sqrt(w)] = O(Lambda A)

whenever w>=A. The exact common-carrier rotations have operator norm one; combining independent complementary root blocks with the shared carrier permits this triangle bound for the complete bank and retained-caller derivatives. Independence across different outer nodes is unnecessary.

No intervening native consumer requires each physical R_j to be O(A): that node is consumed directly by the original terminal F. In the terminal derivative, the potentially O(A) Hessian difference on the common carrier is symmetric. Every nonsymmetric or complementary-root path contains the additional original A factor multiplying R_j, apart from the already bounded w A^2 composition term. Positive summation thus gives aggregate curl O(Lambda A^2). The literal source-zero/origin records, complete root dimensions and the same weighted displacement estimate give the aggregate O(Lambda A^2) energy profile. The reusable total first is O(Lambda A). These are direct graph bounds, not derivatives of an averaged LAW estimate.

For fixed eta=A^beta with `1<beta<3/2`, the worst native powers are `1-beta/2`, `3/2-beta`, and `1-beta/3`; all are positive. Every actual logarithmic factor and numerical radius threshold must still be substituted separately. The requirement u>=A^3 follows for sufficiently small A. **There is no beta=3/2 admission claim:** its normalized self-reserve power is zero. The explicit beta=7/5 powers 3/10, 1/10 and 8/15, and near-RAW grade 37/10, are correct.

### Scalar envelope and scope

Write w=A^alpha and eta=A^beta. On the initial window `0<alpha<=beta<=1`, the unmodified optimum is 16/5 at alpha=4/5, with any beta in [4/5,1]. With weighted closure the same optimum holds for beta in [4/5,3/2). The unchanged prefix and A^4/w lines already prove the upper envelope; taking beta=1 attains it within the admitted range.

If only the leading prefix term is deleted from the displayed ledger, A^3 w and A^4/w instead cap its minimum at 7/2, attained at alpha=1/2 and beta=1. Weighted closure permits beta in [1,3/2) at this same grade. This is a prospective ledger calculation, not a prefix construction or a consequence of this addendum alone. No improved original forcing exponent or combined prefix/skew theorem is audited here.

`check_independent_cutoff.py` passes **1,201 exact rational assertions**, covering conditional variance/share bounds, the complete vertex enumeration of the two-scale exponent ledgers, and strict-window radius powers. Its use of the closed relaxation beta<=3/2 only supplies a scalar upper bound; both claimed optimal grades have admitted witnesses at beta=1. `independent-cutoff-checks.json` records the addendum hash and scope. Neither diagnostic constructs a new producer or executes a native compiler.
