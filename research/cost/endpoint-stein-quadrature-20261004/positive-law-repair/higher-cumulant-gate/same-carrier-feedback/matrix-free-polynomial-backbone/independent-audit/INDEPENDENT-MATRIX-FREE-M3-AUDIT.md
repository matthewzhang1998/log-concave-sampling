# Independent audit: matrix-free Gaussian-backbone m3 consumer

2026-10-05.

## Verdict and frozen evidence

PASS for the stated source-qualified canonical m3 target comparison and the guarded positive fixed-order mean-law component, relative to the already pinned gradient/near-gradient mean compilers. The finite source is genuinely original-VALUE-only, uses four VALUES per terminal occurrence, and is exact on every anchored admissible matrix quadratic without recovering or querying the analytical comparison matrix.

The analytical comparison requires all five stated centered-Gaussian residual norms. It does not provide general C2 closure. Independently, the smooth D=1 source g_A(x)=A(x+sin x)/2 gives an actual order-A^2 target-bias obstruction to the unrestricted interpretation; the complete proof is in `GENERAL-C2-BIAS-OBSTRUCTION.md`.

Audited source: `../MATRIX-FREE-GAUSSIAN-BACKBONE-M3.md`.

Final source SHA256:

    772a0b758f7a4294ec9127629bec5bdd78b1f2ba167f32912ed5e7495e3ba06b

The auditor did not edit the author theorem. Source changes proposed during review concerned the anchored caller constant, nested numerical VALUE floors, and the distinction between adding previously excluded sources and claiming set inclusion between all sufficient certificates. The final source includes the correct distinctions.

Independent diagnostic script: `check_matrix_free_m3_independent.py`.

Recorded output: `independent_matrix_free_m3_checks.json`.

The checker passes 1,540 assertions. It does not import the author's checker. Its maximum complete directional chain-rule discrepancy is 3.611e-11 and maximum quadratic value discrepancy is 2.289e-16. It detects a genuinely nonsymmetric G derivative on a smooth convex-gradient fixture, with observed skew norm 0.002923; thus the tests would reject the erroneous scalar-backbone symmetry shortcut. It also verifies the exact first-Hermite obstruction formula by separate quadrature and encloses its sign and a decimal lower bound using exact rational Taylor bounds.

The numerical tests are diagnostics for the raw source, not an implementation or reproof of the whole imported finite mean compiler. The analytic derivations below discharge the new ports; original compiler requirements remain explicitly imposed.

## 1. The shared two-root covariance is exact

The independent covariance derivation starts from the stationary OU Laplace kernel

    L(a,b)=[1/(a+b)][1/(a+1)+1/(b+1)].

Differentiating its exponentially integrable scalar kernels at a=b=1 gives

    Var(H1)=1/2, Cov(H1,H2)=3/8, Var(H2)=3/8.

The regressions on X0=x are x/2 and x/4. Subtracting their products yields

    Cov((H1,H2)|X0)= [[1/4,1/4],[1/4,5/16]] tensor I.

The scalar root matrix [[1/2,0],[1/2,1/4]] has exactly this product with its transpose. Therefore

    (H1,H2)|x =_law (x/2+N/2, x/4+N/2+M/4).

The same N must occur in both entries. The extra independent M supplies the remaining covariance. This is an exact conditional Gaussian-law factorization, not Gaussianization of a nonlinear force.

For an analytical symmetric K in [0,AI], the backbone argument is x-Kv+K^2w, with conditional mean (I-K/2+K^2/4)x and covariance K^2/4-K^3/2+5K^4/16. Orientation is retained because all powers are actions of the same K. None is implemented as a new oracle: on the quadratic source g=K Id, g(v) and g(g(w)) implement those actions through the original VALUES.

The outer same-G,N,M sharing also matters. The raw quadratic source has covariance

    beta^2 K^2(I-K/2+K^2/4)^2
      +(1/4)K^4(I-K)^2+(1/16)K^6.

This is its own raw covariance. It is not identified with a full canonical outer-history covariance. Its conditional mean is exactly (1/2)K(I-K/2+K^2/4)Z because the positive rule has first moment 1/2. The finite mean consumer needs that own mean, not an unsupported raw-law identity.

## 2. The integrated five-residual comparison is valid

Let r=g-K Id, k=||K||, and let L be a valid Lipschitz bound on r. Since -AI<=Dr<=AI, L=A is always available. The sharper value max(lambda_max(K),A-lambda_min(K)) is also valid; no such improvement is necessary.

On the genuine future history,

    Delta1,t=int e^(-u)r(X_(t+u))du,
    ||Delta1,t||2<=e0.

The same-history Y_t=X_t-KH1,t has covariance I-K+K^2/2 exactly, so

    ||r(X_t-F1,t)||2<=eY+L e0.

Expanding the actual second substitution preserves its genuine future ancestry and yields

    F2=KH1-K^2H2+R2,
    ||R2||2<=eY+(k+L)e0.

The H2 kernel t e^(-t) comes from Fubini on the actual double history integral. There is no independent inner resampling. Applying A-Lipschitz g and conditional Jensen gives the history target error A[eY+(k+L)e0].

For the executed matrix-free graph, the exact comparison identity is

    g(g(w))-K^2w=[g(g(w))-g(Kw)]+r(Kw).

It implies |T_g-(x-Kv+K^2w)|<=|r(v)|+A|r(w)|+|r(Kw)|. Unconditionally under original standard x,N,M,

    v~N(0,I/2), w~N(0,3I/8), Kw~N(0,3K^2/8).

Minkowski does not require these variables to be independent. Hence the additional finite action-replacement error is A[ev+A ew+eKw]. This independently confirms the sharper A coefficient on ew, rather than a needless k+L coefficient.

Adding both errors and applying R1 contraction proves exactly B_res in the author theorem. It is integrated over the standard endpoint; there is no uniform conditional version or differentiated target-bias estimate. K and all five expectations are analytical qualifications, absent from the executing graph. The theorem does not claim that finitely many black-box VALUES can automatically certify them.

## 3. Positive outer rule and value energy

At a fixed retained Z, every x_i=t_iZ+c_iG has the correct OU conditional marginal. Shared roots across nodes do not change linearity of expectation. Thus the exact source own mean is Q psi_*(Z).

At an independent standard endpoint,

    ||T_g||2<=sqrt(D)+A sqrt(D/2)+A^2 sqrt(3D/8).

Therefore ||psi_*||2<=A C_T sqrt(D), and the imported Hermite-operator quadrature certificate gives delta A C_T sqrt(D). No pointwise quadrature claim is made.

The positive rule has beta=sum w_i sqrt(1-t_i^2)<=sqrt(3)/2 by concavity and the exact mean t=1/2. All c_i are strictly positive. The baseline B is genuinely a gradient in G with Hessian sum w_i c_i Dg(x_i), bounded between zero and A beta I. Its analytical potential is never an executing query.

The literal displacement estimate is

    |E|<=A^2 sum w_i |v_i|+A^3 sum w_i |w_i^arg|.

Given Z=z, the conditional means are t_i z/2 and t_i z/4 and conditional variances are (c_i^2+1)/4 and (c_i^2+5)/16. Gaussian moment bounds and the first moment 1/2 give the author's explicit amplitude (5.4). This bound is valid without Hessian continuity or any smallness of A sqrt(D).

## 4. Full first, curl, caller and zero ports

Write H=Dg(T), V=Dg(v), J=Dg(g(w)), W=Dg(w), and H0=Dg(x). All are symmetric in [0,AI]. The complete input derivatives are

    D_x f=H(I-V/2+JW/4),
    D_N f=H(-V/2+JW/2),
    D_M f=HJW/4.

The baseline subtraction changes only the x derivative, subtracting H0. Products retain their actual order. The checker also implements the reverse sweep independently, with WJH in the adjoint path, and verifies it against the transpose of the full Jacobian.

Let eta=A/2+A^2/4, q=A/2+A^2/2, r0=A^2/4. Since ||H-H0||<=A, the raw bounds are

    ||D_G E||<=A beta(1+eta),
    ||D_N E||<=A q, ||D_M E||<=A r0,
    ||D_(G,N,M)E||<=A sqrt(beta^2(1+eta)^2+q^2+r0^2),
    ||D_Z E||<=A(1+eta)/2.

The large leading G-block difference is symmetric. The remaining G-block antisymmetry has norm at most 2 beta A eta. The (N,M) off-diagonal row has norm at most A sqrt(q^2+r0^2). Block addition proves the complete lifted curl bound

    ||D(P_G^*E)-D(P_G^*E)^*||
        <=2 beta A eta+A sqrt(q^2+r0^2).

This is O(A^2), but its G block is generally NOT zero. No local Hessians have been commuted. With half variance shares, A<=1/2, beta<=sqrt(3)/2, normalized first<=2A and curl<=3A^2. These inequalities follow by monotonicity in A,beta after dividing by A and A^2 and evaluation at the endpoint; the checker separately tests them on 101 values. Declaring ell_E=2A, a_seed=2A therefore safely bounds the curl by ell_E a_seed.

The source origins B0(z),E0(z) are their literal G=N=M=0 executions. They obey the stated A|z|/2 and A^2(1/4+A/8)|z| bounds, but generally do not vanish. The anchored correction's caller derivative is bounded by A(1+eta), twice the raw one; the baseline's is bounded by A. Restoring the origins retains their live caller paths. At total z=G=N=M=0 all nested descendants vanish by g(0)=0. At the zero source the complete imported mean service supplies its known Gaussian row, without an energy division.

## 5. Guarded mean-law completion and executable costs

Relative to the existing imported mean contracts, the new ports suffice. The B and E branches use independent COMPLETE banks conditional on the same endpoint and exterior labels. The near-gradient square lift has dimension 3D and projection onto the G coordinates. It retains the recorded-coisometry Gaussian completion; it does not append G,N,M or any force as an observer after completion.

Use rho_B=sqrt(2)A beta, ell_E=2A, a_seed=2A and padding mu=A. The actual native gradient radius, ell_E<=1/4, curl, active-dimension, finite-clock/filter, precision, finite-mode and caller guards are additional hypotheses. A<=1/8 discharges only the displayed ell_E guard. No broader smallness claim is inferred from it.

The imported near-gradient error bracket ell_E(a_seed+mu)+ell_E^3(1+mu^(-1/2)) is O(A^2). Its anchored energy is O(A^2)(|z|+sqrt(D)); the completed branch thus has O(A^4)(|z|+sqrt(D)) error plus absolute floors. The order-at-least-four gradient branch has its separate stated allowance. Conditional product coupling is applied only after both complete returns, yielding the unit-buffer own-mean law. Equal-covariance Gaussian mean distances then combine the residual and quadrature errors.

A raw F occurrence evaluates g(v), g(w), g(g(w)), and g(T): exactly four original VALUES per node. E adds g(x), for five. B needs one. Every changed argument replays the complete nested graph; g(w) is an original VALUE ancestor of g(g(w)). The safe expanded count is therefore

    Q_certificate_setup+Q_captured
       +N_out N_B+5N_out N_E+Q_known/numerical/replay.

The caller origins themselves cost N_out and 5N_out before complete-key reuse. N_B,N_E include every native occurrence and replay. Raw F/E own 3D private roots, B owns D; all compiler roots, fills, filters, clocks, banks and numerical work are additional charged entries. Independent complete banks may not share an old cached random ancestor. First/adjoint sweeps use original HVPs at the recorded original VALUE sites and retain the full caller graph. They never differentiate an HVP. Discarded primals pay replay.

If original VALUE error is bounded by nu at every actually requested site, propagation through the nested descendants gives

    raw F error <=nu+2A nu+A^2 nu=(1+A)^2nu,
    raw E error <=[(1+A)^2+1]nu.

The independently captured anchored E-E0 doubles this safe absolute allowance. The nodewise version nu_f+A nu_v+A nu_b+A^2 nu_w is correct. Original node/root/clock/filter/mode/replay floors remain separate and are never divided by a measured energy or residual. Literal numerical zeros require coherent exact-key reuse or the known source anchor.

The inverse-A original VALUE exponent remains zero at the declared fixed grade, relative to the same native occurrence bounds. This does not establish growing-order complexity or dimension-independent arithmetic.

## 6. Nonlinear class and the general obstruction

The sample-free radial-tail addition to an arbitrary K with spectrum in [A/3,2A/3] is admissible: the added radial Hessian lies in [0,(A/3)I], leaving Dg in [0,AI]. Its residual about K is (A/3)-Lipschitz and bounded by (A/3)(|x|-R)_+.

Every Gaussian linear map in the five norms is a contraction for A<=1/2. A radial envelope therefore bounds all five norms by sqrt(2)A^3/3 at R=sqrt(D)+sqrt(8 log(1/A)); this uses a pointwise envelope, not false monotonicity of arbitrary residual norms under Gaussian contraction. Inserting k+L<=A gives precisely (sqrt(2)/3)(3+2A)A^4. No density ratio exponential in dimension appears. With distinct eigenvalues of K, no scalar slope can have bounded residual at infinity. The included anisotropic quadratics also violate the old scalar target at order A^2 sqrt(D), while the present graph is exact.

These examples add sources excluded by the old scalar-backbone certificate. The present five-norm sufficient class and the earlier two-norm Gaussian-RMS sufficient class are not asserted to be set-ordered on every source.

The independent obstruction file proves, for the original smooth admissible source g_A=A(x+sin x)/2,

    ||R1 psi_*-m3||2>=c A^2-3.04A^3,
    c=0.009194688258588945... .

The coefficient is an exact first-Hermite moment, not a fitted numerical rate. The positive outer quadrature at delta<=A^3 and the admitted own-mean law error contribute only O(A^4), so they do not turn this example into general order-four closure. This establishes a genuine obstruction to that interpretation while preserving the source-qualified theorem.

## 7. Imported provenance and final scope

The following unchanged sources were read or pinned through their existing independent audit:

- Gaussian-RMS scalar predecessor: SHA256 `18640f9c6c567acba4fe830c7b0d64394924887ee4335a7a49f57f817435dccf`.
- Bounded scalar finite-consumer theorem: SHA256 `798d20b87124e7069e465ce4e7bcf92c5af3734d86c8a62aeb288b5fec8eaef3`.
- Its independent finite-consumer audit: SHA256 `49d135ba6692c54012acd9d0d94f603c0015fca98b5d1be1164457ff45d1213a`.

That audit pins the LOW30 value-mean and fixed-radius gradient compiler, the positive Hermite rule, complete independent-branch composition, the recorded-coisometry adapter, origins, source-zero row and full replay ledger. The matrix-free source changes the raw graph, adds one private D-root, and supplies the new chain/curl/residual bounds audited here. The present checker does not silently substitute for those finite-compiler proofs.

No blocker remains for this scoped component. Still absent are general C2 m3 closure, automatic residual certification, arbitrary conditional rescaling, a strong conditional mean/covariance statistic, any law retaining private raw roots, the reverse-OU endpoint join, or an all-order complexity recurrence. Under conditional normalization the actual five translated residual norms and original-g caller implementation must be re-established; the original bulk certificate alone is insufficient.
