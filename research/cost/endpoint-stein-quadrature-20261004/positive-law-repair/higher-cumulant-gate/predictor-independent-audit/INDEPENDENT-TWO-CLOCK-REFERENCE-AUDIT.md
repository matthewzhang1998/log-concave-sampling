# Independent bounded audit: two-clock current and stationary C2 reference

2026-10-04. **PASS at the pinned scope.**

Audited source: `../EXACT-TWO-CLOCK-CURRENT-AND-UNIFORM-C2-REFERENCE.md`

SHA256: `0a966fbd6076dff94d2b928ec1a4a05aa6d33c9bc8222eff051b221ee9881477`

## Verdict

The following claims are sound:

1. The displayed two-clock formula is the exact second-order Gaussian tilt current for every admissible C2 gradient potential, with both rank-one and rank-two terms and the literal nested Gaussian correlations.
2. The stationary resolvent gives a **uniform analytical** two-substitution comparison `W2(P2,mu_U)<=A^3 sqrt(D)` using only the global Lipschitz and anchored monotonicity bounds. It does not require a Hessian modulus.
3. The finite Markov-clock VALUE graph has the stated original query count, root dimension, zero, absolute numerical propagation, and dimension-free first/energy/caller ports.
4. The finite quadratic covariance screen is correct and is a genuine additional floor.

**No finite order-three law theorem is proved.** The comparison of the finite graph to P2 in the required one-energy law/current norm remains open. This audit neither substitutes ordinary path-grid strong error for that comparison nor resumes that discretization route.

The audit supplies 52 exact symbolic Wick checks of the full current and 41 independent finite-port/quadratic-screen assertions. The numerical checks are supplementary; the proofs below justify the scoped claims.

## 1. Reusable exact bilinear current lemma

Let Z,G,H be independent standard Gaussian vectors. For r,t in [0,1], put

    X=rZ+sqrt(1-r^2)G,
    B=tX+sqrt(1-t^2)H.

For an anchored C2 potential U0 with bounded Hessian and a smooth compactly supported phi,

    (1/2) Cov(phi(Z),(U0(Z)-EU0)^2)
      =integral dr integral dt E[
          Dg0(X)g0(B) dot grad phi(Z)
          +r Sym(g0(X) tensor g0(B)):Hess phi(Z)].

To verify, apply the Gaussian covariance identity first to phi and `(U0-EU0)^2/2`, obtaining

    integral dr E[(U0(X)-EU0) g0(X) dot grad phi(Z)].

Represent `Z=rX+sqrt(1-r^2)G'`, where G' is standard and independent of X. A second Gaussian covariance identity, in X with G' held independent, differentiates

    g0(X) dot grad phi(rX+sqrt(1-r^2)G').

Its two derivatives are `Dg0(X)grad phi(Z)` and `r Hess phi(Z)g0(X)`. Contracting with `g0(tX+sqrt(1-t^2)H)` gives exactly the claimed formula. The second identity introduces no extra t factor.

The covariance of the triple `(Z,X,B)` is

    [[1,r,rt],[r,1,t],[rt,t,1]] tensor I.

The fresh H is independent, but X inside B is the same actual X that the outer force reads. Independent replacement of that ancestor changes the current. The rank-two term is the full symmetrized product, not only its conditional covariance; the conditional means also contribute.

Only g0 and Dg0 appear. Gaussian domination, bounded Dg0, and smooth test derivatives justify the covariance identities directly, or through a smooth approximation. No derivative of Dg0 is used. The Gaussian tilt expansion is legitimate because U0 has quadratic growth and its small-amplitude exponential moments are controlled. This coefficient identity alone is not a uniform cubic remainder proof.

`check_bilinear_identity.py` independently uses Gaussian Wick recursion to check potential Hermite pairs `(1,2),(2,3),(2,5),(3,4),(3,3),(4,4),(4,5),(5,5)` against carrier ranks through eight. All 52 exact symbolic tests pass. Its JSON records the separate r,t polynomials for the two current terms and their exact integrals, so the all-bilinear identity is reusable without relying only on the sinusoidal third-cumulant fixture.

## 2. Markov genealogy and the continuous feedback expression

A stationary Gaussian OU field indexed by 0<r<=1 has

    Cov(X_r,X_s)=min(r,s)/max(r,s), X_1=Z.

For s=rt, its triple `(Z,X_r,X_s)` has the preceding covariance matrix. The half-square of `h=integral_0^1 g0(X_r)dr` splits into the two ordered triangular domains. The change of variables s=rt supplies `ds=r dt`, producing precisely the r-weighted rank-two term. The amplitude-two rank-one term is the nested feedback in

    P2=Z-integral dr g(X_r-integral dt g(X_(rt))).

This validates the continuous genealogy. It does not identify that OU field with the old common-two-root field, whose covariance is `rs+sqrt(1-r^2)sqrt(1-s^2)`.

## 3. Uniform C2 stationary resolvent proof

Let X and Y be stationary solutions driven by the same two-sided Brownian motion:

    dX=-X dt+sqrt(2)dB,
    dY=-(Y+g(Y))dt+sqrt(2)dB.

The total Y drift is globally Lipschitz and strongly monotone. Coupling finite-past solutions and passing to their stationary L2 limit therefore provides the joint stationary coupling used in the source. Its Y marginal is the target density proportional to `exp(-|y|^2/2-U(y))`.

Anchoring and monotonicity imply `y dot g(y)>=0`. Target integration by parts gives

    D=E[Y dot(Y+g(Y))]>=E|Y|^2.

Subtracting the two SDEs and taking the stationary limit of variation of constants gives

    Y_t=X_t-integral_0^infinity e^(-s)g(Y_(t-s))ds.

The boundary term has L2 norm tending to zero because of its factor e^-T and the uniform stationary second moments. Minkowski and `|g(y)|<=A|y|` give

    ||X_t-Y_t||2<=A sqrt(D).

Define P0=X and the next analytical substitution by replacing Y in the integral with Pj. Then

    ||P_(j+1),t-Y_t||2
       <=A integral_0^infinity e^(-s)||P_j,(t-s)-Y_(t-s)||2 ds
       <=A sup_s ||P_j,s-Y_s||2.

All Pj inherit stationarity from X, so the supremum is finite and the recurrence gives

    ||P_j,t-Y_t||2<=A^(j+1)sqrt(D).

For j=2, changing variables `r=e^-s`, `t=e^-u` gives exactly the continuous feedback expression of Section 2. The coupling is an admissible law comparison, hence

    W2(Law(P2),mu_U)<=A^3 sqrt(D).

This argument is uniform over the Hessian sandwich and never expands the Hessian. It is materially stronger than merely matching formal amplitude coefficients, but it remains an **analytical reference** with continuum Gaussian and force integrals. It carries no finite cost statement.

For a conditional standardized source with alpha=A s^2, the same calculation yields normalized alpha^3 sqrt(D), then physical `s alpha^3 sqrt(D)=A^3 s^7 sqrt(D)`. Multiplication by the exact reverse-OU coefficient `Delta/(t s^2)` gives `(Delta/t)A^3 s^5 sqrt(D)`. Finite-mode, tilt, and numerical restoration floors remain additional terms, as stated in the source.

## 4. Literal finite graph ports

The repaired source expressly requires `0<r_i,t_j<=1`; thus every clock has a finite logarithmic time and the Gram ratio is defined. Draw the finite set of distinct times

    T={1} union {r_i} union {r_i t_j}

with its known PSD scalar Gram, tensor I_D. Sorted-time successive OU transitions realize exactly this law in the stated Gaussian real-arithmetic convention. The independently constructed transition factor and a Cholesky factor yield the same Gram to numerical accuracy. Numerical covariance-factor precision is separately billed.

Execute

    I_i=sum_j v_j g(X_(r_i t_j)),
    Y_Q=X_1-sum_i w_i g(X_(r_i)-I_i).

There are at most `n_r n_t+n_r` original VALUES after captured mode and anchor work. Repeated times/sites can save work only under identical complete semantic keys; all feedback reads the actual recorded I_i. Requested first/adjoint sweeps use original HVPs at those VALUE sites. A discarded primal graph requires its complete replay.

If each X_t=L_t R for the complete independent Gaussian root R, `||L_t||op=1`. For a root perturbation dR, the correction's derivative has two pieces. The outer direct piece is bounded by

    sum_i w_i A ||L_(r_i)|| <= A,

and the inner feedback piece by

    sum_(i,j) w_i v_j A^2 ||L_(r_i t_j)|| <= A^2.

Thus `Lip(Y_Q-X_1)<=A(1+A)` without any factor involving the number of clocks. The weighted-stack proof in the source is equivalent; its scalar Gram has trace one. Every X_t is marginally standard, so Minkowski yields the stated `C_p A(1+A)sqrt(D)` correction energy.

At all complete roots zero, every source input and output is zero under recorded-anchor reuse. The raw output is exactly zero. In an original conditional problem this corresponds to the captured finite mode after the physical readout, not to a newly obtained exact mode.

For `f(y)=s[g(x+s y)-g(x)]` and a captured mode derivative `||D_a x||<=2`, the direct fixed-input force caller bound is `CAs`. Inner weighted forces preserve that bound; their outer feedback multiplies it by at most `1+alpha`. Physical readout contributes s. The private residual first is `s alpha(1+alpha)=As^3(1+As^2)`. Actual physical caller and mode/anchor dependencies remain in the graph. No derivative of a law-error estimate is taken.

The complete Gaussian dimension is `D |T|`, not 2D. Fresh quadratures and Gaussian covariance roots remain outside the retained caller; after a law comparison the private nodes are not additional retained observers. Positive weights propagate absolute force error with factor at most `1+A`; coefficient, Gram, clock and source versions must be frozen before first/adjoint sweeps.

The independent executable check uses 5 outer and 4 inner interior Gauss nodes, 26 distinct Gaussian times, and at most 25 original VALUES. It verifies the Markov Gram, row and weighted-stack norms, nonlinear recorded-first finite differences, exact zeros, and the quadratic screen. All 41 assertions pass. These tests check ports, not finite quadrature accuracy in law.

## 5. Finite quadratic screen and the remaining gate

For `g(x)=Bx`, the finite graph is

    Y_Q=Z-B H_Q+B^2 J_Q.

When both scalar first quadrature moments are 1/2, the two cross covariances are `Cov(Z,H_Q)=I/2` and `Cov(Z,J_Q)=I/4`. Therefore the B^2 coefficient is `S_Q+1/2`, with

    S_Q=sum_(i,l)w_i w_l min(r_i,r_l)/max(r_i,r_l).

For the continuum field, direct integration gives S=1/2. For the independent 5-node outer fixture, `S_Q=.5267409643853054`, so the degree-two covariance defect is `.026740964385305377`. It does not disappear because all individual Gaussian marginals are correct. The full finite polynomial is also checked independently.

As an additional analytical-only screen, the continuous P2 covariance for a scalar quadratic with coefficient b is

    1-b+b^2-(3/4)b^3+(3/8)b^4.

This follows from `Var(H)=1/2`, `Cov(H,J)=3/8`, and `Var(J)=3/8`; it is consistent with the cubic resolvent bound and differs from a generic finite quadrature polynomial.

The unproved step remains a finite original-VALUE law/current approximation to P2 at order `A^3 sqrt(D)` with public-logarithmic clocks, plus all numerical/restoration floors. Finite moment fitting, the raw first bounds, and the analytical resolvent contraction do not imply this step. Covariance-only correction also does not automatically retain the full rank-two conditional product or the higher descendants. The source correctly leaves this comparison as the active gate.
