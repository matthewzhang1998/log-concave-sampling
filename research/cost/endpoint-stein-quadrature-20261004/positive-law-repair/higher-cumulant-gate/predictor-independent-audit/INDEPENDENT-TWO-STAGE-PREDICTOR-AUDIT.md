# Independent audit: the two-stage positive VALUE predictor

2026-10-04. **PASS at the pinned scope.** Independent algebra, executable checks, and complete proof-file cross-check.

Audited source: `../TWO-STAGE-VALUE-REPAIR-OF-THE-THIRD-CUMULANT-WORD.md`

SHA256: `8c271f17971c1f24cb42e27e1ace4179a7cf19716d24f32755a53ac892b57d98`

## Result and exact scope

The candidate passes the following bounded claims:

- It is an explicit positive pushforward of the same two Gaussian roots, executing `2n+1` anchored original-force VALUES.
- Its finite known coefficients match the full matrix-quadratic covariance through degree two and the **third-cumulant** quadratic-times-odd amplitude word through degree two.
- Its captured-u private-v chord source has an actual order-alpha-squared curl, exact executable origin, order-alpha-squared VALUE energy, and an order-alpha first bound under the general C2 assumption.
- The finite quadrature, finite coefficient arithmetic, anchor, mode, and replay floors remain part of the result. Fixed finite quadrature does not imply exact first-order target matching.
- A fixed C2 shape permits the displayed degree-two coefficient expansion with `o(A^2)` remainder. A uniform `O(A^3)` claim over the full C2 class is not licensed.

The candidate **does not** match the full quadratic-times-odd second-order current. Its fifth cumulant already fails at degree two on the same smooth sinusoidal family. The next-row formula and a two-variable Gaussian kernel identity are recorded below. Neither is a completed positive current compiler.

Independent checks: 130 assertions in `check_predictor_independent.py`, plus 52 exact symbolic assertions in `check_bilinear_identity.py`. The independent code neither imports nor runs the author's checker.

## 1. Executed graph and finite calibration

Let the frozen positive quadrature have mass one and exact first moment `m1=1/2`; put `q_i=sqrt(1-r_i^2)`,

    K(u,v)=sum_i w_i f(r_i u+q_i v),
    K*(u,v)=sum_i w_i f(r_i[u-aK-bf(u)]+q_i[v-cK]),
    Y=u-K*.

All K values and roots in this formula are the actual shared records. Negative a, b, or c would be deterministic arithmetic, not signed probability; here a is positive and b,c are negative.

Write `M_j=sum w r^j`, `N_j=sum w r^j q`, `beta=N_0`, `C=1/4+beta^2`, and

    Kq=M2+2 beta N1,
    Ju=M4+2 beta N3, Jv=N3+2 beta(M2-M4), J0=2M4,
    Lu=3M2, Lv=2(N1+beta M2), L0=1+4M2.

The calibration is exactly

    [1/2 beta 1] [a,c,b]^T = 1-C,
    [Ju  Jv   J0][a,c,b]^T = 1/3,
    [Lu  Lv   L0][a,c,b]^T = 2-2Kq.

For the 241-node checker version the independent solve gives

    a= 1.3809857645648347,
    c=-0.32011162098909196,
    b=-0.3059280782101049,
    determinant=-0.12976040799766486.

In the continuum limit, with `p=pi`, the exact determinant is

    -(5p^2-7p-4)/180,

and the coefficients are

    a=(2p^3+5p^2-6p-28)/(2(5p^2-7p-4)),
    c=-(7p^3-13p^2-20p+4)/(4(5p^2-7p-4)),
    b=(p^2+4)(p^2-7p+8)/(8(5p^2-7p-4)).

Their high-precision values are respectively

    1.38098576500332354245635079841,
    -.320111621161164597928534206107,
    -.305928078227588260032793038581.

The tiny distinction between finite and limiting coefficients is real. The finite determinant gap verifies the checker version; an assertion uniform over a family of rules must retain its moment-accuracy guard and the resulting perturbation bound away from zero.

## 2. Entire matrix-quadratic covariance

For `f(x)=Bx`, with B symmetric, put `d=a/2+c beta` and `e=d+b`. Then

    Y=(I-B/2+eB^2/2)u+(-beta B+d beta B^2)v.

Its exact covariance is

    I-B+(C+e)B^2-(e/2+2d beta^2)B^3
       +(e^2/4+d^2 beta^2)B^4.

The first calibration says `C+e=1`, so this matches the entire matrix coefficient of `(I+B)^(-1)` through degree two, including all old shared-v cross-node products. The check uses rotated positive matrices, not just diagonal fixtures. For small norm B, the actual covariance stays uniformly positive and the Gaussian root comparison gives the claimed cubic quadratic-family error. This is not a nonlinear-law theorem.

## 3. Third-cumulant matching, without a hidden higher-derivative assumption

For `U_A=A[q x^2+epsilon psi(x)]`, with psi odd and `psi'(0)=0`, let `g=psi'` and define

    c2=E[H2(N) g(N)], c4=E[H4(N) g(N)], N~N(0,1).

Assume the original anchored force is in the stated C2 sandwich class; the Gaussian integrals then exist. The executed graph, for each fixed shape, has

    Y=u-Ah+A^2 T+o_Lp(A^2),
    T=a h_u h+c h_v h+b h_u f0(u).

An independently obtained centered-cumulant expansion gives the quadratic-times-odd coefficient

    3 q epsilon [(aJu+cJv+bJ0)c4
             +(aLu+cLv+bL0+2Kq)c2].

The target coefficient is

    q epsilon(c4+6c2).

For an explicit regression check, write X=r u+q_i v and L=u+2 beta v. Then

    E[H2(u) g(X)]=r^2 c2,
    E[H2(u) L g'(X)]
        =r^2(r+2 beta q_i)c4+2r c2,
    E[u L g(X)]=E g(N)+r(r+2 beta q_i)c2.

The first two formulas give the three feedback contractions. The last gives the rank-two term; its `E g(N)` contribution cancels the centered-cumulant mean/variance correction exactly, leaving `6q epsilon Kq c2`.

The last two calibration rows equate them exactly. This proof uses Gaussian regression and integration by parts with g and g'; higher formal derivatives can be transferred to Hermite polynomials. Thus the extension to every admissible odd perturbation in this **one third-cumulant word** does not require a Fourier completeness claim or a C5 hypothesis.

For `psi(x)=sin(kx)-kx`, `c2=-k^3 exp(-k^2/2)` and `c4=k^5 exp(-k^2/2)`, recovering the author's all-frequency formula.

The first-order third cumulant of the finite rule is `-3 A epsilon M2 c2`; the true coefficient is `-A epsilon c2`. This floor has not been eliminated by the degree-two calibration.

## 4. Original-VALUE, caller, origin, and replay ports

The raw graph performs n VALUES for K, one original VALUE for f(u), and n changed VALUES for K*. It uses no potential, density, expectation, covariance, Hessian-valued, or HVP-valued producer. For the conditional original force

    f(y)=s[g(x+s y)-g(x)], alpha=A s^2,

the safe original count is

    M+1+(2n+1)=M+2n+2,

including the captured mode work M and anchor g(x). The extra f(u) leaf is not the anchor unless its actual input is zero. Identical captured caller-only leaves can be shared; accidental equal outputs do not authorize sharing unequal inputs.

The two changed root coordinates are

    u*=u-aK-bf(u), v*=v-cK.

The pointwise shift and Lipschitz force give

    ||K*-K|| <= alpha[sqrt(a^2+c^2)||K||+|b| ||f(u)||].

Consequently the chord E=K*-K has `Lp` size at most `C_p alpha^2 sqrt(D)` for the complete roots and at most `C_p alpha^2(|u|+sqrt(D))` with u captured. The general first bound is `O(alpha)`, not `O(alpha^2)`.

With symmetric first blocks `Ku*,Kv*,Ku0,Kv0`, direct chain rule gives

    Dv K*=Kv*-(aKu*+cKv*)Kv0,
    Du K*=Ku*(I-aKu0-b Df(u))-cKv*Ku0.

The latter includes the indispensable extra f(u) dependency. The independent checker compares both derivatives to finite differences on a genuinely noncommuting nonlinear gradient, rather than checking only a copied identity. The private-v curl is exactly

    DvK*-(DvK*)^T=-(aKu*+cKv*)Kv0+Kv0(aKu*+cKv*),

and is bounded by

    2 (|a|/2+|c| beta) beta alpha^2.

Subtracting K leaves this same curl for E. With u retained, the actual origin E(u,0) must be executed as a caller-only two-packet record, including its f(u) leaf. The source `E(u,v)-E(u,0)` is exactly zero at private v=0, but its captured first also differentiates the origin record. The origin is not generically zero and is not free. A private-v mean replacement integrates that old v and the old K(u,v); it does not preserve either as a separately readable carrier.

At complete root zero, K=f(0)=K*=E=0. Under approximate original force evaluations this is literal only with the recorded-anchor reuse specified by the graph. The physical output is the captured finite mode x, not a newly asserted exact mode.

If `||D_a x||<=2`, then `||D_a f||<=4As` and the extra root feedback changes the K* caller bound by only `1+C alpha`. Thus `||D_a K*||<=C As`, `||D_a Q||<=C`, and the physical private first is `O(s alpha)`. These absolute first bounds do not differentiate a small law error or its small VALUE energy.

Requested first/adjoint sweeps use original HVPs at every recorded VALUE site and differentiate through K and f(u). No saved HVP is differentiated. Discarded primal records require full replay. An imported consumer's N changed occurrences therefore incur the complete occurrence cost, safely `M+1+(2n+1)N`, plus any additionally executed captured origin and consumer-specific tape. The same complete source cannot be billed as one new VALUE per occurrence.

## 5. Finite floors and the C2 remainder boundary

The checker uses exact Gaussian panels only as a mathematical quadrature description. Its final interval of length `h=2^-20` is replaced by a midpoint. Hence

    M1=1/2 exactly,
    M2-1/3=-h^3/12=-1/13835058055282163712,

before floating arithmetic. This nonzero floor is smaller than double precision here; the fact that a printed m2 equals 1/3 is not a proof that it is zero. Other moments, beta, and solved coefficients have their own finite versions. The exact degree-two finite calibration does not erase the first-order finite-law floor.

A proof must freeze the quadrature, known coefficient solve, target tolerance, original numerical precision, and imported clocks before differentiation. Leaf errors propagate through positive weight sums and bounded feedback constants. The determinant gap controls coefficient-solving errors. No numerical floor can be divided by the small actual chord energy. If a vanishing-error asymptotic is stated, the finite quadrature and numerical tolerances must be scheduled accordingly.

For a fixed C2 shape, differentiability of f0 plus domination by its global Lipschitz bound yields `o_Lp(A^2)` in the graph expansion. C2 alone supplies no uniform modulus of continuity for Df0 and hence no uniform cubic remainder over the class. Bounded second derivative of f0, available for each fixed sinusoidal fixture, does give `O_Lp(A^3)`. Neither pointwise coefficient cancellation nor finite-amplitude numerical residuals upgrade this distinction.

## 6. Exact next odd-rank ledger: the fifth cumulant still fails

For every odd n>=3, define

    Jn=a(M_(n+1)+2 beta N_n)
        +c(N_n+2 beta(M_(n-1)-M_(n+1)))+2b M_(n+1),
    Ln=a n M_(n-1)
        +c(2 beta M_(n-1)+(n-1)N_(n-2))
        +b(1+2(n-1)M_(n-1))
        +(n-1)(M_(n-1)+2 beta N_(n-2)).

For the sinusoidal fixture, put `sigma_n=(-1)^((n-1)/2)` and `tau=exp(-k^2/2)`. The predictor's degree-two quadratic-times-odd cumulant coefficient is

    sigma_n n q epsilon tau [k^n Ln-k^(n+2) Jn],

whereas the target is

    sigma_n q epsilon tau [2n k^n-k^(n+2)].

Matching rank n for every k therefore requires `Jn=1/n` and `Ln=2`. Rank three is exactly the calibration already enforced. The continuum rank-five values are instead

    J5=.22202903945726732, L5=1.9520385728793879.

For `q=1/4, epsilon=1/10, k=1`, this gives

    [A^2](kappa5(Y)-kappa5(X_A))=-.00530642048338454.

Independent literal-map differences divided by A^2 are approximately

    A=.04:  -.00558181381,
    A=.02:  -.00543876275,
    A=.01:  -.00537107366,
    A=.005: -.00533834389.

The exact cumulant-polynomial checker verifies ranks 3,5,7 for four frequency/shape choices. The fixture is globally strongly convex for small positive A and has all moments. This obstruction prevents interpreting the successful rank-three word as the whole second-order current. A fresh Gaussian bridge retains it with factor `t^5`; smoothing does not make the word identically zero.

## 7. Two-variable target-current identity

Let `X_r=rZ+sqrt(1-r^2)G` and `X_(r,t)=tX_r+sqrt(1-t^2)H`, with the displayed fresh independent roots. For a sufficiently integrable test phi and potential U1 with g=grad U1, the exact amplitude-two target coefficient is

    C2(phi)=integral_0^1 dr integral_0^1 dt E[
        Dg(X_r)g(X_(r,t)) dot grad phi(Z)
        + r Sym(g(X_r) tensor g(X_(r,t))):Hess phi(Z)].

Proof: start from

    C2(phi)=1/2 Cov(phi,(U1-EU1)^2)
           =integral dr E[grad phi(Z) dot g(X_r)(U1(X_r)-EU1)].

Represent `Z=rX_r+sqrt(1-r^2)G'` and apply the Gaussian covariance identity in X_r, with G' held independent. Differentiating `g(X_r) dot grad phi(Z)` gives exactly the Dg term and the r-weighted Hessian term. There is no additional t factor.

The exact symbolic Wick audit uses covariance matrix

    Cov(Z,X_r,X_(r,t)) = [[1,r,rt],[r,1,t],[rt,t,1]],

and tests Hermite interactions `(1,2),(2,3),(2,5),(3,4),(3,3),(4,4),(4,5),(5,5)` against carrier ranks through eight. All 52 checks pass. `bilinear_identity_checks.json` records each nonzero first-kernel and second-kernel polynomial and its exact integrated target, giving a reproducible two-variable moment ledger.

This is an analytical source identity. Its Dg term is not permission to insert an HVP-valued producer, and the rank-two current is not a global mean source. A positive original-VALUE implementation still must retain the correct common-root carrier, include both current terms, handle the two quadrature clocks and origins, and meet the actual consumer's radius/share/padding/numerical guards.

## 8. Imported conditional mean and covariance-action ports

The proof's conditional-u mean statement (8) is valid with its imported-domain conditions, and is narrower than a carrier-preserving replacement of the original (u,v) graph. To make the normalization concrete, take fixed positive shares `v_K=v_E=1/2`, use the executable anchored sources

    K0(v)=K(u,v)-K(u,0),
    E0(v)=E(u,v)-E(u,0),

and restore the two literal origins after their completed mean programs. Both are square D-dimensional sources in private v; u stays outside both independent banks. The full-gradient source K0 has radius at most `alpha beta` and anchored energy at most `C alpha sqrt(D)`. The source E0 has energy at most `C alpha^2(|u|+sqrt(D))` and the private curl established above.

Because both Kv* and Kv0 are positive semidefinite and at most `alpha beta I`, their difference has norm at most `alpha beta`. Thus sufficient normalized radii are

    rho_K=sqrt(2) alpha beta,
    ell_E<=sqrt(2) alpha beta[1+alpha(|a|/2+|c|beta)].

Require the imported actual gradient-compiler threshold `rho_K<=r_*(k,...)`, the near-gradient gap guard `ell_E<=1/4`, padding `mu=alpha` with `0<alpha<=1`, and all inherited finite clock/filter/numerical guards. The alpha=0 case is the separate source-zero branch and does not evaluate inverse-alpha expressions.

The already audited imported near-gradient mean bound is literally

    Lambda e_E {alpha(alpha+mu)+alpha^3(1+mu^(-1/2))}
       + absolute floors.

With `mu=alpha` this is `Lambda alpha^2 e_E`, hence at most `Lambda alpha^4(|u|+sqrt(D))`. An order-k gradient bank has error at most `Lambda alpha^k sqrt(D)`; `k>=4` fits the same profile. Independent completed-output convolution gives target covariance `v_K I+v_E I=I`. This proves the scoped mean-law comparison. It does not provide a strong estimate of the conditional mean, preserve the old private v, or preserve the actual random old K(u,v).

For the buffer allocation in (9), the conditional Gaussian reference variance is exactly

    t^2 I+(1-2t^2)I=(1-t^2)I.

It is positive for the stated `0<t<=1/2`. The raw predictor bridge instead has the additional conditional covariance `t^2 Cov_v(K*)` and higher descendants. The proof correctly leaves these as debts.

The covariance estimate (10) has a dimension-safe one-energy proof. Let `D=E-E_v E`, `Kc=K-E_v K`, and `K*c=K*-E_v K*`. Then

    Cov(K*)-Cov(K)=E[K*c D^T]+E[D Kc^T].

For any centered vectors A,D,

    ||E[A D^T]||HS <= ||Cov A||op^(1/2) ||D||2.

Gaussian Poincare and the actual private-v first give `||Cov K||op^(1/2)+||Cov K*||op^(1/2)<=C alpha`. The marked energy gives `||D||2<=C alpha^2(|u|+sqrt(D))`. Combining them proves exactly

    ||Cov_v(K*)-Cov_v(K)||HS<=C alpha^3(|u|+sqrt(D)),

without a second hidden dimension factor.

The legal covariance-action source is the executable anchored full gradient `K(u,v)-K(u,0)`, whose covariance equals `Cov_v K`. The word “centered” in a covariance theorem can refer to its analytical L2 energy; it must not be read as an original VALUE instruction to subtract the unknown mean. The imported forward action has mean calibration `Lambda ell e mu`, energy `Lambda ell e`, and first `Lambda ell^2/sqrt(mu)`, with its own full response roots, finite covariance clocks, buffer and numerical guards. Here `ell<=alpha beta` and `e<=C alpha sqrt(D)` by Poincare. An exact covariance oracle is not being substituted for that program.

This audit imports these already reviewed finite mean/covariance consumer theorems rather than re-proving their internal circuits. It checks the new source ports and their substitution at the actual dimension, normalization, origin, error profile and retained caller. It does not endorse a completed third-order reverse-law join.
