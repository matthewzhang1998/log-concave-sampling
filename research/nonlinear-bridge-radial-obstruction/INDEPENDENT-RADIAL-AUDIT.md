# Independent audit of the radial obstruction to the sealed stage-two graph

2026-10-05. This audit does not modify either sealed bridge packet. It analyzes the actual coherent radial source, the genuine nested OU-history target, every retained scalar-row correlation, and the finite outer rule's first Gaussian chaos.

## Verdict

The two-threshold radial source below gives a genuine obstruction. For the sealed finite VALUE graph, with its moment-exact single rule, its actual positive linear-pair rule, and quadratic-pair error `deltaQ`, the audited bound is

    ||Q_out psi_g - m3(g)||_L2(gamma_D)
      >= [c_* - 8 A - deltaQ/4] A^2 sqrt(D) - 7 A,

    c_* = (17/sqrt(2)-11)/192
        = 0.005316746250892225... .

It holds for `0<A<=1/36` and the smoothed source described below. In particular, `A=1/n`, `D=n^6`, `deltaQ<=A^2`, and `n>=10000` give a lower bound `A^2 sqrt(D)/250`. The stronger dimension-scaling choice `D=(3000n)^2` gives `A^2 sqrt(D)/500` for the same range of n.

This is a lower bound for this particular repaired mean graph, not an oracle lower bound for every possible VALUE method. Transfer to an implemented completed law additionally requires its own-mean completion and numerical floors to be `o(A^2 sqrt(D))`.

## 1. A coherent C2-potential source

Set

    alpha=beta=1/2,
    c1=alpha/2=1/4,
    cv=alpha(1-1/sqrt(2)),
    tau=cv/2,
    rho=1-A tau,
    eta=A^2.

The sharp radial amplitude is

    phi(r)=alpha(r-1/2)_+ + beta(r-rho)_+.

For the admissible source, replace both hinges by the same C1 function

    H_eta(t)=0                           if t<=-eta,
             (t+eta)^2/(4 eta)          if -eta<t<eta,
             t                           if t>=eta.

Write the resulting amplitude as `phi_eta`. Define

    g(x)=A sqrt(D) phi_eta(|x|/sqrt(D)) x/|x|,
    g(0)=0,
    U(x)=A D int_0^(|x|/sqrt(D)) phi_eta(s) ds.

Both thresholds exceed eta. Thus phi_eta vanishes near zero, is C1, and satisfies `0<=phi_eta'<=1` and `0<=phi_eta(r)<=r`. Consequently U is C2, g is its C1 gradient, and its radial/tangential Hessian eigenvalues are respectively

    A phi_eta'(r),   A phi_eta(r)/r,

both in `[0,A]`. No derivative bound beyond the stated Hessian interval is assumed or needed.

The hinge perturbation is at most eta/4, so the source-uniform difference from the sharp source is

    delta <= A eta sqrt(D)/4 = A^3 sqrt(D)/4.

The sharp source is used only in the proof, then this coherent-source perturbation is restored.

## 2. Uniform radial reduction, including growing finite node lists

Let `Y=sum_j q_j G_j` be any Gaussian scalar row, let `q=||q_vector||`, and let phi be either of the amplitudes above. The source's deterministic-row surrogate is

    g_row(Y)=A phi(q) Y/q,

with its zero-row value defined as zero. The function

    r -> phi(r)-[phi(q)/q] r

vanishes at q and has derivative in `[-1,1]`. Therefore

    ||g(Y)-g_row(Y)||_2
      <= A || |Y|-q sqrt(D) ||_2
      <= sqrt(2) A q.

The last inequality follows from `E(|G|-sqrt(D))^2<=2`. This argument does not assume that different retained rows are independent.

Every surrogate original-VALUE argument in the sealed graph has row norm below 2 for `A<=1/36`. For example, with `k=sqrt(3/8)` and `d0=1+1/sqrt(2)`,

    ||S||_row <= 1+A/sqrt(2)+A^2 k,
    ||d(U)||_row <= A d0,
    ||e(U,V)||_row <= A^2(1+k).

The largest terminal bound needed is

    ||S +/- d(U) +/- d(V)||_row
      <= 1+A(2+3/sqrt(2))+A^2 k < 2.

The scalar-row radial map itself is A-Lipschitz, so the same source-error propagation used for the sealed absolute VALUE floor applies in L2. Its coefficient is

    K_VALUE=7/2+14A+(11/2)A^2 <4.

Using the local bound `2 sqrt(2) A` at every surrogate query gives full graph reduction error below `8 sqrt(2) A<12A`. Positive normalized node weights and Minkowski make this independent of all three node counts, even when those counts grow with D or A. The same statement holds if the quadratic-pair sum is replaced by its exact positive clock integral.

### Genuine-history reduction

Let

    H1,a=int exp(-h) X_(a+h) dh,
    H1=H1,0,
    H2=int a exp(-a) X_a da.

The genuine retained histories have exactly the laws of v and w:

    v=x/2+N/2,
    w=x/4+N/2+M/4.

In particular, the reduction does not replace the actual history by independent samples. Put

    k0=A phi(1),
    q2=sqrt(1-k0+k0^2/2),
    ell=A phi(q2)/q2.

Then the actual first and second nested forces reduce to

    F1,a ~= k0 H1,a,
    F2 ~= ell(H1-k0 H2).

The first force has L2 replacement error at most `sqrt(2)A`. The second force has error at most `sqrt(2)A(1+A)`. Both Gaussian surrogate terminal radii are below 1. For the final row, this follows directly from

    ||x-ell v+ell k0 w||_row^2
      =1-ell+ell k0/2+ell^2/2
          -(3/4)ell^2 k0+(3/8)ell^2 k0^2 <1,

using `0<=ell,k0<=A`. Thus the complete genuine terminal target has reduction error at most

    sqrt(2) A [1+A(1+A)] <2A.

The graph and target together cost less than 14A before the outer first-chaos factor 1/2. This proves the `7A` term in the verdict.

## 3. Actual sharp-source descendants

At the stationary, v, and w radii, the sharp source has

    g_row(x)=A cU x,   cU=c1+beta A tau,
    g_row(v)=A cv v,
    g_row(g_row(w))=0.

The last identity is important: the inner g(w) row lies strictly below the lower threshold. Hence the resummed shift is retained literally and equals `S=x-A cv v` in the surrogate; it has not been replaced by a different construction.

For every actual one-time or pair row U,

    d(U)=A(cv v-c1 U)-A^2 beta tau U.

The scalar x coefficients of U and V are exactly `r=exp(-a)` and `s=exp(-(a+h))`. Their remaining coefficients, joint covariance, and correlations with v,w are the sealed genuine conditional rows. No independence assertion across finite quadrature nodes is used.

The linear-pair correction is negligible at the grade under examination:

    ||e(U,V)||_row <= A^2(1+sqrt(3/8)),
    |C_g(S,e(U,V))|_row <2 A^3.

This estimate is uniform in the actual pair times and requires no linear-pair quadrature approximation.

## 4. First-chaos expansion with explicit uniform remainder

Every S, single-stencil, and quadratic-stencil terminal row has form

    Y=x+A y+A^2 z,   ||y||,||z||<=1.

For the largest quadratic terminal these bounds follow from

    ||y|| <=3cv/sqrt(2)+2c1 <0.811,
    ||z|| <=2 beta tau <0.074.

If `n=<x,Y>` and `q=||Y||`, elementary norm identities give

    |q-1-A y_x|<=2A^2,
    |n/q-1|<=A^2.

On these rows the lower hinge is active. Also

    phi(q)<=1/4+A(1+A)+beta A tau <0.3.

Consequently the first-chaos scalar coefficient of the sharp terminal source is

    <x,g_row(Y)>
      = A c1+A^2[alpha y_x+beta(y_x+tau)_+] + R,
    |R|<=3A^3.

The term 3 is conservative: the displayed estimates give less than 2.3.

For the genuine target, write `ell=A c2`. Its inner and final terminal radii lie below the upper threshold; `c2=alpha(1-1/(2q2))`, and `|c2-c1|<0.064A`. The target row has the same form with

    y=-c1 v,
    z=-(c2-c1)v/A+c2 cU w,
    ||z||<0.084.

Its upper hinge is therefore absent at leading order. Its first-chaos scalar coefficient is

    A c1-A^2 alpha c1/2+O(A^3),

with the same explicit `3A^3` remainder.

The lower-branch part of the graph has exactly this A^2 coefficient. Its base contributes `-alpha cv/2`; its single rule contributes

    alpha E[cv/2-c1 r]=alpha(cv-c1)/2,

using only its exact first exponential moment `E r=1/2`. The quadratic combination kills the constant and affine lower-branch terms. Thus there is no unaccounted lower-branch A^2 error.

## 5. The nonzero upper-gate coefficient

Set

    u=cv/(2c1)=1-1/sqrt(2),
    a(r)=cv/2-c1 r=c1(u-r).

S is at the upper gate to first order. The leading upper-gate contribution of a single centered stencil is

    beta A^2 a(r)/2.

It is affine in r. Therefore the actual finite single rule, which need not have symmetric nodes, gives exactly

    beta c1 A^2 (u/2-1/4).

For the quadratic stencil, the leading contribution is

    (beta A^2/8)[|a(r)+a(s)|-|a(r)-a(s)|].

The ordered `Exp(2)(a) x Exp(1)(h)` clock law is the ordering of two independent Exp(1) times. Since the displayed integrand is symmetric, r and s may be integrated as independent Uniform(0,1) variables for this calculation. This is only a clock identity; it does not change the retained Gaussian genealogy.

For `0<=2u<=1`,

    E|r-s|=1/3,
    E|2u-r-s|=1-2u+(8/3)u^3.

Hence the quadratic contribution is

    beta c1 A^2 (1/12-u/4+u^3/3).

The total upper-gate discrepancy is

    beta c1 A^2 (u/4-1/6+u^3/3)
      = [(11-17/sqrt(2))/96] A^2
      = -0.01063349250178445... A^2.

The target has no compensating upper-gate term. This establishes the signed witness using the actual coherent source.

The base, single, and quadratic terminal absolute masses total 5/2. Their expansion error is at most `(15/2)A^3`; the target costs `3A^3`, and the linear-pair correction costs less than `2A^3`. Together this is less than `13A^3` in the inner scalar coefficient, uniformly over the allowed clock lists.

## 6. Finite pairs, smoothing, and the outer resolvent

Only the quadratic pair contribution must be compared with an exact clock integral. The sealed pair theorem applies directly to the real VALUE integrand, without a derivative qualification, and gives

    error <= (deltaQ/2) A^2 sqrt(D).

One can apply it first to the smooth source, then use coherent strong perturbation to pass to the sharp proof source. There is no appeal to the one-time operator rule for the nonlinear single stencil.

For two anchored A-Lipschitz sources at uniform distance delta, the sealed graph's strong coherent difference is at most `K_VALUE delta`. The genuine target difference is at most `(1+A+A^2)delta`. Their sum is below `5delta` for `A<=1/36`. This includes every nested changed source argument. With the eta above, the smoothing cost is at most `(5/4)A^3 sqrt(D)` before the outer factor 1/2.

Let `W=Q_out psi_g-m3(g)`. For any vector-valued Gaussian L2 function f,

    E <Z,R1 f(Z)> = (1/2) E <X,f(X)>.

The actual finite outer rule has the same exact first moment, and therefore the same identity. Thus

    ||W||_2 >= |E<Z,W(Z)>|/sqrt(D).

The negative inner witness survives with exactly half its coefficient. No estimate at an outer endpoint and no outer quadrature norm error is necessary. The accumulated A^3 coefficient is at most

    13/2+5/8=7.125<8.

The concentration cost is below 7A, and the finite quadratic-pair cost is `(deltaQ/4)A^2 sqrt(D)`. These are precisely the constants stated in the verdict.

## 7. Consequences and boundaries

For `A=1/n`, `D=(3000n)^2`, `n>=10000`, and `deltaQ<=A^2`, division by `A^2 sqrt(D)` leaves

    c_* -8/n -1/(4n^2)-7/3000 >1/500.

Consequently, a source-uniform monomial upper bound for this graph of the form

    C A^p D^(1/2+gamma),   p>2,

with C independent of A and D, requires

    gamma >= (p-2)/2.

In particular an `A^(8/3)` grade requires at least an additional `D^(1/3)` factor. The present smoothing certificate pays `D^(2/3)` beyond sqrt(D), so this audit does not prove its dimensional exponent sharp. It does prove that eliminating all dimensional loss at any improved power of A is impossible for the unchanged stencil over the unrestricted source class.

The first proposed radial family with a globally linear low branch, `A[alpha x+beta(|x|-R)_+x/|x|]`, does not produce this witness. In that case `cv=c1`, the leading single displacement is centered symmetrically, and the expected quadratic gate contribution cancels. The lower radial threshold is the feature that creates the coherent discrepancy `cv!=c1`.

## 8. Independent numerical check

`check_independent_radial_audit.py` executes the actual deterministic-row limit using the sealed conditional one-time row and the full two-by-two conditional pair covariance root. It evaluates every original nested source call, both centered corrections, the polarized correction, and the genuine-history reduced target. It verifies row variances and the true U,V covariance before evaluation.

`independent_radial_audit_checks.json` records the output. At A=0.0001 the sharp-source inner coefficient is approximately -0.0106325453 and the eta=A^2 smoothed coefficient approximately -0.0106201737, converging to the proven -0.0106334925. The separate scalar quadratic integral check differs from its exact value by about 7.2e-8. These numerical checks are diagnostics; the explicit inequalities above establish the result.

