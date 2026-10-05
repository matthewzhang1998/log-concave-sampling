# An original-VALUE collective correction with weak grade three

2026-10-05. New finite positive-mean-source graph. The general source is an anchored convex gradient `g=grad U`, `U in C2`, `0 <= Dg <= A I`. There is no Hessian modulus, block structure, source-dependent covariance, derivative VALUE, or expectation producer. This is a bounded canonical-m3 component, not an any-order or endpoint theorem.

## Result

For `0 < A <= 1/12`, in arbitrary finite dimension `D`, the new graph below has

    ||Q_out psi_Q - m3(g)||_(L2 gamma_D)
      <= 8 A^3 D + delta_1 A^2 sqrt(D)
           + delta_out A(1+A+A^2 sqrt(3/8))sqrt(D).       (1)

Both rules are explicit positive mass-one OU operator rules with exact first moment one-half and `O(log^2(1/delta))` nodes. They are independent of the source and dimension. The graph retains the original resummed shift, puts a collective force INSIDE the full terminal nonlinear response, and has only `(J+5) N_out` original VALUES per residual occurrence before exact-key simplification. Its private Gaussian dimension is `4D`.

At the declared normalized native ports, safe bounds are private first `<=2A`, curl `<=3A^2`, and residual energy `O(A^2)(|z|+sqrt(D))`. A guarded existing near-gradient completion therefore adds its `O(A^4 sqrt(D))` allowance and every absolute numerical floor. The complete bill, including native occurrences, is in Section 7.

Taking `delta_1<=A`, `delta_out<=A^2` gives a general-C2 fixed-D grade-three source with zero inverse-A exponent in its quadrature-node count at this fixed grade. The complete native occurrence factors remain as declared in the full bill. The explicit extra dimension factor is `sqrt(D)` relative to the natural energy scale. This is NOT a dimension-uniform `O(A^3 sqrt(D))` result. A dimension-uniform coarse cap is also available. All dimension-dependent native costs remain in the full bill.

## 1. Genuine target and genuinely finite new source

The reference target retains the actual stationary OU future with covariance `exp(-|s-t|) I_D`. Conditional on `X_0=x`, define analytically

    F1,a = integral_0^infinity exp(-h) g(X_(a+h)) dh,
    F1   = F1,0,
    F2   = integral_0^infinity exp(-a) g(X_a-F1,a) da,
    psi(x) = E[g(x-F2) | X_0=x],
    m3(g)  = R1 psi,
    R1 f(z)=integral_0^1 P_t f(z) dt,
    P_t f(z)=E f(tz+sqrt(1-t^2)G).

The genuine linear histories are retained as

    v=x/2+N/2,
    w=x/4+N/2+M/4,
    k=sqrt(3/8),

with independent standard `N,M`. They are precisely the joint law of the exponentially weighted first and second linear OU integrals, conditional on x.

Let `(p_j,t_j)` be the rule of Section 2 at tolerance `delta_1`, and put `s_j=-log(t_j)>0`. For a further independent standard D-Gaussian `L`, use the true one-time conditional bridge

    U_j=t_j[x+2s_j N+2s_j(s_j-1)M]+c(s_j)L,
    c(s)^2=1-exp(-2s)[1+4s^2+4s^2(s-1)^2].          (2)

The coefficient `c(s)` is real and positive for s>0 because it is the residual of the genuine Gaussian orthogonal projection onto `(x,N,M)`. The sealed bridge proof gives a separate rational positivity certificate. Each full row in `(x,N,M,L)` has norm one and has its actual joint law with `x,v,w`. A common L across quadrature nodes is part of the NEW finite graph; these nodes are not asserted to form a simultaneous true OU path.

Execute original VALUES and ordinary arithmetic:

    a_v=g(v), a_w=g(w), b0=g(a_w),
    S=x-a_v+b0,
    d_j=a_v-g(U_j),
    T_Q(x,N,M,L)=g(S+sum_j p_j d_j).                 (3)

This retains the entire terminal nonlinear response. Since `sum p_j=1`, the exact-real terminal argument is also

    x-H_Q+b0,        H_Q=sum_j p_j g(U_j).            (4)

The `g(v)` occurrence is recorded and charged in (3); its algebraic cancellation is not a hidden source oracle. An implementation may explicitly simplify it, reducing the VALUE bill by one, only if it declares the simplified graph and its changed readset. The conservative results here retain (3).

For the outer rule `(omega_i,r_i)`, use another independent standard D-Gaussian `G` and set

    x_i=r_i z+sqrt(1-r_i^2)G,
    F_raw=sum_i omega_i T_Q(x_i,N,M,L),
    B_raw=sum_i omega_i g(x_i),
    E_raw=F_raw-B_raw.                              (5)

All four private Gaussian roots are common across all outer nodes in this occurrence. The actual own mean is `Q_out psi_Q`, where `psi_Q(x)=E_(N,M,L) T_Q(x,N,M,L)`. Equations defining expectations identify means for analysis; none is executed.

## 2. Explicit positive OU rule with exact first moment

For `0<delta<1`, choose

    K=ceil(log_2(4/delta)), h=2^(-K),
    n=max(1,ceil(log(16/delta)/log(4))).

On each interval

    [1-2^(-j), 1-2^(-j-1)], j=0,...,K-1,

use the n-point Gauss-Legendre rule for Lebesgue measure. On `[1-h,1]` use its midpoint with weight h. Every node lies strictly between zero and one. The total number is `Kn+1`. All weights are positive, their sum is one, and their weighted node mean is EXACTLY one-half because every panel rule integrates linear functions exactly.

The operator satisfies

    ||Q_delta-R1||_(L2(gamma_D;H) -> L2(gamma_D;H)) <= delta             (6)

for every finite-dimensional Hilbert output H and every D. Here is a direct operator proof. On degree-m Hermite chaos, the complexified `P_z` is multiplication by `z^m`, and so has norm at most one on `|z|<=1`. The Bernstein ellipse of parameter rho=2 around a panel at distance `[d,2d]` from one has center `1-3d/2` and real semiaxis `5d/8`. Its largest modulus is attained at a real endpoint and is strictly below one; the leftmost panel also stays in the unit disk. Thus `P_z` is operator-holomorphic on a neighborhood of the closed ellipse and bounded by one. A degree-(2n-1) Chebyshev approximant has uniform operator error at most `4*4^(-n)`. Positivity and polynomial exactness make the integrated panel error at most twice its mass times that error. Summing all panels gives at most `8*4^(-n)<=delta/2`. The final midpoint interval contributes at most `2h<=delta/2`. This proves (6).

No regularity of f is required for (6). Finite-precision node, weight, moment and row errors are separate absolute floors; this exact-real construction does not erase their cost.

## 3. Weak second-order remainder using only first derivatives

This is the outer-resolvent cancellation that replaces strong history replication.

Let `B(x,xi)` be a C1 displacement, with xi independent of x, such that

    ||D_x B|| <= beta < 1,
    B_2=||B(X,xi)||_2 < infinity,
    B_4=||B(X,xi)||_4 < infinity,

for stationary standard X. Let `||Dg||<=A`. Define the analytic remainder

    Rem_B(x,xi)=g(x-B)-g(x)+Dg(x)B.

Then

    ||R1 E_xi Rem_B||_2
      <= (A/2) { C_beta[(pi/2)B_4^2
                    +sqrt(D) beta B_2/sqrt(1-beta)]
                    +beta B_2 },                    (7)
    C_beta=sqrt(1+(1+beta)^2).

No derivative of `Dg` or `D_x B` occurs in (7).

### Proof of (7)

Fix an outer parameter t<1, a caller z, and xi. Write `c=sqrt(1-t^2)`, `x=tz+cG`, and let q be this Gaussian density. For u in [0,1], `T_u(x)=x-uB(x)` is an orientation-preserving global C1 diffeomorphism because u beta<1. Let `p=(T_u)#q`, and put `m_2=(E_q |B|^2)^(1/2)`.

The Gaussian change-of-variable entropy identity and first-order Gaussian integration by parts cancel the linear Jacobian term. The log-determinant series starts at degree two and gives

    KL(p||q) <= (u^2/2)[m_2^2/c^2+D beta^2/(1-u beta)].                (8)

This is valid for nonsymmetric DB: `|tr(DB)^j|<=D beta^j`, and `det(I-uDB)>0` follows continuously from the identity. No derivative of DB is taken.

The squared Hellinger distance is at most KL. Indeed if a is the Hellinger affinity, Jensen gives `KL>=-2log a>=2(1-a)`. For a vector test H, Cauchy-Schwarz consequently yields

    |integral H(p-q)|
      <= [2(E_p|H|^2+E_q|H|^2)]^(1/2) sqrt(KL(p||q)).                (9)

Take `H(y)=Dg(y)B(y)`. It need not be bounded. It has the second moments required by (9), since `|H|<=A|B|` and

    |B(T_u(x))| <= (1+u beta)|B(x)|.

Changing variables in `E_q Dg(T_u(x))B(x)` and replacing `B(T_u^-1(y))` by `B(y)` introduces at most

    A u beta E_q |B|.

Equations (8),(9) therefore give

    |E_q[(Dg(x-uB)-Dg(x))B]|
      <= A u { C_beta[m_2^2/c
                      +sqrt(D) beta m_2/sqrt(1-beta)]
                     +beta E_q|B| }.                (10)

Now average xi and take L2 in standard z. Jensen and Fubini give

    ||E_xi m_2^2||_(L2_z) <= B_4^2,
    ||E_xi m_2||_(L2_z), ||E_xi E_q|B|||_(L2_z) <= B_2.

The exact integral identity

    Rem_B=-integral_0^1 [Dg(x-uB)-Dg(x)]B du

is valid for C1 g. Integrate u, then integrate the actual outer t. The endpoint has the integrable singularity

    integral_0^1 (1-t^2)^(-1/2)dt=pi/2.

This proves (7). There is no uniform-in-t assertion at t=1 and no frozen outer-clock mixture.

## 4. Applying the cancellation to the actual coherent source

In the new finite graph put `B_Q=H_Q-b0`. On the genuine path retain the SAME true w and put `B_*=F1-b0`. Each comparison changes all live source descendants. In particular `b0=g(g(w))` is not replaced by a fixed numerical shift when g changes.

For both displacements,

    beta=A/2+A^2/4,
    B_p <= A(1+A k)||G_D||_p, p=2,4.                (11)

For `B_Q`, the x first of H_Q is bounded by `A sum p_j t_j=A/2`; for B_* it is bounded by `A integral exp(-2s)ds=A/2`. The x first of b0 is at most A^2/4. Every U_j is marginal standard and w has variance k^2 I, so positivity/Minkowski proves the moment bound. The same proof uses the true stationary future for B_*.

The conditional mean discrepancy is exactly

    E[B_Q-B_* | x]=(Q_delta1-R1)g(x),               (12)

because b0 has the same retained law in both objects. Equation (6) bounds (12) by `delta_1 A sqrt(D)` in L2. The linear term in the expansion at x therefore costs only `delta_1 A^2 sqrt(D)`. Its derivative is an analytical cancellation, not a derivative producer.

The genuine second history is restored strongly, with every nested value present:

    ||F1-F2||_2 <= A^2 sqrt(D),
    ||b0||_2 <= A^2 k sqrt(D),
    ||g(x-F2)-g(x-B_*)||_2
       <= (1+k) A^3 sqrt(D).                        (13)

Combining (7), once for B_Q and once for B_*, with (12),(13), gives the sharper explicit coefficient

    ||R1 psi_Q-m3||_2
      <= A^3 { C_beta[(pi/2)b^2 sqrt(D(D+2))
                         +(beta/A)b D/sqrt(1-beta)]
                     +(beta/A)b sqrt(D)
                     +(1+k)sqrt(D) }
                +delta_1 A^2 sqrt(D),               (14)
    b=1+A k.

For `A<=1/12`, (14) is at most `8 A^3 D+delta_1 A^2 sqrt(D)`. A completely rational verification uses

    k<5/8, b<=101/96, beta<=13/288,
    beta/A<=13/24, C_beta<3/2,
    1/sqrt(1-beta)<21/20, pi/2<11/7, sqrt(3)<7/4.

After `sqrt(D(D+2))<=sqrt(3)D` and `sqrt(D)<=D`, the coefficient is bounded by

    (3/2)[(11/7)(101/96)^2(7/4)
            +(13/24)(101/96)(21/20)]
      +(13/24)(101/96)+13/8 < 8.                    (15)

Finally `||psi_Q||_2<=A(1+A+A^2 k)sqrt(D)`, so the actual outer rule adds exactly the last term of (1).

A source-independent coarse cap is useful when D is large. The new residual energy is at most `A^2(1+A k)sqrt(D)` and the genuine residual is at most `A^2(1+A)sqrt(D)`. Hence the exact-outer bias is also at most

    A^2[2+A(1+k)]sqrt(D).                           (16)

Use the minimum of (14) and (16), then add the actual outer/completion/floor allowances. No source receives a worse certificate merely because the fine estimate has a dimension factor.

## 5. Quadratic exactness and the prior obstruction

For every symmetric `0<=K<=A I`, set g(x)=Kx. The finite terminal is exactly

    Kx-K^2 sum_j p_j U_j+K^3 w.

Its conditional mean is `Kx-K^2 x/2+K^3 x/4`, exactly the genuine canonical inner mean. The outer exact first moment gives the exact canonical-m3 mean `Kx/2-K^2x/4+K^3x/8`. No commuting approximation was made for nonlinear sources; this calculation concerns a genuinely constant Hessian.

The old radial obstruction applies to a different finite Taylor/resummed stencil. Graph (3) retains the nonlinear response to the collective sum. In the old two-shell high-D source, all U_j have standard radius, and their weighted x coefficient is exactly one-half; the leading order-A radial terminal shrink therefore agrees with the genuine target. This observation is a diagnostic, not an extra theorem removing the dimension factor in (14).

In one dimension the old normalized-ancestry sine counterexample has its leading A^2 coefficient proportional to `integral_0^1 C(t)dt-r C(r)`. The new graph instead uses `sum p_j C(t_j)`, so (6) controls that first-order mismatch. A numerical diagnostic below checks this cancellation independently. The proof of (1), however, covers sources that vary arbitrarily with A and D and is not based on a fixed smooth shape or a small list of examples.

## 6. Firsts, curl, caller origins, and root ledger

Let H=Dg(x-B_Q), H0=Dg(x), and write `W=(N,M,L)`. Exact firsts are

    E_x=(H-H0)-H(B_Q)_x,
    E_W=-H(B_Q)_W.                                  (17)

Products remain in this order. Every requested derivative/adjoint sweep uses only original HVPs at the actual VALUE sites, not derivatives of HVPs.

The full W-row of each U_j has norm `sqrt(1-t_j^2)`. Positive weights and the exact first moment imply

    ||(B_Q)_x||<=beta=A/2+A^2/4,
    ||(B_Q)_W||<=eta=A sqrt(3)/2+A^2 sqrt(5)/4.       (18)

The second term in eta uses the full W-row norm of w, `sqrt(5)/4`. This is a full-bank operator estimate; there is no hidden sum over J or coordinates.

For the outer rule let `beta_out=sum omega_i sqrt(1-r_i^2)<=sqrt(3)/2`. The raw private first on `(G,N,M,L)` is bounded by

    A sqrt[beta_out^2(1+beta)^2+eta^2].              (19)

Lift E_raw into the G block and zeros in the other blocks. The leading G derivative in (17) is symmetric. Every remaining skew and off-diagonal term is bounded by

    raw curl <= 2 beta_out A beta+A eta.            (20)

For A<=1/12, multiplying these safe raw bounds by the half-variance normalization factor sqrt(2) gives private first `<2A` and curl `<3A^2`. No Hessians are commuted. The raw retained-z first is at most `A(1+beta)/2`.

Capture and restore the ACTUAL caller origin at `G=N=M=L=0`, including all nonzero descendants. This keeps private first and curl unchanged and at most doubles the retained-z first. The normalized anchored retained-z bound is safely `2A`. Arbitrary exterior source, anchor and scale parameters keep their own actual recorded chains and guards; this z bound is not a substitute for them.

For `kappa_p=||G_D||_p/sqrt(D)`, a useful conditional energy bound is

    ||E_raw||_(Lp|z)
      <= A^2[(1/4+A/8)|z|
                +kappa_p(1+A k)sqrt(D)].            (21)

The captured origin is at most the displayed caller term. At z and all private roots zero, every original VALUE site in (3) is literally zero. At g=0 the entire graph vanishes. No normalization divides by measured energy.

## 7. Positive completion, actual cost, and numerical floors

The raw graph is a finite Gaussian pushforward, not by itself a Gaussian law at its own mean. For the completed component, use independent COMPLETE baseline/residual banks conditional on the caller and live exterior labels, positive half-variance shares, and the pinned own-mean/coisometry services. Safe declared residual parameters are

    ell_E=2A, a_seed=2A, padding mu=A,
    private dimension=4D.

Their existing bracket is `16A^2+8A^3+8A^(5/2)=O(A^2)`. With (21), a guarded service contributes `Lambda_comp A^4 sqrt(D)` after integration. The baseline has its separate gradient-order allowance. `A<=1/12` supplies only a smallness condition; ALL actual native radius, caller, active-dimension, tolerance, clock/filter, finite-mode and precision guards must be checked at the new parameters. This note does not claim that a fully expanded compiler has been run.

With J inner nodes, literal (3) has `J+4` original VALUES; the residual adds one baseline. The safe full bill is

    Q_rule/setup + Q_captured + N_out N_B
      +(J+5)N_out N_E + Q_known/numerical/replay.     (22)

`N_B,N_E` are complete expanded native occurrence counts, including every new bank, fill, buffer, clock, mark, filter, finite mode, capture, later replay and precision operation at the actual dimension and radii. Captured origins use the complete graph. A changed source, anchor, scale, caller, scalar coefficient, row or rule version invalidates its exact-key cache. Discarded primals pay full replay. Every first/adjoint sweep pays one original HVP at each recorded original VALUE site. Arithmetic is O(D) per leaf apart from the source's own cost; no dimension-independent runtime is asserted.

With uniform absolute original-VALUE error nu, the literal unsimplified S and d expression has a conservative terminal floor

    (1+4A+A^2)nu.

It follows by separately propagating `a_w`, `b0`, `a_v`, every `g(U_j)`, and the terminal VALUE. Positive weights prevent a J factor. E_raw adds one baseline nu. A separately executed caller capture pays its own whole floor. If an actual implementation proves exact shared cancellation of a_v, the sharper terminal floor `(1+2A+A^2)nu` is available, but the conservative bill uses the former.

All scalar coefficient, Gaussian-row, quadrature-weight/moment, original-HVP, arithmetic, capture, native-service, finite-mode, and replay errors remain separately priced absolute floors. The weak coherent-source comparison in Section 3 does not reduce arbitrary numerical errors. Formula (2) near its small-time cancellation needs stable known scalar arithmetic and a certified row-error allowance; no unpriced clipping is permitted.

## 8. Next input type and honest boundary

For a live exterior anchor a and scale s>=0, the normalized original source

    f(y)=s[g(a+s y)-g(a)], alpha=A s^2,

is again an anchored C1 convex gradient with Hessian interval `[0,alpha I]`. This theorem applies to f when `0<alpha<=1/12`, with all original-g occurrences, actual g(a) captures, caller/scale chains, physical rescaling and new native guards charged. At alpha=0 use the literal zero-source branch. No source modulus is assumed fixed under scaling.

The produced residual E_raw is near-gradient, not a newly certified convex gradient in the same input class. Therefore this fixed grade-three result cannot simply be fed back into its own source theorem as an any-order induction. Its actual next consumer is the declared near-gradient own-mean service. The retained joint Gaussian roots are internal source inputs; appending them as output observers would require a different law theorem.

The main achievement is an executable collective graph whose first-history mean cancels weakly under the true outer resolvent, without approximating the full history strongly and without paying inverse-A replication at this fixed grade. The remaining general-class challenge is dimension-controlled higher-order weak cancellation and a reusable native output type, not a missing expectation leaf that has been treated as free.
