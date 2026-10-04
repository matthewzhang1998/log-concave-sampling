# A two-stage VALUE repair of the exposed third-cumulant word

2026-10-04. Constructive bounded continuation of the independently audited same-root chord separator. Independent review pending. This matches a specific nonlinear cumulant row, not the whole second-order law current.

## 1. Executed positive predictor

Let f be an anchored original gradient VALUE source, f(0)=0, 0<=Df<=alpha I. At the same two Gaussian roots u,v and positive finite quadrature (w_i,r_i), q_i=sqrt(1-r_i^2), execute

    K=sum_i w_i f(r_i u+q_i v),
    F0=f(u),
    u*=u-a K-b F0, v*=v-c K,
    K*=sum_i w_i f(r_i u*+q_i v*),
    Y=u-K*.                                             (1)

The coefficients a,c,b below are fixed known scalars computed from the finite rule. Some are negative; (1) is still an ordinary positive pushforward. Every force is an original VALUE at its literal site. There is no Hessian-valued leaf, tensor coefficient query, conditional-mean oracle, or signed probability.

Set

    beta=sum w_i q_i, C=1/4+beta^2,
    m2=sum w_i r_i^2, m4=sum w_i r_i^4,
    m11=sum w_i r_i q_i, m31=sum w_i r_i^3 q_i,
    K0=m2+2beta m11,
    Ju=m4+2beta m31, Jv=m31+2beta(m2-m4), J0=2m4,
    Lu=3m2, Lv=2(m11+beta m2), L0=1+4m2.

Solve the known three-by-three system

    [ 1/2   beta   1  ] [a]   [1-C       ]
    [ Ju    Jv     J0 ] [c] = [1/3       ]                (2)
    [ Lu    Lv     L0 ] [b]   [2-2K0      ].

This small scalar coefficient solve is known arithmetic. Enforce a fixed deterministic determinant gap, e.g. |det|>=.1, and |a|+|b|+|c|<=3, before invoking this version. The fine admitted dyadic rules satisfy these guards; no source is queried to verify them. At the continuum rule the coefficients are

    a=1.380985765..., c=-0.320111621..., b=-0.305928078...,
    det=-0.1297604079... .

## 2. What is matched exactly

### Full matrix quadratic covariance word

For f(x)=Bx, 0<=B<=alpha I, put e=a/2+c beta. The actual map (1) is

    Y=[I-B/2+(e+b)B^2/2]u+[-beta B+e beta B^2]v.

Its covariance through degree two is

    I-B+[C+e+b]B^2=I-B+B^2,

by the first row of (2). Its full law is therefore O(alpha^3 sqrt(D))-close to the exact quadratic target, with a fixed covariance gap for small alpha. This is a full matrix identity including the common-v descendants.

### Complete third-cumulant row for quadratic plus odd perturbation

Take, first, the smooth scalar family

    U_A(x)=A[q x^2+epsilon(sin(kx)-kx)],
    f_A(x)=A[2qx+epsilon k(cos(kx)-1)],                  (3)

with q>0 and epsilon k^2 small enough for the stipulated Hessian sandwich. Write h=K/A, h_u=D_u h, h_v=D_v h, f0=f_A/A. The actual graph expands as

    Y=u-Ah+A^2[a h_u h+c h_v h+b h_u f0(u)]
                    +O_Lp(A^3),                        (4)

for this fixed smooth family. Let tau=exp(-k^2/2). The three centered feedback contractions are

    E[(u^2-1)h_u h]=q epsilon tau(Ju k^5-Lu k^3),
    E[(u^2-1)h_v h]=q epsilon tau(Jv k^5-Lv k^3),
    E[(u^2-1)h_u f0(u)]=q epsilon tau(J0 k^5-L0 k^3).

The unchanged rank-two/centering terms contribute -6q epsilon tau K0 k^3. Consequently

    [A^2]kappa_3(Y)
       =3q epsilon tau[(aJu+cJv+bJ0)k^5
                  -(aLu+cLv+bL0+2K0)k^3]
       =q epsilon tau(k^5-6k^3),                        (5)

which is exactly the target coefficient. The first-order coefficient is 3epsilon k^3 tau m2, versus target epsilon k^3 tau. Keep its actual finite quadrature floor; it is not exactly zero merely because floating point rounds m2 to 1/3. With delta_Q<=A^3 it is O(A^4).

The extension to any admissible C2 odd potential psi does not require a Fourier-completeness argument. Put g=psi', c2=E[H2(N)g(N)], c4=E[H4(N)g(N)]. The predictor's quadratic-times-odd **third-cumulant row** is

    3q epsilon[(aJu+cJv+bJ0)c4
                    +(aLu+cLv+bL0+2K0)c2],

and the target row is q epsilon(c4+6c2). Gaussian integration by parts transfers the Hermite derivatives onto the Gaussian; only the bounded first g' is required. Thus rows two and three of (2) match these coefficients for general admissible odd psi, not only a finite frequency list. This does not say every Hermite rank or the full quadratic-times-odd current is matched.

For fixed smooth sine or finite sine combinations, the residual third cumulant is O(A^3)+the finite quadrature floor. Under C2 alone, amplitude differentiability and dominated convergence give cancellation through A^2 with an o(A^2) remainder for each fixed potential, not a uniform O(A^3) theorem over the C2 class. No Hessian-continuity rate has been invented.

## 3. Literal conditional-v mean port

There is a valid caller-retained mean construction at this particular source interface. Freeze u, every outer caller and all finite versions before the private v bank. Write

    E(u,v)=K*(u,v)-K(u,v).

For each K occurrence define the symmetric first blocks Ku=D_uK, Kv=D_vK. Their norms obey ||Ku||<=alpha/2 and ||Kv||<=alpha beta. At the shifted point use Ku*,Kv*. Direct differentiation gives

    D_v K*=Kv*-(aKu*+cKv*)Kv,
    Curl_v E=-(aKu*+cKv*)Kv+Kv(aKu*+cKv*),
    ||Curl_v E||<=2(|a|/2+|c|beta)beta alpha^2.           (6)

The difference Kv*-Kv can have norm O(alpha), but is symmetric. Thus E has private first O(alpha), small curl O(alpha^2), and the actual marked energy

    |E|<=alpha[sqrt(a^2+c^2)|K|+|b||f(u)|],
    ||E(u,v)||_(Lp(v))<=C_p alpha^2(|u|+sqrt(D)).         (7)

Its fixed-u origin E(u,0) need not vanish. It is an executable caller-only record: evaluate the two packets at v=0 and the single site f(u), with complete keys. Charge those 2n+1 VALUES if they were not already captured. Its anchored energy is bounded by the displayed profile together with the literal origin, and its u first is O(alpha), not O(alpha^2). At u=0 and v=0 all values are exactly zero under original-anchor reuse.

K(u,v) is a genuine gradient in v, since D_vK is symmetric. Split its known origin K(u,0) off exactly. Use fixed positive shares vK=vE=1/2 and run independent COMPLETE mean banks for K and E, conditional on the retained u. The normalized gradient radius is rhoK=sqrt(2)alpha beta; impose its literal rhoK<=r_*(k,...) guard. The normalized E radius is at most ellE=sqrt(2)alpha beta[1+alpha(|a|/2+|c|beta)]; impose ellE<=1/4 and all inherited padding, finite-clock, caller and numerical guards. These are actual-radius conditions, not consequences of alpha<=1/2 alone. The former gradient bank can have arbitrary fixed target order; the E bank with padding mu=alpha has its literal O(alpha^2 e_E) law error. Their sum yields

    W2(Law(M(u)|u),N(E_v K*(u,v),I))
       <=Lambda alpha^4(|u|+sqrt(D))+absolute floors,    (8)

for a gradient-bank order at least four, and ordinary small-alpha/public-log guards. This is a Gaussian **mean law**, not a strongly accurate mean, and its private banks have been integrated. The whole captured u is legitimately retained because it was outside both banks from the outset.

For the r=0 reverse bridge with 0<t<=1/2, one positive buffer allocation is

    t u-t M(u)+sqrt(1-2t^2)N.                           (9)

Conditionally on u, its Gaussian reference has mean t[u-E_vK*] and variance 1-t^2. This is an executed allowed use of the complete output M and an independent fresh reserve. It is **not** yet the same law as t[u-K*(u,v)]+sqrt(1-t^2)N: the latter includes conditional covariance t^2 Cov_v(K*|u) and all further private-v descendants. They remain explicit debts.

The leading covariance can be approximated by the genuine-gradient source K:

    ||Cov_v(K*)-Cov_v(K)||HS
       <=C alpha^3(|u|+sqrt(D)),                        (10)

using the operator covariance bound O(alpha^2) for both K and K* and the one-marked-energy estimate (7). An admitted full-gradient forward covariance action on the executable anchored source K(u,v)-K(u,0), whose covariance is unchanged, can therefore supply this leading covariance at the corresponding one-energy scale, with all its own positive reserve, padding, caller, and replay costs. Analytical centered energy is used only in its estimate; E K is not queried. This statement identifies a legal covariance source; it does not silently Gaussianize K* or remove its higher conditional cumulants. A complete positive join at an asserted third-order W2 rate still needs its actual conditional current estimate.

## 4. Every offspring and complete finite cost

One raw occurrence of (1) costs 2n+1 private original VALUES, plus the common captured conditional mode and g anchor. There are no new Gaussian roots beyond u,v. Reusing f(u) is allowed only under the identical source/caller/version key. K* reuses the exact old K and f(u); it does not replay them independently.

At fixed u for (8), the complete caller-only origin subgraph costs up to another 2n+1 sites, shared only under that exact retained-u key. Every changed private source input inside a mean or covariance compiler executes a fresh complete raw graph. If N_K and N_E are the actual complete mean-source occurrence counts, a safe VALUE bill is

    Q_captured + n N_K +(2n+1)N_E + Q_known/numerical,

where Q_captured includes all original-mode, conditional-mode, K(u,0), E(u,0), and f(u) records actually needed. Covariance actions add their full occurrence counts. The total private tape is the complete bank/clock/fill tape of those consumers, not merely the two roots of (1).

The actual full private first of K* is bounded by

    alpha[1+(sqrt(a^2+c^2)+|b|)alpha].

The chord E has first O(alpha) and energy O(alpha^2 sqrt(D)) when u is again Gaussian. Under the original conditional scaling f(y)=s[g(x+s y)-g(x)], alpha=As^2, every fixed caller-a first is O(As) before the final state readout s, with the captured finite-mode derivative retained. The two stage shifts only multiply that by known O(1+alpha) factors. Original physical-caller and numerical profiles propagate through those actual same-site records. No derivative of a law error is used.

First/adjoint sweeps use original HVPs only at the 2n+1 recorded VALUE sites and through the old K/F0 feedback. A discarded primal graph incurs its full replay. All original sites, caller origins, coefficients and quadrature versions are frozen before differentiation. Numerical errors propagate with absolute quadrature mass one and known bounded predictor coefficients. No inverse actual energy is used. The mean consumers retain their separate inverse shares, finite filters, and absolute precision floors.

## 5. The next higher row and the full-current target

The construction has repaired the specific third-cumulant row that broke the single chord. It has not removed every second-order descendant. For odd n>=3, define

    M_j=sum w r^j, N_j=sum w r^j q, sign_n=(-1)^((n-1)/2),
    J_n=a(M_(n+1)+2beta N_n)
           +c(N_n+2beta[M_(n-1)-M_(n+1)])+2bM_(n+1),
    L_n=a nM_(n-1)+c[2beta M_(n-1)+(n-1)N_(n-2)]
           +b[1+2(n-1)M_(n-1)]
           +(n-1)[M_(n-1)+2beta N_(n-2)].

The predictor's quadratic-times-odd order-A^2 cumulant coefficient is

    sign_n n q epsilon tau[k^n L_n-k^(n+2)J_n],

whereas the true target is

    sign_n q epsilon tau[2n k^n-k^(n+2)].

Hence full matching of this family would require L_n=2 and J_n=1/n for every odd n, not just n=3. In the continuum rule the present n=5 values are J_5=0.222029039457..., L_5=1.952038572879..., giving a still nonzero A^2 defect -0.005306420483... at q=.25,epsilon=.1,k=1. This is an explicit offspring ledger, not a claim that fitting more finite rows proves a full-law order.

The companion two-clock current note identifies the exact rank-one/rank-two target for **all** bilinear interactions, with literal nested Gaussian correlation. That is the next finite positive compiler target. Formal row matching alone does not supply its uniform C2 remainder or its finite clock norm estimate.

## 6. Diagnostics

check_two_stage_word_repair.py reports 213 passing assertions: two quadratic coefficients q, four nonlinear frequencies, four small amplitudes, a mixed-frequency odd perturbation, non-diagonal matrix VALUE identities, conditional-v curl with noncommuting Hessians, and literal source zeros. The 241-node finite rule gives (a,c,b)=(1.380985764564835,-0.32011162098909224,-0.3059280782101048). These finite tests support the displayed analytical word identities; they do not establish a general full-law order-three theorem.
