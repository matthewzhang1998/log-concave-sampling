# A radial shell refutes the unshifted single-history covariance formula

2026-10-05. Genuine same-g counterexample, with a true continuous OU source and the actual finite common-H cheap source. Independent exact-text review requested. This is distinct from the earlier port-only trace counterexample.

## Result

The proposed uniform expansion

    R1(j_H-j_Q)
      ?= (1/2)R1[(C_H-C_Q):D2g(x)] + O(Lambda A^4 sqrt(D))       (1)

is FALSE in the admitted Hessian class, even for C-infinity g, positive accurate finite clocks, and the literal same original g in every ancestor and covariance vertex. Here I_H is the true conditional OU integral, I_Q is the finite shared-root packet, j_a(x)=E g(x-I_a), and C_a(x)=Cov(I_a|x).

We construct A=D^(-1/2) and 0<=Dg<=A I such that the L2(gamma_D) discrepancy in (1) is at least c_* A^2, whereas A^4 sqrt(D)=A^3. Every fixed public-log multiplier leaves the separation intact.

The force has a coherent radial mean displacement of order one. The high-dimensional covariance trace probes the radial shell at that shifted radius; evaluating D2g at x loses the shift. A centered evaluation at x-E I_H matches the leading term for THIS family. No general centered/resummed fourth-order theorem or finite consumer is claimed.

## 1. A literal smooth convex-gradient family

Let D>=16, R=sqrt(D), A=1/R, and c=d=1/4. Define a normalized smooth bump

    eta(s)=exp(-1/(1-s^2/4)) for |s|<2, and eta(s)=0 otherwise,
    Z_eta=int_R eta(s) ds,
    psi(s)=Z_eta^(-1) int_(-infinity)^s eta(u)du.

Thus psi is nondecreasing and C-infinity, zero for s<=-2 and one for s>=2. Its derivative is an even, nonzero bump, strictly decreasing on (0,2). All derivatives have fixed bounds independent of D. Moreover

    ||psi'||infinity <= e^(-1)/(2e^(-4/3))=(1/2)e^(1/3)<1,

because eta>=e^(-4/3) on [-1,1]. For x!=0 put r=|x|, n=x/r, P=I-nn*, and

    f_D(x)=psi(r-R)n,
    g_D(x)=c A x+d A f_D(x).                                      (2)

Set f_D(0)=0. It vanishes on the entire ball r<=R-2, so there is no origin singularity. The function is the gradient of

    U_D(x)=(c A/2)|x|^2+d A int_0^|x| psi(u-R)du.

Its Hessian is

    Dg_D=c A I+d A[psi'(r-R)nn*+(psi(r-R)/r)P].                    (3)

It is positive semidefinite. On the support of the radial perturbation r>=R-2>=2, so its radial eigenvalue is at most (c+d)A=A/2 and its tangential eigenvalue at most (c+d/2)A=3A/8. Consequently g_D is smooth, g_D(0)=0, and 0<=Dg_D<=A I globally.

## 2. True history and finite cheap clocks

Use a positive clock rule (v_j,tau_j), all tau_j in (0,1), with

    sum_j v_j=1, sum_j v_j tau_j=1/2,
    ||Q-R1||_(L2(gamma)->L2(gamma))<=delta, delta<=A^2,
    (Q F)(x)=sum_j v_j P_(tau_j) F(x).

The admitted Hermite-multiplier quadratures have these properties, with their stated polylogarithmic node counts. In particular

    |sum_j v_j tau_j^2-1/3|<=delta.                                (4)

No quadrature of a true Markov path is executed. Define beta_Q=sum_j v_j sqrt(1-tau_j^2). Positivity and sqrt(1-u)>=1-u give

    1>=beta_Q>=1-sum_j v_j tau_j^2>=2/3-delta>=7/12,
    beta_Q^2-1/4>=13/144.                                         (5)

The last inequalities hold since delta<=A^2<=1/16<1/12. Thus a fixed covariance gap is already forced by the second-moment accuracy; no continuum limit for the cheap clocks is needed.

At the same fixed endpoint x, the two ORIGINAL-g sources have exact decompositions

    I_a=c A[x/2+sigma_a G_a]+d A K_a, a in {H,Q},
    sigma_H=1/2, sigma_Q=beta_Q,                                  (6)

where G_a is standard Gaussian independent of x and |K_a|<=1 pathwise. For H,

    K_H=int_0^infinity e^(-t) f_D(X_t^x)dt,

and G_H is the normalized linear OU integral of its genuine Brownian history. For Q,

    K_Q=sum_j v_j f_D(tau_j x+sqrt(1-tau_j^2)H), G_Q=H.

G_a and K_a are generally correlated; no decoupling between them is assumed. The identity sigma_H^2=1/4 is the exact conditional covariance of int e^(-t)X_t^x dt. Also

    E K_H=(R1 f_D)(x), E K_Q=(Q f_D)(x),
    ||E K_H-E K_Q||2<=delta ||f_D||2<=delta.                      (7)

This is a mean-operator error on a bounded original VALUE field, not a strong path approximation. The linear clock means agree exactly. Therefore the total mean-source mismatch is at most d A delta in L2.

## 3. Uniform radial estimates

Let X~gamma_D, S_D=|X|-R, and y0=(1-cA/2)X. For any sigma in [0,1], take G standard independent of X and write

    y_sigma=y0-c A sigma G.

All constants below depend only on the fixed c,d,psi. For every fixed p, ||S_D||p is bounded uniformly in D. The following estimates hold uniformly in sigma:

    ||Df_D(y_sigma)-Df_D(y0)||_(Lp;op)<=C_p A,                   (8)
    E_G f_D(y_sigma)
       =f_D(y0)+(c^2 sigma^2/2)A psi'(S_D-c/2)n(X)
                                      +O_(L2(X))(A^2).         (9)

Here and below the arbitrary definition of n at zero is irrelevant.

For completeness, these are radial small-increment estimates, not a small Euclidean displacement assumption. On |X|>=R/2 and |w|<=|y0|/2, with w=-c A sigma G, rho=|y0| and n0=y0/rho,

    |y0+w|-rho=n0.w+(|w|^2-(n0.w)^2)/(2rho)
                                      +O(|w|^3/rho^2),
    n(y0+w)=n0+P0 w/rho+O(|w|^2/rho^2).

The Gaussian moments obey ||n0.w||p=O(A), ||w||p=O(1), and rho is comparable to R. Thus the radial increment and direction change have Lp size O(A), even though |w| itself has order-one energy. The complements of these events have exponentially small Gaussian probability and are harmless under the bounded derivative and polynomial moment majorants.

Formula (3) without c A I yields (8), since psi' is Lipschitz and the direction changes by O(A). Taking expectations in the radius expansion gives

    E_G[|y0+w|-rho]=c^2 A^2 sigma^2(D-1)/(2rho)+O_(L2(X))(A^2),
    E_G[(|y0+w|-rho)^2]=O_(L2(X))(A^2).

The expected first direction increment is zero, its remainder is O(A^2), and the radius/direction cross term is O(A^2). Taylor expansion of the scalar psi therefore gives (9), using

    rho-R=S_D-c/2-(cA/2)S_D,
    A^2(D-1)/rho=A+O_(L2(X))(A^2).

All uses of Taylor in this counterexample are justified by its FIXED smooth bump. No higher-derivative assumption is imposed on the general class being refuted.

## 4. Expansion of the actual force means

The true terminal argument is y_(sigma_a)-d A K_a. Since |K_a|<=1 and the bilinear second derivative of f_D is globally bounded by an absolute constant,

    f_D(y_(sigma_a)-d A K_a)
      =f_D(y_(sigma_a))-d A Df_D(y_(sigma_a))K_a+O(A^2).

After expectation and multiplication by the outer d A, replacing Df_D(y_(sigma_a)) by Df_D(y0) costs O_L2(A^3), by (8). This remains valid with the genuine correlation between G_a and K_a. The resulting mean term is -d^2 A^2 Df_D(y0)E K_a.

The linear part of the outer g has mean difference

    -c A(E I_H-E I_Q),

whose L2 norm is at most c d A^2 delta by (7). The difference between the two displayed K mean terms also costs at most C A^2 delta. Equations (7)-(9) therefore imply, directly for the true history versus FINITE Q,

    j_H-j_Q
      =(d c^2/2)A^2(1/4-beta_Q^2)
                         psi'(S_D-c/2)n(X)
                              +O_L2(A^3+A^2 delta).             (10)

The error bound is uniform over all the admitted Q satisfying the stated properties. In particular delta<=A^2 makes it O(A^3).

## 5. The actual covariance contraction at x

From (6), conditional on x,

    C_a=c^2 A^2 sigma_a^2 I+R_a,
    ||R_a||HS<=C A^2,                                            (11)

uniformly in x,D,Q. Indeed the mixed coefficient is Cov(G_a,K_a). Conditional Gaussian first-chaos Bessel, applied to each physical component of K_a, gives

    ||Cov(G_a,K_a)||HS^2<=E|K_a-E K_a|^2<=1.

This uses only that G_a is a standard Gaussian vector on the same probability space; its independence from K_a is neither needed nor true. Also ||Cov(K_a)||HS<=tr Cov(K_a)<=1. Expanding the covariance in (6) proves (11).

The exact radial derivative formula, for any matrix B, is

    D2 f_D(x):B
      =psi'' n(n* B n)
        +a_r[n tr(PB)+P(B+B*)n],
    a_r=psi'/r-psi/r^2,                                         (12)

where psi and its derivatives are evaluated at r-R. Consequently

    ||D2 f_D(x)||_(HS->vector)
        <=|psi''|+(sqrt(D)+2)|a_r|<=C,
    ||D2 g_D(x)||_(HS->vector)<=C A.                             (13)

For r<R-2 everything vanishes; on the transition shell r is comparable to R; outside it psi'=psi''=0 and the remaining bound is even smaller. These are global estimates for this smooth family. Thus the covariance remainders in (11), contracted with D2g, cost only O(A^3), not a spurious dimension factor.

Setting B=I in (12) gives

    Delta g_D=d A[psi''+(D-1)(psi'/r-psi/r^2)]n
              =d psi'(S_D)n+O_L2(A).                            (14)

Combining (11)-(14) proves

    (1/2)(C_H-C_Q):D2g_D(x)
      =(d c^2/2)A^2(1/4-beta_Q^2)psi'(S_D)n(X)
                                             +O_L2(A^3).        (15)

This uses the ACTUAL covariances of the same original-g sources. It is not a covariance substituted from an unrelated quadratic model; the exact remainder has been controlled.

## 6. A linear R1 witness proves the separation

Put

    h(s)=psi'(s-c/2)-psi'(s),
    k_Q=(d c^2/2)(1/4-beta_Q^2).

Equations (10) and (15) give the pre-resolvent discrepancy

    k_Q A^2 h(S_D)n(X)+O_L2(A^3+A^2 delta).                       (16)

By (5), k_Q is negative with absolute value bounded below by a positive constant independent of D and Q.

The Gaussian shell central limit theorem gives S_D -> S in distribution, S~N(0,1/2), while |X|/R->1 in L2. Since psi' is a nonzero symmetric strictly unimodal bump, its convolution with a nondegenerate centered Gaussian has a strict maximum at zero. Equivalently, integrate its positive layer-cake intervals and use strict decrease of Gaussian interval mass after a nonzero translation. Hence

    E h(S)=E psi'(S-c/2)-E psi'(S)<0.                            (17)

Take the vector L2 test T_D(X)=A X=X/R; its norm is exactly one. The Gaussian Mehler resolvent R1 is self-adjoint and R1 T_D=T_D/2. Thus

    <R1[h(S_D)n],T_D>
       =(1/2)E[h(S_D)|X|/R] -> (1/2)E h(S),                    (18)

which is nonzero. This avoids any assumption about conditional convergence of the full radial history, or a joint shell/path limit.

R1 is an L2 contraction, so (16)-(18) and delta<=A^2 prove, for all sufficiently large D and every such Q,

    ||R1(j_H-j_Q)-(1/2)R1[(C_H-C_Q):D2g_D]||2
                                      >=c_* A^2.                (19)

Since A^4 sqrt(D)=A^3, (19) refutes (1). The lower-to-claimed ratio is c_*/A and diverges faster than every fixed polynomial in log D and log(1/A). The allowed mean-quadrature floors are smaller and cannot absorb the gap.

## 7. Exact scope: coherent centering survives

This is a GENUINE same-g separator. Both I sources, the outer force, and their covariance coefficients use exactly (2), with the literal Markov or shared-root genealogy. It refutes the unshifted candidate in Section 4 of `P3-MEAN-RAW-SOURCE-AND-SECOND-DECOUPLING-GATE.md` insofar as that candidate is specialized to the single-I comparison, and companion equation (13) with the asserted uniform A^4 remainder. It does not use or revive an ordinary Markov path-grid algorithm.

For this family,

    mu_H(x)=c A x/2+d A E K_H,
    |d A E K_H|<=d A.

Therefore x-mu_H has shell coordinate S_D-c/2+O_L2(A) and the same leading direction. Repeating (14)-(15) at x-mu_H changes psi'(S_D) to psi'(S_D-c/2). Consequently the CENTERED covariance evaluation matches (10) to O_L2(A^3+A^2 delta) for this family. The exact centered/resummed current identities remain valid and are not refuted.

This last observation is not a proof that the centered covariance formula has a uniform fourth-order remainder for every admissible g. Such a theorem still needs the actual one-energy contracted-trace cancellation and its own native original-VALUE consumer. The radial counterexample closes only the unshifted-candidate question.
