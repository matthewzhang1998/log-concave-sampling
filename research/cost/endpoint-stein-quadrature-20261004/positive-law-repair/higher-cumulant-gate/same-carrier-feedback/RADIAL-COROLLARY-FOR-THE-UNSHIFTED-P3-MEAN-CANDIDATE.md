# Radial corollary: the unshifted P3 mean candidate also loses the coherent shift

2026-10-05. Corollary candidate for independent review. It concerns only equation (9) of `P3-MEAN-RAW-SOURCE-AND-SECOND-DECOUPLING-GATE.md`; the literal F3_Q and a correctly resummed correction remain viable targets.

## Statement

Use exactly the same smooth radial gradient family, constants, bump and finite clock guards as the independently PASS source

    RADIAL-SHELL-REFUTES-UNSHIFTED-SINGLE-HISTORY-COVARIANCE.md
    SHA256 2807f7307646cc16c79049319b574259874b09e731a13f686696bdd60a05b035.

Thus A=D^(−1/2), c=d=1/4, g_D(x)=A[c x+d psi(|x|−sqrt(D))n_x], and all occurrences use that same original g_D.

Let F2_H(x) be the TRUE second-substitution force conditioned on X1=x. Let F2_Q(x;H,J) be the actual finite two-level common-root force used INSIDE the three-root raw source:

    y_j=tau_j x+c_j H,
    L_j=sum_k u_k g_D(sigma_k y_j+e_k J),
    F2_Q=sum_j v_j g_D(y_j−L_j).

Require the original inner and middle rules to have their actual admitted Hermite multiplier tolerances at most A2 and exact first moments 1/2. Define

    psi2_H(x)=E g_D(x−F2_H), psi2_Q(x)=E g_D(x−F2_Q),
    C_H(x)=Cov(F1_H|x), C_Q(x)=Cov(F1_Q|x),
    F1_Q=sum_j v_j g_D(tau_j x+c_j H).

Then the candidate

    m3=R1 psi2_H
        ?=R1 psi2_Q+(1/2)R1[(C_H−C_Q):D2g_D(x)]
                              +O(Lambda A4 sqrt(D))            (1)

has discrepancy at least c_* A2 for all sufficiently large D. The proposed allowance is Lambda A3. No assertion is made about a centered/resummed terminal current.

## 1. The second substitution still has the same leading Gaussian field

Write g_D=c A id+d A f_D, with |f_D|<=1. Positivity of all clock weights gives the exact decompositions

    F2_a=c A[x/2+sigma_a G_a]+A K_a, a=H,Q,
    sigma_H=1/2, sigma_Q=beta_Q=sum_j v_j c_j.          (2)

For the true history, with I_t the first force displacement,

    K_H=−c integral_0^infinity e^(−t) I_t dt
                      +d integral_0^infinity e^(−t) f_D(X_t−I_t)dt.

The finite K_Q is the identical algebraic decomposition of the displayed finite VALUE graph. The first term in (2) uses only the original linear Gaussian sites before the inner force shift. It is standard Gaussian conditional on x, and may be correlated with K_a.

For every fixed p,

    ||K_a||_(Lp(all original Gaussian roots, x~gamma_D))<=C_p,  (3)

uniformly in D and in the positive rules. Indeed ||I_t||p<=C_p A sqrt(D)=C_p, the corresponding finite averages have the same Minkowski bound, and |f_D|<=1. The actual conditional profile is at most C_p(1+A|x|); no pointwise boundedness of K_H is asserted.

The independently audited old nested-mean theorem gives

    ||E F2_H−E F2_Q||_(L2(gamma_D))<=Lambda A3 sqrt(D)
                                                     =Lambda A2.

The conditional means of the leading Gaussian terms in (2) agree exactly, hence

    ||E K_H−E K_Q||2<=Lambda A.                       (4)

This uses the old source's actual mean guarantee, not a joint coupling or a false equality of the two histories.

## 2. Apply the uniform radial expansion with Lp rather than bounded K

The canonical radial proof gives, for y0=(1−cA/2)x and y_sigma=y0−cA sigma G,

    E f_D(y_sigma)=f_D(y0)
       +(c2 sigma2/2)A psi'(S_D−c/2)n_x+O_L2(A2),
    ||Df_D(y_sigma)−Df_D(y0)||_(L4;op)<=C A,
    ||D2 f_D||op<=C.                                (5)

The terminal inner argument is y_sigma−A K_a. By Taylor and (3),

    E f_D(y_sigma−A K_a)
       =E f_D(y_sigma)−A Df_D(y0)E K_a+O_L2(A2).     (6)

To justify the correlation term explicitly, use Holder:

    ||[Df_D(y_sigma)−Df_D(y0)]K_a||2
      <=||Df_D(y_sigma)−Df_D(y0)||4 ||K_a||4<=C A.

The second-order Taylor remainder costs A2||K_a||4². No independence of G_a and K_a is needed.

Now (4) makes the difference of the linear K terms in (6) O_L2(Lambda A2). The linear c A part of the OUTER g likewise costs c A2||E K_H−E K_Q||2=O(Lambda A3). Therefore

    psi2_H−psi2_Q
       =k_D A2 psi'(S_D−c/2)n_x+O_L2(Lambda A3),
    k_D=(d c2/2)(1/4−beta_Q2).                       (7)

The covariance current in (1) is the SAME first-displacement covariance current as in the canonical proof; it gives

    (1/2)(C_H−C_Q):D2g_D(x)
       =k_D A2 psi'(S_D)n_x+O_L2(A3).               (8)

## 3. The discrepancy survives the outer mean and its finite quadrature

Subtract (8) from (7). The canonical Ax/R1 witness applies without alteration:

    ||R1[psi2_H−psi2_Q−(C_H−C_Q):D2g_D/2]||2
                        >=c_* A2−C Lambda A3.       (9)

Here beta_Q>=2/3−delta gives a fixed nonzero k_D, and A times any fixed public-log polynomial tends to zero along D=A^(−2). Thus (9) is at least a fixed multiple of A2 for all sufficiently large D.

If the raw F3_Q program also uses a finite top rule with actual multiplier error delta_out<=A3, then

    ||(Q_out−R1)psi2_Q||2
       <=delta_out ||psi2_Q||2
       <=C delta_out A sqrt(D)<=C A3.               (10)

This finite outer floor is smaller than the separator. The raw source's guarded OWN-mean compiler also has only the stated Lambda A4 sqrt(D)=Lambda A3 error. Even an ideal exact implementation of the UNshifted counterterm in (1) would therefore target the wrong mean at order A2.

## Scope

The proof retains every nonlinear ancestor inside K_a and uses its actual positive-weight moment bounds and the independently admitted mean comparison. It does not require a small coherent displacement or a strong shared-root law coupling.

This corollary rejects the particular unshifted P3 repair in (1). It does not refute the F3_Q raw source itself, the exact centered/resummed current, an original-VALUE shifted-current constructor, or a full order-four sampler obtained by a different valid method. The correct m3 remains the same-Gaussian-carrier conditional mean of g(x−F2_H).
