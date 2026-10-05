# Independent audit: independent-node source, covariance, and Gaussian cuts

2026-10-05. **PASS for this bounded new-source/covariance/conditioning test.**
The resummed same-carrier m3 current and a finite current-to-law compiler remain
OPEN. This audit neither supplies those objects nor reuses any old common-H
covariance receipt for the new source.

## Frozen input and assumptions

Input: `../INDEPENDENT-NODE-MEAN-SOURCE-COVARIANCE-AND-CUT-TEST.md`.
SHA256: `c463ea3ce5f2d7584dca346c00bbec345592e07e90481f10bd708048de56476c`.
The pin was checked before and after the work; the author input was not edited.

The note is read under the surrounding problem's standing assumptions:
g = grad U, g(0) = 0 (or the identical saved zero anchor),
0 <= Dg <= A I, and the small-A regime 0 < A <= 1/2. The zero and moment
statements require that anchor assumption; the Hessian bound alone would not
remove an arbitrary constant g(0). No third derivative bound is used below.
All clock rules are positive, interior, of mass one and first moment 1/2.
Statements invoking a completed OWN-mean compiler remain subject to its
inherited active-dimension, radius, clock, floor, and caller guards.

## 1. Inner source and changed covariance

At fixed x, write K_j = Dg(tau_j x + c_j H_j). Direct first differentiation gives

    L = D_x I_ind = sum_j v_j tau_j K_j,
    M = D_H I_ind = [v_1 c_1 K_1 | ... | v_N c_N K_N].

Thus 0 <= L <= A I/2 and ||M||op <= A beta_2, with
beta_2^2 = sum_j v_j^2 c_j^2. In fact beta_2^2 <= 3/4: use v_j^2 <= v_j,
then Jensen on sum v_j tau_j^2 >= 1/4. The note's weaker beta_2 <= 1 is safe.
For T = Dg(x-I_ind), the executed terminal firsts are T(I-L) and -TM.
This proves the note's stated J_ind ports without differentiating a Hessian.

Positive-weight Minkowski gives

    ||I_ind||_(Lp|x)
      <= A [|x|/2 + C_p sqrt(D) sum_j v_j c_j].

This is an energy estimate using each node's D-dimensional marginal. It does
not turn the actual ND independent input tape into a D-dimensional tape.
The conditional mean is sum_j v_j P_tau_j g. Independence gives precisely
sum_j v_j^2 Cov(g(tau_j x+c_j H_j)); there are no cross-node terms.
Gaussian Poincare applied to each scalar projection proves
0 <= C_ind <= A^2 beta_2^2 I and its stated HS bound.

For g(x)=Bx, the three coefficients are beta_2^2 B^2, (sum v_j c_j)^2 B^2,
and B^2/4 for independent, common-root, and true continuous Markov histories.
The last factor follows because the double integral of min(t,s)/max(t,s)
over the unit square is 1/2, while conditioning on x subtracts
(integral t dt)^2 = 1/4. Hence the covariance replacement really changes.
A resummed correction would have to use C_true-C_ind and its own live record.
No full-private-gradient property of I_ind is needed or inferred.

The earlier THIRD-order outer-mean decoupling argument applies with the
displayed actual private first and energy: its constants are dimension-free
in the number of Gaussian roots. Its inner conditional mean is unchanged.
The outer source mean is still a one-variable Mehler function, so the genuine
Hermite-multiplier error is charged in exactly the earlier way. This proves
only the inherited third-order comparison, not the missing fourth-order bias.

## 2. Outer source: first, energy, curl, origins, and counts

For outer node i write T_i = Dg(x_i-I_i), K_i = Dg(x_i), and use L_i,M_i above.
The correction E=F_ind-B has exact first blocks

    D_G E = sum_i w_i d_i [(T_i-K_i)-T_i L_i],
    D_H E = -sum_i w_i T_i M_i,
    D_Z E = sum_i w_i r_i [(T_i-K_i)-T_i L_i].

T_i-K_i is symmetric, with operator norm at most A. For
beta_out = sum_i w_i d_i, these identities yield the explicit safe bounds

    first(E) <= A(beta_out+A),
    ||D_Z E||op <= A/2+A^2/4,
    ||Curl(P_G*E)||op <= A^2(beta_out+beta_2).

Here Curl means J-J*. For the first bound, the remaining G and H blocks have
combined norm at most A^2 sqrt(beta_out^2/4+beta_2^2) <= A^2. For the curl
bound, the G-G skew comes only from T_i L_i and is bounded by A^2 beta_out;
the full G-H/H-G off-diagonal block is bounded by A^2 beta_2. This explicitly
checks the full (N+1)D lifted matrix, rather than only its G-G corner.
The baseline is a genuine G-gradient. Consequently F_ind also has O(A)
private first and captured-Z first. These are actual first bounds, not an
invalid differentiation of a small energy estimate.

Pointwise |g(x-I)-g(x)| <= A|I| and positive-weight Minkowski give
||E||_(Lp|Z) <= C_p A^2(|Z|+sqrt(D)). Shared H_j across outer nodes causes
no problem for this estimate or the one-variable outer mean identity.
The entire shared bank must still be regenerated for another complete
source occurrence or changed source arguments.

The literal counts and dimensions are correct:

    raw F: N_out(N+1) original VALUES, private dimension (N+1)D;
    raw B: N_out original VALUES, private dimension D;
    raw E: N_out(N+2) original VALUES, private dimension (N+1)D.

These are safe full occurrence counts, not a count after averaging inner
roots or suppressing replay. The bill Q_captured+N_B Q_B+N_E Q_E must retain
all actual consumer roots and all numerical/replay work. Increasing N can
change active-dimension guards even when its growth is only polylogarithmic.

At Z=G=H=0 all values are literally zero under the saved g(0) anchor. At fixed
nonzero Z, the G=H=0 programs are executed origins, not artificial zeros.
One also has |E(Z,0,0)| <= A^2|Z|/4. An inner origin with x_i depending on G
is private; it cannot be cached as a caller-only constant. The original g
VALUES are the only leaves. Requested first/adjoint sweeps may use HVPs at
the recorded leaves, through the complete graph; no HVP is a producer or
is itself differentiated.

If each inner VALUE has an absolute error floor eps_in and the terminal has
eps_term, then a raw terminal has floor at most A eps_in+eps_term. The outer
positive sum does not amplify it. E additionally pays its baseline floor.
These absolute floors do not improve with cancellation or realized energy.
All response, inverse-heat, and covariance-encoding multipliers would need
their own precision charge. The source note correctly leaves those charges
as obligations rather than claiming them paid.

## 3. Gaussian cut and full-record conditioning

At fixed Z=z, all original G_0,H_j are independent standards and
x=rz+qG_0, Y_j=tau_j x+c_jH_j. The direction coefficients
(1/q, -1_{j in S} tau_j/c_j) give

    D_A x = I,
    D_A Y_j = 0 for j in S, and tau_j I otherwise.

Their Gaussian score S_S has variance
sigma_S^2 I = [q^-2+sum_(j in S) tau_j^2/c_j^2] I and zero covariance with
all marked Y_j. Joint Gaussianity proves independence of that whole marked
vector, hence of its measurable g values and Jacobians. The unmarked
derivative sum v_j tau_j Dg(Y_j) is real and has norm at most A/2. Therefore
the changed terminal force and all unmarked descendants must be replayed.

For all sites, an explicit posterior calculation gives

    mu = sigma_all^-2 [q^-2 rz + sum_j tau_j Y_j/c_j^2],
    x-mu = sigma_all^-2 S_all,
    Var(x | Z,Y_all) = sigma_all^-2 I.

Thus exactly D linear Gaussian innovation dimensions remain if x is unread.
The all-site direction preserves I_ind. It still changes the terminal input
x-I_ind, so it does not preserve the terminal force record. Its precision
contains no v_j. For a node tau=1-epsilon, its contribution is asymptotic to
1/(2 epsilon). Polylogarithmic node count alone cannot pay this score norm
or establish a public-log current/response estimate. This is a cost warning,
not an impossibility theorem for a future weighted cancellation construction.

If x and Y_all are captured, H_j=(Y_j-tau_j x)/c_j and G_0=(x-rz)/q are
determined. If only Y_all and the terminal VALUE t=g(x-I_ind) are captured,
a strictly monotone g is injective on its image and

    x = I_ind(Y_all) + g^{-1}(t).

So that full record also determines the original roots. If the inner record
contains g(Y_j) instead of Y_j, injectivity first recovers each Y_j. An
invertible positive linear g already gives this obstruction within the
admitted class. Noninjective examples cannot justify a uniform fresh-Gaussian
claim. These statements concern the same original marked roots; they do not
rule out introducing a separate, explicitly paid conditional construction.

Before the marked descendants are evaluated, conditioning on x and the
unmarked sites leaves the marked H_j independent standard Gaussians. This
legitimate earlier cut is different from conditioning on the completed
force/terminal record. A native response would have to propagate that earlier
cut through all its descendants and cannot append integrated roots later as
observers. The note's stopping boundary is therefore correct.

## 4. Independent finite checks and verdict

`check_independent_node_audit.py` was written independently of the author's
checker. `independent_checks.json` records **3,174 passing assertions** over
72 multivariate Gaussian configurations, 72 smooth nonlinear-gradient
configurations, and four endpoint-precision stress values. These cover exact
covariance, Schur complement and rank identities; marked/unmarked motion;
actual private dimensions; noncommuting Hessian first/curl formulas and a
directional finite-difference check; full VALUE counts, energies, and origins.
Finite diagnostics support, rather than replace, the proof above.

The author's checker was copied into this audit directory and rerun without
editing the original checker or its JSON. The copy reproduced **960 passing
assertions on 96 Gaussian configurations**. Its output is
`independent_node_cut_checks.json`; the copy is `author_check_reproduced.py`.
Those 960 assertions are not counted as new independent assertions.

**Final decision:** admit the bounded independent-node mean-source,
changed-covariance, and cut/conditioning test. Do not promote it to an m3
current construction, a current-to-law compiler, or a fourth-order endpoint.
No correction to the frozen author pin is required for this scoped verdict.
