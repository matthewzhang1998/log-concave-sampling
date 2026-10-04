# Positive square-clock covariance service from original gradient VALUES

2026-10-04. New bounded constructive gate. The operator identities, quadrature, original-VALUE source, and final action-to-reserve specialization have an independent audit against the already admitted finite VALUE action. This is a cubic covariance service, not yet an all-order endpoint sampler.

## Result first

The continuous stationary-OU first displacement `H=integral_0^1 g(X_r)dr`, conditional on `X_1=Z`, has the exact covariance identity

    Cov(H|Z)=2 integral_0^1 q P_q[(Dv)^2](Z)dq,
    v=integral_0^1 P_r g dr.                            (S)

Its finite positive square-clock approximation is

    B_Q=sum_i w_i r_i P_(r_i) Dg,
    C_Q(Z)=2 sum_j v_j q_j P_(q_j)[B_Q^2](Z).            (Q)

For positive rules of mass one and exact first moment one half, and uniform Hermite moment error at most delta,

    0<=C_Q(Z)<=A^2 I/4,
    ||C_Q-Cov(H|Z)||_(L2(Z);HS)
        <=(3/2) delta A^2 sqrt(D).                     (E)

Thus covariance quadrature alone has arbitrary fixed heat order with polylogarithmically many original-VALUE source sites and no dimension-dependent tolerance. It avoids the artificial self-diagonal of a finite weighted H product.

Neither Dg nor B_Q is a producer oracle. The literal square-action source is

    f_Q,X(u)=sum_i (w_i r_i/c_i)
                 [g(r_i X+c_i u)-g(r_i X)],
    c_i=sqrt(1-r_i^2).                                  (V)

It is a D-dimensional genuine gradient with actual radius A/2, uniform one-energy A sqrt(D)/2, zero at u=0, and bounded captured-X first. Its analytical mean Jacobian is exactly B_Q(X). The admitted zero-clock VALUE square action consumes (V), producing the target B_Q(X)^2 without an HVP-valued leaf. Averaging complete actions over positive q clocks and an owned Gaussian X=qZ+sqrt(1-q^2)G, followed by a single positive Gaussian reserve, gives covariance `v0 I+C_Q(Z)` at error `Lambda A^3 sqrt(D)` plus the separate clock/filter/numerical floors.

## 1. Exact square identity

Let `F=Dg`, only for analysis. The ordered conditional covariance is

    C(Z)=2 integral_0^1 r dr integral_0^1 d tau
      Sym{P_r[g(P_tau g)^T]-(P_rg)(P_(r tau)g)^T}(Z).

The OU covariance identity yields

    P_r[g(P_tau g)^T]-(P_rg)(P_(r tau)g)^T
      =2 tau integral_r^1 s P_(r/s)
          [(P_s F)(P_(s tau)F)^T] ds.

Set r=sq and exchange the integrals to obtain

    C=4 integral_0^1 q P_q[
           integral_0^1 s^3 ds integral_0^1 tau d tau
                Sym((P_sF)(P_(s tau)F)^T)] dq.

On the other hand `B=Dv=integral_0^1 s P_s F ds` is symmetric, and splitting its square into two ordered triangles gives

    B^2=2 integral_0^1 s^3 ds integral_0^1 tau d tau
                 Sym((P_sF)(P_(s tau)F)^T).

This proves (S). All covariance centering remains coupled through these identities. For g(x)=Bx, the answer is B^2/4. An independent polynomial check g(x)=x^2 gives `C(Z)=2(1+Z^2)/9` from both sides.

Because `0<=F<=A I`,

    0<=B,B_Q<=A I/2,
    ||B||_(L2;HS), ||B_Q||_(L2;HS)<=A sqrt(D)/2.         (2)

The square in (S) and (Q) is pointwise positive semidefinite; no signed tensor integral has to be treated as a covariance oracle.

## 2. Finite positive quadrature and its exact error

Use the previously admitted positive dyadic Gauss rule in distance `1-r`, with K dyadic panels, m Gauss nodes on each, and one terminal midpoint. It has

    sum w_i=1, sum w_i r_i=1/2,
    sup_(ell>=0)|sum w_i r_i^ell-1/(ell+1)|<=delta,
    N=O(log^2(1/delta)).                                (3)

Use an independent known rule `(v_j,q_j)` satisfying the same identities; it may be the identical finite version. Shift the chaos degree by one to get

    ||B_Q-B||_(L2;HS)<=delta ||Dg||_(L2;HS)
                         <=delta A sqrt(D).            (4)

The noncommuting matrix identity

    B_Q^2-B^2=B_Q(B_Q-B)+(B_Q-B)B

and (2) give

    ||B_Q^2-B^2||_(L2;HS)<=delta A^2 sqrt(D).            (5)

The outer positive operator `2 sum_j v_j q_j P_(q_j)` has L2 norm at most one, because its total weight is one. Its difference from `2 integral qP_q` has L2 operator norm at most 2delta. Therefore

    ||C_Q-C||<=||B_Q^2-B^2||+2delta||B^2||
             <=(3/2)delta A^2 sqrt(D).                  (6)

This proof uses a single Hilbert-Schmidt energy, not an L4 product of two dimension-sized vector fields. It is valid for the original C2 Hessian class. Dg appears only inside this analytical certificate. The actual g-only program is below.

Both clocks preserve the exact linear/quadratic calibration: when g(x)=Bx, `B_Q=B/2`, hence `C_Q=B^2/4` exactly. Finite coefficient encoding errors receive their numerical budget rather than being called exact floating-point identities.

## 3. Literal genuine-gradient VALUE source

At any captured X, the source (V) uses only g at `r_i X+c_i u` and at the executable same-X anchor `r_i X`. Its analytical potential is

    sum_i (w_i r_i/c_i^2)
       [U(r_i X+c_i u)-c_i g(r_i X) dot u],

but no potential value is queried. Its actual derivative is

    D_u f_Q,X(u)=sum_i w_i r_i Dg(r_i X+c_i u).

Consequently

    0<=D_u f_Q,X<=A I/2,
    |f_Q,X(u)|<=A|u|/2,
    ||f_Q,X||_Lp<=C_p A sqrt(D)/2,
    f_Q,X(0)=0.                                        (7)

These are actual bounds uniform in X. The anchored caller derivative is

    D_X f_Q,X=sum_i (w_i r_i^2/c_i)
                 [Dg(r_i X+c_i u)-Dg(r_i X)],

so

    ||D_X f_Q,X||<=A L_Q,
    L_Q=sum_i w_i r_i^2/c_i<=C.                         (8)

For the last inequality, on a dyadic panel `1-r in [a,2a]`, use `r^2/c<=1/sqrt(a)` and panel weight a. Summing sqrt(a) is a geometric series; the final midpoint contributes at most `sqrt(2h0)`. Thus L_Q is bounded uniformly in K and m (e.g. a conservative constant 4 suffices for the stated rule). There is no endpoint inverse-heat factor in the aggregated source's actual caller first.

Likewise `sum w_i r_i/c_i<=C`, so raw numerical VALUE errors propagate through a bounded absolute coefficient sum. Individual endpoint coefficients and finite square-action response widths remain recorded when absolute leaf tolerances are chosen.

The analytical mean Jacobian is

    E_u D_u f_Q,X(u)=sum_i w_i r_i P_(r_i)Dg(X)=B_Q(X).  (9)

A raw occurrence has N shifted original VALUES and up to N captured same-X anchor VALUES. Identical anchors are shared only under the identical source/caller/version key. Under a changed X, all anchors are new original sites. At zero u, reuse the recorded anchor, so exact zero is literal even for a frozen numerical primitive.

## 4. Aggregate a finite family of actual square actions

Use the admitted LOW30 `C_sq(f;p)` action with fixed `s0 in (0,1/4)` and padding `0<mu<=1`. For a genuine full-gradient source its analytical mean matrix is

    (s0 E Df)^2.

Its finite VALUE program consists of the recorded first-coefficient filter and nested signed VALUE differences. It is pointwise odd in p. Its calibration, energy and complete first are

    mean error <=Lambda ell e mu+e_num,
    action energy <=Lambda ell e+e_num,
    private/p first <=Lambda ell^2/sqrt(mu),
    captured parameter first <=Lambda ell L_X/sqrt(mu).
                                                               (10)

For (V), `ell<=A/2`, `e<=C A sqrt(D)`, and `L_X<=C A`. Enforce the imported small-radius and finite-filter guards at these actual values.

For every finite q_j draw an owned independent coarse root G_j and put

    X_j=q_j Z+sqrt(1-q_j^2)G_j,
    a_j=2v_j q_j, sum a_j=1.

Run independent COMPLETE action-private banks for the sources `f_Q,X_j`, but give all actions the SAME incoming standard p. Execute

    D_Q(p)=sum_j a_j C_sq(f_Q,X_j;p)/s0^2.              (11)

This is an ordinary vector-valued VALUE program. It is not a random matrix multiplication instruction. The common p is intentional: its conditional mean is one matrix times p. Each entire private bank, including its owned G_j, is integrated before the final comparison. Endpoint Z remains captured throughout.

Conditioned on Z, (9)-(11) give

    ||E_private D_Q(p)-C_Q(Z)p||_(L2_p)
        <=Lambda A^2 mu sqrt(D)+e_num,
    ||D_Q||_(L2_(p,private))<=Lambda A^2 sqrt(D)+e_num,
    Lip_(p,all private) D_Q<=Lambda A^2/sqrt(mu),
    ||D_Z D_Q||<=Lambda A^2/sqrt(mu).                   (12)

The first bound follows after averaging the owned X_j; the energy uses the positive sum of actual source energies; and the first bounds include the path `G_j ->X_j->f_Q,X_j` using (8) and the last line of (10). Weighted summation has total mass one, so there is no hidden number-of-clocks or ambient-dimension power. The actual full tape and complete source occurrence count are still retained. The aggregated action is used only in this terminal buffered reserve; it is not asserted to be a genuine-gradient source for later re-entry.

In particular the averaging in (11) is not a strong estimate of B_Q, its square, or the conditional covariance. The program is used only through its completed positive law with the buffer below.

## 5. One positive reserve, exact centering, and the error bill

Fix an available variance share `v0>0` bounded away from zero, choose fixed positive eta,zeta with `eta^2+zeta^2=v0`, and draw an independent standard buffer z. Execute

    R=eta p+D_Q(p)/(2eta)+zeta z.                       (13)

Every sample in (13) is from a positive Gaussian-root pushforward. No signed probability or pointwise square root of a sampled covariance is used. Since (11) is pointwise odd in p, R has exactly zero mean conditional on Z (up to separately assigned numerical symmetry errors). Its source-zero carrier is the actual `eta p+zeta z`.

The admitted private-action buffer estimate from the existing curvature reserve applies with the entire action-private tape, including all owned G_j. It is the Gaussian-Lipschitz one-energy estimate: conditional on p,Z, interpolation from the private action to its conditional mean costs its actual private first times its conditional energy, against the independent z. Square-integrating the conditional energy and using (12) gives

    W2(R, eta p+E_private D_Q(p)/(2eta)+zeta z)
       <=Lambda A^4 mu^(-1/2) sqrt(D)+floors.           (14)

Mean calibration then costs `Lambda A^2 mu sqrt(D)`. The exact linear reference is

    (eta I+C_Q(Z)/(2eta))p+zeta z.

Its covariance is

    v0 I+C_Q(Z)+C_Q(Z)^2/(4eta^2).

The final positive quadratic term is not erased: the fixed-gap Gaussian-root inequality, (2) and (Q), price it by `C A^4 sqrt(D)`. Thus, uniformly in the retained Z at the service level,

    W2(Law(R|Z),N(0,v0 I+C_Q(Z)))
       <=Lambda [A^2 mu+A^4(1+mu^-1/2)]sqrt(D)
                    +e_filter+e_num.                  (15)

Take `mu=A` and enforce the imported actual small-A/public-log guards. This gives `Lambda A^3 sqrt(D)` plus floors. Using (6) and the same fixed-gap Gaussian-root inequality restores the continuous covariance in L2 of the standard endpoint:

    ||W2(Law(R|Z),N(0,v0 I+Cov(H|Z)))||_(L2(Z))
       <=Lambda A^3 sqrt(D)+C delta A^2 sqrt(D)
                    +e_filter+e_num.                  (16)

The service error before this clock restoration is uniform in Z; the quadrature restoration is explicitly an L2(standard Gaussian endpoint) bound. It is not silently assigned to arbitrary deterministic endpoints.

The residual/private/captured-Z first of (13) is at most `Lambda A^2/sqrt(mu)`, in addition to its fixed known Gaussian carrier row. With mu=A this is `Lambda A^(3/2)`, a conservative acceptable absolute first. Original physical caller and finite-mode profiles must be restored by the actual source chain rule; no derivative of (16) is asserted.

## 6. Full finite count, quadratic ledger, and limits

Let N_r,N_q be the inner/outer positive clock counts and N_sq the ACTUAL complete raw-source occurrence count of the finite square action, including its first-coefficient filter and all response replays. A safe original g-VALUE bound is

    2 N_r N_q N_sq + Q_captured + Q_numerical.          (17)

For the printed LOW30 zero-clock action itself, a K_f-node first filter uses 2K_f raw f calls; the outer two signed responses use four more. Thus N_sq=2K_f+4 before any separately required restoration/replay. If each bank captures its N_r identical same-X anchors once, the explicit basic count is

    N_r N_q (2K_f+5)

original g-VALUES, plus external mode/anchor/numerical restoration. The basic Gaussian tape is `D(4N_q+2)`: each bank owns its G_j,z1,z0,z2; all banks share p, and the terminal buffer z is independent. The first-filter nodes share the recorded z1 as prescribed by the action. With `delta=A^J` for fixed J and K_f=O(log(1/A)) at the declared padding/filter floor, the displayed basic count is a public-log polynomial (at most the indicated log-fifth pattern), not an inverse-heat count.

Identical same-X anchors may reduce the count only when genuinely captured within one complete action bank. Discarded primal records incur their full replay. Each requested first/adjoint is an original HVP at a recorded VALUE site; no HVP is a producer leaf. The private tape contains every q-owned root, action root/filter root, p, and reserve z. All versions, padding, counts and numerical tolerances are frozen before caller differentiation.

For fixed target clock order, N_r,N_q are public-log counts. The square action has the already admitted fixed-depth polylog count. There is no inverse-heat replication bank and no hidden algebraic dimension count beyond the original D-dimensional vector operations. Precision budgets retain all recorded response widths and coefficient encodings.

For g(x)=Bx, (V) is exactly `Bu/2` for every X and finite clock rule. Every literal action response is linear and all private source noise cancels from its square coefficient. Equation (11) has exact target `B^2p/4`. The final reserve covariance is `v0 I+B^2/4+B^4/(64eta^2)`; its B-squared coefficient is exactly correct and the displayed degree-four term is paid. This is full matrix calibration, not scalar matching.

The purpose is to replace the bad finite product-clock covariance by the correct continuous covariance at cubic law grade, with an arbitrary-order finite clock certificate. It does not by itself approximate the nonlinear rank-one feedback, prove a final endpoint/full-law join without its required reserve, or create an arbitrary-order covariance action from the bounded action (10). Those are separate gates. In particular the earlier finite-path diagonal separator remains valid for direct nested positive-clock predictors; (13) escapes its scope by using a completed VALUE square action and independent Gaussian reserve.

## Source pins

- Exact OU two-clock/current reference: `../EXACT-TWO-CLOCK-CURRENT-AND-UNIFORM-C2-REFERENCE.md`.
- Existing finite one-clock proof: `../../../ENDPOINT-STEIN-PACKET-AND-POLYLOG-QUADRATURE.md` relative to this folder's ancestry (the named file in endpoint-stein-quadrature-20261004 is authoritative).
- LOW30 original `t30:lem:actions`, especially the zero-clock square variant, first-coefficient filter, and actual captured-parameter derivative bound.
- Admitted action-to-positive-reserve proof: `../../reverse-ou-curvature/THIRD-ORDER-VALUE-CURVATURE-RESERVE.md`, Section 5, with its independent audit.
- Independent checks and review: `independent-audit/` in this clock-compression folder.
