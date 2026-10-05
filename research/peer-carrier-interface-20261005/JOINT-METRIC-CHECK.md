# Why the retained carrier metric matters

5 October 2026. Independent mathematical check of the carrier interface. This note establishes a distinction in metrics; it does not assert that the peer intended the weaker metric.

## Sufficient composition

Let E be the unchanged external record, P a fresh standard Gaussian independent of E, and p_i=A_i^T P fixed known projections. Suppose

    Law((Y_i^a)_i | E,P) = product_i K_i^a(E,p_i),   a in {actual,ideal}.

Define epsilon_i^2=E_(E,P) W2^2(K_i^actual(E,p_i),K_i^ideal(E,p_i)). Couple E and P identically and take the product of measurable conditional optimal couplings (or arbitrarily close measurable couplings). Each marginal has its required conditional product law. Therefore

    ||sum_i a_i (Y_i^actual-Y_i^ideal)||_L2
       <= sum_i |a_i| epsilon_i.

The complete output tuple instead has the Euclidean-product bound sqrt(sum_i epsilon_i^2). Any unchanged C(E,P) cancels exactly. A nonlinear caller requires its actual coordinate sensitivity bounds. Projection overlap does not change the argument.

This is a constrained joint transport distance: the coupling must retain E and P identically almost surely. It is not ordinary joint Wasserstein distance, which permits movement of those coordinates.

## Counterexample for ordinary joint Wasserstein

Let P=(P_1,P_2) be standard two-dimensional Gaussian, n a positive integer, theta_n=n^(-1/2), and

    f_n(P)=sin(n P_1),
    g_n(P)=sin(n (cos(theta_n) P_1 + sin(theta_n) P_2)).

The descriptors are fixed before P is drawn. There are no private tapes, so conditional packet independence holds trivially. The same identity projection P_i=P is legal for every packet.

A Gaussian rotation Q can make g_n(Q)=f_n(P) identically while

    E|P-Q|^2 = 4(1-cos(theta_n)).

Consequently ordinary joint transport obeys

    W2(Law(P,f_n(P)),Law(P,g_n(P)))
       <= 2 sqrt(1-cos(theta_n)) -> 0.

Now leave a second packet f_n(P) unchanged and replace the first by g_n(P). The sum caller is just addition and has unit sensitivity to each packet. Its actual and replaced outputs are 2f_n(P) and f_n(P)+g_n(P), respectively. Both have mean zero. Gaussian characteristic functions give

    Var(f_n) = (1-exp(-2 n^2))/2,
    Cov(f_n,g_n) = [exp(-n^2(1-cos(theta_n)))
                    - exp(-n^2(1+cos(theta_n)))]/2.

Thus Var(2f_n)->2 and Var(f_n+g_n)->1. The L2 triangle inequality in any coupling implies

    liminf W2(Law(2f_n),Law(f_n+g_n)) >= sqrt(2)-1.

Hence small ordinary joint packet error does not compose even for a sum caller, with deterministic descriptors, no hidden private dependencies, and a completely erased carrier in the final output. An environment's carrier sensitivity or a carrier-fixed transport condition is genuinely necessary. The example is not a no-go theorem under additional root-smallness or global-Lipschitz bounds; such bounds must be stated and charged.

## Application to the sealed geometry

Take peer P=H/sqrt(v), and retain the identity projection for each packet. Our original source-zero bank is

    P_h=(L_h^T/sqrt(v)) P+Pi_h V_h.

V_h and the remaining native tapes stay private. A carrier-fixed packet comparison must integrate them conditionally on the same P. The sealed ordinary marginal pre-keep floor does not imply this comparison. The finite-dimensional taming proof instead allows H to move, proves global Lipschitz bounds for every mixed actual/tamed environment, and pays that movement explicitly.

The source papers remain unchanged. These interface conclusions were independently checked against COMMON-CARRIER-AMALGAMATION §4 and FINITE-DIMENSION-CARRIER-HYBRID §§1–4. This check does not re-prove the imported native compiler or peer theorem.
