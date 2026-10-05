# Current-carrying rank-two resummation with a genuine retained-mark alternative

2026-10-05. Constructive finite-source lemma. This advances the mark-consumption boundary of the existing two-replica covariance transport; it does not construct the still-missing true-history base/singleton source or prove the next mean induction.

## 1. Executive result and scope

A covariance transport can carry a small-energy history mark without declaring its first derivative small. There are two different, explicitly distinguished contracts:

1. A genuine joint coupling retains the original mark unchanged and compares the transported endpoint with its averaged Gaussian reference. This has a first-order covariance-fluctuation cost, rather than the smaller unmarked marginal cost.
2. An exact positive-likelihood representation transfers every linear marked observable to one fixed Gaussian carrier. Its new mark is the old coherent mark times an explicit positive rank-two Gaussian likelihood. It preserves the current exactly, not the distribution of the old mark. The likelihood and its contractions are executable in linear vector work and a two-dimensional eigensolve. Its L2 amplification is independent of ambient dimension under a numerical rank-two guard.

The second contract admits finite Hermite current truncations of any chosen degree. Their original-gradient call count is unchanged; only scalar/two-dimensional polynomial work increases. The truncations are current estimators carried by a positive joint pushforward, not positive density approximants or unweighted samplers of the old law.

The true-history mark is not available as an executing leaf. A supplied finite original-VALUE mark graph is required for execution. The construction therefore supplies an actual missing local join, not a falsely complete history theorem.

## 2. Finite input graph and ownership

Retain a caller Y and genuinely unread labels. Let Omega be a finite source bank. It produces two capped vectors F,F' in R^d, together with a finite mark M in R^m. All of their original-gradient ancestors are executed and recorded. The mark may share any source roots and saved ancestors with F,F'; it is not made independent artificially.

Assume the two vectors satisfy |F|,|F'|<=R. Set

    H=(F F'^T+F' F^T)/2,
    B=H/(2v), K=I+B,
    C(Y)=E[H|Y].

Assume C(Y) is PSD, as it is when F,F' are conditionally independent nuisance replicas of one source given their common outer bank. A NEW standard d-Gaussian G is independent of the ENTIRE source bank Omega conditional on Y. The mark M and H are Omega-measurable given Y; no independence of M and H is needed. All retained coordinates for ordinary product-space W2 have finite second moments. The numerical guard is

    h=R^2/(2v)<=1/8.                                      (2.1)

Since ||H||op<=R^2, K is positive definite, rank(B)<=2, and its eigenvalues lie in [1-h,1+h]. Execute

    Z=sqrt(v) K G.                                        (2.2)

This is precisely one positive Gaussian-root pushforward. It uses the same H as the established pair producer and keeps M on its actual source bank. Its conditional covariance is vI+C+E[H^2|Y]/(4v), including the true source/mark cross ancestry.

## 3. A genuine retained-mark joint-law comparison

Define Zstar=(vI+C(Y))^(1/2)G, with the SAME G in the proof. The square root is an analytical comparison only. Conditional on Y, Zstar is independent of Omega and therefore of M. Retain M unchanged under the coupling. Then

    W2(Law(M,Z|Y), Law(M,Zstar|Y))
      <= ||H-C||_(L2(Omega;HS)|Y)/(2sqrt(v))
         + ||C^2||HS/(8v^(3/2)).                         (3.1)

Here W2 on the product space uses its ordinary Euclidean cost; the mark displacement is exactly zero. The estimate also holds with any extra retained source label in place of M, since that coordinate never moves. It is not the order-H-squared unmarked marginal comparison.

Proof: insert sqrt(v)(I+C/(2v))G. The first difference has squared mean norm E||H-C||HS^2/(4v). For every scalar c>=0,

    0<=sqrt(v)+c/(2sqrt(v))-sqrt(v+c)<=c^2/(8v^(3/2)).

Apply spectral calculus to C and use the common G. Minkowski proves (3.1). No differentiation of the mark or of a law error occurs.

For physically independent blocks of dimension at most b, with ell-Lipschitz conditionally centered capped replica sources, C_j<=ell^2 I. If each block cap is R_j<=ell(sqrt(d_j)+sqrt(2L)), then the integrated version of (3.1) is at most

    C ell^2 (sqrt(b)+L) sqrt(D)/sqrt(v)
      + ell^4 sqrt(D)/(8v^(3/2)).                       (3.2)

The first term follows from E||H_j-C_j||HS^2<=E||H_j||HS^2<=R_j^4 and sum_j(d_j+2L)^2<=C(b+L^2)D. The displayed bound is deliberately loose. All source roots remain visible to the retained mark; none is silently consumed.

For the existing pair source ell<=C A^2, this is O(A^4(sqrt(b)+L)sqrt(D)/sqrt(v)). It is a valid stronger joint contract with a larger error than the previous marginal A^8 contract. After an A-Lipschitz readout it costs O(A^5(sqrt(b)+L)sqrt(D)/sqrt(v)); that arithmetic fact does not remove any other induction debt.

## 4. Exact Gaussian likelihood current transfer

Let X=sqrt(v)G with the new G independent of Omega. Define

    L_H(X)=det(K)^(-1)
           exp{ X^T(I-K^(-2))X/(2v) }.                 (4.1)

The determinant and nontrivial quadratic form live only on span(F,F'). If that span has rank zero or one, use its actual rank. Formula (4.1) has no conditional moment, derivative, covariance-square-root, or unknown density oracle.

Execute the ordinary positive joint pushforward

    (X, Msharp, Wsharp)=(sqrt(v)G, M L_H(X), L_H(X)).    (4.2)

For every integrable scalar/vector test phi and almost every Y,

    E[M phi(Z)|Y]=E[Msharp phi(X)|Y],
    E[phi(Z)|Y]=E[Wsharp phi(X)|Y].                     (4.3)

The same identity holds with grad phi, Hess phi, or any tensor-valued test when the displayed expectations exist. Thus M=-q Delta carries the exact linear current -q E[Delta grad phi(Z)|Y] onto the fixed Gaussian carrier. It retains all correlations between Delta and the original H. The entire source genealogy appears inside Msharp; it is not replaced by an independent mark bank.

Proof: conditional on (Y,Omega), (4.1) is the ordinary density ratio of N(0,vK^2) with respect to N(0,vI). Change variables x=sqrt(v)K g. Then integrate Omega. This proves more generally the exact positive joint-measure identity

    E[f(Omega,Z)|Y]=E[L_H(X) f(Omega,X)|Y].             (4.4)

Equation (4.4) is a weighted representation. In the unweighted executing law (4.2), the marginal of X is N(0,vI). It is NOT the old transported law, and Msharp does NOT have the original mark distribution. No rejection sampler, importance-weight removal, or conversion of a LAW return to RAW noise is being claimed.

## 5. Exact energy calculation; no false small first

Let E=K^2-I=2B+B^2 and e=2h+h^2. Under (2.1), e<=17/64<1. Direct Gaussian integration gives

    E_G[L_H(X)]=1,
    E_G[L_H(X)^2]=det(K^2(2I-K^2))^(-1/2)
                =det(I-E^2)^(-1/2).                   (5.1)

The determinant has at most two nonunit factors, so

    E_G[L_H^2]<= (1-e^2)^(-1),
    E|Msharp|^2 <= (1-e^2)^(-1) E|M|^2,               (5.2)
    E|Msharp-M|^2 <= [(1-e^2)^(-1)-1] E|M|^2.          (5.3)

In (5.3), M is evaluated on the same Omega and independent X. Equations (5.2)-(5.3) are uniform in all source correlations and retained Y, and also hold conditionally whenever the mark energy is finite. They require no bound on the derivative of M. In particular a genuine small-energy mark O(A^k sqrt(D)) stays that small under this local current transfer. It is never normalized as though its first were O(A^k).

For any p>=1 with (p-1)e<1, the exact formula is

    E_G[L_H^p]
       =det(I+E)^(-(p-1)/2)
          det(I-(p-1)E)^(-1/2).                       (5.4)

This makes all used higher-moment thresholds explicit. Source dependence on Y remains captured. There is no globally bounded first: differentiating (4.1) produces Gaussian-polynomial factors times an exponential. Consequently this primitive is not admitted to a bounded-first native port without a new clipping/derivative proof.

## 6. Finite current polynomials of any selected order

Write He_n for the standard Gaussian Wick tensor. The likelihood has the L2-convergent orthogonal expansion

    L_H(sqrt(v)G)=sum_(k>=0) T_k(E,G),
    T_k(E,G)=(E^(tensor k):He_(2k)(G))/(2^k k!).        (6.1)

Only the rank-at-most-two range of E is read. Define L_[p]=sum_(k=0)^p T_k and execute M_[p]=M L_[p]. This is a finite mark graph over the positive Gaussian/source bank, with possibly signed mark coordinates. It is not a positive density approximation when p>=1.

The exact norm generating identity is

    sum_(k>=0) t^k ||T_k||_(L2 G)^2
        =det(I-t E^2)^(-1/2).                         (6.2)

Taking t=1/(2e^2) when e>0 yields

    ||L_H-L_[p]||_(L2 G)
        <=2^(1/2) (sqrt(2)e)^(p+1).                   (6.3)

For e=0 the likelihood is exactly one. Because M is fixed when conditioning on Omega,

    ||Msharp-M_[p]||_2
       <=sqrt(2)(sqrt(2)e)^(p+1)||M||_2.              (6.4)

Therefore every test with uniform norm <=1 has a current-observable error bounded by the right-hand side of (6.4), including gradient tests with ||grad phi||infinity<=1. The scalar density representation has the same bound with M=1. This is a relative-energy certificate; no mark-first estimate or cumulant matching is used.

For actual eigenvalues lambda_1,lambda_2 of E and independent standard coordinates z_1,z_2 of G in their eigenspace,

    L_[p]=sum_(a+b<=p)
       lambda_1^a lambda_2^b He_(2a)(z_1) He_(2b)(z_2)
       /(2^(a+b) a! b!).                              (6.5)

Hermites are generated by He_(n+1)(z)=z He_n(z)-n He_(n-1)(z). Compute all factors through degree 2p and the triangular sum. Cost is O(p^2) scalar operations after O(d) projections. The original-VALUE and Gaussian-root counts are exactly unchanged from the underlying H,M graphs.

## 7. Literal implementation and costs

Given F,F', form an orthonormal basis U of their span by a rank-revealing two-vector QR. On this at-most-two-dimensional span form

    Bsmall=U^T H U/(2v), Ksmall=I+Bsmall.

For G compute z=U^T G, solve y=Ksmall^(-1) z, and evaluate

    log L=-log det Ksmall+(|z|^2-|y|^2)/2.             (7.1)

Then L=exp(log L). No D-by-D matrix is materialized. Basis orientation is immaterial. At degeneracy use the lower actual rank or continuous invariant formulas. In particular do not divide by a tiny Gram determinant in an unguarded implementation.

If the fully expanded underlying graph costs Q_H VALUES and d_H roots and the same-bank finite mark adds Q_M new VALUES and d_M new roots, the exact current packet costs

    Q=Q_H+Q_M,
    roots=d_H+d_M+d,
    arithmetic=underlying work+O(d+m),                (7.2)

where shared exact-key original calls and roots are counted once. The finite p current uses O(p^2+d+m) additional arithmetic and no additional source calls. For the existing fourteen-VALUE pair graph, Q_H=14 and its prior total 17d roots already includes G; hence the new packet is 14+Q_M VALUES and 17d+d_M Gaussian coordinate roots, with no extra carrier, plus the original scalar clock randomness and any separately requested inactive-slot buffer. Unknown true histories may not be inserted into Q_M at unit cost.

As a completely finite local instantiation, use the five-VALUE coherent mark Delta=U-V from BOUNDED-SHIFTED-CURRENT.md on affine rows drawn from the same saved source bank. This gives at most nineteen original VALUES in total (fewer only with exact-key reuse), and its proved energy ab A^3||x0||_2 is multiplied by at most (1-e^2)^(-1/2). This is the actual local mark, not the genuine integrated history Delta3.

Changing Y, a clock, smoothing/secent scale, or any affine row replays all affected ancestors in H and M. Changing a source Gaussian argument replays every dependent original VALUE. No old completed-law carrier is subtracted.

## 8. Numerical execution and thresholds

The exact-real formulas require only finite VALUES, Gaussian roots, vector arithmetic, a size-two spectral solve, log, and exp. A finite-precision implementation pays a Gaussian tail before bounding the exponential. With rank r<=2, on |U^T G|<=B_g,

    |log L|<=r[-log(1-h)]+c_h B_g^2,
    c_h=max_{|b|<=h}|1-(1+b)^(-2)|/2.

The derivative with respect to Ksmall is bounded by a fixed polynomial in (1-h)^(-1), B_g and L on this event. The Gaussian tail is paid in the L2 norm using a moment (5.4) strictly above two and Holder; choose its p so (p-1)e<1. Under h<=1/8, p=3 is allowed. The required B_g is logarithmic in the absolute tolerance and the used source/mark moment profile. Gram near-degeneracy is handled by a stable spectral/QR method or by computing the same matrix function from the rank-two update without a small denominator. An eigenvalue/rank threshold contributes its explicit operator perturbation, not an uncharged deletion.

For finite polynomial currents, fixed-degree Hermite recurrences have only polynomial Gaussian profiles. Evaluate on a declared root cap and pay Gaussian moment tails. In either case choose original-VALUE, cap, clock, solve and scalar-function floors against every actual downstream amplification. The substantive errors (3.1) and (6.4) are not numerical floors.

## 9. What this does and does not resolve

The exact local covariance graph no longer has to discard a coherent mark entirely. One can either keep its full law and pay (3.1), or carry its linear current on a fixed Gaussian carrier with (4.3) and controlled energy (5.2). Arbitrary selected current orders in (6.4) do not add source queries. These are repeatable local contracts whenever the next node really has a fresh independent carrier and a qualified finite source/mark graph; at a changed carrier every dependent mark ancestor must be replayed.

It would be invalid to infer an arbitrary-depth history invariant solely from these identities. The missing base/singleton fields still contain the whole true nonlinear ancestor. No source graph for those fields has been constructed here. The fourteen-VALUE pair supplier still has its original bounded-block actual-history comparison and does not acquire a dimension-uniform history comparison from this lemma. The likelihood graph also lacks a bounded native first port, and weighted current identities do not become unweighted sampling laws without another construction.

The proposed full m4 task is therefore unresolved. The concrete new deliverable is an executable marked covariance/current join with explicit actual joint-law scope, positivity, costs and guards, rather than another assertion that the old marginal theorem somehow preserves marks.

## 10. Bounded native continuation

RETAINED-CARRIER-NATIVE-CURRENT-PORT.md retunes the pair secant, clips the first Hermite current correction, proves its complete/caller first and generic full-bank curl, and applies the imported zero-baseline native mean service while retaining G as a caller. The exact likelihood remains unbounded; that continuation does not submit it directly to a native compiler.

The likelihood energy formula is a SINGLE rank-at-most-two block statement. Multiplying likelihoods from many physical blocks generally increases the rank and can amplify the second moment with the number of blocks. There is no dimension-free global product-likelihood claim. The native continuation acts on block-local current fields; arbitrary globally coupled callbacks need a separate joint join.
