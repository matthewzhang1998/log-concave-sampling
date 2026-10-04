# Independent review: dual-completed ambient lift and common-noise defect

2026-10-04. Reviewed `FENCHEL-INVERSE-LIFT-AND-LITERAL-K3-GRAPH-CURRENT.md` at final SHA-256 `7357f4e314c2d8024baea2a6cb1fa8021762728a6456a951f1831e7c5b78cbe4` (verified 17:40 UTC). The main positive claim and the retained obstruction are correct. This review is complementary to `FENCHEL-GAP-INVERSE-AND-CONSTRAINED-LAW-AUDIT.md`.

## Verified claims

- Adding omega[V0*(p0)-V0*(p1)] to the corresponding signed primal terms makes the graph p blocks ra M*gt(t0) and -ra M*gt(t1). Their sum is exactly a M*E. The old unmarked r M*d term cancels.
- At finite inverse depth N, the common residual is exactly r omega[e_N(p0)-e_N(p1)]. The same finite circuit/cached arithmetic makes it vanish literally when p0=p1. A small absolute VALUE error does not give a q^N-small derivative or relative-E error under the supplied Hessian bounds.
- For g0=gt=m id and M=I, E=ra²m³ epsilon Z. The opposite-mode expression (16) is exact. Its L2 norm is at least sqrt(2)ra m sqrt(d) when a,c,s>0. Thus its ratio to native E energy is at least sqrt(2)/(a m² epsilon), and it is nonzero at epsilon=0.
- The literal W=(S,U,Z) density with product conditional Fenchel fibers has score

      -W + sum_i (D_W x_i)*[p_i-m_tau(x_i)]/tau.

  This includes the actual nonlinear derivative of x1. Without conditional normalization, replace m_tau(x_i) by g0(x_i), and the W marginal is reweighted by the product of the two normalizers.
- Independent product fibers fail at epsilon=0: their K3 energy covariance is 2r²a²m³tau I, despite the zero native field.
- The shared-noise graph p_i=g0(x_i)+sigma V produces exactly a common terminal displacement a sigma M V. It is the existing common terminal heat operation, not an additional distinct smoothing mechanism. The estimate |E_sigma-E|<=2ra sigma ||M||||V|| is valid for Lip(gt)<=1.
- Since the signed Fenchel gaps vanish identically along the exact feature graph, their full pullback is zero. The old natural terminal pullback and its HVP-valued ancestor companion therefore remain.

## Two finite-budget corrections applied and verified in the pinned version

1. The unnormalized common-sum readout has norm sqrt(2). Therefore the residual in (15) is bounded by sqrt(2) times the displayed full-source error (11), or directly by

       (5/2)r omega q^N || |p0|+|p1| ||_Lp.

   The normalized common coordinate (p0+p1)/sqrt(2) has the full-source bound without that factor.

2. For a literal total inverse tolerance delta_h, truncation and primitive VALUE noise must share the budget. A concrete sufficient choice is

       N >= max(0,ceil(log(5 P_p/delta_h)/log 5)),
       delta_g <= delta_h/5.

   These give truncation <=delta_h/2 and accumulated call error <=delta_h/2. Choosing truncation <=delta_h and allowing positive additional noise gives a constant multiple of delta_h instead, not the literal stated tolerance. Finite arithmetic errors also require a share of the budget.

## Sharp elementary bound for the inverse-constraint difference

For a fixed shared displacement delta, put

    F_delta(x)=h(g0(x)+delta)-x.

It obeys

    |F_delta(x)| <= (5/2)|delta|,
    DF_delta(x)=B^{-1}(A-B),
    A=Dg0(x), B=Dg0(h(g0(x)+delta)).

Since A,B lie between .4I and .6I, ||A-B||<=.2, so

    Lip F_delta <= .5,
    |F_delta(x0)-F_delta(x1)|
      <= min{5|delta|, .5|x0-x1|}.

These are actual first identities and bounds. No derivative of a Hessian, and no HVP-valued VALUE, appears. Multiplication by r omega gives the exact corresponding bounds for the extra common companion (29).

## Explicit obstruction to a uniform product of the two small factors

For any k>0 use the anchored scalar primitive

    g_k(x)=x/2+sin(kx)/(10k),
    x0=0, x1=pi/(2k),
    delta=(pi/4+.1)/k.

Let C=pi/4+.1. Then g_k(x1)=delta, so

    F_delta(0)=x1.

Define y by y/2+.1 sin y=2C=pi/2+.2. Its derivative lies in [.4,.6], hence

    1/3 <= y-pi <= 1/2.

The other inverse value is h(g_k(x1)+delta)=y/k, and consequently

    |F_delta(0)-F_delta(x1)|=(y-pi)/k
      lies in [1/(3k),1/(2k)],
    |delta| |x1|=C pi/(2k²).

There is therefore no uniform constant multiplying |delta||x0-x1| that controls this defect over the original Hessian-bounded class. Both quantities tend to zero; the defect is order 1/k while their product is order 1/k². This supports retaining (29) as a separate actual field rather than assigning it an unproved product mark.

## Scope

The corrected common mode is a real constructive improvement. Neither this review nor the independent audit rejects an exact inverse reference with coherent finite VALUE restoration, provided the receiving theorem actually supports that restoration. The remaining baseline and constrained-law current are substantive and are correctly retained in the main artifact. No all-rank closure is established here.

## Verified extension: centered actual-energy lower bound and inverse floor

The final main artifact at the SHA pinned above adds a correct centered-energy estimate. It specifically requires M=I and BOTH g0 and gt to have Hessians in [.4I,.6I]. For positive r,a,epsilon, literal differentiation of the existing query graph gives

    D_Z E = ra² epsilon Ht(t1) H0(x1) Ht(S+epsilon Z).

No HVP is used as an executed VALUE in this analytical identity. Write each of the three ordered factors as .5I+K_i, with ||K_i||<=.1. Without commuting any factors, the seven nonconstant terms have summed norm at most

    3(.5²)(.1)+3(.5)(.1²)+.1³=.091.

Thus the product is .125I+R, ||R||<=.091, and its symmetric part is at least .034I.

Here is the full conditional-chaos argument. Fix S,U, let f(Z)=E(S,U,Z), and let

    B(S,U)=E_Z[f(Z)Zᵀ]=E_Z[D_Z f(Z)].

The equality is coordinatewise Gaussian integration by parts and uses only the bounded first just displayed. Therefore sym B >= .034 ra² epsilon I. Since the coordinate functions Z_j are orthonormal in Gaussian L2, Bessel's inequality applied to each component f_i gives

    E_Z|f-E_Z f|² >= ||B||_F²
      >= sum_i B_ii² >= d(.034 ra² epsilon)².

The law of total variance then gives the claimed centered bound

    e_E=||E-E E||_2 >= .034 ra² epsilon sqrt(d).

This is stronger for the stated purpose than only bounding the uncentered norm; no noncancellation assumption is hidden in it. The separate pointwise strong-monotonicity bound |E|>=.4³ ra² epsilon |Z| is also correct.

For the original graph, anchoring and Lip(g0),Lip(gt)<=.6 give

    ||p0||_2<=.6 sqrt(d),
    ||p1||_2<=.6(1+.6a epsilon)sqrt(d).

For a epsilon<=1 their sum is at most 1.56 sqrt(d). Consequently the common inverse residual satisfies the concrete bound

    ||r omega[e_N(p0)-e_N(p1)]||_2
      <=3.9r omega q^N sqrt(d)
      <=(3.9/.034) omega q^N/(a² epsilon) e_E.

The common-noise changed inputs remain O(sqrt(d)) when sigma is bounded, with their actual larger profile retained. Similarly (28) and ||V||_2=sqrt(d) give

    ||E_sigma-E||_2 <=(2/.034) sigma/(a epsilon) e_E.

These ratios are for epsilon>0. At epsilon=0 the literal common finite-inverse discrepancy cancels by identical-query arithmetic instead of dividing by zero energy. To make the inverse residual at most rho e_E, it is enough, for the displayed original-profile constants, to choose

    q^N <= (.034/3.9) rho a² epsilon/omega,
    N >= max(0,ceil(log(3.9 omega/(.034 rho a² epsilon))/log 5)),

and to allocate provider/arithmetic errors to the same residual budget. At any fixed polynomial grade, with polynomially parameterized a and epsilon, this is logarithmic iteration work and therefore contributes b_extra=0 to the inverse-heat copy exponent. It supplies no small derivative of the residual, no opposite-baseline cancellation, and no constrained-current theorem.

An additional 400 random noncommuting matrix-triplet diagnostics passed; the proof is the ordered expansion and conditional-chaos argument above, not the diagnostics.

Final scope clarification verified: L_F is a signed saddle reference, not a jointly convex or asserted normalizable sampling potential. Its S,d Hessian cross block is nonzero while the d,d block is zero, which precludes joint convexity. This caveat changes no formulas or previous verification.
