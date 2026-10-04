# Two precise remaining-root and calibration boundaries for common auxiliary heat

New supplement to COMMON-AUXILIARY-TERMINAL-HEAT-SCREEN.md, held b03b0c2f. Its program, assumptions and six-VALUE cost are unchanged.

## 1. The coherent-drift root field keeps its original orientation grade

Let m_sigma(W)=E_V E_sigma(W,V). For the actual held coherent-drift family f8363e7d (and its uniform finite-K extension a0ee7ab3), the original centered energy is comparable to kappa0=ra, and

    ||m_sigma-E||_2 <= C a sigma e,
    ||O(m_sigma)-O(E)||HS <= C a sigma kappa0 e.

The full ambient or conditional-gradient program need not be invoked to establish these inequalities; they follow from the literal common terminal heat and the output-swap Hilbert identity in b03.

The same held source has O(E)_00>=c kappa0 e. Hence, for sufficiently small a sigma,

    O(m_sigma)_00 >= (c/2) kappa0 e.                         (1)

This is a same-grade remaining ORIGINAL-root field. An arbitrarily accurate conditional gradient mean in V has not eliminated the original orientation; it has left (1) in its W-dependent conditional mean. This conclusion does not exclude a later joint W-current correction. It identifies the exact field that correction must own.

The extension to every declared finite K>=3 is uniform: b03's pointwise mixed-difference argument uses only the final common terminal displacement, while a0ee7ab3 provides the uniform old energy and orientation lower bounds. No union bound or new finite-K derivative claim is used.

## 2. Bounded Hessians alone do not give small relative VALUE restoration

There is a separate actual scalar native K3 fixture showing why the Hessian modulus in the positive drift screen cannot be silently dropped. It is NOT a counterexample to the weaker kappa e coefficient tolerance.

Choose c=s=1/sqrt(2), M=1, r=a, epsilon=a^.9 and SAME primitive at both original nodes

    g(y)=m y+(d/k)sin(ky),  m=1/2, d=1/10.

Its Hessian lies exactly in [.4,.6]. Let sigma=a^p with 0<p<.9, and choose a fixed q satisfying

    1+p < q < 1.9,   k=a^(-q).

The source uses the exact S,U,Z Gaussian records and complete feedback:

    x=cS+sU,
    h=a[g(S)-g(S+epsilon Z)],
    t0=S+a g(x),   t1=S+a g(x+h),
    E=r[g(t0)-g(t1)].

Use the same common terminal heat as b03. Its conditional mean is EXACTLY

    m_sigma=r{m(t0-t1)+(d/k)exp[-(k a sigma)^2/2]
                                         [sin(kt0)-sin(kt1)]}.       (2)

Since k epsilon tends to infinity while k a epsilon tends to zero, the original feedback and ancestor finite differences give, in L2 after division by r a² epsilon,

    E = m Z [m+d cos(kx)] [m+d cos(kt0)] + o_L2(1).                 (3)

Here and in the next display the left side denotes the normalized field. For clarity, the errors are bounded by constants times

    (k epsilon)^(-1) + k a epsilon + k a² epsilon,

using bounded Gaussian moments and the explicit sine Taylor remainder. Every displayed error tends to zero in the stated q window.

The terminal phase satisfies kt0=kS+a m kx+a d sin(kx). The final term tends to zero. The two remaining linear Gaussian phases (kx,kS+a m kx) become independently uniform modulo 2pi, jointly with the independent Z: their covariance has smallest eigenvalue comparable to k². Thus ordinary bounded trigonometric moment calculations give

    e/(r a² epsilon) -> m(m²+d²/2),

where e is the ACTUAL centered energy. Since k a sigma tends to infinity, (2) gives

    m_sigma/(r a² epsilon)
                     =m² Z[m+d cos(kx)]+o_L2(1).

Consequently

    ||m_sigma-E||_2/e -> d/sqrt(2m²+d²)
                              =0.1400280084... .                  (4)

The source and the heat use the SAME bounded-Hessian gradient, the same feedback, and no independently resampled original root. Therefore a uniform claim ||m_sigma-E||_2=o(e) cannot follow for this whole source class merely from sigma tending to zero.

There is no conflict with section 1 or with the desired weaker orientation port. In this fixture e is comparable to r a² epsilon, while kappa0=ra. The Hilbert orientation bound ||O(F)||HS<=2||F_centered||_2², together with the pointwise mark preservation, already puts the entire possible orientation error at O(e²)=O(a epsilon kappa0 e). This is smaller than the original kappa0 e scale. Formula (4) only blocks a stronger relative-VALUE restoration premise; an appropriate SAME-E weak/current theorem may still succeed.

## 3. Costs and scope

The additional parameter k changes the allowed original potential, not the number of source queries. The source remains six original VALUES, known scalar rows and the same one auxiliary Gaussian. No query-count or precision uniformity in k is asserted beyond the supplied original-oracle model; finite encoding and arithmetic retain their actual costs. The fixture uses a smooth potential with uniformly bounded Hessian but no uniform third-derivative bound, exactly the original first-only scope.

Section 1 is a definite remaining-root obligation for the naive conditional mean route. Section 2 explains the correct weaker calibration target. Neither section substitutes a joint Gaussian law on the original source variables or supplies the unresolved matrix current.
