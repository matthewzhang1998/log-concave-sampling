# Independent audit: same-root chord and nonlinear reverse-bridge gate

2026-10-04. Independent analytical review and independently implemented diagnostics. Author files were neither edited nor imported by the checker.

## Decision

**PASS** for the bounded construction and nonlinear separator, with the precise qualifications below. The author incorporated every requested correction before the final source pin recorded at the end.

- The direct conditional-packet/affine-bridge composition is a valid positive full-law fallback, at the stated second-order normalized error.
- The literal same-root local VALUE chord removes the complete matrix-quadratic second-order covariance discrepancy.
- The scalar third-cumulant obstruction survives exact affine moment matching and the fixed reverse bridge at r=0, t=1/2. It gives a genuine normalized Omega(A^2) W2 lower bound.
- The whole-root gradient lift has the claimed small full curl. Freezing u does not destroy every small-curl certificate: the projected v-source also has quadratic curl, as detailed below. The remaining consumer limitation is the unavailable joint retention of the old private-root current.
- No all-order impossibility, new general third-order transition, sublinear complete-cost recurrence, or retained-carrier consumer is proved.

## Independent deliverables

`check_chord_independent.py` is a self-contained checker using SciPy Gaussian quadrature, direct tilt-series extraction, literal VALUE evaluation, and recorded-site first replay. It never imports or runs the author checker. Its output is `chord_independent_checks.json`.

The current run passes **366 assertions**. It includes:

1. Independent first/second target cumulant coefficients from normalized partition-function series.
2. Independent direct centered-moment expansion of the actual chord, for three positive quadrature versions.
3. Literal finite-A output laws with favorable exact mean/variance repair.
4. Fixed reverse-bridge cumulants and bounded fourth-moment diagnostics.
5. Full non-diagonal matrix covariance coefficients and literal matrix VALUE identities.
6. Original-gradient moving-anchor VALUE/JVP/VJP replay and exact zeros, with maximum finite-difference error 2.522e-11.
7. Whole-root lifted curl and the separately valid fixed-u projected curl, including noncommuting Hessians.
8. A bounded sine-test witness for first-order Gaussianization failure.

The analytical arguments below, not quadrature numerics, establish the asymptotic claims.

## 1. Independent target coefficient calculation

Write U_A=A U_1, with

    U_1(x)=q x^2+epsilon(sin(kx)-kx),
    q=1/4, epsilon=1/10, k=1,
    tau=exp(-k^2/2).

The potential is anchored at zero and has Hessian between .4A and .6A. Consequently .2A x^2 <= U_A(x) <= .3A x^2. Small signed-A derivatives of Gaussian integrals are dominated by an integrable Gaussian times a polynomial, so the following Taylor calculations have genuine remainder bounds.

The log moment-generating function is

    log E_gamma exp(hZ-A U_1(Z))
      =h^2/2-A E_(N(h,1))U_1
           +(A^2/2)Var_(N(h,1))(U_1)+O(A^3).

At odd Hermite rank three, q^2 and epsilon^2 terms vanish by parity. The mixed term is controlled by

    Cov_(N(h,1))(Z^2,sin(kZ)-kZ)
      =tau[2hk cos(kh)-k^2 sin(kh)]-2kh.

Taking three h-derivatives yields

    kappa_3(X_A)
      =epsilon k^3 tau A
       +q epsilon tau(k^5-6k^3) A^2+O(A^3).

For the given fixture these coefficients are 0.06065306597126335 and -0.07581633246407918. Independent normalized tilt-series extraction agrees to floating-point quadrature accuracy.

## 2. Independent actual-chord coefficient calculation

Let h=sum_i w_i f_0(r_i Z+s_i G), with f_0=U_1', and let h_z be its actual Z-first. Both node families use the same G. The actual VALUE map has

    Y=Z-Ah+lambda A^2 h_z h+O_Lp(A^3).

Here the remainder is uniform over positive unit-mass rules: f_0 has bounded first and second derivatives, and h has a quadrature-uniform linear-growth envelope.

Set a_1=-E h and b_1=-2E Zh. Expanding the centered third moment directly gives

    [A] kappa_3(Y)=-3E[(Z^2-1)h],
    [A^2] kappa_3(Y)
      =3lambda E[(Z^2-1)h_z h]+3E[Z h^2]-3a_1 b_1.

This expansion retains the mean correction; omitting it would produce the wrong coefficient. Direct Gaussian integration gives

    E[(Z^2-1)h_z h]
      =q epsilon tau[k^5(m4+2beta m31)-3k^3 m2],
    E[Z h^2]
      =2q epsilon k[tau(1-k^2(m2+2beta m11))-1],
    E h=epsilon k(tau-1), E Zh=q.

It follows that

    [A] kappa_3(Y)=3epsilon k^3 tau m2,
    [A^2] kappa_3(Y)
      =3q epsilon tau[lambda(k^5J-3k^3m2)-2k^3K],

with K=m2+2beta m11 and J=m4+2beta m31. The epsilon^2 contributions vanish by parity, rather than being dropped as a small-epsilon approximation.

In the continuum limit,

    beta=pi/4, lambda=3/2-pi^2/8,
    K=(1+pi/2)/3, J=1/5+pi/15.

The resulting exact second-order defect is

    delta_3=q epsilon exp(-1/2)
       [3lambda(J-1)-2(1+pi/2)+5]
      =-0.009301002921086040033035066414838218609891922676302...

The independent checker obtains the author's finite-rule coefficient -0.009301002921166715, agreeing with the author result to roundoff. A positive three-point Gaussian tail variant gives -0.009301002921091053; an independent global 32-point Gaussian rule gives -0.00930141452095723. These are corroborating diagnostics, not substitutions for the exact limit.

### Finite-quadrature qualification

For the fixed author's dyadic rule with a final midpoint panel of width h=2^-24,

    m2=1/3-h^3/12

in exact arithmetic. Its first-coefficient discrepancy is therefore

    -(epsilon k^3 tau/4) h^3
      =-3.2109465768530444e-24.

This nonzero number is below IEEE double resolution in m2. Thus a fixed finite-rule difference is generally a tiny linear-in-A quadrature floor plus its A^2 defect and O(A^3) remainder; it is not literally delta_3 A^2+O(A^3) with no floor.

The proof's growing-rule statement is valid when quadrature accuracy is tied to A so that the linear floor is o(A^2), for example its stated delta_Q<=A^3. Coefficients then converge to the displayed delta_3 and the clean conclusion is delta_3 A^2+o(A^2). A positive quadrature exactly integrating r^2 also removes this particular linear cumulant floor. The independent checker reports this qualification explicitly.

## 3. Moment matching and the W2 lower bound

Both actual and target variances are 1-2qA+O(A^2), because the quadrature has exact first moment 1/2. Hence their standard-deviation ratio is 1+O(A^2). Translation does not affect third cumulants. Since the actual third cumulant is O(A), granting the exact affine repair changes it only by O(A^3).

Uniform centered fourth moments follow analytically as follows:

- The actual chord is zero at zero and has a quadrature-uniform linear-growth bound in (Z,G).
- Gaussian fourth moments are finite uniformly in A.
- The target quadratic bounds above give uniform fourth moments and a normalization bounded away from zero.
- The exact affine repair factors converge to one.

For any coupling of the two common-mean laws, write their centered variables as X and Y. If both fourth moments are at most M,

    |E X^3-E Y^3|
      <= ||X-Y||_2 ||X^2+XY+Y^2||_2
      <= [sqrt(E X^4)+(E X^4 E Y^4)^(1/4)
                            +sqrt(E Y^4)] ||X-Y||_2
      <=3sqrt(M)||X-Y||_2.

Taking the infimum over couplings yields the claimed positive W2/A^2 lower limit. This requires uniform fourth moments; W2 alone would not control a third moment.

At r=0, the exact bridge is tX+sqrt(1-t^2)N. Independent Gaussian convolution multiplies third cumulants by t^3, preserves equal first two moments, and preserves a uniform fourth-moment bound. The fixed t=1/2 bridge therefore retains one eighth of the nonzero cubic witness. No deconvolution or identification of an unobserved coupling Gaussian is used.

## 4. Matrix calibration, VALUE execution, and positivity

For f(x)=Bx, every node uses the same G, and

    K=BZ/2+beta BG,
    E=-(B/2)K.

Thus the literal map has coefficients

    C_Z=I-B/2+lambda B^2/4,
    C_G=-beta B+lambda beta B^2/2.

Its covariance is C_Z C_Z^*+C_G C_G^*, with B^2 coefficient

    1/4+beta^2+lambda/2=1.

The remaining B^3 and B^4 coefficients in the proof are correct. The error relative to (I+B)^-1 is O(alpha^3), and the fixed covariance gap converts that to O(alpha^3 sqrt(D)) Gaussian W2. This is a full matrix calculation; it is not merely a trace test. The signed chord arithmetic is an ordinary deterministic transformation of Gaussian roots and hence defines a positive probability law.

There are two complete node banks. The second evaluates original gradients at genuinely shifted points depending on the entire first K. It does not execute a Hessian-vector product as a VALUE source. In a requested first or adjoint, only original first matrices at the stored VALUE sites occur. Differentiating such first matrices is unnecessary.

## 5. Whole-root lift and the corrected fixed-caller distinction

With R_i=(r_i I,s_i I), P=(I,0), and Pi=P^*P, set

    J(W)=sum_i (w_i/r_i)R_i^* f(R_iW),
    E_J(W)=J(W-P^*K(W))-J(W).

Then J is a genuine full gradient, PJ=K, and PE_J=E exactly. Writing B0=DJ(W), B1=DJ(W-P^*K(W)),

    D E_J=B1-B0-B1 Pi B0,
    D E_J-(D E_J)^*=-B1 Pi B0+B0 Pi B1.

Although ||Bj|| can cost alpha L_Q, the product factors have

    ||Bj P^*||=||P Bj||=||DK||<=alpha.

The full skew is therefore at most 2alpha^2, without an L_Q penalty. The stated first bound 2alpha L_Q+alpha^2 is valid. The stated energy bound with L_Q is conservative: integrating DJ along the P^*K displacement actually gives ||E_J||_p<=alpha||K||_p.

Any imported mean consumer must apply its small-radius, covariance-gap, numerical, clock, and complete-replay guards to the actual lifted radius 2alpha L_Q+alpha^2, rather than assuming alpha<=1/2 suffices. Lifted captured-caller and numerical ports also retain their actual L_Q coefficient bill.

A review correction is important when u is frozen as caller. Let K_u0,K_v0 and K_u1,K_v1 denote packet firsts at the original and shifted root points. All of these are symmetric because each is a positive scalar combination of symmetric local force derivatives. Consequently

    D_v E=K_v1-K_v0-K_u1 K_v0,
    D_v E-(D_v E)^*=-K_u1 K_v0+K_v0 K_u1,
    ||D_v E-(D_v E)^*||<=alpha^2 beta_Q.

So the projected D-dimensional v-source retains small curl too. Its full first can still be O(alpha), and its origin E(u,0) and energy are caller-dependent. An anchored conditional mean implementation must retain and price that caller-only origin graph.

Neither certificate supplies the missing joint observation. A global completed mean law integrates out u,v. A completed conditional mean law at fixed u integrates out its private v. The chord's later expression still reads the same K(u,v). Reading an independent replacement K, or appending the old private v/K to a marginal mean-law comparison, is not licensed by that comparison. This law-only retention restriction, not a false loss of projected curl, is the valid consumer boundary.

## 6. Conditional caller, zero, query, and replay ports

For f(y)=s[g(x_M+sy)-g(x_M)] and alpha=As^2, the proof correctly retains every occurrence of the captured mode and anchor. The first of the mode in a=rz has norm at most 2. Hence D_a f and D_a K are O(As), and the K1 feedback introduces only K_u1 with norm at most alpha/2. This gives the stated absolute caller bounds and

    ||D_(u,v) Q_chord-sP||=O(s alpha)=O(As^3).

These are actual graph derivatives, not derivatives of law-error estimates. Physical original-caller derivatives carry their imported sqrt(A) factor; full lifted output readouts additionally carry L_Q where applicable.

The complete conditional VALUE count is M+1+2n_alpha, before a separately requested terminal force or downstream completed consumer. The original global-mode ledger remains an imported cost, not a newly free service. A repeated complete source occurrence pays its own two private node banks; only the identical captured caller-only mode/anchor graph is reusable. The bridge adds one fresh D-Gaussian root when its coefficient is positive, for raw dimension 3D, otherwise 2D at t=1. Completed mean consumers have their own larger complete tapes.

At all private zeros, same-site reuse of the recorded anchor makes f, K, and E exactly zero. The physical packet returns its captured mode. Requested first/adjoint sweeps visit the stored original VALUE sites, including the K-to-K1 feedback; discarded primals incur replay. The independent checker verifies these facts for a moving physical anchor and noncommuting Hessians, in addition to the analytic chain-rule review.

## 7. Direct affine fallback and scope

For the existing conditional packet Q,

    W2(Law(Q|z),nu_(r,z))<=C A^2 s^5 sqrt(D)+eta_Q.

Use the same exposed caller and a fresh Gaussian N in the affine bridge. Scaling the conditional coupling by B_rt=Delta/(t s^2) gives exactly

    W2<=C (Delta/t)A^2 s^3 sqrt(D)+B_rt eta_Q.

The complete actual non-Gaussian law passes through this composition. It does not claim that a W2 bound by itself bounds every higher cumulant, nor that all actual cumulants equal the target's. The zero reserve at t=1 causes no division problem. The mean-only or mean/covariance Gaussian reference has the separately demonstrated first-order bounded sine-test discrepancy.

The scope is appropriately narrow: this specific quadratic-calibrated local chord fails the next nonlinear full-law order even with ideal affine repair. Its failure says nothing prohibiting a different positive carrier-retained current compiler. No recurrence improving the general full-law grade has been established.

## Final source disposition

**PASS.** The final source incorporates all review corrections: growing-quadrature qualification, corrected continuum decimal, actual lifted consumer guards and caller/precision bills, and the valid fixed-u projected small-curl identity with the actual private-root retention limitation.

Audited source:

- `../SAME-ROOT-CHORD-AND-FIRST-NONLINEAR-BRIDGE-GATE.md`
- SHA-256: `fc248e0a5b145cd325e5dae47c9a2728a0b28b6ebc8a9cdb85344044125bbb5b`

This verdict covers the positive affine fallback, literal chord, matrix calibration, scalar nonlinear separator, and stated VALUE/caller/zero/replay ports. It does not promote either marginal mean-law consumer to an unproved joint retained-carrier law.
