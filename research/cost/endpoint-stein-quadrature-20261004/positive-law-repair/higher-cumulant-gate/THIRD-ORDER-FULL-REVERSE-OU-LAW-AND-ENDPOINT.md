# A third-order full reverse-OU law and finite endpoint

2026-10-04. Integrated constructive candidate, independent audit in progress. This is a fixed general-C2 order gain with the imported public-log guards, not an all-order recurrence.

## Result first

Let g=grad U, g(0)=0, 0<=Dg<=A I, with A sufficiently small for the explicitly imported fixed-order mean/square-action guards. There is a finite positive original-gradient VALUE program whose full endpoint law obeys

    W2(output,mu_U)<=Lambda(A) A^3 sqrt(D)
                         +restored mode/numerical floors,        (1)

where mu_U is proportional to exp(-|x|^2/2-U(x)) and Lambda is a fixed public-log polynomial. The complete original VALUE count is also a fixed public-log polynomial. Every higher conditional cumulant and shared-root descendant is included in a full-law comparison, rather than declared absent after moment matching.

The notation Lambda(A) suppresses the other declared public logarithms, including actual active dimensions and precision/caller parameters where the imported ports contain them. All their literal guards are retained; no algebraic dimension-versus-heat restriction is introduced or silently assumed.

In the original physical posterior normalization, multiplication by sqrt(A) gives grade A^(7/2) up to the displayed public logarithms, improving the previous A^(5/2) endpoint grade. Original finite-mode and numerical restoration terms retain their actual physical weights and caller profiles. No strong mean oracle, Hessian-valued producer, signed probability, or tensor oracle is introduced.

The fixed power three is the new claim. No general c(P)/P->0 theorem or recursive higher-order source admission follows from (1).

## 1. Audited component interfaces and labels

The construction uses these separately proved interfaces, with their actual scopes:

- `EXACT-TWO-CLOCK-CURRENT-AND-UNIFORM-C2-REFERENCE.md`: stationary-OU second substitution P2=Z-F_cont has W2<=A^3 sqrt(D) to mu_U, uniformly under C2. It is an analytical reference only.
- `THIRD-ORDER-NESTED-FORCE-MEAN-WITH-C2-DECOUPLING.md`: a finite shared-root nested VALUE source, completed by gradient/near-gradient mean banks, gives M(Z) with conditional target N(m_cont(Z),I), integrated law error C A^3 sqrt(D)+Lambda A^4 sqrt(D)+floors, where m_cont=E[F_cont|Z]. Z is a fresh standard Gaussian carrier exposed before all its private banks.
- `clock-compression/POSITIVE-SQUARE-CLOCK-COVARIANCE-SERVICE.md`: a positive VALUE square-clock action targets

      C(Z)=Cov(H_cont|Z)=2 integral_0^1 q P_q[(Dv)^2](Z)dq,
      H_cont=integral_0^1 g(X_r)dr,
      v=integral_0^1 P_r g dr.

  Its positive finite quadrature has L2(Z;HS) error <=1.5 delta A^2 sqrt(D). Its completed positive Gaussian reserve has error Lambda A^3 sqrt(D) at padding mu=A, with a fixed positive baseline variance. The covariance action includes all owned X=qZ+sqrt(1-q^2)G roots in its actual first/energy/parameter port.
- `CUBIC-GAUSSIANIZATION-AND-POSITIVE-BUFFERED-LAW-RETURN.md`: for any centered Gaussian-Lipschitz source of first L and energy e, convolution with N(0,Q), Q>=qI, differs from its matching Gaussian by at most sqrt(2)L^2e/(3q). This is a full-law statement, with no second derivative of the source.

Each use below is conditional on the same previously exposed original caller and fixed source versions. M and the covariance reserve use independent COMPLETE private banks after Z is captured. No old force sample, private clock or coupling Gaussian is appended as an observer after it was integrated.

## 2. Full buffered comparison for one anchored source

Fix 0<h<=1/2. The analytical buffered second-substitution law is

    T_ref=h[Z-F_cont]+sqrt(1-h^2)N.                       (2)

Conditional on Z, H_cont and F_cont are Gaussian-functionals of the same true conditional Markov path. Let E_cont=F_cont-H_cont. Positive integral weights, Lip(g)<=A and the actual conditional Gaussian row norms give

    Lip_private H_cont<=A,
    Lip_private F_cont<=A(1+A),
    ||F_cont-E[F_cont|Z]||_(L2(private))<=C A sqrt(D),
    ||E_cont||_(L2(private)|Z)<=C A^2(|Z|+sqrt(D)).        (3)

These statements may be obtained on finite positive Gaussian path approximations and passed to their L2 limit. Their constants do not grow with the number of private Gaussian rows. The same positive-weight/row-norm argument as the finite reference-port proof controls the first; the moment bound uses one physical D-dimensional Gaussian energy per query and Minkowski, not the size of the white-noise record.

Gaussian Poincare gives covariance operator bounds C A^2 for both H_cont and F_cont. The one-marked-energy covariance comparison therefore yields

    ||Cov(F_cont|Z)-C(Z)||HS
       <=C A^3(|Z|+sqrt(D)).                             (4)

The cubic Gaussianization lemma, applied after scaling by h, compares (2), conditional on Z, with

    N(h[Z-m_cont(Z)],(1-h^2)I+h^2 Cov(F_cont|Z))

at error C h^3 A^3 sqrt(D)/(1-h^2). All higher private-path cumulants are handled by this estimate. Restoring (4) costs C h^2 A^3(|Z|+sqrt(D)) by the fixed-gap Gaussian-root inequality. Thus

    ||W2(Law(T_ref|Z),
         N(h[Z-m_cont(Z)],(1-h^2)I+h^2 C(Z)))||_(L2(Z))
       <=C h^2 A^3 sqrt(D).                             (5)

This is not obtained by claiming equal moments imply equal laws.

## 3. The finite positive program and variance budget

Set d=1-2h^2>=1/2 and eta=zeta=sqrt(d/2). Run the finite mean service M(Z). Independently run the positive square-clock aggregate D_Q(p) from the covariance service, keeping its SAME incoming p across all its complete action banks. Execute

    R=eta p+[h^2/(2eta)]D_Q(p)+zeta zeta_root,
    T=h Z-h M(Z)+R.                                    (6)

Here zeta_root is a fresh standard Gaussian independent of the whole action and mean tapes. The scalar zeta is its displayed coefficient. The actual action D_Q is the original-VALUE square program, not a matrix oracle.

Conditional on Z, R targets N(0,dI+h^2 C(Z)) with error

    Lambda[h^2 A^3+h^4 A^(7/2)+h^4 A^4]sqrt(D)
                     +C h^2 delta A^2 sqrt(D)+floors
       <=Lambda h^2 A^3 sqrt(D)+floors                 (7)

when delta<=A^2 and A<=1. The source radius and padding are still those of the original covariance service; h^2 is applied after the full action. Its positive quadratic covariance term is paid in (7), not erased.

The completed M target has unit covariance. Its independent sum with R therefore has exactly the reference covariance

    h^2I+dI+h^2 C(Z)=(1-h^2)I+h^2 C(Z).

The reference mean is h[Z-m_cont(Z)]. The actual output is not asserted to have these exact moments; it has the component law comparisons. Their conditional product coupling, then integration over the fresh standard Z, proves that (6) differs from the Gaussian reference in (5) by

    C h A^3 sqrt(D)+Lambda(h A^4+h^2 A^3)sqrt(D)+floors.

Finally (5) and the analytical reference W2(P2,mu_U)<=A^3 sqrt(D) give

    W2(Law(T),Law(h X+sqrt(1-h^2)N))
       <=Lambda h A^3 sqrt(D)+floors,
    X~mu_U, N independent.                              (8)

Every operation in (6) is a positive Gaussian-root pushforward. Negative deterministic coefficients such as -hM are allowed arithmetic; no signed distribution is sampled.

## 4. Exact conditional reverse-OU bridge

Return to an exposed caller z at 0<=r<t<=1. Let

    s^2=1-r^2, Delta=t^2-r^2,
    A_rt=r(1-t^2)/(t s^2), B_rt=Delta/(t s^2),
    v0=Delta/t^2, h=sqrt(Delta)/s.

The exact conditional endpoint X has law

    nu_(r,z)(dx) proportional to exp(-|x-rz|^2/(2s^2)-U(x))dx,

and the exact bridge is A_rt z+B_rt X+sqrt(v)N, v=(1-t^2)Delta/(t^2s^2).

Compute the finite conditional mode x=x_M and saved g(x) as in the existing conditional packet. The anchored normalized force is

    f(y)=s[g(x+s y)-g(x)], alpha=As^2.

Invoke the construction (6) for this original-VALUE source f, with alpha in place of A, and return

    K_new(z)=A_rt z+B_rt x+sqrt(v0)T.                  (9)

If Delta<=s^2/4, h<=1/2, so its reserve is uniformly positive. The exact affine geometry satisfies B_rt s=sqrt(v0)h and v=v0(1-h^2). Equation (8) therefore gives the full conditional bridge bound

    W2(Law(K_new(z)),K_exact(r,t;z))
       <=Lambda (Delta/t) A^3 s^5 sqrt(D)
                       +B_rt|R_M|+restored floors,     (10)

where R_M=x_M-rz+s^2g(x_M) is the actual finite-mode residual. The factor Lambda retains the actual fixed-order public-log source guards. All intrinsic estimates above are uniform over the original captured caller after anchoring; their L2(Z) quadrature estimates apply to the NEW private standard carrier Z of this call, not to the external z distribution.

The zero-buffer case t=1 is not covered by (9). It will be handled by a separately priced terminal call below.

## 5. Exact-kernel contraction and finite endpoint schedule

The exact conditional nu_(r,z) satisfies

    W2(nu_(r,z),nu_(r,z'))<=r|z-z'|.                   (11)

For example, use the synchronous strong-convexity comparison of potentials with quadratic center a=rz; the force shift is (a-a')/s^2 and the strong-convexity constant is at least 1/s^2. The same-noise affine bridge then has Lipschitz constant

    A_rt+B_rt r=r/t.                                  (12)

This is the **exact reference kernel** contraction. It is not substituted for an actual finite program's same-seed first.

Choose rho=3/4 and r_j^2=1-rho^j, j>=0. At stage j, s_(j-1)^2=rho^(j-1), Delta_j=s_(j-1)^2/4, and h=1/2. Run (9) along these stages. Let J be the first integer with

    s_J^5<=A, equivalently rho^J<=A^(2/5).              (13)

At the final stage r_J->1, use the already audited original conditional positive Stein packet directly, with no small-buffer Gaussianization. Its physical conditional error is

    C A^2 s_J^5 sqrt(D)+mode/numerical floors
       <=C A^3 sqrt(D)+floors.                         (14)

If E_j is the W2 error after stage j, the kernel comparison gives E_j<=(r_(j-1)/r_j)E_(j-1)+epsilon_j. The terminal contraction is r_J. Therefore the intrinsic errors obey

    E_final<=sum_(j=1)^J r_j epsilon_j+epsilon_terminal
       <=Lambda A^3 sqrt(D) sum_j Delta_j s_(j-1)^5
                         +C A^3 sqrt(D)
       <=C Lambda A^3 sqrt(D).                        (15)

The fresh private Z in every local mean/covariance service is standard by construction, even when its external entering caller has a non-Gaussian approximate law. This is why the L2 standard-carrier interfaces compose here.

## 6. Finite modes, floors, and original counts

The conditional-mode fixed-point iteration contracts by alpha=As^2. Taking a fixed count M>=3 gives

    |R_M|<=s^2 alpha^M |g(rz)|
             <=A^(M+1)s^(2M+2)r|z|.

After (10) and the contraction weights, its total expected contribution is bounded by C A^4 times the entering caller moment for M=3. The exact OU marginals have second moment <=D. The error recursion, with these small linear caller terms retained, gives a discrete Gronwall bound ||z_j||2<=C sqrt(D) and a total mode floor O(A^4 sqrt(D)). One may instead retain the displayed mode floors symbolically; no uniform pointwise bound for an unbounded external caller is claimed.

All numerical/known-coefficient/filter floors are assigned at their actual source readouts and summed with the same reference-kernel weights. Preserve the global finite-mode linear-tilt restoration from the original physical posterior construction. These are absolute floors with their original caller profiles, never divided by a realized source energy. They may be selected below the displayed A^3 law allowance using the existing finite precision model; no exact potential value or exact mode is used.

The schedule has J=O(log(1/A)). Before the terminal call,

    rho A^(2/5)<s_J^2<=A^(2/5),

so every invoked local alpha=As^2 has log(1/alpha)=O(log(1/A)). There is no hidden exponentially small endpoint clock. All fixed-order mean, square-action, and dyadic quadrature counts therefore remain public-log polynomials uniformly through the schedule.

For one mean call, with actual complete source occurrence counts N_B,N_E,

    Q_mean<=Q_captured+n_out N_B+n_out(n_in+2)N_E
                         +Q_known/numerical.

For the covariance action, use its independently checked complete count; its basic original g count is N_r N_q(2K_f+5), before separately requested restoration/replay. Each stage also has its actual conditional mode/anchor and known Gaussian arithmetic. The terminal packet costs M+1+n_alpha original VALUES, plus any separately requested canonical terminal force. Summing over J preserves zero inverse-A query exponent at this fixed grade.

No random source ancestor is shared across independent complete mean/covariance banks or across changed callers. Cache only identical captured caller-only modes, anchors, and origins under the same finite version. The same p is deliberately shared within one square-clock covariance aggregate; its owned X and response banks retain the independence required by that action. Complete private Gaussian dimensions include every mean/filter/action root and stage carrier, not just the schematic Z,p,zeta_root.

## 7. Actual first/caller/zero scope

All producers in (9) use original g VALUES. Requested first/adjoint sweeps use original HVPs only at those recorded VALUE sites, through the complete nested ancestors and conditional modes. Discarded primals pay their full replay. No saved HVP is differentiated.

The raw nested mean source has actual first O(alpha), small curl O(alpha^2), and captured-carrier first O(alpha), with its recorded nonzero fixed-carrier origin. Its completed mean carries its admitted known unit Gaussian carrier and finite residual/caller firsts. The square-clock reserve has residual/private/captured-carrier first Lambda h^2 alpha^(3/2), besides its fixed Gaussian carrier, and its original physical-caller profiles follow its checked parameter port. These are actual source bounds, independent of the law couplings.

For an exterior normalized conditional center a, f's actual caller first is O(As), with the finite-mode derivative retained. The final s,B_rt,sqrt(v0) readouts are applied before composing caller/first bounds. The complete finite endpoint has the literal finite chain-rule/JVP/VJP program, with its known source-dependent first majorants and public-log stage count. The exact contraction (12) is used only to propagate W2 errors, never as a claim about those actual derivatives. No new retained/protected/proxy or strong paired-source theorem is asserted by this endpoint law result.

More explicitly, the actual incoming-state first at a positive r obeys the safe bound

    (r/t)[1+C Lambda A Delta
                    +C Lambda A^(3/2)Delta^(3/2)].     (16)

The finite mode has first at most r/(1-alpha), contributing relative O(A Delta) after its B_rt readout. The mean's captured-source first is O(Lambda Asr), giving relative O(Lambda A Delta). The square action's captured-source first is O(Lambda sqrt(alpha)Asr); after its h^2 sqrt(v0) readout this gives the final term of (16). At r=0 the incoming caller is unread and that first is zero. Over the schedule, sum Delta<=1 and sum Delta^(3/2)<=1, so the excess product in (16) is bounded by exp(C Lambda A+C Lambda A^(3/2)), at most a fixed constant under the actual small-A guards. This is an actual finite-graph bound, separate from (12).

The source-zero Gaussian stage rows have incoming coefficient r/t and fresh variance v0=1-(r/t)^2. Their complete endpoint composition therefore has one known unit coisometric carrier. The source-dependent private residual firsts are O(Lambda alpha) in each normalized buffered stage. Their final weighted sum is bounded by C Lambda A, because along this schedule sum_j sqrt(Delta_j) s_(j-1)^2 is a fixed geometric sum. In addition to these direct fresh-private terms, the source-dependent incoming-state Jacobian corrections in (16) contribute O(Lambda A sum_j Delta_j)=O(Lambda A); the same bounded product/discrete Gronwall argument closes their propagation through earlier roots. The terminal conditional packet has the corresponding original O(alpha) residual first. Thus the complete standardized endpoint has private first bounded by 1+C Lambda A relative to its literal known carrier, without pretending that a reference coupling noise is one of those roots. Under original physical normalization this supplies the usual sqrt(A) private scale; original physical caller profiles follow their recorded chain rule and finite-mode restoration.

At all owned zero roots, each raw anchored source follows its exact saved-anchor zero trajectory. Fixed nonzero carrier origins are computed and restored explicitly inside the mean services. The source-zero carriers in every completed mean and reserve remain their actual Gaussian rows. At total input/private zero the anchored global state is zero (or the saved physical mode after restoring normalization), not an exact uncomputed mode. Every schedule, selector, quadrature, padding, variance share and numerical version is frozen before differentiation.

## 8. Scope of the gain

The construction uses the exact continuous current only as a comparison reference, realizes its nonlinear mean by the finite C2 decoupling source, and realizes its conditional covariance by a positive square-clock action. The cubic Gaussianization theorem prices all remaining higher private cumulants. It does not resume an ordinary unsmoothed path grid and does not inherit that grid's diagonal clock error.

The bounded conclusion is a first full nonlinear order gain: standardized Lambda A^3 sqrt(D), physical Lambda A^(7/2) sqrt(D), with polynomial-logarithmic original VALUE work at the declared fixed grade. Public logs are explicit; without further padding optimization this is not asserted as a log-free C A^3 bound. The higher-order recurrence, and any eventual-sublinear c(P)/P claim, remain open.
