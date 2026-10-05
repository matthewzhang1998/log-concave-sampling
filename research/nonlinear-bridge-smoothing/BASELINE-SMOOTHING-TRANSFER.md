# Baseline-preserving Gaussian smoothing has an extra factor of A

2026-10-05. Sections 1–6 prove baseline-preserving target transfer under the original Hessian interval, with a genuine OU history throughout. Sections 7–8 strengthen this to weak stability of the actual finite residual stencil: analytical mollification then gives grade 8/3 for the ORIGINAL unsmoothed stage-two graph, with no new producer leaf or VALUE occurrence.

## 1. Result

Let g,h:R^D→R^D be anchored C1 maps satisfying

    g(0)=h(0)=0,   0≤Dg,Dh≤A I,   0<A≤1/2,
    δ=sup_x |g(x)−h(x)|<∞.

For each source q, use the stationary OU history X_s, started at X_0=x, and define

    F1,q,a = ∫_0^∞ e^(−r) q(X_(a+r)) dr,
    F2,q   = ∫_0^∞ e^(−a) q(X_a−F1,q,a) da,
    ψ2,q(x)=E[q(x−F2,q) | X_0=x].

The outer resolvent convention is

    R1 f(z)=∫_0^1 E_G f(tz+sqrt(1−t²)G) dt,
    m3(q)=R1 ψ2,q,
    R_q=m3(q)−R1 q.

Put

    L=A/2+A²/4,
    M(A)=(1+A)π/2 + L/[A sqrt(1−L)],
    C_D(A)=(1+A)+M(A)sqrt(D).

Then

    ||R_g−R_h||_(L2 γ_D) ≤ C_D(A) A δ.                 (1)

In particular, uniformly for A≤1/2,

    C_D(A) ≤ 3/2 + [3π/4+5/(2sqrt(11))]sqrt(D)
            < 5sqrt(D),
    ||R_g−R_h||_2 ≤ 5A sqrt(D) ||g−h||_∞.             (2)

This uses no modulus of continuity of Dg or Dh, no bound on D²g, no commuting-Hessian reduction, and no replacement of the genuine nonlinear history by its conditional mean. Symmetry and positivity are available from the source hypothesis but the transfer proof itself only needs the stated A-Lipschitz and anchoring properties.

For the anchored Gaussian smoothing

    h_ε(x)=E_Y[g(x+εY)−g(εY)],   Y∼N(0,I_D),

we have 0≤Dh_ε≤A I, h_ε(0)=0 and

    ||g−h_ε||_∞ ≤ 2A ε E|Y| ≤ 2A ε sqrt(D).

Consequently

    ||R_g−R_hε||_2 ≤ 2 C_D(A) A² ε sqrt(D)
                    ≤ 10 A² ε D.                    (3)

The ambient-dimension loss in (3) is explicit. This is not a dimension-uniform O(A² ε sqrt(D)) assertion for arbitrary high-dimensional sources.

If g acts in a fixed orthogonal decomposition into blocks of dimensions d_j≤b, apply the argument separately to every block. Isotropic Gaussian smoothing preserves the decomposition, even if its basis is not supplied to an executor. Then

    ||R_g−R_hε||_2
      ≤ 2 C_b(A) A² ε sqrt(D)
      ≤ 10 A² ε sqrt(bD).                            (4)

Thus bounded blocks yield a uniform-in-D extra-A smoothing gain with the natural sqrt(D) scaling. For a general source, taking the single block b=D gives the valid fixed-D result (3).

## 2. A Gaussian near-identity transport lemma using only first derivatives

Let Y∼N(0,I_d), let B:R^d→R^d be C1 with ||DB||≤L<1, and put T(y)=y−B(y). The map T is a global orientation-preserving C1 diffeomorphism: the inverse equation y=z+B(y) is a contraction, and every I−θDB is invertible for 0≤θ≤1.

Let ν=T_#γ_d. Change of variables gives its relative entropy

    H(ν|γ_d)
      = E[−Y·B(Y)+|B(Y)|²/2−log det(I−DB(Y))].

Stein's identity is valid because B has at most linear growth and bounded first derivative:

    E[Y·B(Y)]=E tr DB(Y).

For every real matrix K of operator norm at most L<1, no symmetry being needed,

    −tr K−log det(I−K)
       = ∑_(k≥2) tr(K^k)/k
       ≤ d ∑_(k≥2) L^k/k
       ≤ d L²/[2(1−L)].

The log determinant is the real branch obtained continuously from I; det(I−K)>0. Therefore

    2H(ν|γ_d) ≤ E|B(Y)|² + d L²/(1−L).

Pinsker's inequality, in its L1-density convention, implies

    ||ν−γ_d||_L1 ≤ sqrt(E|B(Y)|²+d L²/(1−L)).         (5)

For any vector-valued bounded f with sup|f|≤δ, this gives

    |E f(T(Y))−E f(Y)|
       ≤ δ sqrt(E|B(Y)|²+d L²/(1−L)).                (6)

The important cancellation is between Stein's divergence and the linear log-Jacobian term. It avoids differentiating the transport Jacobian, so no uncontrolled D²B is introduced. A direct integration by parts after multiplying by (I−DB)^−1 would generally introduce precisely such an uncontrolled derivative and is not the proof used here.

The same statement applies to an input N(m,c²I), c>0, after taking

    B(y)=F(m+cy)/c.

Then ||DB||=||DF||≤L and the displacement term in (6) is E|F(m+cY)|²/c².

## 3. Genuine-history displacement and stability estimates

Realize the conditional OU path as

    X_s=e^(−s)x+ξ_s,

where the entire Gaussian process ξ is independent of x. All nested F1 and F2 values use this same actual future path, including its true cross-time covariance. For almost every ξ, differentiation under the absolutely convergent integrals is justified by bounded Dg:

    ||D_x F1,g,a|| ≤ A e^(−a)/2,
    ||D_x F2,g||
      ≤ ∫ e^(−a) A[e^(−a)+A e^(−a)/2] da
      = A/2+A²/4 = L <1.                             (7)

For stationary X, Minkowski and anchoring yield

    ||F1,g,a||_2 ≤ A sqrt(D),
    ||F2,g||_2 ≤ A(1+A)sqrt(D).                      (8)

The corresponding pathwise cross-source estimates are

    |F1,g,a−F1,h,a|≤δ,
    |F2,g−F2,h|≤(1+A)δ.                             (9)

Equation (9) preserves the nested source dependence: compare g(X_a−F1,g,a) with h(X_a−F1,h,a), using the source-value difference δ and the A-Lipschitz shift difference separately.

## 4. The exact nonlinear decomposition

Set f=g−h. Before taking any conditional expectations, the difference of the inner nonlinear remainders is exactly

    g(x−F2,g)−g(x)−h(x−F2,h)+h(x)
      = f(x−F2,g)−f(x)
        +h(x−F2,g)−h(x−F2,h).                       (10)

The last term has pointwise norm at most A(1+A)δ by (9). This is already the desired extra A.

For the first term fix 0≤t<1, z, and the complete future noise ξ. Put c=sqrt(1−t²), x=tz+cG. The map x↦x−F2,g(x,ξ) has derivative defect at most L by (7), so (6) gives

    |E_G[f(x−F2,g)−f(x)]|
       ≤ δ sqrt(E_G|F2,g(x,ξ)|²/c²+D L²/(1−L)).

Average ξ and take L2 in z∼γ_D. Conditional Jensen and (8), with x again standard Gaussian, give

    ||E_(G,ξ)[f(x−F2,g)−f(x)]||_(L2_z)
       ≤ δ sqrt(D)
              sqrt(A²(1+A)²/c²+L²/(1−L)).           (11)

The outer endpoint t→1 has not been discarded. Minkowski and the integrable singularity

    ∫_0^1 dt/sqrt(1−t²)=π/2

imply from (11)

    ||R1 E_ξ[f(x−F2,g)−f(x)]||_2
       ≤ Aδ sqrt(D)[(1+A)π/2+L/(A sqrt(1−L))].       (12)

The value at t=1 is immaterial to the Lebesgue integral. No false uniform-in-t estimate is asserted. Adding the last term of (10) proves (1).

## 5. What the transfer permits for stage two

The baseline-preserving analytical target is

    T_ε = R1 g + [m3(h_ε)−R1 h_ε].

Its error to the original target is exactly R_g−R_hε and hence satisfies (3) or (4). By contrast, replacing the full target m3(g) by m3(h_ε) also pays the unsuppressed baseline difference R1(g−h_ε), of size O(Aε sqrt(D)). Retaining R1g is therefore essential to this argument.

Gaussian smoothing supplies derivative bounds with B_ε=O(ε^−1) and C_ε=O(ε^−2) for the stage-two theorem. For example, subtract A I/2 from the Hessian inside the centered Gaussian derivative kernels to obtain the valid bounds

    Lip(Dh_ε) ≤ A/[sqrt(2π) ε],
    Lip(D²h_ε) ≤ A/[sqrt(2) ε²].                     (13)

For the second inequality, if u,v are unit vectors, the Gaussian kernel

    (u·Y)(v·Y)−u·v

has variance 1+(u·v)²≤2. These estimates are in the required induced Euclidean multilinear norms and need no source derivative beyond bounded Dg.

On blocks of size at most b, the stage-two shifted-bias theorem therefore contributes

    O_b(A⁴[ε^−1+ε^−2]sqrt(D)).

Combining this with (4), for 0<ε≤1, gives the analytical target-bias bound

    O_b([A²ε+A⁴ε^−2]sqrt(D)).                       (14)

Choosing ε=A^(2/3) gives O_b(A^(8/3)sqrt(D)), a genuine quantitative improvement over A² under only the original Hessian interval within the bounded-block class. For arbitrary fixed D, the same exponent holds with an explicitly dimension-dependent constant by taking b=D. The lower ε^−1 term is A^(10/3); an O(A⁴) own-mean completion allowance is smaller still, subject to all original guards.

Equation (14) is a target-transfer and stage-two-bias conclusion, not a completed original-VALUE implementation theorem. An exact h_ε VALUE is an expectation, and cannot simply be installed as a producer leaf. Any finite realization or approximation of smoothing must separately preserve or re-establish the source ports, anchor/caller dependencies, positive own-mean completion, numerical floors, and full occurrence/replay charges. In particular, no claim of zero inverse-A original-VALUE exponent follows from (14) alone. This proof removes the mathematical target-transfer obstruction; it does not remove that implementation obligation.

## 6. Checks and scope boundaries

- For an anisotropic linear source g(x)=Kx, h_ε=g exactly. The transfer error is identically zero, as it should be; the upper bounds need not be sharp there.
- The genuine path is used only in analysis. There is no continuous-history producer, strong-mean target substitute, source-dependent covariance reconstruction, or differentiation of an HVP.
- The argument explicitly handles the outer R1 region arbitrarily close to t=1 through an integrable 1/sqrt(1−t²) bound.
- No source-independent high-dimensional bound better than the displayed dimension factors is proved here.
- No sealed stage-one or stage-two file was edited. The stage-two theorem and its shifted derivation are inputs; this note proves the new baseline-preserving transfer lemma.

## 7. Stronger extension: weak stability of the actual stage-two residual stencil

The implementation obligation at the end of Section 5 can be avoided entirely: smoothing can be used only in the proof while executing the original unsmoothed stage-two graph. The missing ingredient is the following weak stability statement for that graph, which we now prove.

Use exactly the finite inner stage-two rules, with positive normalized weights and their actual one-time/pair Gaussian rows, fixed identically for both sources. Write E_q(x,W) for its raw residual F_stage2,q(x,W)−q(x), and

    S_q^out=R1 E_W E_q.

This notation denotes the exact outer mean of the finite inner graph, not the completed output and not the finite outer quadrature. For q=g,h as in Section 1 and 0<A≤1/36, put

    k=sqrt(3/8), s0=1/sqrt(2)+Ak,
    d0=1+1/sqrt(2), e0=1+k,
    L*=7A/2+A²/4,
    HΣ=(7/2)s0+2d0+A e0
       =2+11/(2sqrt(2))+A[1+(9/2)k],
    C_st,d(A)=14+11A/2
               +sqrt(d)[(π/2)HΣ
                    +(7/2)(L*/A)/sqrt(1−L*)].

Then, in dimension d,

    ||S_g^out−S_h^out||_2 ≤ C_st,d(A) Aδ
                            ≤ 37A sqrt(d) δ.        (15)

The constants and proof are uniform in all finite positive node counts. No derivative regularity beyond the Hessian interval is required.

### 7.1 Every terminal argument is a near-identity transport

Condition on all private inner roots W=(N,M,L1,L2). The actual finite row arguments U,V are affine in x, their x coefficients have absolute value at most one, and each has standard Gaussian marginal when x and the roots are standard Gaussian. Set

    b=q(q(w)), S=x−q(v)+b,
    d(U)=q(v)−q(U),
    e(U,V)=q(U)−q(U−q(V))−b.

Since v_x=I/2 and w_x=I/4,

    ||DS−I||≤A/2+A²/4,
    ||Dd(U)||≤3A/2.

The exact ordered derivative of e is

    De=[Dq(U)−Dq(U−q(V))]U_x
          +Dq(U−q(V))Dq(V)V_x
          −Dq(q(w))Dq(w)/4.

The two Hessians in the bracket both lie in [0,A I], so their difference has norm at most A. No factors are commuted. Therefore

    ||De||≤A+5A²/4.

The terminal arguments are S, S±d, S±e, and S±d(U)±d(V). Their derivative defects from I are bounded, respectively, by

    A/2+A²/4,
    2A+A²/4,
    3A/2+3A²/2,
    7A/2+A²/4.

Thus all are bounded by L* when A≤1/36 (indeed this comparison holds for A≤1/2), and L*<1. The transport lemma applies to every terminal argument separately.

The corresponding integrated L2 displacement bounds are

    ||S−x||_2 ≤A s0 sqrt(d),
    ||d(U)||_2 ≤A d0 sqrt(d),
    ||e(U,V)||_2 ≤A² e0 sqrt(d).

The four terminal families have total absolute coefficient masses 1,1,1,1/2. Consequently their absolute-weighted displacement constants sum to

    s0 +(s0+d0)+(s0+A e0)+(s0+2d0)/2 = HΣ.

Only the actual nodewise marginals enter these estimates. Shared quadrature roots are not treated as a genuine multi-time history.

### 7.2 Baseline cancellation supplies the weak gain

Expand the finite stencil into terminal evaluations. Excluding the residual baseline −q(x), its signed terminal coefficients c_i have

    ∑c_i=1,   ∑|c_i|=7/2.

For f=g−h, and T_i,q denoting each source-dependent terminal argument, the exact residual difference is

    E_g−E_h
      =∑c_i[f(T_i,g)−f(x)]
         +∑c_i[h(T_i,g)−h(T_i,h)].                 (16)

There is no uncancelled f(x) term. The first sum is bounded, after the exact outer R1, using the transport lemma separately for each terminal map. Their displacement sum is HΣ and each derivative defect is at most L*. Hence its L2 norm is at most

    Aδ sqrt(d)[(π/2)HΣ
                       +(7/2)(L*/A)/sqrt(1−L*)].   (17)

The same integrable outer singularity as in Section 4 is retained; no estimate is asserted at a fixed outer node t=1.

For the second sum, pathwise source replacement gives

    |S_g−S_h|≤(2+A)δ,
    |d_g−d_h|≤2δ,
    |e_g−e_h|≤(3+2A)δ.

Apply the A-Lipschitz bound of the terminal h and sum actual absolute coefficients. The result is

    Aδ[(7/2)(2+A)+2+(3+2A)+2]
      =Aδ[14+11A/2].                              (18)

Equations (17) and (18) prove the first bound in (15). All coefficients are increasing in A. At A=1/36, the dimension-one constant is less than 36.487; since d≥1, C_st,d(A)≤37sqrt(d). This proves the second bound.

## 8. Corollary: the original unsmoothed graph has quantitative grade 8/3

Assume only the original anchored Hessian interval and block separability into dimensions d_j≤b. No uniform modulus of the block Hessians is assumed. Execute the ORIGINAL stage-two graph on g itself, without adding any smoothing node, convolution oracle, grid, randomization, derivative query, or VALUE occurrence. Retain the ordinary source-independent positive rules already required by that packet.

For analysis alone, introduce h_ε. The exact telescoping identity is

    S_g^out−R_g
      =(S_g^out−S_hε^out)
        +(S_hε^out−R_hε)
        +(R_hε−R_g).                               (19)

The first and last terms are controlled by (15) and (2). The middle term is exactly the original stage-two shifted-bias estimate applied to the smooth analytical source h_ε. Combining blockwise and using δ_j≤2Aεsqrt(d_j) gives

    ||S_g^out−R_g||_2
      ≤ [84sqrt(b) A²ε
            +K_(b,Bε,Cε,A) A⁴
            +δ1 A²+sqrt(11/8)δL A³+(δQ/2)A²]
           sqrt(D),                               (20)

where valid declared smoothing constants from (13) are

    Bε=1/[sqrt(2π)ε],   Cε=1/[sqrt(2)ε²].

The 84sqrt(b) constant is the conservative sum 2(37+5)sqrt(b). The function K is exactly the stage-two packet's shifted-bias constant; all its original coefficients remain intact.

Now take ε=A^(2/3). The terms in (20) arising from the Bε and Cε parts of K are respectively O_b(A^(10/3)) and O_b(A^(8/3)). The weak transfer term is also O_b(A^(8/3)). Thus the exact-outer mean of the original finite inner graph has

    ||R1 E_W F_stage2,g−m3(g)||_2
       ≤ O_b(A^(8/3)sqrt(D))
            + the displayed finite-rule errors.    (21)

This is an outer-mean theorem. It does not claim the stronger inner L2 stage-two bias under unrestricted Hessian regularity.

Restore the actual finite outer rule. Its universal L2 operator error applies to the actual unsmoothed mean, whose norm is bounded by A C_F2 sqrt(D) using the original source-energy proof. Therefore the additional error is exactly the already-budgeted

    δ_out A C_F2 sqrt(D).

For example, the original choices δ1≤A², δL≤A, δQ≤A², δ_out≤A³ make every rule error O(A⁴sqrt(D)), smaller than A^(8/3). All finite node counts and all original VALUE occurrence counts are unchanged.

Finally, the actual executed source is still the original g graph, whose normalized first/curl/caller/energy ports never required B or C. Its same guarded positive own-mean completion therefore adds its existing O(Λ_comp A⁴sqrt(D)) allowance and restored absolute floors, with all native guards, replay/capture obligations, independent completed banks and positivity requirements unchanged. Subject to those same guards, this gives a completed canonical-m3 mean-law component with grade 8/3 on the unrestricted-Hessian-regularity bounded-block class.

For an arbitrary fixed dimension D, take b=D: the exponent 8/3 remains valid, with the explicit dimension losses in (20) and K. This does not prove dimension-uniform high-dimensional grade 8/3 without a structural qualification. It also does not prove an induction to arbitrary grades or the full fourth-order endpoint join.

The practical conclusion differs from the hypothetical construction in Section 5: no smoothing implementation is needed. The virtual smoothed source is eliminated by (19), and the ORIGINAL finite graph gains a stronger theorem. The zero inverse-A original-VALUE exponent is retained at this fixed grade because the executed graph, rules, and completion are unchanged; no hidden Gaussian-convolution calls remain.

### 8.1 Optimizing the displayed dimension dependence

It is not necessary to restrict ε≤1 in (20) or in the derivative estimates (13). For b≥1, choose instead

    ε=A^(2/3)b^(1/6).

Use sqrt(b+2)≤sqrt(3b) and sqrt((b+2)(b+4))≤sqrt(15)b. Define

    P(A)=[d0 e0+(A/2)e0²+1/2]sqrt(3)/sqrt(2π),
    Q(A)=[(d0+A e0)³+4d0³+A³e0³]sqrt(15)/(6sqrt(2)).

The target-bias part of (20), excluding rule errors, is then bounded by

    [84+P(A)A^(2/3)b^(−1/3)+Q(A)]
       A^(8/3)b^(2/3)sqrt(D)
      ≤100 A^(8/3)b^(2/3)sqrt(D).                   (22)

For the last convenient constant, A≤1/36 implies P(A)<2.28 and Q(A)<11.54, so even replacing A^(2/3)b^(−1/3) by one gives a total below 98. The stated 100 leaves slack. Alternatively retain the exact displayed P,Q expression without decimal estimates.

For arbitrary dimension with b=D, (22) is

    100 A^(8/3)D^(7/6).

Relative to an O(A²sqrt(D)) benchmark, its power improvement has factor (AD)^(2/3), apart from constants. Thus the general-D certificate is useful in an AD-small regime; it is not a dimension-free improvement. The bounded-block result has the same comparison factor (Ab)^(2/3), uniformly in ambient D.

The bound C_st,d≤37sqrt(d) used above can be checked without floating-point arithmetic. For A≤1/36, use π<22/7, sqrt(2)>7/5, k<5/8, L*<1/10, and sqrt(1−L*)>18/19. Then

    C_st,d(A)/sqrt(d)
      ≤14+11/72
          +(11/7)(83/14+61/576)
          +(7/2)(505/144)(19/18)
      <37.
