# What the shrinking-buffer mean join actually exports

2026-10-05. A finite-depth reentry theorem and a certified-grade ceiling for the unmodified LAW-only join. This is not a new accuracy exponent, an endpoint sampler, or an impossibility theorem for improved constructions.

## 1. Named inputs and outputs

Let X_t be the stationary Gaussian OU history and define analytical force substitutions

    F_0,t = 0,
    F_k,t = integral_0^infinity exp(-s) g(X_(t+s)-F_(k-1,t+s)) ds,
    m_k(x) = E[F_k,0 | X_0=x].

The original source satisfies g=grad U, g(0)=0, and 0<=Dg<=A I. All executed leaves remain original-g VALUES. Histories below exist only in analysis and are justified by finite Gaussian-bank approximation before passing to the limit.

The actual reusable input at depth k is a FINITE RAW GRAPH S_k(y,Omega), together with:

- Its complete Gaussian bank and exact source-zero and caller-origin records.
- Its own-mean target certificate ||E S_k-m_k||_(L2 gamma)<=e_k.
- A common-known-carrier baseline B_k, genuine gradient in one D-dimensional private root, with first <=A.
- E_k=S_k-B_k with full private AND retained-caller first <=Lambda_k A, square-lift curl <=Lambda_k A^2, and centered Lp energy <=Lambda_(k,p) A^2(|y|+sqrt(D)).
- Its expanded original-VALUE count Q_k, baseline count Q_B,k, private dimension d_k, all ancestor/first/adjoint replay paths, and absolute floors.

Lambda_k denotes declared public-log factors. Native compiler guards must hold at these actual constants and dimensions. The source need not be a convex gradient as a whole; its literal split is the input to raw:split.

The optional completed unit-variance LAW N(m_k(y),I) is NOT the reusable source above. For a new output buffer u, rerun the compiler on S_k/sqrt(u), with its complete graph and tapes. Scaling an old completed unit LAW would scale its target mean, and subtracting its known carrier would not preserve its LAW certificate.

At k=2 the existing F_Q supplies this input with numerical O(A) private first and e_2<=C A^3 sqrt(D). At k=3 the pinned shrinking-buffer packet supplies it with

    e_3 <= C A^(16/5)sqrt(D)+Lambda_3 A^(19/5)sqrt(D)+floors.

## 2. Genuine depth propagation and covariance debt

Coherent stationary contraction gives, uniformly in finite k,

    ||F_k||_2 <= A/(1-A) sqrt(D),
    ||F_k-F_(k-1)||_2 <= A^k sqrt(D),
    ||F_k-F_2||_2 <= A^3/(1-A) sqrt(D), k>=2.

The complete Gaussian-bank first obeys

    L_k <= A(1+L_(k-1)) <= A/(1-A).

It does NOT shrink when e_k improves. Also

    ||m_k-m_2||_2 <= A^3/(1-A) sqrt(D).

The existing Cov(F2|y) service remains usable at every such depth. The required conditional covariance restoration is

    ||Cov(F_k|y)-Cov(F_2|y)||_(L2_y;HS) <= C A^4 sqrt(D).   (1)

Proof: couple the true future histories coherently at retained y. Put Delta=F_k-F_2. The Gaussian covariance identity, with the Riesz derivative applied to Delta, gives

    ||Cov(Delta,F_j | y)||HS <= L_j ||Delta-E[Delta|y]||_2.

The covariance difference is Cov(Delta,F_k)+Cov(F_2,Delta). Integrate in y and use the previous coherent energy and first bounds. This proof requires no small first for Delta and introduces no extra dimension factor. The conditional future source banks in the executed services remain fresh given y; the coherent coupling is used only in the target comparison.

Thus higher force depth does not currently create a new order-three covariance target: the first unpaid change to the existing F2 covariance is order A^4. Dividing by the actual positive Gaussian gap and applying the outer A-Lipschitz consumer costs C A^5/sqrt(u) sqrt(D). This is already dominated by the current variable mean/Gram bill A^5/u^(3/2).

### A useful adjacent third-cumulant restoration lemma

There is also a dimension-free stability bound for the third centered moment, which can be used by a genuine skew producer. Let a,b be centered maps of a complete Gaussian bank with operator firsts L_a,L_b, and let Delta be any centered vector source with energy e_Delta. For a fixed matrix M, Gaussian Poincare and Cov(a)<=L_a^2 I, Cov(b)<=L_b^2 I give

    Var(a^T M b)
      <= E|Da^T M b+Db^T M^T a|^2
      <=4 L_a^2 L_b^2 ||M||HS^2.

Thus the operator f -> E[f a b^T], restricted to centered scalar L2 inputs, has norm <=2 L_a L_b. Apply it separately to the components of Delta and sum squares:

    ||E[Delta tensor a tensor b]||HS <=2 L_a L_b e_Delta.   (1a)

No first of Delta is needed. Telescoping the three slots of the centered third moment and applying (1a) proves

    ||kappa3(F_k|y)-kappa3(F_2|y)||_(L2_y;HS)
       <= C A^2 ||F_k-F_2||_2 <= C A^5 sqrt(D).

Similarly the old first force H=F1 satisfies

    ||kappa3(F_2|y)-kappa3(H|y)||_(L2_y;HS)
       <= C A^4 sqrt(D).

These are target-restoration statements for actual coherently coupled forces. They do not execute a third tensor or establish a skew LAW producer. They resolve a possible extra depth debt if such a producer is independently admitted.

## 3. An actual finite reentry program

Take w=A^(4/5), q=1-w. Use exactly the pinned disintegration and positive outer quadrature, with bulk t<=q and near endpoint t>q. On the bulk, actual conditional variance v is in [w,2w].

At a bulk node:

1. Rerun the mean compiler on -q S_k(y,Omega) at baseline variance v/2, integrating its COMPLETE fresh private bank.
2. Run the existing whole Gram LAW at variance v/(4q^2), then multiply its complete output by q.
3. Run the rebalanced mixed-K LAW at variance v/(4q^2), then multiply its complete output by q.
4. Add those three complete returns to the known conditional center a(z,y), and apply F(u)=g(u-wg(u)).

The conditional Gaussian references add to the right mean -q m_k(y), with mean discrepancy at most q e_k, and covariance vI+q^2 Cov(F_k|y), up to the explicitly paid (1). No covariance matrix, force mean, or analytical Gaussian root is an executing leaf.

At EVERY depth, retain the original F_Q as the fresh RAW near-endpoint branch. It has numerical O(A) complete-bank first. Its mean differs from m_k by C A^3 sqrt(D), using m_k-m_2 above. Thus it meets the old endpoint estimate without inheriting the new Lambda_k private-first factor. Replacing this branch by S_k is permitted but inferior: its fresh-RAW comparison would instead pay Lambda_k^2 A^3 sqrt(w).

The true short-prefix estimate is uniform in k: replacing its integrand by g(X_0) costs C[A w^(3/2)+A^2 w]sqrt(D) in the force. After the terminal g and coherent commutation the cost is C[A^2 w^(3/2)+A^3 w]sqrt(D). The mean input e_k changes none of this.

Align the known outer Gaussian carrier across positive outer nodes by the same scalar-block rotations as in the pinned packet. The resulting S_(k+1) again has a genuine-gradient B_(k+1), near-gradient E_(k+1), and the declared actual first, curl and energy ports. This is a direct graph statement, not a derivative of a LAW estimate.

## 4. Error recurrence and the unchanged certified ceiling

For w in the supported source-closure range A<=w<=1/2 and sufficiently small actual native radii, the construction above yields

    e_(k+1) <= A e_k
      + C sqrt(D)[A^2 w^(3/2)+A^3 w+A^3 sqrt(w)+A^4/w+A^4 w+A^4]
      + Lambda_(k+1) sqrt(D)[A^5/w^(3/2)+A^6/w^(7/6)]
      + absolute floors.                                  (2)

The A^4 w is the old F_Q mean mismatch on the endpoint interval; leaving it as A^4 would also be safe. The A^4 term includes finite target quadrature/restoration. All constants in the numerical C terms are uniform in finite k when A<=1/2. All new compiler/replay factors are retained in Lambda_(k+1).

Set w=A^alpha. The existing forcing exponents are

    2+3alpha/2, 3+alpha, 3+alpha/2, 4-alpha,
    4+alpha, 4, 5-3alpha/2, 6-7alpha/6.

On 0<alpha<=1 the maximum of their minimum is 16/5, uniquely achieved at alpha=4/5 by balancing the short-prefix and covariance-Gaussianization terms. In particular the next input e_k is multiplied by A, but the forcing term is unchanged. At fixed depth the recurrence has the schematic form

    e_(k+1) <= A e_k+B_(k+1),
    e_k <= A^(k-2)e_2+sum_(j=3)^k A^(k-j) B_j.

Each B_j has the same certified leading grade 16/5. Calling 16/5 a ceiling means the BEST GRADE CERTIFIED BY THIS LEDGER; it is not a lower bound on the actual graph's error and does not rule out cancellation or a redesigned source.

The infinite-history fixed-point mean m_infinity satisfies

    ||m_infinity-m_3||_2 <= A^4/(1-A)sqrt(D).

Therefore S_3 already approximates that limiting mean at the same certified grade. Reentry does not buy a better limiting-mean exponent, and the limiting mean is not by itself the endpoint distribution.

## 5. Why the actual raw mean graph cannot substitute for a covariance-matched force

There is a sharp same-g one-dimensional separator. Let g(x)=A x, and use the common-carrier source S_k=B_k+E_k. If its outer nodes are (omega_j,t_j),

    B_k(z,G)=A[a_Q z+b_Q G],
    b_Q=sum_j omega_j sqrt(1-t_j^2).

The positive quadrature's Hermite-degree-two certificate gives sum omega_j t_j^2<=1/3+epsilon. Since sqrt(1-t^2)>=1-t^2,

    b_Q>=2/3-epsilon.

For epsilon<=1/12 this is >=7/12. The remainder energy certificate gives

    E_z Var(S_k|z) >= A^2 b_Q^2-C Lambda_k A^3.

For the genuine force, the first-order term F1=A integral_0^infinity exp(-s)X_s ds has conditional variance A^2/4. Indeed its unconditional variance is A^2/2 and its conditional mean is A z/2. Coherent contraction gives

    Var(F_k|z)=A^2/4+O(A^3),

constant in z in this linear example. Hence, in the declared small-parameter window,

    ||Var(S_k|z)-Var(F_k|z)||_(L2_z) >= c A^2.          (3)

Already (7/12)^2-1/4=13/144>0 supplies a fixed positive leading gap. As the actual positive rule tends to the uniform integral, b_Q tends to pi/4 and the leading gap is pi^2/16-1/4. The lower bound above does not need that limiting argument.

Thus even an excellent conditional mean certificate and complete first/curl/energy ports can coexist with an actual order-A^2 covariance mismatch. Inserting S_k as a RAW covariance-matched tail would reintroduce the order-three weak debt; the shrinking-buffer LAW join avoids this by using S_k solely as the mean input and paying the WHOLE independent covariance services.

A completed mean LAW also cannot fix this issue by carrier subtraction. For example, with m=0, both G and cos(A)G+sin(A)H have exactly law N(0,I), while their residuals after subtracting the same known G have respectively zero and 2(1-cos(A)) covariance. Both residual firsts are O(A), and both vanish at A=0. A marginal LAW does not specify the needed residual coupling.

## 6. Minimal additional port that changes the forcing term

An actual covariance match is already present in the bulk. The remaining A^4/w term is the third centered Stein-field defect of the true conditional force. Its one-energy bound is O(A^3 sqrt(D)); the terminal A-Lipschitz factor and two Gaussian integrations by parts give A^4/v. The improved m_k certificate shrinks neither that field nor the complete-bank first.

One sufficient next port is a genuine finite positive conditional skew/third-cumulant correction, usable at the SAME retained (z,y) with buffer u comparable to w, which lowers this consumer error to O(A^5/u^(3/2)sqrt(D)), preserves the necessary shifted-current correlations, and has actual baseline/residual/caller first and full replay costs that close the preceding source class. A tensor expectation alone or a terminal mean LAW is insufficient.

Conditionally on THAT proved new producer, replacing only A^4/w by A^5/w^(3/2) in (2) permits w=A and grade 7/2. The prefix, near-RAW endpoint, and existing variable mean/Gram service terms all then have exponent 7/2. This is a prospective exponent ledger, NOT a claimed rank-three construction. The independent rank-three worker is checking an existing native skew packet against this exact contract.

If instead only the short-prefix term is improved, the unchanged near-RAW and Gaussianization bills cap this ledger at 10/3, balancing A^3 sqrt(w) and A^4/w at w=A^(2/3), provided the improved prefix reaches that grade.

Beyond either improvement, an arbitrary-order conclusion still needs explicit further moment/cumulant consumers, shrinking-buffer service refinements, and a count/guard recurrence. Generic C2 smoothing alone does not assign a new exponent: its actual source, bias, derivative/first constants, and all VALUE replays must enter the same ledger.

## 7. Complete reentry costs and guards

For the new variable-buffer mean service on S_k, let N_B,k(u), N_E,k(u) be the literal native occurrence counts, INCLUDING their internal mean/pair/filter/response/origin/restoration recurrences. Then

    Q_M,k(u) <= Q_cap,k(u)
                +N_B,k(u) Q_B,k
                +N_E,k(u)(Q_k+Q_B,k)+Q_ext,k(u).

All origins are caller-specific complete records. A changed source argument requires a full new S_k replay; no source-private root, covariance return, or discarded primal is cached across arguments.

If Q_H,j and Q_K,j are the fully expanded pinned service costs at bulk node j,

    Q_(k+1) <= sum_bulk[Q_M,k(u_j)+Q_H,j+Q_K,j+2]
               +N_near[Q_FQ+2]+Q_ext.

This is preferable to replaying Q_k on the endpoint branch. With d_M,k,j,d_H,j,d_K,j the COMPLETE respective private dimensions, after common-carrier alignment,

    d_(k+1)=D+sum_bulk[d_M,k,j+d_H,j+d_K,j]+3D N_near.

Each d_M,k,j includes a fresh full d_k tape at every reentered S_k occurrence. In the near term F_Q has 2D roots, plus its fresh conditional-Y root after sharing the outer carrier, giving 3D per node. All outer and native node counts and dimensions divided by D are public-log polynomials at any FIXED depth. The polynomial degree/constants and Lambda_k can grow with k; no depth-uniform or sublinear-order count is asserted.

At the actual normalization, require

    r_k,rho_k=Lambda_k A/sqrt(u) <= the imported order-four radius,
    Lambda_k r_k^2/sqrt(A) <= its imported self-reserve radius,
    mixed-K radius <= C Lambda A/u^(1/3) <= its imported order-five radius,
    u>=A^3 and all declared positive root gaps.

For w=A^alpha, the substantive mean/source radii scale as A^(1-alpha/2) and A^(3/2-alpha). The actual unscaled residual first is bounded by Lambda[A+A^(3/2)/sqrt(u)]. Keeping it O(Lambda A) requires alpha<=1; native admission alone could allow more, but that would no longer close the same source port. At alpha=4/5 the radius powers are 3/5 and 7/10; at alpha=1 both are 1/2. Mixed-K is respectively 11/15 and 2/3. Actual Lambda_k, dimensions and counts must be substituted anew at every depth.

Primitive precisions are assigned over the complete frozen graph, with each local absolute floor divided by its actual downstream norm/profile and the total sum bounded by the requested absolute budget. Include all retained caller, mode, shield, finite clock, scalar root and rotation factors. These finite inverse scales affect logarithmic precision, not uncharged query replication. An inherited target floor in e_k appears in (2) multiplied by A; new replay floors must still be paid independently.

## 8. Honest comparison

The imported general-C2 F_Q result already supplied order A^3 conditional force mean, including m3 at that grade by target comparison. The pinned join genuinely improves canonical m3 to grade16/5. The present note adds reusable finite-depth typing, covariance restoration, and an explicit stable recurrence, but NO new accuracy exponent. The existing fixed-order native constructions remain the actual implementation engine; merely reentering them does not establish an any-order force-mean or endpoint accuracy theorem.

## Source pins

- /workspace/shared/law-only-variance-join-20261005/MANIFEST.json, SHA256 673bdfcc6dc54c00ea380c942039add5ecda9f32686e6af0398243923dc6f3a7.
- /workspace/shared/law-only-variance-join-20261005/LAW-ONLY-SHRINKING-BUFFER-JOIN.md, SHA256 40891db38221e3b02f2e657a4ef93c3e4fd06c69b29cc03e7db9c5761fe2fc0a.
- /workspace/shared/collective-nonlinear-redesign-20261005/outer-resolvent/REENTRY-AND-PRIOR-RESULT-COMPARISON.md.
- /workspace/shared/recovery-20261004/cost/endpoint-stein-quadrature-20261004/positive-law-repair/higher-cumulant-gate/THIRD-ORDER-NESTED-FORCE-MEAN-WITH-C2-DECOUPLING.md.
