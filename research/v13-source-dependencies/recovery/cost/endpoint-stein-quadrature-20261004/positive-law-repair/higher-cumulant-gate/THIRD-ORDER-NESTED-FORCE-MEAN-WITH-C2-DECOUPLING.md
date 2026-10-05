# Third-order nested-force mean with uniform C2 decoupling

2026-10-04. Candidate finite mean branch for the exact two-clock reference. Independent audit requested. This is a conditional Gaussian mean-law source, not a strong mean statistic. Its covariance branch is a separate continuous-current construction.

## 1. Target and finite result

Assume g=grad U, g(0)=0, 0<=Dg<=A I, 0<A<=1/2. Let X_r be the analytical Gaussian Markov OU path anchored at X_1=Z, with covariance min(r,s)/max(r,s). Its exact second-substitution force is

    F_cont=int_0^1 g(X_r-I_r)dr,
    I_r=int_0^1 g(X_(r tau))d tau,
    m_cont(Z)=E[F_cont|Z].                              (1)

The exact-reference note proves W2(Z-F_cont,mu_U)<=A^3 sqrt(D) by stationary resolvent contraction. This statement does not execute the path.

Choose the already admitted positive quadrature Q_in=(v_j,tau_j) and Q_out=(w_i,r_i), each with uniform Hermite-moment error delta<=A^2. Their node counts are O(log^2(1/A)), all nodes lie strictly inside (0,1), and both first moments are exactly 1/2. Define the actual finite inner packet

    H_Q(x,H)=sum_j v_j g(tau_j x+sqrt(1-tau_j^2)H).

At exposed Z, draw fresh independent D-roots G,H, using the same G at every outer node and the same H at every inner occurrence. Execute

    x_i=r_i Z+sqrt(1-r_i^2)G,
    F_Q(Z,G,H)=sum_i w_i g(x_i-H_Q(x_i,H)).             (2)

The candidate intrinsic conditional-mean comparison is

    ||E_(G,H)F_Q-m_cont||_(L2(Z~gamma))
       <=C[A^3+delta A+delta A^2]sqrt(D).                (3)

Thus delta<=A^2 gives a uniform C2 third-order force-mean approximation with only n_out(n_in+1) private original g VALUES. The bound is obtained by a Gaussian decoupling estimate integrated over the outer heat; it is not a pointwise Hessian Taylor remainder.

The actual finite source also admits a caller-retained positive mean compiler. Under its stated normalized-radius/clock/padding guards, the resulting complete M(Z) targets N(E F_Q|Z,I) with error Lambda A^4(|Z|+sqrt(D)) plus floors. Combining this with (3) yields an integrated conditional W2 mean-law error Lambda A^3 sqrt(D). Every source-private G,H/filter/clock bank is integrated; Z and all earlier external callers are retained.

## 2. Literal conditional Markov ancestry

Fix 0<r<1 and expose Z=z. Put q=sqrt(1-r^2) and X_r=rz+qG. For each finite set of tau nodes, Gaussian Markov conditioning gives

    X_(r tau)=tau X_r+L_tau V,
    ||L_tau||=sqrt(1-tau^2),
    Cov(L_tau V,L_sigma V)
        =[min(tau,sigma)/max(tau,sigma)-tau sigma]I.      (4)

V is independent of Z and G. Indeed the residual X_(r tau)-tau X_r has zero covariance with both X_r and Z. This is the same inner ancestor X_r, not an independent force site. The inner residuals remain correlated with each other according to (4).

The continuous I_r is the L2 limit of finite positive averages of these nodes. On any interval away from tau=0 the path is L2-continuous; the omitted interval near zero costs its length times a uniform moment bound. Therefore the finite approximation argument does not posit an undefined Gaussian point at tau=0.

For either the continuous I_r or a finite positive inner Markov average, as a function J(x,V),

    Lip_x J<=A, and for the exact integral Lip_x J<=A/2,
    Lip_V J<=A,
    ||J(rz+qG,V)||p<=C_p A(|z|+sqrt(D)).                (5)

The finite approximations can be chosen with exact first moment 1/2; the weaker A bound is already enough for the analytic decoupling lemma when A<=1/2. The private first follows from the operator norms of the actual Gaussian rows L_tau and positive weights, with no root-count factor. The amplitude moment follows from each point's marginal mean tau r z and covariance at most I.

The cheap packet H_Q has the same bounds, D_x H_Q between 0 and A I/2, and exact conditional mean

    E_H H_Q(x,H)=v_Q(x),
    v(x)=int_0^1 P_tau g(x)d tau,
    ||v_Q-v||_(L2(gamma))<=delta A sqrt(D).              (6)

It has a different private-root law from the continuous inner path. That difference is priced next, rather than erased.

## 3. Uniform decoupling at an exposed outer caller

Apply the previously independently audited Gaussian conditional-decoupling lemma to the normalized near-identity map

    G-J(rz+qG,V)/q.

Its derivative with respect to G of the displacement J/q is D_xJ, bounded by A<=1/2. Its other-private first is at most A/q, its displacement energy is at most C A(|z|+sqrt(D))/q, and its G-Jacobian HS energy is at most A sqrt(D). The lemma requires the G-first to be below one; it does not require the other-private first or displacement energy to be small.

Combining the lemma's marginal-to-flow estimate and its flow-to-Euler estimate, and multiplying the state back by q, gives

    W2(Law(X_r-J(X_r,V)|Z=z),
                    Law(X_r-E_V J(X_r,V)|Z=z))
       <=C A^2(|z|+sqrt(D))/q.                          (7)

For the continuous J=I_r, apply this to finite Gaussian averages first and take their L2 limit. All constants and relevant first/energy bounds are uniform in the size of that private Gaussian bank. This is an analytical limit, not an executed unbounded graph.

Equation (7) is applied separately to the continuous Markov inner path, whose mean is v, and to H_Q, whose mean is v_Q. These are marginal couplings at the same exposed z; neither keeps the integrated private V/H, and no strong same-tape conditional-mean approximation is asserted.

The outer g is A-Lipschitz. Therefore the difference of its two conditional expected values costs at most A times (7). After squaring/integrating z~gamma and then using Minkowski in r,

    int_0^1 C A^3 sqrt(D)/sqrt(1-r^2) dr
        =(C pi/2)A^3 sqrt(D).                           (8)

The singular factor is integrable. No compiler is executed on the inflated field J/q, so no inverse-q source replication or radius guard has been hidden in (8). The analytical r=1 endpoint has measure zero and is handled by the integral limit.

To compare the two deterministic shifts, use Lip(g)<=A and (6). Unconditionally in Z,G, X_r is standard Gaussian at every r. Conditional Jensen hence gives

    ||E_G[g(X_r-v_Q(X_r))-g(X_r-v(X_r))|Z]||_(L2(Z))
       <=delta A^2 sqrt(D).                             (9)

Thus replacing the continuous inner Markov integral by the cheap shared-root packet changes the integrated force mean by at most C A^3 sqrt(D)+delta A^2 sqrt(D). This comparison uses the actual Gaussian caller law, not a false uniform-in-x quadrature estimate.

## 4. Outer quadrature remains a one-variable Mehler operator

Define only in analysis

    Psi_Q(x)=E_H g(x-H_Q(x,H)).

For each outer r,

    E_(G,H)[g(X_r-H_Q(X_r,H))|Z]=P_r Psi_Q(Z).

No r-dependent source or private source version remains inside Psi_Q. Since x,H are independent standard Gaussians and every inner point has standard Gaussian marginal,

    ||Psi_Q||_(L2(gamma))<=A(1+A)sqrt(D).                (10)

The admitted uniform Hermite multiplier quadrature therefore gives

    ||sum_i w_i P_ri Psi_Q-int_0^1 P_r Psi_Q dr||2
       <=delta A(1+A)sqrt(D).                           (11)

Combining (8)-(11) proves (3). In particular there is no assertion of pointwise clock analyticity of g, and no new two-dimensional tensor quadrature theorem is smuggled into the mean branch. The inner and outer quadratures are used only at their actual one-variable Gaussian conditional means.

## 5. Original-VALUE gradient plus near-gradient split

On the literal shared G,H record of (2), set

    B(Z,G)=sum_i w_i g(x_i), E=F_Q-B.

B is a genuine gradient in G: D_GB=sum_i w_i q_i Dg(x_i) is symmetric positive semidefinite. Its private first is at most A and centered energy at most C A sqrt(D). Its caller origin B(Z,0) is computed with n_out original VALUES.

For one outer node let B1=Dg(x_i-H_Q), B0=Dg(x_i), Jx=D_xH_Q, JH=D_HH_Q. Then

    D_G E_i=q_i[(B1-B0)-B1 Jx],
    D_H E_i=-B1 JH.                                    (12)

Both B1-B0 and Jx are symmetric. The first difference may be O(A), but contributes no G-G curl. Under the coisometry P_G=(I,0), the square lift P_G*E therefore has

    Lip E<=C A,
    ||Curl(P_G*E)||op<=C A^2,
    ||E||_(Lp(G,H)|Z)<=C_p A^2(|Z|+sqrt(D)).             (13)

The off-diagonal curl blocks are the actual D_HE_i, already O(A^2). All old inner packet values and their common H remain inside each E occurrence. No independent inner force is substituted.

The exact private origin E(Z,0,0) is executable and obeys |E(Z,0,0)|<=A^2|Z|/4 using both quadrature first moments. Subtract that recorded origin before invoking the near-gradient compiler; its anchored energy has the same profile as (13). The actual Z first of B and E is O(A), not the derivative of their smaller energy. All original physical callers and finite-mode anchors keep their actual first paths.

Use fixed shares vB=vE=1/2. Put beta_out=sum w_i q_i and beta_in=sum v_j sqrt(1-tau_j^2). The baseline's normalized gradient radius is sqrt(2)A beta_out and must satisfy the literal fixed-order gradient-mean radius guard. A sufficient explicit normalized E first bound is sqrt(2)A(beta_out+A), and its curl is at most sqrt(2)A^2(beta_out+beta_in). Impose sqrt(2)A(beta_out+A)<=1/4, padding mu=A, and every inherited finite-clock, caller and precision guard. Under these actual-radius conditions, independent COMPLETE mean banks give

    W2(Law(M(Z)|Z),N(E F_Q|Z,I))
       <=Lambda A^4(|Z|+sqrt(D))+absolute floors,        (14)

for baseline mean order at least four. The deterministic origins are added back exactly. Applying (3) and integrating Z~gamma gives the stated third-order Gaussian force-mean law. This law retains Z but integrates every bank/private G,H/filter/clock root. It does not permit reusing those roots as a later observer.

## 6. Complete original work, zeros, and restored versions

A raw B occurrence costs n_out private original VALUES. A raw E occurrence costs n_out(n_in+2): n_in inner nodes, one shifted terminal value, and one baseline per outer node. H_Q uses the same H at every inner node and for every outer node of that occurrence. Exact aliases may be reused only on that complete record. The source dimensions are D for B and 2D for E; the completed consumers have all their own replay/filter/clock dimensions.

Capture B(Z,0), E(Z,0,0), all original/conditional modes and g anchors under their complete finite keys before private sampling. A safe bill is

    Q_captured+n_out N_B+n_out(n_in+2)N_E
                    +known/numerical work.              (15)

Every changed raw compiler argument is a new full private graph, including its n_out n_in original ancestors. Requested first/adjoint sweeps use original HVPs at these recorded VALUE points, and a discarded primal pays its replay. No HVP is a producer leaf and no saved HVP is differentiated.

All raw zeros at Z=G=H=0 are literal under saved-anchor reuse. At fixed nonzero Z the private origins are the separately recorded caller-only values above, not zero by a Gaussian moment argument. The completed source-zero Gaussian carriers remain those of the admitted mean compilers, never noises from a law coupling.

The only coefficient singularity 1/q appeared in the analytical proof (7), where it was integrated in (8). The executed source (2) has no inverse-q multiplier. Its VALUE precision floors propagate through positive weight sums one, the original source first A, and the finite mean consumers' actual readouts/filters. All quadratures, node counts, source versions, origins, modes, shares, padding and numerical tolerances are fixed before differentiation. Finite physical mode/conditional mode tilts and numerical floors are restored exactly as in the existing conditional velocity service, at their actual original caller laws.

## 7. Join scope

This note supplies the previously missing nonlinear **mean branch** for the exact continuous two-clock reference, with a uniform C2 A^3 error and no inverse-heat query exponent. A complete buffered transition additionally needs the separately proved continuous conditional covariance compression, its positive VALUE reserve, and the cubic Gaussianization lemma. Those branches must use independent complete banks conditional on the same retained Z and all exterior callers.

No all-order recurrence is inferred. If the covariance branch is admitted, relative reverse-OU increments h=sqrt(Delta)/s<=1/2 can use the already specified positive reserve split and target an A^3 full buffered bridge. The zero-buffer terminal step remains separately priced. These are stated as integration obligations, not silently supplied by a Gaussian mean law.
