# A matrix-free original-VALUE Gaussian backbone for canonical m3

2026-10-05. Source-qualified, fixed-order component theorem. Independent audit and a separate analytical scope obstruction accompany this packet. The general C2 m3 gate, fourth-order endpoint join, and all-order recurrence are not claimed.

## 0. Result

There is no need to reconstruct a D by D Hessian or execute its matrix square root to resum the anisotropic second-substitution Gaussian backbone. Two **known scalar Gaussian rows** factor the entire conditional history covariance. Three original gradient VALUES replace the linear actions, and a fourth retains the terminal nonlinearity:

    v=x/2+N/2,
    w=x/4+N/2+M/4,
    T_g(x,N,M)=x-g(v)+g(g(w)),
    f_g(x,N,M)=g(T_g(x,N,M)).                         (0.1)

N,M are independent standard D-Gaussians. The N in v and w is the SAME root. This graph is defined for every original anchored C2 convex-gradient source

    g=grad U, g(0)=0, 0<=Dg<=A I, 0<A<=1/2.

Its useful approximation theorem is source-qualified. Let K be any fixed symmetric analytical witness with 0<=K<=AI, r(x)=g(x)-Kx, k=||K||, and let L be a certified Lipschitz constant of r. One may always take L=A. Neither K, Kx, a Hessian action, nor an expectation is a producer leaf. The algorithm does not need the entries, eigenbasis, or action of K.

Define the following FIVE explicit centered-Gaussian residual norms, with each occurrence of Z standard:

    e0 = ||r(Z)||2,
    eY = ||r((I-K+K^2/2)^(1/2) Z)||2,
    ev = ||r(Z/sqrt(2))||2,
    ew = ||r(sqrt(3/8) Z)||2,
    eKw= ||r(sqrt(3/8) K Z)||2.

The square root above specifies an analytical Gaussian law; it is not computed. Put

    B_res=A[eY+(k+L)e0+ev+A ew+eKw].                 (0.2)

For the genuine-history canonical target m3=R1 psi2, the finite positive original-VALUE mean program described below obeys

    ||W2(Law(M_*(Z)|Z),N(m3(Z),I))||_(L2 gamma_Z)
       <= B_res + C_T delta A sqrt(D)
          + Lambda A^4 sqrt(D) + restored absolute floors,   (0.3)

where C_T=1+A/sqrt(2)+A^2 sqrt(3/8), delta is the positive outer Hermite-rule error, and all native fixed-order gradient/near-gradient mean-compiler guards remain imposed. No smallness of A sqrt(D) is required. The result is integrated over the original standard Gaussian Z, not uniform in Z.

With delta<=A^3 and B_res<=C_res A^4 sqrt(D), this is a guarded order-four canonical m3 mean-law component. Its original VALUE count is

    Q_certificate_setup+Q_captured
      +N_out N_B+5 N_out N_E+Q_known/numerical/replay,        (0.4)

at actual normalized compiler radii. N_out=O(log^2(1/delta)); N_B,N_E are fully expanded native occurrence counts, including every later replay. Raw F and E own 3D private Gaussian coordinates; B owns D. No D-query matrix recovery or hidden Hessian oracle appears.

For EVERY matrix quadratic g(x)=Kx, B_res=0 and the source mean is exactly canonical m3. This includes the balanced A/3,2A/3 quadratic defeating every scalar backbone by A^2(2-A)sqrt(D)/(36 sqrt(10)). A genuinely nonlinear anisotropic class with an unbounded residual and a sample-free O(A^4) certificate is given in Section 7.

## 1. Genuine OU-history comparison with an analytical matrix witness

On the stationary standard OU history, Cov(X_t,X_u)=exp(-|t-u|)I, define at every genuine future time t

    H1,t=int_0^infinity e^(-u) X_(t+u) du,
    F1,t=int_0^infinity e^(-u) g(X_(t+u)) du,
    Delta1,t=F1,t-K H1,t.

Since r is globally Lipschitz and r(0)=0, these are square integrable. They are comparison variables only, never continuous-path VALUE leaves. Stationarity and Minkowski give ||Delta1,t||2<=e0. The coupled Gaussian

    Y_t=X_t-KH1,t

has exact marginal covariance Sigma_Y=I-K+K^2/2: Cov(X_t,H1,t)=I/2 and Cov(H1,t)=I/2. Hence

    ||r(X_t-F1,t)||2 <= eY+L e0.

This does not assert independence of Y_t, Delta1,t or the retained endpoint. The actual second substitution is

    F2=int_0^infinity e^(-t) g(X_t-F1,t) dt
      =K H1-K^2 H2+R2,
    H1=int e^(-t) X_t dt,
    H2=int t e^(-t) X_t dt,
    R2=-K int e^(-t) Delta1,t dt
          +int e^(-t) r(X_t-F1,t) dt.

The convolution producing H2 uses the same genuine future ancestry. Fubini is justified by linear growth and Gaussian moments. Consequently

    ||R2||2 <= eY+(k+L)e0.                            (1.1)

Both Dg and K lie in [0,AI], so -AI<=Dr<=AI and L=A is always valid, without any Hessian modulus. A sharper certified L may be used.

## 2. Factor the joint conditional history, not a source covariance matrix

The exact conditional Gaussian block at X0=x is

    E[H1|x]=x/2, E[H2|x]=x/4,
    Cov((H1,H2)|x)=[[1/4,1/4],[1/4,5/16]] tensor I.

Its known scalar factorization is precisely

    (H1,H2)|X0=x =_law (v,w)
    v=x/2+N/2, w=x/4+N/2+M/4.                         (2.1)

Thus the complete matrix backbone at the terminal argument has law

    S_K=x-Kv+K^2w.                                    (2.2)

Its conditional mean and covariance are, only for analysis,

    E[S_K|x]=(I-K/2+K^2/4)x,
    Cov(S_K|x)=K^2/4-K^3/2+5K^4/16.                  (2.3)

No nonscalar matrix square root is executed. The two correlated scalar rows in (2.1) already encode every orientation and cross term when K actions are available through a truly quadratic original g. Keeping a single N in both v and w is indispensable.

Let psi_K(x)=E[g(S_K)|x]. Applying A-Lipschitz g on the genuine coupling and then conditional Jensen gives

    ||psi2-psi_K||_(L2 gamma) <= A[eY+(k+L)e0].         (2.4)

Only the LAW of the entire Gaussian pair was replaced by fresh N,M. This does not preserve a previously exposed old-history observer.

## 3. Replace linear actions by original VALUES and price the residual

Set psi_*(x)=E_(N,M) f_g(x,N,M) for the literal graph (0.1). Couple it to (2.2) using the same x,N,M. The elementary identity

    g(g(w))-K^2w
       =[g(g(w))-g(Kw)]+r(Kw)

gives

    |T_g-S_K| <= |r(v)|+A|r(w)|+|r(Kw)|.

The term g(Kw) is used only in this proof, never as an additional executing leaf. When x is standard Gaussian independently of N,M,

    v~N(0,I/2), w~N(0,3I/8), Kw~N(0,3K^2/8).

Their correlation does not invalidate Minkowski. Therefore

    ||psi_*-psi_K||2 <= A[ev+A ew+eKw].                (3.1)

Combining (2.4), (3.1), and the L2(gamma) contraction of R1 proves

    ||R1 psi_*-m3||2 <= B_res.                         (3.2)

This is a residual-qualified comparison of the original finite nonlinear graph, not the false assertion g(g(w))=K^2w for arbitrary sources. There is no expected-action oracle, empirical matrix approximation, Hessian-vector producer, or covariance reconstruction. All five norms are source qualifications, not callable producer expectations.

## 4. Positive outer quadrature, own mean, and exact quadratic test

Freeze the previously audited positive interior rule (w_i,t_i), of mass one, exact first moment 1/2, and Hermite operator error delta. Let c_i=sqrt(1-t_i^2), beta=sum_i w_i c_i. At retained standard Gaussian Z draw independent D-roots G,N,M; each is shared across ALL nodes of its own occurrence. Execute

    x_i=t_i Z+c_i G,
    v_i=x_i/2+N/2,
    w_i^arg=x_i/4+N/2+M/4,
    T_i=x_i-g(v_i)+g(g(w_i^arg)),
    F_*=sum_i w_i g(T_i),
    B=sum_i w_i g(x_i), E=F_*-B.                       (4.1)

The superscript on w_i^arg distinguishes the vector argument from the scalar positive weight. F_* has exactly 4N_out original VALUES; E has at most 5N_out before complete-key sharing. Its exact own mean is Q psi_*(Z), because each x_i has the required OU conditional marginal. No independence across nodes is needed for this identity.

At independent standard x,N,M,

    ||T_g||2 <= [1+A/sqrt(2)+A^2 sqrt(3/8)]sqrt(D).

Anchoring and the A-Lipschitz property imply ||psi_*||2<=A C_T sqrt(D), and hence

    ||Q psi_*-R1 psi_*||2 <= delta A C_T sqrt(D).       (4.2)

For g(x)=Kx, put L_K=I-K/2+K^2/4. The ACTUAL shared-root raw program is exactly

    F_*=(1/2)K L_K Z + beta K L_K G
              +(-K^2/2+K^3/2)N+(K^3/4)M.             (4.3)

Thus its conditional mean is the exact canonical m3(Z)=(1/2)K L_K Z, and its actual raw covariance is

    beta^2 K^2 L_K^2 + (1/4)K^4(I-K)^2+(1/16)K^6.     (4.4)

Every matrix here is a polynomial in the SAME K; simultaneous diagonalization is legitimate in this quadratic calculation only. In particular beta^2, not sum_i w_i^2 c_i^2, records all shared-G cross-node products. The covariance in (4.4) is the raw own-source covariance, not the covariance of the full canonical outer history. The mean consumer needs the exact own mean, not that unsupported covariance identification.

As a separate quadratic specialization, the exact unit-buffer law could also be executed by

    H+(1/2)g(Z)-(1/4)g(g(Z))+(1/8)g(g(g(Z))),

with shared nested VALUES and independent H~N(0,I). It uses three VALUES and is exactly N(m3(Z),I). The general source theorem uses (4.1), preserving the full terminal nonlinear response to the Gaussian backbone; it does not extend this special quadratic polynomial identity to arbitrary g.

## 5. Complete source ports: C2 firsts, curl, amplitudes, callers and zeros

Let, at each actual node,

    H_i=Dg(T_i), V_i=Dg(v_i),
    J_i=Dg(g(w_i^arg)), W_i=Dg(w_i^arg), H_i^0=Dg(x_i).

All are symmetric and in [0,AI]. They occur in the analysis and in requested first/adjoint sweeps at their recorded VALUE sites; the producer never returns a derivative action. Define

    P_i=I-V_i/2+J_i W_i/4,
    Q_i=-V_i/2+J_i W_i/2,
    R_i=J_i W_i/4.

The actual first paths are

    D_G E=sum_i w_i c_i(H_i P_i-H_i^0),
    D_N E=sum_i w_i H_i Q_i,
    D_M E=sum_i w_i H_i R_i,
    D_Z E=sum_i w_i t_i(H_i P_i-H_i^0).                (5.1)

No two Hessians at different points have been commuted. Let

    eta=A/2+A^2/4, q=A/2+A^2/2, r0=A^2/4.

Since H_i-H_i^0 is a symmetric difference of matrices in [0,AI], it has norm at most A. Equations (5.1) give the explicit raw bounds

    G-first(E)<=A beta(1+eta),
    N-first(E)<=A q,
    M-first(E)<=A r0,
    full private first(E)<=A sqrt(beta^2(1+eta)^2+q^2+r0^2),
    caller-Z first(E)<=A(1+eta)/2.                    (5.2)

The baseline is a true gradient in G with private first at most A beta; its analytical potential sum_i (w_i/c_i)U(t_iZ+c_iG) need not be queried.

For the square lift P_G^*E=(E,0,0), use curl to mean the operator norm of its Jacobian minus its transpose. The leading G block sum_i w_i c_i(H_i-H_i^0) is symmetric. The remainder has skew norm at most 2 beta A eta. The N,M off-diagonal row has norm at most A sqrt(q^2+r0^2). Hence

    curl(P_G^*E)<=2 beta A eta+A sqrt(q^2+r0^2)=O(A^2). (5.3)

This replaces the scalar backbone's exact G-block symmetry with a legitimate small-commutator bound. No Hessian continuity or third derivative is required. For A<=1/2 and beta<=sqrt(3)/2, a safe normalized bound under half variance shares is private first <=2A and curl<=3A^2. Thus declaring ell_E=2A and a_seed=2A safely gives curl<=ell_E a_seed, without division by a measured first or residual norm.

For any fixed finite p>=2 let kappa_p satisfy ||N(0,I_D)||p<=kappa_p sqrt(D). Given Z=z, each v_i has mean t_i z/2 and covariance (c_i^2+1)I/4; each w_i^arg has mean t_i z/4 and covariance (c_i^2+5)I/16. The literal displacement bound gives

    ||E||_(Lp|z)
      <=A^2[(1/4+A/8)|z|
                  +kappa_p(1/sqrt(2)+A sqrt(3/8))sqrt(D)].    (5.4)

No energy estimate is differentiated to obtain the firsts. The actual caller origins are

    B0(z)=sum_i w_i g(t_i z),
    E0(z)=sum_i w_i {g(t_i z-g(t_i z/2)+g(g(t_i z/4)))
                                                       -g(t_i z)}.

They cost N_out and 5N_out original VALUES before exact-key sharing and obey

    |B0(z)|<=A|z|/2,
    |E0(z)|<=A^2(1/4+A/8)|z|.                         (5.5)

Subtract and restore these ACTUAL executions. Their source-anchor paths retain every dependence on z and exterior labels. The anchored E-E0 keeps (5.2)'s private firsts and (5.3)'s curl; its caller-Z first is at most A(1+eta), twice the raw bound. Its energy is bounded by (5.4)+(5.5). The anchored baseline's caller first is at most A.

At total Z=G=N=M=0 all descendants, including g(w), g(g(w)), T and g(T), are zero coherently from g(0)=0. The captured origins are not declared zero at nonzero Z. At the zero source g=0 the entire raw source vanishes for every root; the completed mean service supplies its literal known unit-Gaussian row. No normalization divides by B_res or an energy that can vanish.

## 6. The finite positive mean consumer and full executable ledger

Use the same pinned fixed-order gradient/near-gradient own-mean compilers and recorded-coisometry adapter as the audited scalar packet, with independent COMPLETE banks conditional on Z and all actual exterior labels. For variance shares v_B=v_E=1/2 the literal form is

    M_B=B0+sqrt(v_B) N_k((B-B0)/sqrt(v_B)),
    M_E=E0+sqrt(v_E) P_G N^v(P_G^*(E-E0)/sqrt(v_E)),
    M_*=M_B+M_E.

Here the square lift now has private dimension 3D, not 2D. Use gradient order k>=4, rho_B=sqrt(2)A beta, ell_E=2A, a_seed=2A, padding mu=A. Impose rho_B<=the actual native gradient radius threshold, ell_E<=1/4, and every native curl, active-dimension, finite-clock/filter, precision, mode, and caller guard. For this safe declared ell_E, A<=1/8 is sufficient only for the ell_E guard; it is not a waiver of the other native guards. Less conservative actual radii may be used only after their guards are explicitly verified.

The imported near-gradient law-error row is

    Lambda e_E {ell_E(a_seed+mu)
                   +ell_E^3(1+mu^(-1/2))}+absolute floors.

Its bracket is O(A^2) and the actual anchored energy is O(A^2)(|z|+sqrt(D)), so this row is O(Lambda A^4)(|z|+sqrt(D)). The gradient order k>=4 branch has its own O(Lambda A^4 sqrt(D)) allowance. The independent completed targets add to N(Q psi_*(Z),I), by conditional product coupling only after their COMPLETE programs return. Gaussian equal-covariance means then differ in W2 by their Euclidean distance; (3.2) and (4.2) establish (0.3).

This imports existing finite consumers, not an oracle for a conditional mean. No raw root, evaluated force, private source record, training root, or old noise is later appended as an observer. No output is read as a strong mean/covariance statistic or stripped of its executed Gaussian buffer.

Every output is an ordinary finite Gaussian pushforward. Vector subtraction and polynomial coefficients are deterministic arithmetic, not signed probability measures. All positive shares, fills and covariance buffers remain those of the imported consumers. The source-zero rows are actual rows in those complete programs, not coupling inventions.

The safe original-VALUE bill is exactly (0.4): B costs N_out, E costs 5N_out per occurrence, and N_B,N_E include every source replay, mean bank, filter, clock, pair and mark in the fully expanded native programs. B0 and E0 are actual caller captures with N_out and 5N_out costs before complete-key sharing; their own repeats and any finite-mode original anchors must be charged too. A raw F_* or E owns 3D Gaussian coordinates G,N,M; a raw B owns D. All compiler Gaussian roots, fills, clocks, arithmetic, coefficient reads and scalar roots are additional ledger entries inherited without deletion. The only new scalar coefficients are 1/2 and 1/4; no source-dependent matrix coefficient or matrix square root is read or formed. Original leaf evaluations have vector-arithmetic cost O(D), apart from the original source's own evaluation cost; this is not a dimension-independent arithmetic claim.

A changed input or numerical version replays the WHOLE VALUE graph. The middle VALUE g(g(w)) owns g(w) as an original VALUE ancestor, not as an unpriced local Hessian action. Requested first/adjoint sweeps use one original HVP at each recorded original VALUE site per sweep, with actual product chain rules (5.1). A discarded primal pays complete replay. No HVP is differentiated. All quadrature/source-certificate/guard/share/schedule/mode/tolerance versions are fixed before such a sweep; actual exterior source parameters and anchors remain live through their recorded caller graph. The analytical witness K is absent from that graph, so no derivative of a fitted K is being stopped.

For nodewise original VALUE errors nu_v,nu_w,nu_b,nu_f at the respective g(v),g(w),g(g(w)),g(T) leaves, with errors measured at their actually requested sites, a valid terminal absolute allowance is

    nu_f+A nu_v+A nu_b+A^2 nu_w.                     (6.1)

This prices displacement of both nested descendants. A uniform nu gives (1+A)^2 nu for F_*, [(1+A)^2+1]nu for E, and twice that for the anchored E-E0 before exact-key reuse. B and B-B0 have nu and 2nu allowances. Errors in x,N,M have the absolute Lipschitz multipliers from (5.1)-(5.2); weight/node/Gaussian-row/arithmetic/clock errors, finite-mode residuals and replay errors retain their existing absolute floors. Never divide any floor by an RMS residual, A-dependent measured norm or vanishing energy. Coherent exact-key zero reuse is necessary for a literal numerical zero.

At fixed grade and with declared/certified source parameters, the original VALUE inverse-A exponent remains zero: multiply the existing public-log occurrence count by 5N_out. Certificate fitting or validation, if performed, adds all its original VALUES, Gaussian roots, arithmetic and rigorous error accounting in Q_certificate_setup; it is not hidden in this exponent claim.

## 7. Sample-free nonlinear anisotropic class beyond scalar backbones

The five residual expectations can be certified without executing K or any expectation. Suppose the source has an analytical witness 0<=K<=AI and known scalar envelope

    |r(x)|<=e_core+L(|x|-R)_+,
    R=sqrt(D)+T, T>=0.

Every Gaussian linear map in Section 0 is an operator contraction: Sigma_Y has eigenvalues 1-u+u^2/2<=1 for u in [0,A], the v,w scales are below one, and sqrt(3/8)||K||<=1. Consequently each norm is at most

    e=e_core+L sqrt(J_D(1,R))
      <=e_core+L sqrt(2) exp(-T^2/4),                 (7.1)
    J_D(1,R)=E[(|Z|-R)_+^2].

The first quantity is an explicit one-dimensional chi integral; the final bound follows from P(|Z|>=sqrt(D)+t)<=exp(-t^2/2) and integrating the survival function. This uses pointwise radial-envelope monotonicity under contractions, not the false general claim that shrinking a Gaussian lowers an arbitrary residual norm. With the always-valid k,L<=A,

    B_res<=A(3+3A)e.                                 (7.2)

The certificate inputs e_core,L,R can be derived from source structure. No matrix entries or unknown Gaussian expectations must be recovered at runtime. Finite black-box VALUES alone are not claimed to discover or certify this envelope for an arbitrary source.

For a concrete nonquadratic anisotropic family, choose any fixed orthogonal orientation and a matrix K with eigenvalues in [A/3,2A/3]. In D>=2 require at least two distinct eigenvalues; taking equal multiplicities of A/3,2A/3 is the exact earlier scalar obstruction when D is even. Let d=A/3,

    T=sqrt(8 log(1/A)), R=sqrt(D)+T,
    h(s)=0                              for s<=0,
         (5/2)s^4-3s^5+s^6              for 0<s<1,
         s-1/2                          for s>=1,
    g(x)=Kx+d h(|x|-R)x/|x|, g(0)=0.                 (7.3)

The nonlinear term vanishes near zero. It is the gradient of d int_0^|x| h(q-R)dq. Its Hessian is PSD, with radial eigenvalue d h'(|x|-R) and tangential eigenvalue d h(|x|-R)/|x|, both between zero and d. Thus U is C2 (indeed smoother here), 0<=Dg<=AI, r is d-Lipschitz, and |r(x)|<=d(|x|-R)_+. The generally anisotropic K does not commute with the radial Hessian; the first/curl proof did not require it to do so.

All five norms are therefore <=sqrt(2) A^3/3. Since k<=2A/3 and L=d=A/3, the SHARPER combined certificate is

    B_res<= (sqrt(2)/3)(3+2A) A^4,                    (7.4)

for every D>=1, hence at most that constant times A^4 sqrt(D). Original certificate setup needs zero VALUE queries: it is the displayed source-algebra and Gaussian-concentration calculation. The executing graph is unchanged when the eigenbasis of K is unknown. Supplying a dense K inside a source implementation may itself affect that original oracle's arithmetic cost; our algorithm neither reads nor reconstructs it.

The residual around K is unbounded. For rho>=R+1,

    g(rho n)=rho(K+dI)n-d(R+1/2)n.

When K has two distinct eigenvalues, sup_x|g(x)-lambda x| is infinite for EVERY scalar lambda. More importantly, the included anisotropic quadratic subfamily r=0 is exactly handled by the new graph, while every scalar backbone has the prior sharp order-A^2 sqrt(D) gap on the balanced A/3,2A/3 source. Thus this adds sources excluded by the old scalar-backbone order-four certificate, not merely a new matrix formula with no executable action implementation.

The new and old sufficient residual certificates need not be ordered on every scalar-near source; (0.1) incurs three additional explicit Gaussian residual terms. The theorem claims the concrete strict gain just demonstrated, not domination of every old bound.

## 8. Scope boundaries and conditional normalization

The graph (0.1) is executable for any anchored C2 convex gradient, but (0.2) is useful only when its actual source has a sufficiently small certified residual about SOME linear witness. Globally Lipschitz r alone generally gives e_j=O(A sqrt(D)), hence only an O(A^2 sqrt(D)) comparison. Nothing here supplies general C2 order-four closure.

For the actual normalized conditional source

    f(y)=s[g(a+s y)-g(a)], alpha=A s^2,

the witness is K_f=s^2K and the residual is

    r_f(y)=s[r(a+s y)-r(a)].

Apply the theorem with its FIVE actual translated Gaussian residual norms: r_f(Z), r_f((I-K_f+K_f^2/2)^(1/2)Z), r_f(Z/sqrt(2)), r_f(sqrt(3/8)Z), and r_f(sqrt(3/8)K_f Z). Original unshifted envelopes do not imply the needed grade uniformly in a,s. The executable f VALUE and its caller captures must retain and price their actual original-g implementations. In particular a caller in (7.3)'s shell need not inherit the small original Gaussian bulk residual.

No automatic data-driven witness selection, small-residual discovery, generic conditional normalization, strong conditional mean, raw covariance closure, retained-private-root law, reverse-OU endpoint join, growing-order recurrence or sublinear all-order exponent is established. This is a complete matrix-free source-qualified fixed-order canonical m3 mean component relative to the pinned finite mean compilers.

## 8.1 An explicit smooth obstruction to removing the residual qualification

The accompanying independent proof `independent-audit/GENERAL-C2-BIAS-OBSTRUCTION.md` supplies a quantitative boundary for this PARTICULAR executable graph. In D=1 let

    g_A(x)=A(x+sin x)/2, 0<A<=1/2.

This is an anchored smooth convex gradient with derivative in [0,A]. For m_*=R1 psi_* the actual infinite-outer-rule target mismatch obeys

    ||m_*-m3||2 >= c A^2-3.04 A^3,
    c=0.009194688258588945... >0.                     (8.1)

The constant is half the absolute first-Hermite coefficient of

    q(x)=h'(x)[R1 h(x)-E h(x/2+N/2)], h=(x+sin x)/2,
    <q,Z>=(1/4)[(1/2)e^(-1/2)+2e^(-1)-1/2
                 -(3/2)e^(-2)-(1/4)e^(-1/4)-(3/4)e^(-5/4)].

The proof retains the actual history and bounds the full Taylor remainder in L2 by 3.04 A^3; it is not a formal expansion or Monte Carlo estimate. Since R1 Z=Z/2, (8.1) is a direct quantitative lower bound. For A<=c/6.08 it is at least (c/2)A^2. A finite delta<=A^3 outer rule and the guarded own-mean consumer subtract only their order-four allowances from this lower bound.

The source of the failure is specific: nonlinear evaluation of the averaged Gaussian history, E h(H1)|X0, need not equal the average of its nonlinear evaluations, E int e^(-t)h(X_t)dt|X0. Linear sources have no such debt. This obstructs unrestricted closure by (0.1), while leaving its proved anisotropic residual-qualified class intact and leaving other nonlinear corrections open.

## 9. Provenance

Incoming scalar Gaussian-RMS manifest:

    ../gaussian-rms-backbone/MANIFEST.json
    SHA256 321cd1922f7b8a32b3ea76f138ca3b463a4cfb016c9b7d79056a6b822041acd0.

The pinned scalar theorem is SHA256 18640f9c6c567acba4fe830c7b0d64394924887ee4335a7a49f57f817435dccf, and the unchanged bounded-consumer source is SHA256 798d20b87124e7069e465ce4e7bcf92c5af3734d86c8a62aeb288b5fec8eaef3. Its independent audit is SHA256 49d135ba6692c54012acd9d0d94f603c0015fca98b5d1be1164457ff45d1213a and pins the actual LOW30 gradient/near-gradient compilers, complete branch composition, coisometry adapter, caller origins, scalar positive quadrature and original VALUE replays.

The earlier exact quadratic VALUE amplifier in ../../../POSITIVE-SECOND-ORDER-LAW-AND-QUADRATIC-AMPLIFIER.md motivates preserving the original polynomial action graph. Unlike that special quadratic proof, the present theorem separately prices nonlinear action-replacement error with (0.2), and proves actual noncommuting source ports. Its new proof does not use the older amplifier's nonlinear extension, which remains unproved.
