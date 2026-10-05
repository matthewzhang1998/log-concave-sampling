# A one-energy law bound for Gaussian covariance mixtures

2026-10-04. New analytical consumer, independent review requested. This supplies the covariance-mixture step of the orientation-free completed-kernel route. It does not itself supply the finite coefficient-mean compiler, the mixed K-current, or the same-carrier m3 mean.

## Result

Let X be a standard Gaussian private tape, and let C(X) be a symmetric PSD D-by-D matrix. Regard the symmetric matrices as a Hilbert space with Frobenius norm. Assume

    Lip_HS(C)<=L_C,
    e_C=||C−E C||_(L2;HS)<infinity.

For fixed nu>0, let Y conditional on X be N(0,nu I+C(X)), and put Sigma=nu I+E C. Then

    W2(Law(Y),N(0,Sigma))
        <= L_C e_C /(sqrt(2) nu^(3/2)).               (1)

The right side is one Hilbert-Schmidt energy. Covariance matching by itself is not enough; the actual Gaussian-to-HS Lipschitz bound is a hypothesis. A coherent high-dimensional random scalar variance times I does not satisfy the same dimension-free L_C hypothesis.

Applied to the finite j coefficient kernel B_Q from the frozen source below, with its owned coarse Gaussian X banks, this gives O(A4 sqrt(D)) at a fixed positive nu. Every bound here is uniform in the captured standard endpoint Z; no derivative of an output-law error is taken.

## 1. Matrix-space Stein kernel with one energy

Write F=C−E C as a vector in the Hilbert space Sym_D, using any fixed orthonormal basis. Let L be the Gaussian OU generator on the private X tape. Gaussian integration by parts supplies the matrix-space Stein kernel

    tau_F(F)=E[D_X(−L)^(-1)F (D_X F)* | F].          (2)

It satisfies E[F_a phi(F)]=E[sum_b tau_ab partial_b phi(F)] for smooth scalar phi. The bounded operator first of F and the Gaussian spectral-gap resolvent estimate imply

    ||tau_F||_(L2;HS)
       <=L_C ||D_X(−L)^(-1)F||_(L2;HS)
       <=L_C e_C.                                   (3)

No source derivative beyond D_X C is used. Gaussian Sobolev approximation supplies (2)-(3) for Lipschitz C. The conditional expectation in (2) is analytical, not an executed matrix oracle.

## 2. Posterior covariance fluctuation after the Gaussian convolution

Let S=nu I+C, and let phi_S(y) be its centered Gaussian density. Differentiation in a symmetric matrix direction E gives

    D_C phi_S(y)[E]
      =(1/2)phi_S(y) trace[E(S^(-1) y y* S^(-1)−S^(-1))]. (4)

Apply the Stein identity of Section 1 to the scalar function C -> phi_(nu I+C)(y). After division by the positive mixture density,

    E[F | Y=y]
      =(1/2) E[tau_F contracted with
         {S^(-1)y y* S^(-1)−S^(-1)} | Y=y].           (5)

The contraction takes place in Sym_D. For a fixed row T of tau_F, condition on C and set Y=S^(1/2)N with N standard. Its centered quadratic contraction is

    trace[T S^(-1/2)(NN*−I)S^(-1/2)].

Gaussian Wick isometry gives squared L2 norm

    2 ||Sym(S^(-1/2) T S^(-1/2))||HS2
                         <=2 nu^(-2)||T||HS2.       (6)

Sum over all output coordinates of F, use conditional Jensen in (5), and then (3). Thus

    ||E[C−E C | Y]||_(L2(Y);HS)
        <=||tau_F||_(L2;HS)/(sqrt(2)nu)
        <=L_C e_C/(sqrt(2)nu).                       (7)

This is where covariance uncertainty is reduced by a second small factor without a second sqrt(D). Separately bounding Gaussian quadratic entries or using a dimension-sized matrix norm would lose this fact.

The density calculation is justified first with bounded smooth C and smooth test cutoffs. S>=nu I controls the inverse factors, the displayed L2 bounds control the contractions, and Gaussian Sobolev/truncation limits give the stated case. In the finite-j application C is uniformly bounded in operator norm, so no upper-tail covariance issue occurs.

## 3. A self-contained Gaussian Stein-to-W2 estimate

Conditional Gaussian integration by parts shows that Y has the Stein matrix

    T_Y(y)=E[S | Y=y].

Consequently T_Y−Sigma=E[C−E C |Y]. For any centered vector Y with Stein matrix T_Y and any positive definite Sigma, one has

    W2(Law(Y),N(0,Sigma))
        <=|| (T_Y−Sigma) Sigma^(-1/2)||_(L2;HS).     (8)

For completeness, take an independent G~N(0,Sigma) and interpolate

    Y_t=tY+sqrt(1−t2)G, 0<=t<=1.

Stein and Gaussian integration by parts show that the derivative of a smooth test expectation is t E[(T_Y−Sigma):D2 phi(Y_t)]. One more integration by parts in the independent G realizes a continuity velocity

    v_t(Y_t)= t/sqrt(1−t2)
       E[(T_Y−Sigma) Sigma^(-1)G | Y_t].              (9)

Its L2 norm is at most t/sqrt(1−t2) times the right side of (8). The dynamic Wasserstein bound and integral_0^1 t/sqrt(1−t2)dt=1 prove (8). Endpoint limits follow from finite second moments. The coupling Gaussian G is only analytical and is never substituted for an actual program root.

Since Sigma>=nu I, (7)-(8) prove (1).

## 4. Why the actual finite B_Q mixture meets the hypothesis

Use `RESOLVENT-COVARIANCE-STABILITY-AND-FINITE-J-KERNEL.md`, SHA256 359cdfffa1ce4b7467d24ed24d6b8efdedb85a16c5de0b650d3b4ea37aa52cfd. Its j_Q is globally Lipschitz with L_j<=A(1+A/2), and

    B_Q(X)=sum_i w_i r_i P_(r_i) D j_Q(X),
    ||B_Q||op<=L_j/2.

For a fixed direction u, commute its derivative with D j_Q and perform one Gaussian integration by parts:

    D B_Q(X)[u]
      =sum_i (w_i r_i2/c_i)
         E[(D j_Q(r_iX+c_iG)u)G*], c_i=sqrt(1−r_i2).

Gaussian Bessel gives

    ||D B_Q(X)[u]||HS
        <=L_j |u| sum_i w_i r_i2/c_i <=C A|u|.      (10)

The dyadic positive rule has the same bounded coefficient sum already proved in the existing square-clock service. This holds for non-gradient j_Q and needs only its bounded first. It does not differentiate an original Hessian.

Let G_k be the owned independent outer coarse roots, X_k=q_kZ+sqrt(1−q_k2)G_k, with positive weights a_k summing to one. Set only as a conditional covariance target

    C_mix(G;Z)=sum_k a_k B_Q(X_k) B_Q(X_k)*.

The product rule, (10), and Cauchy–Schwarz over the complete G tape imply

    Lip_(G -> HS) C_mix <=C A2,
    0<=C_mix<=C A2 I,
    ||C_mix−E_G C_mix||_(L2(G);HS)<=C A2 sqrt(D).     (11)

All constants are uniform in captured Z. The positive weights prevent any clock-count factor: sum a_k |h_k| <=(sum a_k2)^(1/2)||h||<=||h||. Therefore (1) gives

    W2(Law(N(0,nu I+C_mix)|Z),
          N(0,nu I+E_G C_mix|Z))<=C A4 sqrt(D).      (12)

When a_k=2 v_k q_k, its mean is exactly the finite covariance target C_Q(Z) in the frozen coefficient-source theorem. No covariance mixture is sampled by forming a matrix in the eventual program: the finite completed-kernel construction must produce the mixture through its actual positive Gaussian-root outputs, with the required independent complete banks and variance shares.

## 5. Finite-program use and limits

Suppose a separate admitted original-VALUE mean compiler gives, conditional on captured X and fresh coarse p, an output with law close to N(B_Q(X)p/sqrt(nu),I), at its stated A4 error and actual caller/root ports. Independent complete banks with weights sqrt(nu a_k), each retaining its own p_k and X_k through that comparison, then have target conditional on all X_k

    N(0,nu I+sum_k a_k B_Q(X_k)B_Q(X_k)*).

The p_k integration creates the transpose orientation exactly. Equation (12) now replaces that covariance mixture by its mean at order four. Its Gaussian variance nu was explicitly budgeted in the completed means; no completed output was read as a strong coefficient value.

This paragraph is an interface composition, conditional on the separate finite mean/compiler proof. Its actual original-query count is inherited from that construction, including all first-coefficient filter and inner j_Q replays. This analytical consumer introduces no original g calls. All its owned coarse roots are part of the actual finite program and remain owned until the law comparison integrates them; the Gaussian G in Section 3 is never appended as an observer.

This fills only the mixture-law step. The mixed K-current, same-carrier m3 mean, original mode/numerical floors, and final complete order-four endpoint join remain separate gates.
