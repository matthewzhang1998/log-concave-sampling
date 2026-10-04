# Reconstructed whole-source OU orientation stability

Publication copy: nonmathematical context and/or local paths were sanitized. Original and public SHA-256 values are recorded in `INVENTORY.json`; source/audit pins below identify their historical versions, not these edited bytes.

## Provenance

Reconstructed mathematical exposition dated 2026-10-04. These are NEW bytes, not a claim to recover the old file byte-for-byte. The historical source SHA was4b2b72a4da96f944b1981de748a4db75b22e22458f20fd6035f7b1f16e5494a6; R's historical independent audit SHA was2cb198c0b5db4dcd2ee508e1ff4c84db1ddd0e54df916a3175845d684b546b80. The1080-check diagnostic is being reconstructed separately. The mathematical statement and proof below retain the reviewed normalization and scope.

## Statement

Let f:R^n->R^n be the SAME finite square source, with centered Gaussian energy e, first at most A and curl C=Df-Df* bounded in operator norm by kappa. Define

    O_f=Cov(f)-Sym B_f,
    B_f=integral_0^1 E[J_u²]du,
    J_u(X)=E Df(sqrt(u)X+sqrt(1-u)Z).

Let P_t=exp(-t N) be the Gaussian OU semigroup on the ENTIRE source input array. Then

    ||O_f-O_(P_t f)||HS <= kappa e eta(t),
    ||O_f-O_(P_t f)||op <= A kappa eta(t),
    eta(t)=sup_(k>=1) [1-exp(-2kt)]/sqrt(k)
          <=min(1,sqrt(2t)).                              (1)

This is an analytical target-restoration estimate. It neither executes P_t f as a differentiable source nor permits independent heating of the primitive Hessians of f's graph. All aliases must remain inside one complete common-input substitution.

## Exact Gaussian coefficient convention

For the degree-k vector Hermite coefficient f_(i;a1,...,ak), symmetric in its k input indices, define

    (T_k f)_(i;a1,...,ak)
      =(1/k) sum_l f_(a_l;i,a1,...,omit(a_l),...,ak),
    T_0=0,       Aop=I-E-T.

The Gaussian forward-transpose identity is

    O_f=Sym E[(f-Ef) outer Aop f].                         (2)

It follows by Hermite isometry: the integral defining B_f gives the coefficient k²(k-1)!/k=k! multiplying the averaged output/input swap. T is a bounded self-adjoint contraction on the input-symmetric coefficient space. It commutes with OU.

The degree(k-1) curl coefficient is

    C_(ij;beta)=k[f_(i;j,beta)-f_(j;i,beta)].

Gaussian divergence delta_j=x_j-partial_j creates one Hermite input index. Symmetrizing it with beta gives the exact identity

    delta(Curl f_k)=k Aop_k f_k.                          (3)

The denominator is k, NOT k+1. For k=1 this reads delta(F-F*)=(F-F*)x, as required.

For a vector field U in degree k-1,

    ||delta U||2<=sqrt(k)||U||2.                          (4)

Indeed the input squared norm is(k-1)! times its coefficient norm, while the output is k! times the squared norm of its symmetrization. Symmetrization is an orthogonal contraction. No dimension factor appears.

For each fixed output row u, equations(3)--(4) imply

    ||u* Aop_k f_k||2<=k^(-1/2)||u*C_(k-1)||2.

Orthogonality across degrees therefore yields

    ||u* Aop(I-P_(2t))f||2²
       <=eta(t)² E|u*C|²
       <=eta(t)² kappa²|u|².                             (5)

The elementary inequality1-exp(-x)<=min(1,sqrt(x)) proves the last bound on eta in(1). More generally the proof only needs E[CC*]<=kappa_row² I; the same result holds with kappa_row=sqrt(||E CC*||op).

## One-energy conclusion

OU self-adjointness and its commutation with Aop give

    O_f-O_(P_t f)
      =Sym E[(f-Ef) outer Aop(I-P_(2t))f].                 (6)

The first factor, as a map from Gaussian L2 to the output, has Hilbert--Schmidt norm e. The second has output-row operator norm at most eta(t)kappa by(5). The Hilbert--Schmidt/operator product bound gives the first inequality of(1). Gaussian Poincare gives Cov(f)<=A²I, so using the ordinary output-row norm of the first factor gives the operator bound. Symmetrization does not increase either norm.

Finite Hermite sums prove the identities algebraically. Gaussian Sobolev chaos expansion and L2 closure extend them to the actual Lipschitz source; no derivative of a saved Hessian is used. The argument applies conditionally on fixed exterior labels with the corresponding conditional energy profile. A coisometric physical projection cannot increase the stated norms.

## Native-use boundary

P_t f evaluates every original node, origin, twin bank and repeated query on the SAME substituted Gaussian input, before taking its analytical expectation. It preserves the whole finite source and its energy cancellations. If h=sqrt(2t), equation(1) supplies an O(h kappa e) constant-matrix restoration tolerance. Inside a fixed positive covariance gap this can be consumed by ordinary covariance transport without differentiating that error.

A product of separately averaged primitive Hessians generally has a different query genealogy and is not identified by this theorem. Neither a finite VALUE executor, its complete root/public/caller returns, nor the true-orientation adjoint ordering is supplied by this analytical lemma.
