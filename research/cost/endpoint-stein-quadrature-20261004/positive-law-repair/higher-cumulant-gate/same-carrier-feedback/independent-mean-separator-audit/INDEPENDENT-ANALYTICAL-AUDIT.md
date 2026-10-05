# Independent analytical audit of two mean-gate separators

2026-10-05. Reviewer: independent single-history-current analyst.

## Exact source pins

- `../COHERENT-SHIFT-SEPARATES-UNSHIFTED-M3-FEEDBACK.md`, SHA256 `cb04bb4e1cd3886ed378bb45593ea855b800e1cd1e532916c5f238cd3141e285`.
- `../GRADIENT-PORTS-DO-NOT-CONTROL-A-CONTRACTED-TRACE.md`, SHA256 `e6050915d148b7248e97bd26647f5948f14ef50d47dbd365d9c668b18790faa1`.
- The shared coherent-shift construction was checked against `../SAME-CARRIER-P3-MEAN-AND-FULL-COVARIANCE-GATE.md`, SHA256 `a63c38f14206ec278f7c12a56542968a7d1675af6f550610e4539e500a653eca`.

## Verdict

PASS for both explicitly stated separator claims. Neither result is a counterexample to the actual same-g, single-history covariance expansion. Neither supplies a finite fourth-order mean consumer.

## A. Unshifted small-feedback formula

1. The gradient family is anchored and has Hessian between zero and c(2+epsilon)A I<A I. It remains in the original C2 class; in fact it is smooth.
2. The common-path limit v.H -> d is strengthened to L4 by the fourth-moment law of large numbers for independent coordinate OU functionals. The rank-one remainder is bounded and O(D^(-1/2)). The corresponding v.E2/A limit is L2 with bounded L4 norm by the same coordinate Taylor calculation and finite Gaussian moments.
3. The scalar k is positive: simultaneous sign reversal annihilates the log-cosh contribution, and Gaussian regression gives int_0^1 r dr E[N tanh N]=(1/2)E[N tanh N]. Hence C0<0, while h is strictly negative everywhere.
4. The diagonal-force contribution to v.R_D is O_L2(A^3). In detail H_i=O_L4(A), E2_i=O_L4(A^2), so the scalar b remainder is O_L2(A^3) before multiplying the outer coefficient A and summing with 1/sqrt(D). The result is A^4 sqrt(D)=A^3 in this counterfamily.
5. The rank-one tanh expansion has an O_L2(A^3) quadratic remainder. Its leading A^2 term is exactly c[sech^2(U1)-sech^2(U1-d)]e(U). The L4/L2 bounds above justify multiplying the limiting factors.
6. Conditioning does not assume a limit law for the full force: L2 contraction preserves the approximation. The scalar OU path U is independent of the orthogonal endpoint coordinates, so the limiting conditional function is exactly the stated F(v.x).
7. F is bounded and nonzero. Its second bracket is strictly negative and its first bracket vanishes only at u=d/2. The Gaussian R1 multipliers 1/(n+1) are positive, so R1 F is nonzero. Its scalar subspace is preserved in every D.
8. Projection onto the unit vector v and contraction of the vanishing error prove the full-vector lower bound c_* A^2. This exceeds A^4 sqrt(D)=A^3 by 1/A, beyond every fixed public-log multiplier.

Thus replacing the shifted Hessian by Dg(x) in this SMALL-FEEDBACK current is invalid. This argument does not refute the differently scaled single-I covariance correction.

## B. Gradient/first/one-energy port insufficiency for a trace

1. The inner cutoff cancels the gradient at zero exactly. The outer cutoff gives a globally smooth compactly supported potential. The construction is a genuine gradient, hence curl zero.
2. Under x=sqrt(D)z, the potential equals A^2 D^(3/2) f_D(z). The radial support and cutoff derivative bounds give a dimension-independent operator bound on D^2 f_D. Hence Lip(E_D)<=C A^2 sqrt(D)=C A.
3. On the Gaussian annulus, the gradient is the stated second-chaos vector. Direct orthogonality gives

       ||(E_0)_1||2^2/A^4 =18+2(D-1)=2D+16,
       sum_(i>1)||(E_0)_i||2^2/A^4=4(D-1),
       ||E_0||2^2=A^4(6D+12), E E_0=0.

4. Direct coordinate differentiation gives Delta E_0=A^2(2D+4)e1.
5. The chi-square tails outside [D/2,3D/2] are exponentially small; the cutoff field and the uncut gradient have polynomial growth. Consequently their L2 difference is a fixed polynomial in D times A^2 exp(-cD). This proves the claimed energy estimate after absorbing finitely many small dimensions.
6. Gaussian integration by parts gives E Delta E_D=E[(|X|^2-D)E_D]. Cauchy-Schwarz bounds the difference from the uncut expectation by sqrt(2D)||E_D-E_0||2, still exponentially small. No pointwise fourth derivative estimate is needed for this tail comparison.
7. R1 preserves Gaussian expectations. For W=A^2I the current's L2 norm is at least its mean norm, giving c A^4 D. At A=D^(-1/2) this is c A^2, versus the proposed A^3 allowance.
8. W is indeed a preserved product of two constant symmetric Hessians, but E is independently chosen. If those Hessians were the Hessians of the actual source g(x)=Ax, its nonlinear chord would not equal this E. This separation is decisive: the example proves the listed PORTS INSUFFICIENT, not the actual SAME-g conditional expansion false.

## Consequence for the open mean gate

A generic trace/Hodge estimate derived solely from energy, first, curl, zero, and a preserved operator-bounded Hessian product cannot close the fourth-order mean gate. A valid theorem needs an additional same-source identity or a stronger typed current hypothesis, and then its own finite original-VALUE consumer. The exact centered/resummed single-history identity remains a distinct analytical comparison, not an order-four reduction.
