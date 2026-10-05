# The cheap mean source cannot supply the E2 covariance by polarization

2026-10-04. Companion to `SAME-CARRIER-P3-MEAN-AND-FULL-COVARIANCE-GATE.md`, frozen SHA256 a63c38f14206ec278f7c12a56542968a7d1675af6f550610e4539e500a653eca. Independent review requested. This separates one specific existing source from the required covariance target, not all possible finite law compilers.

## Result

The original cheap shared-root source, even with arbitrarily accurate positive clock quadrature and an exact downstream covariance action, does not generally approximate Cov(E2|Z) at the required A4 sqrt(D) scale. The counterfamily of Section 5 of the companion gives a covariance error Omega(A2) at A=D^(-1/2). Exact polarization cannot fix a source's incorrect genealogy.

## 1. The two scalar limits

Keep the companion's c=1/4, epsilon=1/2, b(t)=t+epsilon log cosh t, d=c epsilon E log cosh N, and h(u)=tanh(u−d)−tanh(u). For the true Markov force chord, it proved

    (v dot E2)/A = c[C_M + J_M] + o_L2(1),
    J_M = integral_0^1 h(U_r)dr,
    Cov(U_r,U_s)=min(r,s)/max(r,s), U1=Z0.            (1.1)

The continuum version of the cheap source instead uses

    U_r^C=r Z0+sqrt(1−r2)G0

with the SAME G0 for every r, and an independent common inner H for every inner node. The identical coordinate law-of-large-numbers calculation gives

    (v dot E_C)/A = c[C_C + J_C] + o_L2(1),
    J_C=integral_0^1 h(U_r^C)dr.                      (1.2)

The constant C_C may differ from C_M; it vanishes under conditional covariance. The inner v-projection still converges to d: every inner point tau X+sqrt(1−tau2)H has standard Gaussian marginal, so its b-mean is epsilon E log cosh N. No assertion that the inner paths have the same joint law is needed for this limit.

In both cases the conditional mean of J given Z0 is exactly

    integral_0^1 P_r h(Z0)dr.                         (1.3)

Thus this separator is invisible to the conditional-mean test used by the admitted cheap source.

## 2. A strict one-dimensional covariance gap

For 0<r<s<1,

    rho_C(r,s)=rs+sqrt(1−r2)sqrt(1−s2) > r/s=rho_M(r,s).

Indeed multiplying the difference by s reduces positivity to s sqrt(1−r2)>r sqrt(1−s2), equivalent to s>r. Both correlations are nonnegative.

Expand h in probabilists' Hermite polynomials, h=sum_(n>=0) a_n He_n with a_n=E[h(N)He_n(N)]/n!. Since h is bounded, this converges in L2 and the Gaussian covariance expansion gives

    Var(J_C)−Var(J_M)
      =sum_(n>=1) a_n2 n! integral_0^1 integral_0^1
                       [rho_C(r,s)^n−rho_M(r,s)^n]drds > 0. (2.1)

All summands are nonnegative. Positivity follows already from n=1 because a1=E h'(N)<0, proved in the companion. Direct integration gives

    integral integral rho_M =1/2,
    integral integral rho_C =1/4+pi2/16,
    Var(J_C)−Var(J_M)
        >= (E h'(N))2 (pi2−4)/16 >0.                 (2.2)

Their conditional means agree by (1.3), so the same positive gap equals

    E Var(J_C|Z0)−E Var(J_M|Z0).                      (2.3)

The full D-dimensional endpoints contain only additional orthogonal Gaussian coordinates, independent of the scalar path in (1.1)-(1.2); those coordinates do not alter the scalar conditional mean. L2 convergence in the companion implies convergence of the integrated conditional variances. Therefore

    E v*[Cov(E_C|Z_D)−Cov(E2|Z_D)]v
       =c2 A2 [Var(J_C)−Var(J_M)] + o(A2).            (2.4)

In particular the L2(Z;HS) norm of the covariance difference is Omega(A2). This exceeds A4 sqrt(D)=A3 by an inverse-A factor, and any fixed public-log allowance remains insufficient.

## 3. Relation to the actual finite quadrature

The continuum cheap source in (1.2) is only an analytical limit of the already defined finite common-G,H VALUE packet; it is not proposed as a producer. For the admitted positive quadratures with error tending to zero, their measures converge weakly to uniform Lebesgue measure by polynomial density and the uniform moment accuracy. For fixed standard Z0,G0, the map r -> h(rZ0+sqrt(1−r2)G0) is continuous on [0,1] and bounded. Hence the finite outer sum converges in L2 to J_C. Its conditional means converge likewise. The coordinate-average constants and inner mean d are preserved as above; the uniformly bounded positive-weight moment estimates make the joint dimension/clock limit valid for any such sequence of increasingly accurate rules. Equivalently one can first choose a fixed sufficiently accurate positive rule and then take D large; a fixed positive covariance gap already remains.

Thus improving this source's clock accuracy cannot change its common-root genealogy into the Markov genealogy. An action that exactly computes its own Cov(E_C) targets the wrong covariance. A constructive repair needs a new source/current that explicitly retains and corrects the genealogy, or a direct native covariance service with a proof of the exact continuous target.
