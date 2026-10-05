# A smooth one-dimensional obstruction to general matrix-free m3 closure

2026-10-05. Independent analytical counterexample for the four-VALUE source in `../MATRIX-FREE-GAUSSIAN-BACKBONE-M3.md`. This does not contradict its residual-qualified theorem.

## Result

In dimension D=1 take, for 0<A<=1/2,

    h(x)=(x+sin x)/2,
    g_A(x)=A h(x),
    U_A(x)=A x^2/4+(A/2)(1-cos x).

This is an original anchored smooth convex-gradient source: g_A(0)=0 and

    g_A'(x)=(A/2)(1+cos x) in [0,A].

Let m3 be the genuine canonical history target. Let m_* be the exact infinite-outer-rule mean of the finite matrix-free source,

    v=x/2+N/2,
    w=x/4+N/2+M/4,
    psi_*(x)=E[A h(x-A h(v)+A h(A h(w)))],
    m_*=R1 psi_*.

Then the actual target discrepancy has the explicit lower bound

    ||m_*-m3||_(L2 gamma)
       >= c A^2-3.04 A^3,
    c=0.009194688258588945... >0.                    (1)

In particular it is Theta(A^2) as A tends to zero, not O(A^4). The source is fixed in form as A varies; dimension, regularity, and convex-gradient admissibility do not deteriorate. Exactness on every matrix quadratic therefore does not imply general nonlinear closure.

For the finite positive outer rule with operator error delta<=A^3, replace the right side of (1) by

    c A^2-3.04 A^3-C_T A^4,
    C_T=1+A/sqrt(2)+A^2 sqrt(3/8).

A completed own-mean consumer with its proved O(A^4) law error likewise cannot remove this order-A^2 mean mismatch: the integrated Wasserstein triangle inequality subtracts only that O(A^4) allowance and its separately restored numerical floors.

## 1. Exact first nonlinear term of the canonical target

Use the stationary standard OU history with covariance exp(-|t-u|). At every genuine future time t define

    J_t=int_0^infinity e^(-u) h(X_(t+u)) du.

The actual first substitution is F1,t=A J_t, and

    F2=A int_0^infinity e^(-t) h(X_t-A J_t) dt.

Because h is anchored and 1-Lipschitz, |h(x)|<=|x|. By stationarity and Minkowski,

    ||J_t||p<=||Z||p,
    ||F2-A J_0||2<=A^2,
    ||F2||4<=A(1+A)3^(1/4).                         (2)

The second inequality compares h(X_t-A J_t) and h(X_t) on the SAME actual history. It does not resample any ancestor.

Also |h''|<=1/2. The scalar Taylor inequality, with its integral remainder, gives

    |A h(x-y)-A h(x)+A h'(x)y|<=(A/4)y^2.

Consequently, in joint L2 over the original stationary history,

    A h(X0-F2)
       =A h(X0)-A^2 h'(X0)J_0+R_can,
    ||R_can||2
       <=A^3[1+(sqrt(3)/4)(1+A)^2].                 (3)

Conditional Jensen cannot increase this error. Conditional on X0=x, the genuine future OU marginal gives

    E[J_0|X0=x]
       =int_0^1 E[h(tx+sqrt(1-t^2)Z)]dt
       =R1 h(x).

Hence

    psi2(x)=A h(x)-A^2 h'(x)R1 h(x)+r_can(x),
    ||r_can||2<=the bound in (3).                    (4)

## 2. Exact first nonlinear term of the finite source

When x,N,M are independent standard Gaussians,

    v~N(0,1/2), w~N(0,3/8).

Write its displacement as

    d_A=-A h(v)+A h(A h(w)).

The same anchor and Lipschitz properties give

    ||d_A+A h(v)||2<=A^2 sqrt(3/8),
    ||d_A||4<=A 3^(1/4)[1/sqrt(2)+A sqrt(3/8)].      (5)

Taylor's inequality now gives

    A h(x+d_A)=A h(x)-A^2 h'(x)h(v)+R_fin,
    ||R_fin||2<=A^3 {sqrt(3/8)
          +(sqrt(3)/4)[1/sqrt(2)+A sqrt(3/8)]^2}.    (6)

Conditioning only x and combining (4),(6) proves

    psi_*(x)-psi2(x)=A^2 q(x)+r_A(x),
    q(x)=h'(x)[R1 h(x)-E h(x/2+N/2)],
    ||r_A||2<=C(A)A^3,                              (7)

where

    C(A)=1+(sqrt(3)/4)(1+A)^2+sqrt(3/8)
            +(sqrt(3)/4)[1/sqrt(2)+A sqrt(3/8)]^2.

Every term is nondecreasing for A>=0. At A=1/2,

    C(1/2)=3.0312523067017936... <3.04.

The L2 contraction R1 therefore gives

    m_*-m3=A^2 R1 q+rho_A, ||rho_A||2<=3.04 A^3.     (8)

All remainders are absolute integrated bounds. This proof neither differentiates a Wasserstein estimate nor uses a formal, unpriced asymptotic expansion.

## 3. Explicit nonzero first Hermite coefficient

The linear parts of h cancel in q. Gaussian characteristic functions give

    q(x)=(1/4)(1+cos x)
         {int_0^1 exp(-(1-t^2)/2) sin(tx)dt
                                      -exp(-1/8)sin(x/2)}.    (9)

A quick nonzero check is

    q'(0)=(1/2)[1-exp(-1/2)-(1/2)exp(-1/8)]<0.

For the quantitative lower bound use E[Z sin(aZ)]=a exp(-a^2/2), together with cos z sin(tz)=[sin((t+1)z)+sin((t-1)z)]/2. Direct integration yields

    <q,Z>_(L2 gamma)
       =(1/4)[(1/2)e^(-1/2)+2e^(-1)-1/2
              -(3/2)e^(-2)-(1/4)e^(-1/4)-(3/4)e^(-5/4)]
       =-0.01838937651717789... .                   (10)

For completeness the integral part before the outer factor 1/4 is

    int_0^1 e^(-(1-t^2)/2)
      {t e^(-t^2/2)
        +[(t+1)e^(-(t+1)^2/2)+(t-1)e^(-(t-1)^2/2)]/2}dt
    =(1/2)e^(-1/2)+2e^(-1)-1/2-(3/2)e^(-2),

and the t=1/2 comparison term is

    (1/4)e^(-1/4)+(3/4)e^(-5/4).

R1 is self-adjoint on L2(gamma) and R1 Z=Z/2. Since ||Z||2=1,

    ||R1 q||2>=|<R1 q,Z>|=|<q,Z>|/2=c.

Combining this with (8) proves (1). In particular for 0<A<=c/6.08 the lower bound is at least (c/2)A^2. The upper bound O(A^2) follows already from (8) and q in L2, so the asymptotic order is genuinely Theta(A^2).

## 4. Meaning and limits

The leading discrepancy has a specific source: E h(H1)|X0 generally differs from the genuine exponential-time average E int e^(-t)h(X_t)dt|X0. These agree for linear h and differ for this smooth nonlinear h. Factoring the exact JOINT Gaussian law of (H1,H2) repairs anisotropic quadratic covariance, but it cannot exchange nonlinear evaluation with a time average.

This example does not refute another matrix-free finite graph, a controlled nonlinear correction, a larger finite mean construction, or any residual-qualified result. For K=A/2 its residual is (A/2)sin x, of size O(A) in the relevant Gaussian norms; it does not satisfy the present order-four residual certificate. It also does not claim that a source-dependent K action is necessary: the positive exact-quadratic graph remains valid without one.

## 5. Independent diagnostic evidence

`check_matrix_free_m3_independent.py` verifies (10) using a separately implemented 120-node Gauss-Hermite rule and 160-node Gauss-Legendre integral, then checks the displayed Taylor constant. The proof of nonzero error and the bound (1) are the analytical calculations above; the numerical quadrature is a diagnostic, not an unqualified rigorous-integration oracle. Exact rational Taylor enclosures additionally verify the sign in (10) and the rigorous decimal lower bound c>0.0091946882585889. The same checker independently tests covariance, original VALUE counts, noncommuting first/curl bounds, complete caller paths, adjoints, and numerical leaf-error propagation of the positive source.
