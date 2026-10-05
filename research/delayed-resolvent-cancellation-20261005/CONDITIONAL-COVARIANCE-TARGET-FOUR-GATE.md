# Conditional covariance debt: a genuine source-level target-four gate

## Result and scope

The proposed coefficient, sign, constants and delayed-quadratic defect are correct. For the genuine anchored C-infinity source

    f(x)=(x+eta sin x)/(1+eta),   eta=1/2,   g(x)=A f(x),

the outer-resolvent first-chaos difference between the common-noise and true-OU history terminal sources is

    E[X R1(psi_common-psi_true)(X)] = -c_eta A^3 + remainder,
    c_(1/2) > 7/1728,
    |remainder| <= A^4/6.                                  (1)

Consequently, for 0<A<=1/100,

    |E[X R1(psi_common-psi_true)(X)]| > A^3/500.              (2)

There is also an explicit four-node finite positive source graph pair satisfying the same lower bound. No limiting quadrature, numerical lower-bound certification or unexecuted expectation oracle is needed for that finite example.

This is a conditional-covariance block lower bound. The comparison here omits the extra nonlinear inner-history term in the full m3 definition. It therefore does NOT establish a lower bound for the whole m3 graph or for every target-four construction. It does establish that matching conditional force means alone cannot make this terminal block accurate to order A^4.

## 1. Admissible original source

Take

    U(x)=A/(1+eta) [x^2/2+eta(1-cos x)].

Then U is C-infinity, g=U', g(0)=0, and

    A(1-eta)/(1+eta) <= g'(x) <= A.

In particular the original source meets the exact anchored convex-gradient class. At eta=1/2, f is 1-Lipschitz, |f(x)|<=|x|, and sup|f'''|=1/3.

## 2. Conditional Gaussian histories

Let X be standard Gaussian. In correlation coordinates r in [0,1], write X_r for the stationary OU future conditioned on X, with

    E[X_r|X]=rX,
    v_r=Var(X_r|X)=1-r^2,
    d_T(r,s)=Cov(X_r,X_s|X)=min(r,s)/max(r,s)-rs.

The value at r=s=0 is immaterial to all integrals. Define the common-noise rows

    Y_r=rX+sqrt(1-r^2)N,
    d_C(r,s)=sqrt(1-r^2)sqrt(1-s^2),

where N is independent standard Gaussian. Define H_T=integral_0^1 f(X_r)dr and H_C=integral_0^1 f(Y_r)dr. They have identical conditional means given X because every one-time conditional Gaussian marginal agrees. Also d_C>=d_T>=0: this follows immediately for r<=s by squaring the inequality sqrt(1-r^2)sqrt(1-s^2)>=r(1-s^2)/s.

Let psi_C(x)=E[g(X-AH_C)|X=x], and define psi_T similarly. No factor A is included inside H_C or H_T.

## 3. Coefficient identity and its sign

For any integrable scalar F, the outer resolvent has the exact first-chaos identity

    E[X R1F(X)]=(1/2)E[XF(X)].

Taylor expansion of A f(X-AH) gives terms A f(X), -A^2 f'(X)H and (A^3/2)f''(X)H^2. The first two terms cancel in the C-versus-T comparison, because the conditional means match. Since f''(X)=-eta sin(X)/(1+eta), the cubic coefficient is negative and

    c_eta = eta/[4(1+eta)] E[X sin(X)(H_C^2-H_T^2)].         (3)

Define

    M(a)=(1/2)[(1+a)exp(-(1+a)^2/2)
                 +(1-a)exp(-(1-a)^2/2)].

Thus M(a)=E[X sin(X)cos(aX)]. For jointly Gaussian rows of conditional means rX,sX, conditional variances v_r,v_s and covariance d,

    E[Y_r Y_s|X]=rsX^2+d,
    E[Y_r sin(Y_s)|X]
       =rX exp(-v_s/2)sin(sX)+d exp(-v_s/2)cos(sX),
    E[sin(Y_r)sin(Y_s)|X]
       =(1/2)exp(-(v_r+v_s)/2)
          [exp(d)cos((r-s)X)-exp(-d)cos((r+s)X)].

Substitution into (3) gives exactly

    c_eta=eta/[4(1+eta)^3] integral integral {
       (d_C-d_T)[exp(-1/2)
          +eta(exp(-v_s/2)M(s)+exp(-v_r/2)M(r))]
       +(eta^2/2)exp(-(v_r+v_s)/2)
          [(exp(d_C)-exp(d_T))M(r-s)
             -(exp(-d_C)-exp(-d_T))M(r+s)] } dr ds.         (4)

All factors of 1/2 in this identity are accounted for: one from the quadratic Taylor coefficient and one from the outer first-chaos multiplier.

## 4. Analytic lower bound

For r in [0,1], M(r)>0. The derivative with respect to covariance d of the conditional sin-product expectation is E[cos(Y_r)cos(Y_s)|X]. After multiplying by X sin(X) and averaging, its absolute value is at most

    E|X sin(X)| <= E|X| <1.

This is the needed rigorous interpretation of the sin-product derivative bound; no unsupported pointwise bound on X sin(X) is used. Consequently the braces in (4) are at least

    (d_C-d_T)[exp(-1/2)-eta^2].                            (5)

The covariance mass has the exact integral

    integral integral(d_C-d_T)drds = pi^2/16-1/4.

Indeed integral sqrt(1-r^2)dr=pi/4, while integral integral min(r,s)/max(r,s)=1/2 and integral integral rs=1/4. At eta=1/2 the prefactor in (4) is 1/27. Using exp(-1/2)>3/5 and pi>3 gives

    c_(1/2) > (1/27)(7/20)(5/16)=7/1728.                  (6)

Numerically the actual coefficient is approximately 0.013275858; this numerical value is only a diagnostic. The proof uses (6).

## 5. Uniform fourth-order remainder

For either history, positive averaging and |f(x)|<=|x| give ||H||4<=3^(1/4). Holder gives E[|X||H|^3]<=3. The third-derivative Taylor remainder of A f(X-AH) therefore has first-chaos absolute expectation at most

    (A^4/6)(1/3) E[|X||H|^3] <= A^4/6.

The outer first-chaos factor reduces this to A^4/12 for each history. Comparing two histories costs at most A^4/6, proving (1). For A<=1/100,

    7/1728-A/6 >=7/1728-1/600 >1/500,

which proves (2). The same proof works for any finite positive mass-one clock average.

## 6. A fully finite four-node original-VALUE certificate

Use four correlations r=(1,3,5,7)/8 and equal weights 1/4. One history uses their exact conditional joint OU Gaussian law, and the other uses the single common N in all four rows. Both sources execute their four original f-values (equivalently original g-values with the known scalar factor) followed by the terminal original g. Their conditional means and first exponential moment agree exactly.

Here is the exact executable original-VALUE readset. For the common-noise source, draw one standard N and set U_j=r_j x+sqrt(1-r_j^2)N. For the true-OU source, draw independent standard Z_1,...,Z_4, set W_4=r_4 x+sqrt(1-r_4^2)Z_4, and recurse in the chronological order j=3,2,1:

    W_j=(r_j/r_(j+1))W_(j+1)
          +sqrt(1-(r_j/r_(j+1))^2)Z_j.

Then execute respectively

    b_C=(1/4)sum_(j=1)^4 g(U_j),   T_C=g(x-b_C),
    b_T=(1/4)sum_(j=1)^4 g(W_j),   T_T=g(x-b_T).

Each graph makes exactly five original g VALUE calls, with every terminal descendant retained. There is no division by A and no executed f or derivative oracle. All four known scalar ratios lie strictly between zero and one. Before an outer root is added, the private Gaussian dimensions are one for the common-noise graph and four for the true-OU graph.

The theorem's exact R1 first-chaos comparison also has a fully finite outer realization: choose the single positive outer node of correlation 1/2, replace x by z/2+sqrt(3)G/2 with a new independent standard G, and execute the same complete five-VALUE graph. Its conditional mean is P_(1/2)psi, whose first-chaos multiplier is exactly the same 1/2 as R1. Thus the lower bound applies to these finite outer graphs as well. This does not assert that P_(1/2) approximates R1 on other chaoses. The original five-VALUE bill is unchanged; one extra Gaussian root is used.

For this finite pair, the total covariance mass difference equals

    I_4=[(sqrt(63)+sqrt(55)+sqrt(39)+sqrt(15))/32]^2-127/420.

The second term is the exact average conditional OU covariance. Simple rational root lower bounds give sqrt(63)>7.9, sqrt(55)>7.4, sqrt(39)>6.2 and sqrt(15)>3.8. Hence

    I_4 > (253/320)^2-127/420
        =693949/2150400 >5/16,

with excess 21949/2150400. Thus (5)-(6) and the remainder proof apply verbatim to this actual finite source pair. The lower bound A^3/500 at A<=1/100 is therefore genuinely finite and source-valid.

This finite graph's quadrature does not approximate the full resolvent to arbitrary accuracy, and no such claim is needed. It exhibits a concrete matched-conditional-mean pair with a cubic covariance debt. Higher-order Gauss quadrature diagnostics separately show convergence to the continuous coefficient.

## 7. Exact delayed-graph quadratic defect

For the uncorrected delayed graph in DELAYED-CONDITIONAL-NOISE-THEOREM.md, let g(x)=Ax, w=sqrt(A), q=1-w, and impose the exact first moment sum a_j r_j=1/2. Before the outer resolvent its own mean is

    [A-(1+w^2)A^2/2+(wq^2/2)A^3]x.

The genuine m3 pre-resolvent mean is [A-A^2/2+A^3/4]x. Therefore the exact outer defect is

    A^3[-3/8+wq^2/4]x.                                   (7)

This is uniformly of cubic order. Indeed max_(0<=w<=1) w(1-w)^2/4=1/27, so the absolute coefficient in (7) is at least 73/216. The optional quadratic-exact correction in the companion delayed-source note removes this particular defect; it does not by itself remove the general nonlinear conditional-covariance debt above.

## 8. Reproducible diagnostics

check_covariance_debt.py computes (4) using positive Gauss-Legendre rules on the square and after the smooth-triangle substitution r=su, includes the exact four-node example, and checks the rational lower-bound arithmetic. covariance_debt_checks.json records the results. At 256 triangle nodes per direction, c_eta is approximately 0.01327585825; the four-node coefficient is approximately 0.01188997892. These are numerical checks of the analytically proved signs and scales, not replacements for the rational lower bounds.
